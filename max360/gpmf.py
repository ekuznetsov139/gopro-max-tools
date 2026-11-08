"""GoPro metadata (GPMF) from a .360 file: scaled sensor tracks, DISP, and the factory calibration.

Pure Python over PyAV, so no C parser has to be built.

The decompression in _decompress follows GPMF_Decompress and the code tables of
GPMF_bitstream.h from GoPro's gpmf-parser (https://github.com/gopro/gpmf-parser),
(C) Copyright 2017 GoPro Inc, licensed under Apache-2.0 or MIT.  This project is not
affiliated with GoPro.
"""
import struct

import av
import numpy as np

_TYPES = {'b': 'i1', 'B': 'u1', 's': '>i2', 'S': '>u2', 'l': '>i4', 'L': '>u4', 'f': '>f4', 'd': '>f8', 'j': '>i8', 'J': '>u8'}


def _klv(buf, off=0, end=None):
    """KLV walk: 4-byte key, 1-byte type, 1-byte struct size, 2-byte big-endian
    repeat, payload padded to 4 bytes.  Yields (key, type, size, repeat, offset)."""
    end = len(buf) if end is None else end
    while off + 8 <= end:
        key = buf[off:off + 4]
        if key == b'\0\0\0\0':
            return
        typ, size, rep = buf[off + 4], buf[off + 5], struct.unpack('>H', buf[off + 6:off + 8])[0]
        yield key.decode('latin1'), chr(typ) if typ else None, size, rep, off + 8
        off += 8 + ((size * rep + 3) & ~3)


def _values(buf, typ, size, rep, off):
    dt = np.dtype(_TYPES[typ])
    return np.frombuffer(buf, dt, count=rep * size // dt.itemsize, offset=off).astype(np.float64).reshape(rep, size // dt.itemsize)


def _metadata_packets(path):
    with av.open(path) as c:
        met = [s for s in c.streams if s.type == 'data' and s.metadata.get('handler_name', '').strip().endswith('MET')]
        if not met:
            raise ValueError(f"{path}: no GoPro metadata track")
        for packet in c.demux(met[0]):
            if packet.dts is not None:
                yield bytes(packet)


def _streams(buf):
    """Each STRM of a metadata packet, as a list of its KLV entries."""
    for k, _, s, r, o in _klv(buf):
        if k == 'DEVC':
            for k2, _, s2, r2, o2 in _klv(buf, o, o + s * r):
                if k2 == 'STRM':
                    yield list(_klv(buf, o2, o2 + s2 * r2))


def read_tracks(path, wanted=('CORI', 'GRAV')):
    """{fourcc: (N, K) float array, scaled} of sensor tracks; row k of CORI and GRAV
    is video frame k."""
    out = {w: [] for w in wanted}
    for buf in _metadata_packets(path):
        for entries in _streams(buf):
            first = {}
            for e in entries:
                first.setdefault(e[0], e)
            for w in wanted:
                if w in first:
                    scale = _values(buf, *first['SCAL'][1:]).ravel() if 'SCAL' in first else 1.0
                    out[w].append(_values(buf, *first[w][1:]) / scale)
    missing = [w for w, v in out.items() if not v]
    if missing:
        raise ValueError(f"{path}: no {missing} in the metadata track")
    return {w: np.concatenate(v) for w, v in out.items()}


# ---- Factory calibration, in moov/udta/GPMF -----------------------------------
def _box(f, start, end, want):
    """Payload range of the first box named `want` among the children in [start, end)."""
    off = start
    while off + 8 <= end:
        f.seek(off); hdr = f.read(8)
        if len(hdr) < 8:
            break
        size, typ = struct.unpack('>I4s', hdr); hs = 8
        if size == 1:
            size = struct.unpack('>Q', f.read(8))[0]; hs = 16
        elif size == 0:
            size = end - off
        if typ == want:
            return off + hs, off + size
        off += size
    return None


def read_calibration(path):
    """The factory calibration of docs/FORMAT.md, from moov/udta/GPMF: {block: {field: value}} with the
    blocks keyed by their DVID fourcc ('GEOC', 'FRNT', 'BACK', ...) or, failing that, by their DVNM name.
    Scalars come back as floats or ints, strings as str, longer arrays as lists."""
    with open(path, 'rb') as f:
        f.seek(0, 2); end = f.tell(); r = (0, end)
        for name in (b'moov', b'udta', b'GPMF'):
            r = _box(f, r[0], r[1], name)
            if r is None:
                raise ValueError(f"{path}: no {name.decode()} box; is this a .360 file?")
        f.seek(r[0]); buf = f.read(r[1] - r[0])
    out = {}
    for k, typ, size, rep, off in _klv(buf):
        if k != 'DEVC' or typ is not None:
            continue
        block, key = {}, None
        for k2, t2, s2, r2, o2 in _klv(buf, off, off + size * rep):
            pay = buf[o2:o2 + s2 * r2]
            if t2 in _TYPES:
                v = _values(buf, t2, s2, r2, o2).ravel(); v = [int(x) if t2 in 'bBsSlLjJ' else float(x) for x in v]; v = v[0] if len(v) == 1 else v
            elif t2 in ('c', 'F'):
                v = pay.decode('latin1', 'replace').rstrip('\0')
            else:
                continue
            if k2 == 'DVID' and t2 == 'F':
                key = v
            elif k2 == 'DVNM' and key is None:
                key = v
            block[k2] = v
        if key is not None:
            out[key] = block
    return out


def baseline(path):
    """Lens 0 relative to lens 1, cm: the negative of the header's SHFX/Y/Z (mm), as docs/FORMAT.md defines B."""
    g = read_calibration(path).get('GEOC')
    if not g or any(k not in g for k in ('SHFX', 'SHFY', 'SHFZ')):
        raise ValueError(f"{path}: no GEOC SHFX/Y/Z in the header")
    return np.array([-g['SHFX'], -g['SHFY'], -g['SHFZ']]) / 10.


# ---- Compressed payloads ('#') ----------------------------------------------
# Magnitude codes, as (length in bits, bits, magnitude); a sign bit follows each but m0.
_MAGNITUDES = [
    (1, 0b0, 0), (2, 0b10, 1), (4, 0b1100, 2), (5, 0b11011, 3), (5, 0b11101, 4), (6, 0b110100, 5), (6, 0b110101, 6),
    (6, 0b111110, 7), (7, 0b1110000, 8), (7, 0b1110011, 9), (7, 0b1111000, 10), (7, 0b1111001, 11), (7, 0b1111011, 12),
    (8, 0b11100100, 13), (8, 0b11100101, 14), (8, 0b11110100, 15), (9, 0b111000101, 16), (9, 0b111000110, 17),
    (9, 0b111101010, 18), (10, 0b1110001000, 19), (10, 0b1110001110, 20), (10, 0b1111010110, 21), (10, 0b1111111100, 22),
    (11, 0b11100010010, 23), (11, 0b11100011111, 24), (11, 0b11110101110, 25), (12, 0b111000100111, 26),
    (12, 0b111000111101, 27), (12, 0b111101011111, 28), (13, 0b1110001001101, 29), (13, 0b1110001111001, 30),
    (13, 0b1111010111101, 31), (14, 0b11100010011000, 32), (14, 0b11100011110000, 33), (14, 0b11110101111000, 34),
    (14, 0b11110101111001, 35), (15, 0b111000100110010, 36), (15, 0b111000100110011, 37), (15, 0b111000111100011, 38),
]
_ZERO_RUNS = [(7, 0b1111110, 16), (8, 0b11111110, 32), (9, 0b111111111, 64), (10, 0b1111111101, 128)]
_ESC, _END = 0b1110001111000100, 0b1110001111000101

_CODES = {(n, bits): ('zeros', count) for n, bits, count in _ZERO_RUNS}
_CODES.update({(n, bits): ('delta', m) for n, bits, m in _MAGNITUDES})
_MAX_CODE = max(n for n, _ in _CODES)


def _decompress_u8(data):
    """One channel of 8-bit unsigned samples from a '#' payload."""
    typ, size, count = chr(data[0]), data[1], struct.unpack('>H', data[2:4])[0]
    if typ != 'B' or size != 1:
        raise ValueError(f"unsupported compressed payload: type {typ!r}, sample size {size}")
    out = np.empty(count, np.uint8)
    last, quant = data[4], data[5]
    out[0], pos = last, 1
    bits = int.from_bytes(data[6:], 'big')      # 16-bit aligned from here, already so: 4 + 2 bytes in
    nbits = 8 * (len(data) - 6)
    at = 0

    def peek(n):
        return (bits >> max(nbits - at - n, 0)) & ((1 << n) - 1) if nbits - at >= n else (bits << (n - (nbits - at))) & ((1 << n) - 1)

    while True:
        word = peek(16)
        if word == _END:
            out[pos:] = last
            return out
        if word == _ESC:
            at += 16
            delta = peek(8)
            delta = delta - 256 if delta >= 128 else delta
            last = (last + delta * quant) & 0xff
            out[pos] = last
            pos += 1
            at += 8
            continue
        for n in range(1, _MAX_CODE + 1):
            code = _CODES.get((n, word >> (16 - n)))
            if code:
                break
        else:
            raise ValueError("invalid code in a compressed payload")
        at += n
        if code[0] == 'zeros':
            out[pos:pos + code[1]] = last
            pos += code[1]
        else:
            if code[1]:
                sign = -1 if peek(1) else 1
                at += 1
                last = (last + sign * code[1] * quant) & 0xff
            out[pos] = last
            pos += 1
        if pos > count:
            raise ValueError("compressed payload overruns its declared size")


def read_disp_rows(path):
    """(R, 640) uint8: the DISP rows exactly as stored, in file order.  Most callers
    want read_disp, which sorts out how rows relate to frames."""
    rows = []
    for buf in _metadata_packets(path):
        for entries in _streams(buf):
            for k, typ, size, rep, off in entries:
                if k == 'DISP':
                    if typ != '#':
                        raise ValueError(f"{path}: DISP has type {typ!r}, expected compressed")
                    rows.append(_decompress_u8(buf[off:off + size * rep]))
    if not rows:
        raise ValueError(f"{path}: no DISP in the metadata track; it exists only in 360 modes")
    return np.stack(rows)


def disp_rows_to_frames(rows, n_frames):
    """(F, 640): one DISP row per video frame, row k for frame k; F can be a little
    under n_frames, as the last frames of a clip may have no row.

    Real-time video stores one row per frame.  Time-lapse stores about two, in repeated
    pairs: frames 2m and 2m + 1 use rows 4m and 4m + 1, and rows 4m + 2 and 4m + 3 are
    copies of them.  See docs/FORMAT.md."""
    if len(rows) >= 1.5 * n_frames:
        picks = [(k // 2) * 4 + (k % 2) for k in range(n_frames)]
        return rows[[p for p in picks if p < len(rows)]]
    return rows[:n_frames]


def read_disp(path):
    """(F, 640) uint8: one DISP row per video frame; see disp_rows_to_frames."""
    with av.open(path) as c:
        n_frames = c.streams.video[0].frames
    return disp_rows_to_frames(read_disp_rows(path), n_frames)
