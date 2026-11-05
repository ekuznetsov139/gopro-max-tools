"""The Python package, with no footage needed.  python3 tests/test_python.py

With MAX360_TEST_CLIP=<a .360>, it also reads that clip."""
import fractions
import os
import struct
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import av
import numpy as np

from max360 import compose_cross, disp_rows_to_frames, gpmf, layout, Orientation, Reader

FAILS = []


def check(name, ok, detail=''):
    print(('PASS ' if ok else 'FAIL ') + name + (f'  {detail}' if detail != '' else ''))
    if not ok:
        FAILS.append(name)


# ---- DISP decompression -------------------------------------------------------
# GPMFCompress() from gpmf-write's output for 20 x 50, then 20 x 53, uint8.
PAYLOAD = bytes.fromhex('420100283201fc36e3c5000000000000')
row = gpmf._decompress_u8(PAYLOAD)
check('a zero run keeps the old value', (row[:20] == 50).all() and (row[20:] == 53).all(), row[18:22].tolist())


def encode(values):
    """The writer's logic, for building test payloads."""
    bits = ''
    run = 0
    mags = {m: f'{b:0{n}b}' for n, b, m in gpmf._MAGNITUDES}
    for a, b in zip(values[:-1], values[1:]):
        d = int(b) - int(a)
        if d == 0:
            run += 1
            continue
        for n, code, c in sorted(gpmf._ZERO_RUNS, key=lambda z: -z[2]):
            while c <= run:
                bits += f'{code:0{n}b}'
                run -= c
        bits += '0' * run
        run = 0
        if abs(d) in mags:
            bits += mags[abs(d)] + ('1' if d < 0 else '0')
        else:
            bits += f'{gpmf._ESC:016b}' + f'{d & 0xff:08b}'
    bits += f'{gpmf._END:016b}'
    bits += '0' * (-len(bits) % 16)
    return b'B\x01' + struct.pack('>H', len(values)) + bytes([values[0], 1]) + int(bits, 2).to_bytes(len(bits) // 8, 'big')


rng = np.random.default_rng(0)
for name, v in (('runs and small steps', np.clip(50 + np.cumsum(rng.integers(-1, 2, 640) * (rng.random(640) < 0.2)), 15, 79)),
                ('long runs', np.repeat([40, 41, 60, 20], 160)),
                ('steps too big for the table, via ESC', np.repeat([10, 100, 20, 120], 160)),
                ('every value changes', (np.arange(640) * 7) % 79)):
    v = v.astype(np.uint8)
    check(f'decompress round trip: {name}', (gpmf._decompress_u8(encode(v)) == v).all())
try:
    gpmf._decompress_u8(b'S\x02' + PAYLOAD[2:])
    check('an unsupported payload type is refused', False)
except ValueError:
    check('an unsupported payload type is refused', True)

# ---- rows to frames -----------------------------------------------------------
rows = np.arange(16)[:, None] * np.ones((16, 640), np.uint8)
# Time-lapse: pairs repeat.  Eight frames, rows 0 1 0 1 4 5 4 5 ... written as distinct numbers here.
got = disp_rows_to_frames(rows, 8)[:, 0].tolist()
check('time-lapse: frames 2m, 2m + 1 use rows 4m, 4m + 1', got == [0, 1, 4, 5, 8, 9, 12, 13], got)
got = disp_rows_to_frames(rows, 16)[:, 0].tolist()
check('real-time: row k for frame k', got == list(range(16)), got[:4])
got = disp_rows_to_frames(rows[:16], 14)[:, 0].tolist()
check('real-time: spare rows at the end are dropped', got == list(range(14)), len(got))
got = disp_rows_to_frames(rows[:13], 8)[:, 0].tolist()
check('time-lapse: a frame with no row is left out', got == [0, 1, 4, 5, 8, 9, 12], got)

# ---- KLV ------------------------------------------------------------------------
klv = b'TSMP' + b'L\x04\x00\x01' + struct.pack('>I', 7) + b'NAME' + b'c\x01\x00\x03' + b'abc\0'
e = list(gpmf._klv(klv))
check('KLV walk, with padding', [x[0] for x in e] == ['TSMP', 'NAME'] and e[1][4] == 20, e[1])
check('KLV values', gpmf._values(klv, *e[0][1:]).ravel().tolist() == [7.0])


# ---- layout ---------------------------------------------------------------------
def track(values):
    """A 4096 x 1344 YUV420 frame whose columns from each [(first col, width, luma)] are flat."""
    f = av.VideoFrame(layout.TRACK_W, layout.TRACK_H, 'yuvj420p')
    y = layout._plane(f.planes[0])
    for col, w, v in values:
        y[:, col:col + w] = v
    layout._plane(f.planes[1]).fill(128)
    layout._plane(f.planes[2]).fill(128)
    f.pts = 0
    f.time_base = fractions.Fraction(1, 30)
    return f


T = layout.TILES
LUMA = {name: 20 + 20 * i for i, name in enumerate(T)}
f0 = track([(T[n][1], T[n][2], LUMA[n]) for n in T if T[n][0] == 0])
f1 = track([(T[n][1], T[n][2], LUMA[n]) for n in T if T[n][0] == 1])
y = layout._plane(compose_cross(f0, f1).planes[0])
check('cross size', y.shape == (layout.CROSS_H, layout.CROSS_W), y.shape)
for name, (tr, col, w, rot, x, yy) in T.items():
    h_, w_ = (w, layout.CENTER) if rot is not None else (layout.TRACK_H, w)
    check(f'{name} lands where the table says', (y[yy:yy + h_, x:x + w_] == LUMA[name]).all())
check('the gaps are black', y[0, 0] == 0 and y[0, 2032] == 0 and y[2719, 5439] == 0)

# ---- orientation helpers ---------------------------------------------------------
from max360.orientation import _min_rotation, _quat_to_matrices

q = rng.normal(size=(50, 4))
q /= np.linalg.norm(q, axis=1, keepdims=True)
Rs = _quat_to_matrices(q)
check('quaternions give rotations', np.allclose(Rs @ Rs.transpose(0, 2, 1), np.eye(3)) and np.allclose(np.linalg.det(Rs), 1))
a, b = np.array([0.0, 1, 0]), np.array([0.6, 0.8, 0])
check('the smallest rotation takes a to b', np.allclose(_min_rotation(a, b) @ a, b))
check('...and leaves the axis alone', np.allclose(_min_rotation(a, b) @ np.array([0, 0, 1.0]), [0, 0, 1]))

# ---- panorama ---------------------------------------------------------------------
try:
    import torch
    from max360.panorama import PanoRenderer
except ImportError:
    print('SKIP panorama: no torch')
    PanoRenderer = None
if PanoRenderer:
    r = PanoRenderer(width=1024, overlap=64, device='cpu')
    check('panorama size', r.height == round(np.radians(120) * 1024 / (2 * np.pi)), r.height)
    # A cross with each tile flat: directions must land on the right tile.
    cross = np.zeros((layout.CROSS_H, layout.CROSS_W, 3), np.uint8)
    for name, (tr, col, w, rot, x, yy) in T.items():
        h_, w_ = (w, layout.CENTER) if rot is not None else (layout.TRACK_H, w)
        cross[yy:yy + h_, x:x + w_] = LUMA[name]
    I = np.eye(3)
    for name, d in (('C0', (0, 0, 1)), ('C1', (0, 0, -1)), ('U0', (0, -1, 0.3)), ('U1', (0, -1, -0.3)), ('D0', (0, 1, 0.3)), ('D1', (0, 1, -0.3)),
                    ('R0', (1, 0, 0.3)), ('R1', (1, 0, -0.3)), ('L0', (-1, 0, 0.3)), ('L1', (-1, 0, -0.3))):
        v = np.array(d, np.float32)[None, None]
        v = v / np.linalg.norm(v)
        check(f'direction {d} reads {name}', r.sample(cross, I, v)[0, 0, 0] == LUMA[name], int(r.sample(cross, I, v)[0, 0, 0]))
    p = r.render(cross, I)
    check('the overlap repeats the first columns', p.shape[1] == 1024 + 64 and (p[:, :64] == p[:, 1024:]).all(), p.shape)
    Ry = np.array([[0.0, 0, 1], [0, 1, 0], [-1, 0, 0]])     # turns the camera: what was at +x is now straight ahead
    a = r.sample(cross, Ry, np.array([[[1.0, 0, 0]]], np.float32))[0, 0, 0]
    check('the camera rotation is applied', a in (LUMA['C1'], LUMA['C0']) and a != r.sample(cross, I, np.array([[[1.0, 0, 0]]], np.float32))[0, 0, 0], int(a))

# ---- a real clip, if given --------------------------------------------------------
clip = os.environ.get('MAX360_TEST_CLIP')
if clip:
    with Reader(clip, threaded=False) as rd:
        it = iter(rd)
        f = next(it)
        check('first frame is k = 0', f.k == 0, f.k)
        check('raw size', f.raw().shape == (layout.RAW_H, layout.RAW_W, 3))
        rd.seek(37)
        f = next(iter(rd))
        check('seek lands on the frame', f.k == 37, f.k)
    d = gpmf.read_disp(clip)
    check('one DISP row per frame', abs(len(d) - rd.frames) <= 2, (len(d), rd.frames))
    o = Orientation(clip)
    R = o.camera_from_world(0)
    check('orientation is a rotation', np.allclose(R @ R.T, np.eye(3), atol=1e-6))

print('FAILED: ' + ', '.join(FAILS) if FAILS else 'all passed')
sys.exit(1 if FAILS else 0)
