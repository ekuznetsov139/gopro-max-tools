"""Where the ten tiles of a GoPro Max frame live.

Each frame is a pair of HEVC tracks, 4096 x 1344:

    track 0    L1 L0 C0 R0 R1
    track 1    D0 D1 C1 U1 U0

The digit is the lens: 0 faces forward, 1 backward.  C tiles are 1344 px wide, the others
688.  "Raw" below means track 0 stacked over track 1, 4096 x 2688, which is the frame the
camera model in cpp/ works on.  "Cross" unfolds each lens around its centre tile,
with black gaps, which is what panorama.py samples:

        U0            U1
    L0  C0  R0  R1    C1  L1
        D0            D1
"""
import av
import cv2 as cv
import numpy as np

TRACK_W, TRACK_H = 4096, 1344
RAW_W, RAW_H = TRACK_W, 2 * TRACK_H
CENTER, SIDE = 1344, 688
CROSS_W, CROSS_H = 5440, 2720

# name: (track, first source column, source width, cv.rotate code or None, cross x, cross y)
TILES = {
    'L0': (0, 688, SIDE, None, 0, 688),
    'C0': (0, 1376, CENTER, None, 688, 688),
    'R0': (0, 2720, SIDE, None, 2032, 688),
    'R1': (0, 3408, SIDE, None, 2720, 688),
    'L1': (0, 0, SIDE, None, 4752, 688),
    'U0': (1, 3408, SIDE, cv.ROTATE_90_CLOCKWISE, 688, 0),
    'D0': (1, 0, SIDE, cv.ROTATE_90_CLOCKWISE, 688, 2032),
    'U1': (1, 2720, SIDE, cv.ROTATE_90_COUNTERCLOCKWISE, 3408, 0),
    'C1': (1, 1376, CENTER, cv.ROTATE_90_COUNTERCLOCKWISE, 3408, 688),
    'D1': (1, 688, SIDE, cv.ROTATE_90_COUNTERCLOCKWISE, 3408, 2032),
}


def _plane(plane):
    """A plane's pixels as a (height, width) uint8 array over the frame's own buffer."""
    return np.frombuffer(plane, np.uint8).reshape(-1, plane.line_size)[:plane.height, :plane.width]


def compose_cross(f0, f1):
    """The cross layout as an av.VideoFrame, from a pair of decoded track frames.  Built in
    YUV, so that colour conversion happens once, over the finished picture."""
    fmt = f0.format.name
    if fmt not in ('yuvj420p', 'yuv420p') or (f0.width, f0.height) != (TRACK_W, TRACK_H):
        raise ValueError(f"unexpected track format {fmt} {f0.width}x{f0.height}")
    out = av.VideoFrame(CROSS_W, CROSS_H, fmt)
    full_range = fmt == 'yuvj420p' or int(f0.color_range) == 2
    for i, sub in enumerate((1, 2, 2)):
        dst = _plane(out.planes[i])
        dst.fill((0 if full_range else 16) if i == 0 else 128)
        src = (_plane(f0.planes[i]), _plane(f1.planes[i]))
        for track, col, w, rot, x, y in TILES.values():
            tile = src[track][:, col // sub:(col + w) // sub]
            if rot is not None:
                tile = cv.rotate(tile, rot)
            h_, w_ = tile.shape
            dst[y // sub:y // sub + h_, x // sub:x // sub + w_] = tile
    out.pts = f0.pts
    out.time_base = f0.time_base
    out.color_range = f0.color_range
    out.colorspace = f0.colorspace
    return out
