"""Camera orientation per frame, levelled against gravity.

From the metadata track: row k of each is video frame k.

  CORI  unit quaternion (w, x, y, z), world -> camera, stored as int16 over 32767 and so
        only about unit length until normalised.  It starts near identity however the
        camera is held, so it says nothing about up, and it drifts against gravity.
  GRAV  unit gravity vector in camera coordinates.  It supplies up, per frame.  Raw it
        jitters by about 0.3 degrees per frame, so it is smoothed, in CORI's world.

Camera: x right, y down, z through the centre of C0.  World: y along gravity, heading held
by CORI, z = C0's direction at frame 0, not north.
"""
import math

import cv2 as cv
import numpy as np

from .gpmf import read_tracks

GRAV_SMOOTH_FRAMES = 61


def _quat_to_matrices(q):
    w, x, y, z = q.T
    return np.stack([1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w),
                     2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w),
                     2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)], axis=1).reshape(-1, 3, 3)


def _min_rotation(a, b):
    """The smallest rotation taking unit a to unit b."""
    v = np.cross(a, b)
    s = float(np.linalg.norm(v))
    if s < 1e-9:
        return np.eye(3)
    K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]]) / s
    theta = math.atan2(s, float(a @ b))
    return np.eye(3) + math.sin(theta) * K + (1 - math.cos(theta)) * K @ K


class Orientation:
    """camera_from_world(k): CORI's rotation for frame k, re-levelled so that world +y is
    gravity.  The correction only tilts, so CORI's heading survives."""

    def __init__(self, path):
        d = read_tracks(path, ('CORI', 'GRAV'))
        self.cori, self.grav = d['CORI'] / np.linalg.norm(d['CORI'], axis=1, keepdims=True), d['GRAV']
        self.n = min(len(self.cori), len(self.grav))
        self._Rc = _quat_to_matrices(self.cori[:self.n])
        gw = np.einsum('nji,nj->ni', self._Rc, self.grav[:self.n])      # gravity in CORI's world
        pad = GRAV_SMOOTH_FRAMES // 2
        gw = cv.blur(np.pad(gw, ((pad, pad), (0, 0)), mode='edge'), (1, GRAV_SMOOTH_FRAMES),
                     borderType=cv.BORDER_REPLICATE)[pad:-pad]
        self._gw = gw / np.linalg.norm(gw, axis=1, keepdims=True)
        y = self._gw[0]
        f = self._Rc[0].T @ np.array([0.0, 0.0, 1.0])
        f = f - (f @ y) * y
        f /= np.linalg.norm(f)
        self._level0 = np.c_[np.cross(y, f), y, f]                      # level world -> CORI's world

    def camera_from_world(self, k):
        k = min(max(int(k), 0), self.n - 1)
        return self._Rc[k] @ _min_rotation(self._gw[0], self._gw[k]) @ self._level0
