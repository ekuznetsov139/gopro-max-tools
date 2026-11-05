"""Decoded frames from a .360: the two tracks paired by PTS."""
import collections

import av
import numpy as np

from .layout import compose_cross


class Frame:
    """One video frame.  `k` is its index in the clip, which is also its row in the
    per-frame metadata (gpmf.read_disp, Orientation)."""

    def __init__(self, k, f0, f1):
        self.k, self.f0, self.f1 = k, f0, f1

    def raw(self, format='bgr24'):
        """4096 x 2688: track 0 over track 1, what the camera model in cpp/ takes."""
        return np.vstack([self.f0.to_ndarray(format=format), self.f1.to_ndarray(format=format)])

    def cross(self):
        """5440 x 2720 BGR: the cross layout, what panorama.PanoRenderer takes."""
        return compose_cross(self.f0, self.f1).to_ndarray(format='bgr24')


class Reader:
    """Iterates a clip's frames from the current position.

    Tracks are paired by PTS: a seek lands each on its own keyframe, so the one that comes
    up early is dropped until they meet, and at the end a frame with no partner is dropped.

    With threaded=True, closing PyAV's threaded HEVC decoder after stopping before the end
    of a clip can occasionally hang.  Write out results before closing if that matters."""

    def __init__(self, path, threaded=True):
        self.path = path
        self._c = av.open(path, options={'noautorotate': '1'})
        self._tracks = list(self._c.streams.video)
        if len(self._tracks) != 2:
            self._c.close()
            raise ValueError(f"{path}: expected 2 video tracks, found {len(self._tracks)}")
        if threaded:
            for s in self._tracks:
                s.thread_type = 'AUTO'
        s = self._tracks[0]
        self.frames = min(t.frames for t in self._tracks)
        self.rate = s.average_rate
        self._reset()

    def _reset(self):
        self._demux = self._c.demux(*self._tracks)
        self._queues = {s.index: collections.deque() for s in self._tracks}
        self._skip_to = 0

    def seek(self, k):
        """The next frame yielded is frame k."""
        s = self._tracks[0]
        self._c.seek(int(round(float(k / (s.average_rate * s.time_base)))), stream=s, backward=True)
        self._reset()
        self._skip_to = k

    def __iter__(self):
        s = self._tracks[0]
        q0, q1 = (self._queues[t.index] for t in self._tracks)
        while True:
            while q0 and q1:
                if q0[0].pts < q1[0].pts:
                    q0.popleft()
                elif q1[0].pts < q0[0].pts:
                    q1.popleft()
                else:
                    f0, f1 = q0.popleft(), q1.popleft()
                    k = int(round(float(f0.pts * s.time_base * s.average_rate)))
                    if k >= self._skip_to:
                        yield Frame(k, f0, f1)
            packet = next(self._demux, None)
            if packet is None:
                return
            self._queues[packet.stream.index].extend(packet.decode())

    def close(self):
        self._c.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
