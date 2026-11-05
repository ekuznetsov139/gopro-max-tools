"""A gravity-levelled equirectangular panorama from the cross layout, on the GPU if there is one.

This shows the picture as the camera stitched it: the stitch warp is NOT undone, so the outer
344 px of each side tile are not true rays, and straight lines bend there by amounts that depend
on depth.  That suits detection.  For measuring, use the camera model in cpp/.

Tiles are equi-angular, u = (4 / pi) atan(t), not gnomonic.  A side tile keeps the centre tile's
scale along their shared edge, and the mid-plane z = 0 between the lenses falls HW = 656.5 px from
that edge, so a pair of tiles shares a 63 px band, blended linearly.  HW is a fit to the picture
as stitched, not a property of the optics.
"""
import math

import numpy as np
import torch

from .layout import CENTER as FACE, SIDE

HW = 656.5


class PanoRenderer:
    """Cross layout BGR + camera_from_world -> panorama BGR.  Rows are equally spaced in latitude;
    `width` columns make a full circle, and `overlap` further columns repeat the first ones, so
    that someone on the seam is whole in one copy.  Every pixel samples the composite once per lens,
    blended by a weight that is 0 or 1 outside the shared band.  Sample points are clamped to their
    own tile, or bilinear filtering would bring in the black gap or a tile the cross does not join on."""

    def __init__(self, width=5376, lat_up=60.0, lat_down=60.0, overlap=0, device=None):
        self.dev = torch.device(device or ('cuda' if torch.cuda.is_available() else 'cpu'))
        self.width, self.overlap = width, overlap
        self.px_per_radian = width / (2 * math.pi)
        self.lat_up = math.radians(lat_up)
        self.height = int(round(math.radians(lat_up + lat_down) * self.px_per_radian))
        lon = (torch.arange(width, dtype=torch.float32, device=self.dev) + 0.5) / width * 2 * math.pi - math.pi
        lat = self.lat_up - (torch.arange(self.height, dtype=torch.float32, device=self.dev) + 0.5) / self.px_per_radian
        la, lo = torch.meshgrid(lat, lon, indexing='ij')
        self._world = torch.stack([la.cos() * lo.sin(), -la.sin(), la.cos() * lo.cos()], dim=2)

    @staticmethod
    def _tile(col, row, x0, y0, w, h):
        return col.clamp(x0, x0 + w - 1), row.clamp(y0, y0 + h - 1)

    def _maps(self, R, world=None):
        """Cross (col, row) for each lens, and lens 0's weight, for world directions (H, W, 3)."""
        # world @ R.T spelled out: on one AMD card (gfx1201, torch 2.13 + ROCm 7.2) (N,3) @ (3,3) and the same einsum
        # were silently wrong for large N, in float32 and float64.
        Rt = torch.as_tensor(R.T, dtype=torch.float32, device=self.dev)
        w = self._world if world is None else world
        d = w[..., 0:1] * Rt[0] + w[..., 1:2] * Rt[1] + w[..., 2:3] * Rt[2]
        x, y, z = d[..., 0], d[..., 1], d[..., 2]
        ax, ay, az = x.abs(), y.abs(), z.abs()
        eps = 1e-9
        ang = lambda t: (4 / math.pi) * torch.atan(t)
        full = lambda t: (ang(t) + 1) / 2 * FACE - 0.5

        is_c = (az >= ax) & (az >= ay)
        is_x = ~is_c & (ax >= ay)
        s = torch.where(is_c, az, torch.where(is_x, ax, ay)).clamp_min(eps)
        e = ang(z / s)
        n0, n1 = (1 - e) * HW, (1 + e) * HW             # px from each lens's C edge
        w0 = torch.where(is_c, (z > 0).to(d.dtype), (((HW - n0) / (2 * (SIDE - HW))) + 0.5).clamp(0, 1))
        along = full(y / s)                             # L, R, C rows
        a0, a1 = full(x / s), full(-x / s)              # U, D columns; C's columns too

        c0, r0, c1, r1 = [torch.zeros_like(x) for _ in range(4)]

        def put(mask, lens, col, row, x0, y0, w, h):
            col, row = self._tile(col, row, x0, y0, w, h)
            if lens == 0:
                c0[mask], r0[mask] = col[mask], row[mask]
            else:
                c1[mask], r1[mask] = col[mask], row[mask]
        for lens in (0, 1):
            m = is_c & ((z > 0) if lens == 0 else (z <= 0))
            put(m, lens, (688 if lens == 0 else 3408) + (a0 if lens == 0 else a1), 688 + along, *((688, 688) if lens == 0 else (3408, 688)), FACE, FACE)
        m = is_x & (x > 0)
        put(m, 0, 2032 + n0 - 0.5, 688 + along, 2032, 688, SIDE, FACE)
        put(m, 1, 2720 + SIDE - n1 - 0.5, 688 + along, 2720, 688, SIDE, FACE)
        m = is_x & (x <= 0)
        put(m, 0, SIDE - n0 - 0.5, 688 + along, 0, 688, SIDE, FACE)
        put(m, 1, 4752 + n1 - 0.5, 688 + along, 4752, 688, SIDE, FACE)
        m = ~is_c & ~is_x & (y < 0)
        put(m, 0, 688 + a0, SIDE - n0 - 0.5, 688, 0, FACE, SIDE)
        put(m, 1, 3408 + a1, SIDE - n1 - 0.5, 3408, 0, FACE, SIDE)
        m = ~is_c & ~is_x & (y >= 0)
        put(m, 0, 688 + a0, 2032 + n0 - 0.5, 688, 2032, FACE, SIDE)
        put(m, 1, 3408 + a1, 2032 + n1 - 0.5, 3408, 2032, FACE, SIDE)
        # On C, the other lens has no tile; its weight is 0 and it samples the same point.
        c1 = torch.where(is_c & (z > 0), c0, c1)
        r1 = torch.where(is_c & (z > 0), r0, r1)
        c0 = torch.where(is_c & (z <= 0), c1, c0)
        r0 = torch.where(is_c & (z <= 0), r1, r0)
        return c0, r0, c1, r1, w0

    def upload(self, cross):
        return torch.from_numpy(np.ascontiguousarray(cross)).to(self.dev).permute(2, 0, 1)[None].to(torch.float32)

    def _sample(self, src, c0, r0, c1, r1, w0):
        H, W = src.shape[2:]

        def at(col, row):
            grid = torch.stack([2 * col / (W - 1) - 1, 2 * row / (H - 1) - 1], dim=2)[None]
            return torch.nn.functional.grid_sample(src, grid, mode='bilinear', padding_mode='border', align_corners=True)[0]
        out = at(c0, r0)
        if ((w0 > 0) & (w0 < 1)).any() or (w0 == 0).any():
            out = out * w0 + at(c1, r1) * (1 - w0)
        return out

    @torch.no_grad()
    def render(self, cross, R):
        """(height, width + overlap, 3) BGR uint8.  R is Orientation.camera_from_world(k)."""
        out = self._sample(self.upload(cross), *self._maps(R))
        out = out.round().clamp(0, 255).to(torch.uint8).permute(1, 2, 0)
        # cat / contiguous on the GPU: .cpu() keeps a permuted tensor's planar layout, which numpy then copies slowly.
        return torch.cat([out, out[:, :self.overlap]], dim=1).cpu().numpy()

    @torch.no_grad()
    def sample(self, cross, R, world):
        """BGR uint8 for world directions (H, W, 3), straight from the tiles, so with one
        interpolation.  A crop cut from the panorama instead is interpolated twice."""
        w = torch.as_tensor(world, dtype=torch.float32, device=self.dev)
        out = self._sample(self.upload(cross), *self._maps(R, world=w))
        return out.round().clamp(0, 255).to(torch.uint8).permute(1, 2, 0).contiguous().cpu().numpy()
