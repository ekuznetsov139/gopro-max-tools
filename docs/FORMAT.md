# GoPro Max: file format, camera model and stitch warp

What a GoPro Max `.360` file contains, how its pixels relate to viewing rays, and how the camera warps them. This comes from working on files from one camera (firmware `H19.03.02.00.00`). GoPro has not documented most of it, and this project is not affiliated with GoPro. See [Limits](#limits) for what hasn't been established.

## The file

A `.360` is an MP4 with two HEVC video tracks, 4096×1344, 4:2:0, about 30 Mb/s each at 29.97 fps. A frame is the pair of pictures with the same PTS. Clips split into chapters at about 4 GB (`GS01xxxx`, `GS02xxxx`, …).

Other tracks: audio (real-time video only), timecode, and a `GoPro MET` metadata track in GPMF format.

### Tiles

Each lens's view is folded into five tiles: a centre tile 1344 px wide (C) and four side tiles 688 px wide (L, R, U, D). The digit is the lens: 0 faces forward, 1 backward.

```
track 0    L1 L0 C0 R0 R1
track 1    D0 D1 C1 U1 U0
```

Each side tile's outer edge is at the stitch, and its inner edge joins its centre tile. The stitch lines in the raw frame are at x = 687.5 and x = 3407.5.

Track 1's tiles are rotated. To unfold each lens into a cross, rotate U0 and D0 clockwise, and U1, C1 and D1 counter-clockwise. `max360/layout.py` has the full table.

"Raw frame" below means track 0 stacked over track 1, 4096×2688.

### Tile projection

Tiles are equi-angular cube faces, not gnomonic: a direction with tangent `t` on a face lands at `u = (4/π)·atan(t)`, so 1344 px span 90°, or about 856 px/rad.

If you take the picture as stitched, a side tile keeps the centre tile's scale along their shared edge, and the plane between the lenses falls about 656.5 px from that edge. That gives each pair of tiles a 63 px shared band. This is a fit to the stitched picture, good enough for the panorama in `max360/panorama.py`, not a property of the optics. For true rays, use the [camera model](#camera-model).

## Capture modes

The mode changes how metadata relates to frames:

| | Time-lapse | Real-time video |
|---|---|---|
| Capture rate (from `STMP` timestamps) | 2 fps at the setting used here | about 30 fps |
| Playback rate | 29.97 fps | 29.97 fps |
| Audio track | no | yes |
| `DISP` rows per frame | about 2, in repeated pairs | about 1 |

To tell them apart, use the presence of audio or the ratio of `DISP` rows to frames. **Don't use the `Rate` tag** that exiftool reports: the same value (`2_1SEC`) appears on both time-lapse and real-time clips.

## Metadata

All of this is in the `GoPro MET` track, read by `max360/gpmf.py` with no C parser. GoPro's [gpmf-parser](https://github.com/gopro/gpmf-parser) documents the container.

### CORI and GRAV

Both have one row per video frame, plus 1–2 spare rows at the end, so row k is frame k. In these camera coordinates x is right, y is down, and z points through the centre of C0.

- **`CORI`** is a quaternion (w, x, y, z), world to camera. It's stored as int16 over 32767, so it's only about unit length until normalised. It starts near identity however the camera is held, so it says nothing about up, and it drifts against gravity by more than 10° over a few minutes.
- **`GRAV`** is a unit gravity vector in camera coordinates. Raw, it jitters by about 0.3° per frame.

`max360/orientation.py` combines them: `GRAV`, smoothed, supplies up, and `CORI` holds the heading. Heading is relative to the first frame, not to north.

### DISP

`DISP` is a 1-D depth map along the stitch: 640 values per row, one for each unit of φ, measured by the camera for its stitch. Stored values range from 15 (very close) to 79 (far). It exists only in 360 modes.

**φ** runs around the stitch. It's 0 straight up and increases clockwise, so 160 is right, 320 is down and 480 is left. One unit is 8.4 px.

**Storage.** Each row is a compressed (`#`) payload: Huffman-coded deltas with zero runs, in the scheme described in gpmf-parser's `GPMF_bitstream.h`. Payloads have 2–6 unexplained bytes after the end code; ignoring them does no harm.

**Rows and frames.** `gpmf.read_disp` returns one row per frame, handling both modes:

- **Real-time:** row k is frame k.
- **Time-lapse:** rows come in repeated pairs, so rows 4m+2 and 4m+3 are copies of rows 4m and 4m+1. Frame 2m uses row 4m, and frame 2m+1 uses row 4m+1.

### A decoder bug

`GPMF_Decompress` in gpmf-parser (checked on `main` at `9a71506`, August 2026) gives wrong `DISP` values. When a zero run and the delta after it share one 16-bit lookup, it applies the delta before writing the run, so samples just before a step take the value from after it.

GoPro's own encoder, `GPMFCompress`, settles which decoding is right: re-encoding a row reproduces the file's bytes from the decoding in `max360/gpmf.py`, but not from gpmf-parser's. `tools/gpmf_bug/` has a reproduction and a patch.

## Stitch warp

The lenses are about 3.8 cm apart, so nearby objects don't line up where their views meet. The camera hides this by warping the outer part of each side tile, frame by frame, according to `DISP`. This happens before encoding and can't be turned off, so **the outer 344 px of each side tile are not true rays**, by amounts that depend on scene depth.

### Strip coordinates

The warp is simplest in strip coordinates, where the four quadrants of a lens's outer tiles unroll into one strip: `s` is the distance from the stitch, and `t` runs along the whole stitch, 0 to 5376, wrapping. From a raw pixel (x, y):

| Region | s | t |
|---|---|---|
| x ≥ 2048, y < 1343.5 | \|x − 3407.5\| | y + 0.5 |
| x < 2048, y ≥ 1343.5 | \|x − 687.5\| | 1344 + y − 1343.5 |
| x < 2048, y < 1343.5 | \|x − 687.5\| | 2·1344 + 1343.5 − y |
| x ≥ 2048, y ≥ 1343.5 | \|x − 3407.5\| | 3·1344 + 2687.5 − y |

Then `φ = t/5376 · 640 + 80`, wrapped to 0–640, and `i = floor(φ)`.

### The warp

This takes a **recorded** pixel to a **model** pixel, which is the one the camera model takes. It applies only to side tiles (raw x < 1376 or x ≥ 2720).

```
C(s) = (688 - s) * (1.30 * 69.12 / 688)           static part, every side-tile pixel
f    = 1 - s / 344

if f < 0:     s' = s + C(s),  t' = t              beyond 344 px: static part only

g = Z(f)                                          falloff: linear interpolation in the table below
v(j) = 0.333 * d[j] + 0.667 * (d[j-1] + d[j+1]) / 2       d = this frame's DISP row; j wraps at 640
D = (1 - (φ - i)) * v(i) + (φ - i) * v(i+1)
a = (1 - D / 80) * g

s' = s + C(s) - 69.12  * a * (c1 + c2 * cos(2π i / 160))
t' = t + σ * 59.875 * a * c3 * cos(π i / 320)     σ = +1 on lens 0, -1 on lens 1, t' wraps at 5376
```

With `c1 = 1.279495`, `c2 = −0.188683`, `c3 = 1.002468`, and `Z` taking these values at f = 0, 1/16, …, 1:

```
0.002021 0.021586 0.051468 0.099344 0.161015 0.233702 0.303061 0.378131 0.449792
0.526414 0.591500 0.665341 0.727919 0.792834 0.860367 0.924303 1.0
```

These constants are fits, from footage of one camera.

- **Size:** at the stitch, the static part moves a pixel by about 90 px. The depth-driven part reaches about 82 px across the stitch and 49 px along it for the closest objects.
- **`DISP` of 80 cancels the depth-driven part.** Stored values top out at 79, so it never quite vanishes.
- **The reverse direction has no closed form.** `StitchWarp::model_to_recorded` in `cpp/` solves it by iteration.
- **Some model pixels have no recorded source.** Inversion usually fails within about 100 px of the stitch, which fits those pixels having been discarded by the camera, though this hasn't been proven. Going the other way, recorded → model → recorded is reliable on nearly all side-tile pixels.

## Camera model

`cpp/maxcam.h` maps a model pixel to the viewing ray of its lens. It's a lattice spline on top of the tile layout. `calibration/model4.txt` is an example, fitted to one camera, and representative of the Max's lens distortion rather than an exact fit to any other unit.

- **Coordinates:** x right, y up, z forward, unlike the metadata's y down. Lens 0 faces +z, lens 1 faces −z. Rays aren't normalised.
- **Baseline:** in a camera frame centred between the lenses, lens 0 sits at +b/2 and lens 1 at −b/2, with `b = (−2.328, −0.0005, 3.038)` cm in the example.
- **Lattice:** three levels with 2, 4 and 8 cells per tile, 902 3-D values in all, shared along tile edges so that each lens is continuous. `calibration/README.md` describes the file.

A recorded pixel becomes a ray in two steps, and looking up a ray reverses them:

```cpp
Vec2 m   = StitchWarp::recorded_to_model(recorded_px, disp);
Vec3 ray = model.pixel_to_ray(m);

Vec2 m2  = model.ray_to_pixel(ray, lens);
Vec2 px  = StitchWarp::model_to_recorded(m2, disp);
```

## Limits

- **The warp on real-time video.** The constants come from time-lapse footage and are assumed to apply to real-time video. That hasn't been measured directly.
- **Other cameras.** The lattice, baseline and warp constants were fitted to one unit, on one firmware version.
- **GoPro MAX2 and later.** Nothing here has been tested on them.
- **Colour.** The camera's colour processing isn't covered.

## Shooting calibration footage

For calibration, use time-lapse or time-warp with white balance away from auto, ISO fixed (maximum = minimum), and the shutter fixed. Time-lapse modes have no shutter setting of their own: they use the one from video mode. Keep the camera steady, as frame-to-frame changes in exposure and colour can't currently be compensated.
