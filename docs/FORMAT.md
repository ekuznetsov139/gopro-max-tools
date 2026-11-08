# GoPro Max: file format, camera geometry and stitch warp

What a GoPro Max `.360` file contains, how its pixels relate to viewing rays, and how the camera warps them. This comes from working on files from one camera (firmware `H19.03.02.00.00`). GoPro has not documented most of it, and this project is not affiliated with GoPro. See [Limits](#limits) for what hasn't been established.

**In short:** the camera has already removed lens distortion using its factory calibration, so the recorded tiles are close to a clean equi-angular cube. The one major thing left to undo is a per-frame warp on the outer part of each side tile.

`cpp/` implements the [cube](#ray-from-a-pixel) and [warp](#the-warp) below. `max360/` doesn't include them yet.

## Hardware

| | MAX | MAX2 (not tested here) |
|---|---|---|
| Sensors | 2× Sony IMX577, 1/2.3", 7.857 mm diagonal | 2× Sony IMX677, 1/2.3", 7.85 mm diagonal |
| Pixels | 12.3 MP, 4056×3040, 1.55 µm | 23.91 MP, reported 4216×4216 capture and 3840×3840 active area per lens |
| Mounting | Flat PCBs, one per lens, about 20 mm behind it, offset sideways | Back to back, on one axis |
| Lens spacing | 31–33 mm along the axis, about 23 mm across | |
| 360 video | 5376×2688, up to 30 fps | 7680×3840 up to 30 fps, 5376×2688 up to 60, 3328×1664 up to 100 |

Each lens has its own one-piece sensor, so the tiles below are a digital remapping of its round image, not a physical layout. On the MAX, each sensor sits past the camera's centre, which is why the two can't share an axis. Teardowns: [MAX](https://gethypoxic.com/blogs/technical/gopro-max-teardown), [MAX2](https://gethypoxic.com/blogs/technical/gopro-max-2-review-teardown).

## The file

A `.360` is an MP4 with two HEVC video tracks, 4096×1344, 4:2:0, about 30 Mb/s each at 29.97 fps. A frame is the pair of pictures with the same PTS. Clips split into chapters at about 4 GB (`GS01xxxx`, `GS02xxxx`, …).

Other tracks: audio (real-time video only), timecode, and a `GoPro MET` metadata track in GPMF format. Factory calibration is in a separate GPMF block, in `moov/udta/GPMF`.

### Tiles

```
track 0    L1 L0 C0 R0 R1
track 1    D0 D1 C1 U1 U0
```

The digit is the lens: 0 faces forward, 1 backward. Each centre tile (C) is a full cube face, 1344 px. Each side tile (L, R, U, D) is 688 px: 656 px from the centre-tile edge to the plane between the lenses, then 32 px of overlap past it. So paired tiles share a 64 px band. The `PMOD` field, `[32, 1408, 1376]`, appears to describe this layout: overlap, 1344 + 2·32, and 1344 + 32.

The stitch lines in the raw frame are at x = 687.5 and x = 3407.5. "Raw frame" means track 0 stacked over track 1, 4096×2688.

Track 1's tiles are rotated. To unfold each lens into a cross, rotate U0 and D0 clockwise, and U1, C1 and D1 counter-clockwise. `max360/layout.py` has the full table.

### Ray from a pixel

Tiles are equi-angular, not gnomonic: a direction with tangent `t` lands at `u = (4/π)·atan(t)`, so a full face spans 90° across 1344 px, or about 856 px/rad. A gnomonic cube is off by up to about 60 px within a face.

Coordinates are x right, y up, z forward. For a raw pixel (x, y) outside the warp, with `v = 1 − 2·(y + 0.5 − top)/1344`, where `top` is 0 on track 0 and 1344 on track 1, and `w = tan(v·π/4)`:

| Raw x | Tile | Lens on track 0 / 1 | n (px from the C edge) | Direction before normalising |
|---|---|---|---|---|
| 0–687 | L1 / D0 | 1 / 0 | x + 0.5 | (−1, w, −q) |
| 688–1375 | L0 / D1 | 0 / 1 | 1376 − (x + 0.5) | (−1, w, q) |
| 1376–2719 | C0 / C1 | 0 / 1 | | (tan(a·π/4), w, 1), with a = 2·(x + 0.5 − 1376)/1344 − 1 |
| 2720–3407 | R0 / U1 | 0 / 1 | (x + 0.5) − 2720 | (1, w, q) |
| 3408–4095 | R1 / U0 | 1 / 0 | 4096 − (x + 0.5) | (1, w, −q) |

Here `q = tan((1 − n/656)·π/4)`. On track 1, swap the first two components and negate the third.

All rays share one frame, with the rear lens's small rotation already removed by the camera. Lens 0 sits at +B/2 and lens 1 at −B/2, where B is the baseline from the file's calibration, the negative of `SHFX/Y/Z`: (−2.336, 0.002, 3.580) cm on the camera tested.

## Factory calibration

Every file carries the camera's own calibration, in `moov/udta/GPMF` (`max360.read_calibration` reads it; `max360.baseline` gives B in cm):

| Block | Field | Meaning | Example |
|---|---|---|---|
| `GEOC` | `CALW`, `CALH` | Sensor size the calibration refers to | 4056, 3040 |
| | `SHFX`, `SHFY`, `SHFZ` | Rear lens relative to front, mm | 23.36, −0.02, −35.80 |
| | `ANGX`, `ANGY`, `ANGZ` | Rear lens rotation, degrees | −0.18, 180.31, −0.10 |
| `FRNT`, `BACK` | `KLNS` | Lens polynomial, `r(θ) = k1·θ + k2·θ³ + k3·θ⁵ + k4·θ⁷ + k5·θ⁹`, θ in radians, r in sensor px | 1082.3, −59.8, 2.8, −2.9, 0.4 |
| | `CTRX`, `CTRY` | Principal point offset, px | −1.59, 2.03 |
| | `MFOV` | Half field of view at the sensor's short edge, degrees | 100 |
| | `SFTR` | Close to 1; meaning unknown | 1.0002 |
| Global | `PRJT`, `PMOD` | Projection (`EACO`) and tile layout | `[32, 1408, 1376]` |

These meanings are inferred, not documented. The polynomial form is supported by `r(100°)` coming to about 1520 px, the sensor's half-height.

`SHFX/Y/Z` and `ANGX/Y/Z` are identical to every digit across files from the camera tested, recorded in different sessions, so they're stored values, not measured per clip. Whether they're calibrated per unit or follow from the design isn't known, and their many digits say nothing about their accuracy. They're probably what the camera used to build the picture, but that doesn't establish them as the true offset between viewpoints, and they give one offset per lens, for all ray angles.

**Don't apply `KLNS` to recorded pixels.** The camera has already used it: applying it again, in either direction, makes the fit far worse. The fields are still useful for two things:

- **The baseline.** `SHFZ` puts the viewpoints about 2 mm outside the lens spacing above each side of the camera. At 90° to the axis a ray stays at one z, which has to lie within the glass where it enters, so stitch rays are out at the convex glass rather than back near the sensors.
- **Native resolution.** The sensor has about 1083 px/rad at the centre, falling to about 530 radially at 90°, compared with 856 in the file. The file therefore discards detail near the centre of each lens, and interpolates near the stitch.

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

`DISP` is a 1-D depth map along the stitch: 640 values per row, one for each unit of φ, measured by the camera for its stitch. It exists only in 360 modes. Stored values range from 15 (very close) to 79.

It behaves roughly as inverse distance: `distance ≈ 15 m / (66.8 − DISP)`, so 66.8 is about infinity. That law was checked at under a metre, where the camera has real stereo to work with. Values above 66.8 occur too, and with one lens covered `DISP` stops tracking distance, though the camera still warps by it.

**φ** runs around the stitch. It's 0 straight up and increases clockwise, so 160 is right, 320 is down and 480 is left. One unit is 8.4 px.

**Storage.** Each row is a compressed (`#`) payload: Huffman-coded deltas with zero runs, in the scheme described in gpmf-parser's `GPMF_bitstream.h`. Payloads have 2–6 unexplained bytes after the end code; ignoring them does no harm.

**Rows and frames.** `gpmf.read_disp` returns one row per frame, handling both modes:

- **Real-time:** row k is frame k.
- **Time-lapse:** rows come in repeated pairs, so rows 4m+2 and 4m+3 are copies of rows 4m and 4m+1. Frame 2m uses row 4m, and frame 2m+1 uses row 4m+1.

### A decoder bug

`GPMF_Decompress` in gpmf-parser (checked on `main` at `9a71506`, August 2026) gives wrong `DISP` values. When a zero run and the delta after it share one 16-bit lookup, it applies the delta before writing the run, so samples just before a step take the value from after it.

GoPro's own encoder, `GPMFCompress`, settles which decoding is right: re-encoding a row reproduces the file's bytes from the decoding in `max360/gpmf.py`, but not from gpmf-parser's. `tools/gpmf_bug/` has a reproduction and a patch.

## Stitch warp

The lenses are several centimetres apart, so nearby objects don't line up where their views meet. The camera hides this by warping the outer part of each side tile, frame by frame, according to `DISP`. This happens before encoding and can't be turned off, so **the outer 344 px of each side tile are not true rays**, by amounts that depend on scene depth.

### Strip coordinates

The warp is simplest in strip coordinates, where the four quadrants of a lens's outer tiles unroll into one strip: `s` is the distance from the stitch, and `t` runs along the whole stitch, 0 to 5376, wrapping. From a raw pixel (x, y):

| Region | s | t |
|---|---|---|
| x ≥ 2048, y < 1343.5 | \|x − 3407.5\| | y + 0.5 |
| x < 2048, y ≥ 1343.5 | \|x − 687.5\| | 1344 + y − 1343.5 |
| x < 2048, y < 1343.5 | \|x − 687.5\| | 2·1344 + 1343.5 − y |
| x ≥ 2048, y ≥ 1343.5 | \|x − 3407.5\| | 3·1344 + 2687.5 − y |

Then `φ = t/5376 · 640 + 80`, wrapped to 0–640, `i = floor(φ)` and `ψ = 2π·i/640`. In the [ray table](#ray-from-a-pixel), `n = 688 − s`.

### The warp

This takes a recorded pixel to its unwarped position on the clean cube. It applies only within 344 px of the stitch (`s < 344`):

```
v(j) = 0.333 * d[j] + 0.667 * (d[j-1] + d[j+1]) / 2       d = this frame's DISP row; j wraps at 640
D    = (1 - (φ - i)) * v(i) + (φ - i) * v(i+1)
a    = (66.77 - D) / 80 * Z(s)                            Z: the falloff below

s' = s - 88.4 * a / max(|cos ψ|, |sin ψ|)
t' = t + σ * 88.4 * ρ * a * cos ψ                         σ = +1 on lens 0, -1 on lens 1; t' wraps at 5376
```

Here `ρ = |SHFX / SHFZ|`, 0.653 on the camera tested. `Z` is a parabola that reaches zero 332 px from the stitch, continued towards the stitch by its tangent from 247 px on; past the tile's outer edge (`s < 0`) it keeps its edge value:

```
Z(s) = c (332 - s)²                              247 <= s < 332,  c = 2.03e-5
Z(s) = c (332 - 247) (332 - 247 + 2 (247 - s))   s < 247
Z(s) = 0                                         s >= 332
```

So `Z` is about 1.00 at the stitch, 0.70 at 86 px, 0.41 at 172 px, 0.11 at 258 px and 0.02 at 301 px.

This is the lenses' parallax, as far as tested:

- **Across the stitch**, the shift comes from the lenses' separation along the axis. The `max` term is the cube's own scale at that point on the stitch, largest at the tile corners.
- **Along the stitch**, it comes from their sideways offset, which lies along the stitch straight up and down, and points along the ray at the sides. Hence `cos ψ` and the ratio ρ. Left free, the fit gives 0.646.
- **Size:** for the closest objects (`DISP` = 15), about 57 px across the stitch at the middle of a face, 81 px at a tile corner, and 37 px along it.

The fitted values here are 66.77, 88.4 and the three constants of `Z`. `Z` was measured from pairs of frames of a still camera in which only `DISP` changed, so poses, markers and lens model drop out; it replaces a 17-value table from the older model, which is in this repository's history. The table is higher by up to about 0.05 between 100 and 220 px, and about twice this `Z` from 250 to 300 px.

## Corrections lattice

`CameraModel` in `cpp/maxcam.h` can hold corrections on top of the cube: a lattice spline over the tile layout, in three levels with 2, 4 and 8 cells per tile, 902 3-D values in all. The tile fractions are first remapped so that level 0's bilinear part is exactly the cube. The values are corrections added to the cube's rays; `calibration/model4.txt` has them all zero, so the ray is the cube's, and that file can also carry a baseline correction and a rotation of lens 1.

At present there is no evidence that any lattice corrections improve the fit.

## Limits

- **A small radial term may remain.** Marker data leaves room for up to about 1 px of smooth radial distortion on top of the cube, and isn't sensitive enough to rule it in or out. It also doesn't constrain the area near the stitch independently of the warp.
- **The warp isn't fully modelled.** Errors within the warped strip are several times those outside it, and they grow towards the tile corners, where the shift appears to increase faster than the `max` term says. `Z` was measured at one range of distances, well under 1 m, and predicts a far scene slightly worse than the old table did; the two lenses may not respond equally. The warp's gain implies about 13.5 m in the distance law above, against about 15 m measured, and that gap is unexplained.
- **The warp on real-time video.** The constants come from time-lapse footage and are assumed to apply to real-time video. Row timing has been checked on real-time video, but the warp's size hasn't.
- **Other cameras.** All of this comes from one unit, on one firmware version. If each camera rectifies with its own calibration, the cube should hold for other units, using their own baseline. That's untested.
- **GoPro MAX2 and later.** Nothing here has been tested on them. With its lenses on one axis, the MAX2's along-stitch term should be close to zero, so its warp may differ.
- **Colour.** The camera's colour processing isn't covered.

## Shooting footage

For geometry (calibration, measuring the warp), exposure and colour settings don't matter; what helps is a still camera and a `DISP` that changes while nothing else does.

For feature matching and 3-D reconstruction across frames, colours must not drift, so use time-lapse or time-warp with white balance away from auto, ISO fixed (maximum = minimum), and the shutter fixed. Time-lapse modes have no shutter setting of their own: they use the one from video mode.
