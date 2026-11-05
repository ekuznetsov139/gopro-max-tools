# gopro-max-tools

Tools for GoPro Max `.360` files: decoded frames, metadata including `DISP`, gravity-levelled orientation, a panorama, and a camera model with the stitch warp undone.

The file format is mostly undocumented, so [docs/FORMAT.md](docs/FORMAT.md) writes down what's been worked out. Not affiliated with GoPro.

**Status:** tested on files from one camera (firmware `H19.03.02.00.00`). `calibration/model4.txt` is an example fitted to that camera: it's representative of the Max's lens distortion, but another unit will differ. See [Limits](docs/FORMAT.md#limits).

## Python

The package imports as `max360`.

```bash
pip install .              # add [panorama] for the panorama, which needs torch
```

```python
from max360 import Reader, Orientation, read_disp
from max360.panorama import PanoRenderer

path = 'clip.360'
disp = read_disp(path)                  # (frames, 640), one row per frame in either capture mode
orientation = Orientation(path)
renderer = PanoRenderer()

with Reader(path) as reader:
    reader.seek(100)
    for frame in reader:
        raw = frame.raw()               # 4096 x 2688 BGR: track 0 over track 1
        cross = frame.cross()           # 5440 x 2720 BGR: each lens unfolded
        R = orientation.camera_from_world(frame.k)
        panorama = renderer.render(cross, R)    # 5376 x 1792 BGR, gravity straight down the rows
        break
```

The panorama shows the picture as the camera stitched it, which suits detection. Its outer strips aren't true rays, so use the C++ model for measuring.

## C++

`cpp/maxcam.h` and `cpp/maxcam.cpp` have no dependencies beyond a C++17 compiler.

```cpp
#include "maxcam.h"
using namespace maxcam;

CameraModel model;
if (!model.load("calibration/model4.txt"))
    return 1;

// disp: this frame's 640 DISP values as floats, from max360.read_disp
Vec2 m = StitchWarp::recorded_to_model(Vec2(600, 700), disp);
Vec3 ray = model.pixel_to_ray(m).normalized();
```

## Tests

```bash
python3 tests/test_python.py          # no footage needed; set MAX360_TEST_CLIP to also read a clip
make -C cpp test
```

`tools/gpmf_bug` reproduces a `DISP` decoding bug in gpmf-parser.

## Licence

Apache-2.0; see `LICENSE`. `max360/gpmf.py` follows GoPro's gpmf-parser, and `tools/gpmf_bug/` builds on gpmf-parser and gpmf-write, © GoPro Inc., under Apache-2.0 or MIT.
