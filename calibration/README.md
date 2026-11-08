# Corrections file

`model4.txt` holds corrections for `cpp/maxcam.h`'s `CameraModel`, added on top of the clean equi-angular cube of [docs/FORMAT.md](../docs/FORMAT.md). The file shipped here is all zero: no corrections, the model is the cube, and the baseline comes from the `.360` file's own header (or the built-in default). At present there is no evidence that any lattice correction improves the fit.

It is a text file of 904 rows, three numbers each:

| rows | content |
|---|---|
| 1..902 | lattice corrections `x y z`, added to the cube's rays: 50 + 178 + 674 values for levels 0, 1, 2 (2, 4 and 8 cells per tile) |
| 903 | baseline correction `x y z`, cm, added to the header's (or default) baseline |
| 904 | rotation of lens 1, degrees, as a rotation vector `x y z` |

Load it with `model.load(clip_path, "calibration/model4.txt")`; either path may be empty.
