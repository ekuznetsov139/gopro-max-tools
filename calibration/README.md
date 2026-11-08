# Camera model file

`model4.txt` is the clean equi-angular cube for `cpp/maxcam.h`: the corners of the cube, and zero in every other lattice value, so it holds no fitted shape. Only its baseline is specific to a camera, so another unit should use its own, from its file's calibration ([docs/FORMAT.md](../docs/FORMAT.md#factory-calibration)).

It is a text file of 925 rows:

| rows | content |
|---|---|
| 1..902 | lattice values `x y z`: 50 + 178 + 674 for levels 0, 1, 2. Corners of the cube at level 0, everything else zero |
| 903..904 | lens axis rotations, quaternion `w x y z`; identity in every fitted file |
| 905 | baseline `x y z`, cm: lens 0 relative to lens 1 |
| 906..925 | left over from older models; not read |
