# Example calibration

`model4.txt` was fitted to one GoPro Max. It is representative of the lens distortion, but another unit will differ.

It is a text file of 925 rows:

| rows | content |
|---|---|
| 1..902 | lattice values `x y z`: 50 + 178 + 674 for levels 0, 1, 2 |
| 903..904 | lens axis rotations, quaternion `w x y z`; identity in every fitted file |
| 905 | baseline `x y z`, cm |
| 906..921 | the warp constants: the first 3 rows' first column are `c1 c2 c3`, and the third column is the falloff table. `cpp/maxcam.cpp` has these built in and does not read them |
| 922..925 | unused scale factors |
