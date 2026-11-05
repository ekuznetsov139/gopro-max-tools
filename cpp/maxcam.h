// Camera model and stitch-warp inversion for the GoPro Max.
//
// Two separate maps, both on the raw 4096x2688 frame (track 0 stacked over track 1):
//
//   CameraModel  maps a "model pixel" to the viewing ray of its lens.
//   StitchWarp   maps between model pixels and the frame as recorded, which the
//                camera has warped per frame from its DISP metadata.
//
// So a pixel of a recorded frame becomes a ray by:
//
//   Vec2 m = StitchWarp::recorded_to_model(recorded_px, disp);
//   Vec3 ray = model.pixel_to_ray(m);
//
// and a ray is looked up in a recorded frame by:
//
//   Vec2 m = model.ray_to_pixel(ray, lens);
//   Vec2 recorded_px = StitchWarp::model_to_recorded(m, disp);
//
// See docs/FORMAT.md.  Coordinates: x right, y up, z forward; pixel centres at
// integers.  Lens 0 faces +z, lens 1 faces -z.  Not affiliated with GoPro.
#pragma once

#include <cmath>
#include <string>
#include <vector>

namespace maxcam {

const int FRAME_W = 4096;       // one track's width
const int TILE_H = 1344;        // one track's height, and the centre tile's width
const int SIDE_W = 688;         // width of each side tile
const int FRAME_H = 2 * TILE_H;
const int DISP_LEN = 640;       // entries in one frame's disparity row

struct Vec2 {
    double x = 0, y = 0;
    Vec2() {}
    Vec2(double x_, double y_) : x(x_), y(y_) {}
    Vec2 operator+(Vec2 o) const { return Vec2(x + o.x, y + o.y); }
    Vec2 operator-(Vec2 o) const { return Vec2(x - o.x, y - o.y); }
    double norm2() const { return x * x + y * y; }
};

struct Vec3 {
    double x = 0, y = 0, z = 0;
    Vec3() {}
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    Vec3 operator+(Vec3 o) const { return Vec3(x + o.x, y + o.y, z + o.z); }
    Vec3 operator-(Vec3 o) const { return Vec3(x - o.x, y - o.y, z - o.z); }
    Vec3 operator*(double s) const { return Vec3(x * s, y * s, z * s); }
    double dot(Vec3 o) const { return x * o.x + y * o.y + z * o.z; }
    double norm() const { return std::sqrt(dot(*this)); }
    Vec3 normalized() const { return *this * (1.0 / norm()); }
};

// Which of the ten tiles a pixel falls in (track * 5 + column, columns left to right).
int zone_index(Vec2 px);
// Which lens a pixel belongs to: 0 or 1.
int lens_index(Vec2 px);

class CameraModel {
public:
    CameraModel();

    // Reads a calibration file (calibration/README.md).  Returns false, with error()
    // set, if the file is missing or malformed.
    bool load(const std::string& path);
    const std::string& error() const { return error_; }

    // The ray of a pixel, not normalized.  gx and gy, when given, receive its
    // derivative with respect to the pixel's x and y.
    Vec3 pixel_to_ray(Vec2 px, Vec3* gx = nullptr, Vec3* gy = nullptr) const;

    // The pixel of a ray on lens `side` (0 or 1; -1 tries lens 0, then lens 1), by
    // Newton iteration.  A pixel with x < -10 means the ray has none on that lens.
    Vec2 ray_to_pixel(Vec3 ray, int side = -1) const;

    // Lens 1's centre relative to lens 0's, in centimetres, as fitted.
    Vec3 baseline() const { return baseline_; }

    // The nominal mapping with no calibration: the tile layout alone.
    static Vec2 box_ray_to_pixel(Vec3 ray, int side);
    static Vec3 box_pixel_to_ray(Vec2 px);

private:
    static const int LEVELS = 3;
    std::vector<Vec3> data_;
    std::vector<int> tiles_[2][LEVELS][5];      // indices into data_
    Vec3 baseline_;
    std::string error_;

    Vec3 spline(int side, int tile, double fx, double fy, Vec3* gx, Vec3* gy) const;
    Vec2 ray_to_pixel_on(Vec3 ray, int side) const;
};

class StitchWarp {
public:
    // `disp` is that frame's DISP_LEN disparities.
    static Vec2 recorded_to_model(Vec2 recorded_px, const float* disp);
    // Solved by iteration.  Near the stitch a model pixel may have no recorded source
    // (the camera discarded it); the result is then the closest attempt, not a solution,
    // so check it against recorded_to_model if that matters.
    static Vec2 model_to_recorded(Vec2 model_px, const float* disp);
};

}  // namespace maxcam
