// Camera model and stitch-warp inversion for the GoPro Max.
//
// The model is the clean equi-angular cube of docs/FORMAT.md:
//
//   CameraModel  maps a "model pixel" to the viewing ray of its lens.  By default it is the cube
//                with a built-in baseline; a .360 file supplies the camera's own baseline, and a
//                corrections file (calibration/README.md) can add a lattice of ray corrections, a
//                baseline change and a rotation of lens 1.  With no corrections the cube is used
//                in closed form both ways.
//   StitchWarp   maps between model pixels and the frame as recorded, which the
//                camera has warped per frame from its DISP metadata.
//
// All on the raw 4096x2688 frame (track 0 stacked over track 1).  A pixel of a recorded frame
// becomes a ray by:
//
//   Vec3 ray = model.recorded_to_ray(recorded_px, disp);
//
// and a ray is looked up in a recorded frame by:
//
//   Vec2 recorded_px = model.ray_to_recorded(ray, lens, disp);
//
// A model pixel can lie past its tile's outer edge, in the next tile's pixel range but still on
// its own lens, so the calls that take one also take the lens.
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

    // Both arguments optional (""): a .360 file, whose header gives the camera's baseline, and a
    // corrections file (calibration/README.md), added on top.  Returns false, with error() set, if
    // a given file is missing or malformed; the model is then back at its defaults.
    bool load(const std::string& clip_path, const std::string& corrections_path = "");
    const std::string& error() const { return error_; }

    // A recorded pixel's ray, normalized; and the recorded pixel of a ray on lens `lens` (0 or 1),
    // with x < -10 if that lens has none.  `disp` is that frame's DISP_LEN disparities.
    Vec3 recorded_to_ray(Vec2 recorded_px, const float* disp) const;
    Vec2 ray_to_recorded(Vec3 ray, int lens, const float* disp) const;

    // The ray of a model pixel, not normalized.  gx and gy, when given, receive its derivative
    // with respect to the pixel's x and y.  Without `lens`, the pixel's own tile decides.
    Vec3 pixel_to_ray(Vec2 px, Vec3* gx = nullptr, Vec3* gy = nullptr) const;
    Vec3 pixel_to_ray(Vec2 px, int lens, Vec3* gx = nullptr, Vec3* gy = nullptr) const;

    // The model pixel of a ray on lens `side` (0 or 1; -1 tries lens 0, then lens 1), by
    // Newton iteration.  A pixel with x < -10 means the ray has none on that lens.
    Vec2 ray_to_pixel(Vec3 ray, int side = -1) const;

    // Lens 0's centre relative to lens 1's, in centimetres.
    Vec3 baseline() const { return baseline_; }
    // The corrections file's rotation of lens 1, degrees, as a rotation vector; zero by default.
    Vec3 lens1_rotation() const { return rot1_deg_; }
    // Whether any lattice correction is non-zero.  If not, rays use the cube in closed form.
    bool has_lattice() const { return !lattice_zero_; }

    // The along-stitch to across-stitch ratio of the warp, from the baseline.
    double warp_ratio() const { return std::fabs(baseline_.x / baseline_.z); }

    // The cube alone, with no corrections.
    static Vec2 cube_ray_to_pixel(Vec3 ray, int side);
    static Vec3 cube_pixel_to_ray(Vec2 px, int lens);

private:
    static const int LEVELS = 3;
    std::vector<Vec3> base_;        // the cube's corners at level 0, zero elsewhere
    std::vector<Vec3> data_;        // base_ plus the corrections file's lattice values
    std::vector<int> tiles_[2][LEVELS][5];      // indices into data_
    bool lattice_zero_ = true;
    Vec3 baseline_, rot1_deg_;
    double rot1_[9];                // lens 1's rotation matrix, identity by default
    std::string error_;

    void reset();
    void set_rotation(Vec3 deg);
    Vec3 rotate1(Vec3 v, bool inverse) const;
    Vec3 spline(int side, int tile, double cx, double cy, Vec3* gx, Vec3* gy) const;
    Vec2 ray_to_pixel_on(Vec3 ray, int side) const;
};

class StitchWarp {
public:
    // `disp` is that frame's DISP_LEN disparities; `rho` is CameraModel::warp_ratio().
    static Vec2 recorded_to_model(Vec2 recorded_px, const float* disp, double rho);
    // Solved by iteration.  Near the stitch a model pixel may have no recorded source
    // (the camera discarded it); the result is then the closest attempt, not a solution,
    // so check it against recorded_to_model if that matters.
    static Vec2 model_to_recorded(Vec2 model_px, int lens, const float* disp, double rho);
};

}  // namespace maxcam
