#include "maxcam.h"

#include <algorithm>
#include <cstdio>

namespace maxcam {

namespace {

const int H = TILE_H, W = FRAME_W, Wz = SIDE_W;

double clampd(double x, double lo, double hi) { return std::min(std::max(x, lo), hi); }

// Which tile block of which lens each frame tile is: {lens, block}.
const int TILE_BLOCKS[10][2] = {
    {1, 4}, {0, 0}, {0, 1}, {0, 2}, {1, 3},
    {0, 4}, {1, 0}, {1, 1}, {1, 2}, {0, 3},
};

// The tile (column) a pixel's x falls in, with that tile's left edge and width.
int locate_column(double x, int& x0, int& w) {
    if (x < Wz) { x0 = 0; w = Wz; return 0; }
    if (x < Wz * 2) { x0 = Wz; w = Wz; return 1; }
    if (x < W - Wz * 2) { x0 = Wz * 2; w = H; return 2; }
    if (x < W - Wz) { x0 = W - Wz * 2; w = Wz; return 3; }
    x0 = W - Wz; w = Wz;
    return 4;
}

// ---- Spline weights ---------------------------------------------------------
inline void edge_effects(double fx, double fy, double w[4]) {
    double t1 = 1.0 - 4.0 * (fx - 0.5) * (fx - 0.5);
    double t2 = 1.0 - 4.0 * (fy - 0.5) * (fy - 0.5);
    w[0] = t1 * (1.0 - fy) * (1.0 - fy);
    w[1] = t1 * fy * fy;
    w[2] = t2 * (1.0 - fx) * (1.0 - fx);
    w[3] = t2 * fx * fx;
}

inline void edge_effects_grad(double fx, double fy, double wgx[4], double wgy[4]) {
    double t1 = 1.0 - 4.0 * (fx - 0.5) * (fx - 0.5);
    double t2 = 1.0 - 4.0 * (fy - 0.5) * (fy - 0.5);
    double t1x = -8.0 * (fx - 0.5), t2y = -8.0 * (fy - 0.5);
    double n1 = (1.0 - fy) * (1.0 - fy), n1y = -2.0 * (1.0 - fy);
    double n2 = fy * fy, n2y = 2.0 * fy;
    double n3 = (1.0 - fx) * (1.0 - fx), n3x = -2.0 * (1.0 - fx);
    double n4 = fx * fx, n4x = 2.0 * fx;
    wgx[0] = t1x * n1; wgx[1] = t1x * n2; wgx[2] = t2 * n3x; wgx[3] = t2 * n4x;
    wgy[0] = t1 * n1y; wgy[1] = t1 * n2y; wgy[2] = t2y * n3; wgy[3] = t2y * n4;
}

inline double bump(double x, double pos, double scale) {
    x = (x - pos) * scale;
    if (x <= -1 || x >= 1)
        return 0;
    x = 1 - x * x;
    return x * x;
}

inline double bump_grad(double x, double pos, double scale) {
    x = (x - pos) * scale;
    if (x <= -1 || x >= 1)
        return 0;
    return 4 * x * (x * x - 1) * scale;
}

inline double edge_correction(double x, double y, bool outer_edge, Vec2* grad) {
    double t = 1.0 - 4.0 * (x - 0.5) * (x - 0.5);
    double wx = t * t;
    if (y >= 0.5) {
        *grad = Vec2(0, 0);
        return 0;
    }
    double wy;
    if (outer_edge)
        wy = (1.0 - y * 2) * (1.0 - y * 2);
    else
        wy = (1.0 - y * y * 4) * (1.0 - y * y * 4);
    grad->x = -2 * t * (4.0 * 2 * (x - 0.5)) * wy;
    if (outer_edge)
        grad->y = -wx * 2. * (1.0 - y * 2) * 2.;
    else
        grad->y = -wx * 2. * (1.0 - y * y * 4) * 2 * y * 4;
    return wx * wy;
}

// ---- The stitch warp, in "stripped" coordinates -----------------------------
// All four quadrants of the outer tiles unrolled into one strip: x is the distance
// from the stitch, y runs along the whole stitch (4 * 1344 long, wrapping).
Vec2 to_stripped(Vec2 px) {
    if (px.x >= 2048.) {
        if (px.y >= 1343.5)
            return Vec2(std::fabs(px.x - 3407.5), 3 * 1344 + 2687.5 - px.y);
        return Vec2(std::fabs(px.x - 3407.5), 0.5 + px.y);
    }
    if (px.y >= 1343.5)
        return Vec2(std::fabs(px.x - 687.5), 1344 + px.y - 1343.5);
    return Vec2(std::fabs(px.x - 687.5), 2 * 1344 + 1343.5 - px.y);
}

Vec2 from_stripped(Vec2 s, int cam) {
    if (s.y < 1344)
        return Vec2(3407.5 + s.x * (cam ? 1 : -1), s.y - 0.5);
    if (s.y < 2688)
        return Vec2(687.5 + s.x * (cam ? 1 : -1), s.y - 0.5);
    if (s.y < 2688 + 1344)
        return Vec2(687.5 + s.x * (cam ? -1 : 1), 1343.5 - (s.y - 2688));
    return Vec2(3407.5 + s.x * (cam ? -1 : 1), 2687.5 - (s.y - 2688 - 1344));
}

// Fitted constants.  x_zone is the warp's falloff with distance into the strip; the
// three weights are the c1, c2, c3 of docs/FORMAT.md.
const float DISP_SCALE_X = 69.12f, DISP_SCALE_Y = 59.875f;
const float WARP_WIDTH = 344.f, MAX_DISP = 80.f;
const double WEIGHT_X0 = 1.279495, WEIGHT_X1 = -0.188683, WEIGHT_Y = 1.002468;
const double INTERP = 0.333;
const float X_ZONE[17] = {
    0.002021f, 0.021586f, 0.051468f, 0.099344f, 0.161015f, 0.233702f, 0.303061f, 0.378131f, 0.449792f,
    0.526414f, 0.591500f, 0.665341f, 0.727919f, 0.792834f, 0.860367f, 0.924303f, 1.0f,
};

inline double interpolate4(const double* p, double dx) {
    double v1 = INTERP * p[1] + (1. - INTERP) * (p[0] * 0.5 + p[2] * 0.5);
    double v2 = INTERP * p[2] + (1. - INTERP) * (p[1] * 0.5 + p[3] * 0.5);
    return v1 * (1. - dx) + v2 * dx;
}

Vec2 warp_stripped(Vec2 s, const float* vdisp, int cam) {
    float phase = (s.y / (4 * 1344)) * 640. + 80.;
    if (phase >= 640.)
        phase -= 640.;
    float fx = 1. - s.x / WARP_WIDTH;
    float compression = (1. - 1.30 * DISP_SCALE_X / 688.);
    compression = (688. - s.x) * (1.0 - compression);
    if (fx < 0.0)
        return s + Vec2(compression, 0);
    int ifx = std::min(15, int(fx * 16));
    float delta_fx = fx * 16 - ifx;
    fx = X_ZONE[ifx] * (1. - delta_fx) + X_ZONE[ifx + 1] * delta_fx;

    int iphase = std::min(std::max(int(phase), 0), 639);
    double ldisp[4];
    for (int i = 0; i < 4; i++)
        ldisp[i] = vdisp[(iphase - 1 + i + DISP_LEN) % DISP_LEN];
    float disp = interpolate4(ldisp, phase - iphase);

    double x_weight = WEIGHT_X0 + WEIGHT_X1 * std::cos((iphase / 160.) * 2 * 3.1415926);
    float y_weight = WEIGHT_Y * std::cos((iphase / 320.) * 3.1415926);
    disp = (1. - disp / MAX_DISP) * fx;

    float dx = -disp * DISP_SCALE_X * x_weight + compression;
    float dy = disp * DISP_SCALE_Y * y_weight * (cam ? -1 : 1);
    s = Vec2(s.x + dx, s.y + dy);
    if (s.y < 0)
        s.y += 1344 * 4;
    if (s.y >= 1344 * 4)
        s.y -= 1344 * 4;
    return s;
}

}  // namespace

// ---- Tiles ------------------------------------------------------------------
int zone_index(Vec2 px) {
    double x = px.x + 0.5, y = px.y + 0.5;
    int x0, w;
    return (y >= H ? 5 : 0) + locate_column(x, x0, w);
}

int lens_index(Vec2 px) {
    return TILE_BLOCKS[zone_index(px)][0];
}

// ---- CameraModel ------------------------------------------------------------
CameraModel::CameraModel() {
    // Three levels of lattice, 2, 4 and 8 cells per tile, over the five tiles of each
    // lens.  Tile blocks 0..2 form a band sharing their edges; blocks 3 and 4 share
    // their edges with that band, so the model is continuous across the whole lens.
    int vars = 0;
    int level_offsets[LEVELS];
    for (int k = 0; k < LEVELS; k++) {
        int d = 2 << k;
        level_offsets[k] = vars;
        vars += 2 * (5 * d * d + 2 * d + 1);
    }
    data_.resize(vars);
    for (int i = 0; i < 2; i++)
        for (int k = 0; k < LEVELS; k++)
            for (int j = 0; j < 5; j++) {
                int d = 2 << k, w = d + 1;
                int base = level_offsets[k] + i * (5 * d * d + 2 * d + 1);
                std::vector<int>& t = tiles_[i][k][j];
                t.resize(w * w);
                for (int y = 0; y < w; y++)
                    for (int x = 0; x < w; x++) {
                        if (j < 3) {
                            t[x + y * w] = base + (x + d * j) + y * (3 * d + 1);
                        } else if (j == 3) {
                            if (x == 0)
                                t[y + (d - x) * w] = tiles_[i][k][0][y];
                            else if (x == d)
                                t[y + (d - x) * w] = tiles_[i][k][2][d - y];
                            else if (y == d)
                                t[y + (d - x) * w] = tiles_[i][k][1][x];
                            else
                                t[y + (d - x) * w] = base + (d + 1) * (3 * d + 1) + y * (d - 1) + (x - 1);
                        } else {
                            if (x == 0)
                                t[y + (d - x) * w] = tiles_[i][k][0][d * (d + 1) + d - y];
                            else if (x == d)
                                t[y + (d - x) * w] = tiles_[i][k][2][d * (d + 1) + y];
                            else if (y == 0)
                                t[y + (d - x) * w] = tiles_[i][k][1][d * (d + 1) + x];
                            else
                                t[y + (d - x) * w] = base + (d + 1) * (3 * d + 1) + d * (d - 1) + (y - 1) * (d - 1) + (x - 1);
                        }
                    }
            }
}

bool CameraModel::load(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "r");
    if (!f) {
        error_ = "cannot open " + path;
        return false;
    }
    bool ok = true;
    for (Vec3& v : data_)
        ok = ok && std::fscanf(f, "%lf %lf %lf\n", &v.x, &v.y, &v.z) == 3;
    // Rotations of the two lens axes.  The fitted files all hold identity; nothing here
    // applies them, so anything else must not be silently ignored.
    for (int i = 0; i < 2 && ok; i++) {
        double q[4];
        ok = std::fscanf(f, "%lf %lf %lf %lf\n", q, q + 1, q + 2, q + 3) == 4;
        if (ok && (q[0] != 1.0 || q[1] != 0.0 || q[2] != 0.0 || q[3] != 0.0)) {
            error_ = path + ": lens axis rotations other than identity are not supported";
            std::fclose(f);
            return false;
        }
    }
    ok = ok && std::fscanf(f, "%lf %lf %lf\n", &baseline_.x, &baseline_.y, &baseline_.z) == 3;
    std::fclose(f);
    if (!ok)
        error_ = path + ": not a calibration file";
    return ok;
}

Vec3 CameraModel::spline(int side, int tile, double fx, double fy, Vec3* pgx, Vec3* pgy) const {
    const std::vector<int>& t0 = tiles_[side][0][tile];
    auto at = [&](const std::vector<int>& t, int i) -> const Vec3& { return data_[t[i]]; };

    // Level 0: bilinear between the tile's corners, and quadratic bumps from its edge
    // midpoints and centre.
    const Vec3 &c0 = at(t0, 0), &c1 = at(t0, 2), &c2 = at(t0, 6), &c3 = at(t0, 8);
    const Vec3 &cx0 = at(t0, 3), &cx1 = at(t0, 5), &cy0 = at(t0, 1), &cy1 = at(t0, 7), &cc = at(t0, 4);
    Vec3 def = (c0 * (1 - fx) + c1 * fx) * (1 - fy) + (c2 * (1 - fx) + c3 * fx) * fy;
    Vec3 gx, gy;
    if (pgx) {
        gx = (c1 - c0) * (1 - fy) + (c3 - c2) * fy;
        gy = (c2 - c0) * (1 - fx) + (c3 - c1) * fx;
    }
    double w[4];
    edge_effects(fx, fy, w);
    double cfx = fx - 0.5, cfy = fy - 0.5;
    double wc = (1.0 - 4.0 * cfx * cfx) * (1.0 - 4.0 * cfy * cfy);
    if (pgx) {
        gx = gx + cc * ((-8.0 * cfx) * (1.0 - 4.0 * cfy * cfy));
        gy = gy + cc * ((1.0 - 4.0 * cfx * cfx) * (-8.0 * cfy));
        double wgx[4], wgy[4];
        edge_effects_grad(fx, fy, wgx, wgy);
        gx = gx + cy0 * wgx[0] + cy1 * wgx[1] + cx0 * wgx[2] + cx1 * wgx[3];
        gy = gy + cy0 * wgy[0] + cy1 * wgy[1] + cx0 * wgy[2] + cx1 * wgy[3];
    }
    def = def + cy0 * w[0] + cy1 * w[1] + cx0 * w[2] + cx1 * w[3] + cc * wc;

    // Levels 1 and 2: compact bumps, local to one block of the finer lattice.
    for (int k = 1; k <= 2; k++) {
        int d = 2 << k, lw = d + 1;
        const std::vector<int>& t = tiles_[side][k][tile];
        int bx = std::min(std::max(2 * (int)std::floor(fx * d / 2), 0), d - 2);
        int by = std::min(std::max(2 * (int)std::floor(fy * d / 2), 0), d - 2);
        double dfx = (fx * d - bx) / 2, dfy = (fy * d - by) / 2;

        const Vec3 &ky0 = at(t, bx + 1 + by * lw), &ky1 = at(t, bx + 1 + (by + 2) * lw);
        const Vec3 &kx0 = at(t, bx + (by + 1) * lw), &kx1 = at(t, bx + 2 + (by + 1) * lw);
        const Vec3 &kc = at(t, bx + 1 + (by + 1) * lw);
        const Vec3 &c00 = at(t, bx + by * lw), &c01 = at(t, bx + 2 + by * lw);
        const Vec3 &c10 = at(t, bx + (by + 2) * lw), &c11 = at(t, bx + 2 + (by + 2) * lw);

        def = def + c00 * (bump(dfx, 0.0, 2.0) * bump(dfy, 0.0, 2.0));
        def = def + c01 * (bump(dfx, 1.0, 2.0) * bump(dfy, 0.0, 2.0));
        def = def + c10 * (bump(dfx, 0.0, 2.0) * bump(dfy, 1.0, 2.0));
        def = def + c11 * (bump(dfx, 1.0, 2.0) * bump(dfy, 1.0, 2.0));
        def = def + kc * (bump(dfx, 0.5, 2.0) * bump(dfy, 0.5, 2.0));
        if (pgx) {
            double s = d / 2.;
            gx = gx + kc * (bump_grad(dfx, 0.5, 2.0) * bump(dfy, 0.5, 2.0) * s);
            gy = gy + kc * (bump(dfx, 0.5, 2.0) * bump_grad(dfy, 0.5, 2.0) * s);
            gx = gx + c00 * (bump_grad(dfx, 0.0, 2.0) * bump(dfy, 0.0, 2.0) * s);
            gy = gy + c00 * (bump(dfx, 0.0, 2.0) * bump_grad(dfy, 0.0, 2.0) * s);
            gx = gx + c01 * (bump_grad(dfx, 1.0, 2.0) * bump(dfy, 0.0, 2.0) * s);
            gy = gy + c01 * (bump(dfx, 1.0, 2.0) * bump_grad(dfy, 0.0, 2.0) * s);
            gx = gx + c10 * (bump_grad(dfx, 0.0, 2.0) * bump(dfy, 1.0, 2.0) * s);
            gy = gy + c10 * (bump(dfx, 0.0, 2.0) * bump_grad(dfy, 1.0, 2.0) * s);
            gx = gx + c11 * (bump_grad(dfx, 1.0, 2.0) * bump(dfy, 1.0, 2.0) * s);
            gy = gy + c11 * (bump(dfx, 1.0, 2.0) * bump_grad(dfy, 1.0, 2.0) * s);
        }
        Vec2 egy0, egy1, egx0, egx1;
        double ecy0 = edge_correction(dfx, dfy, by == 0, &egy0);
        double ecy1 = edge_correction(dfx, 1. - dfy, by + 2 == d, &egy1);
        double ecx0 = edge_correction(dfy, dfx, bx == 0, &egx0);
        double ecx1 = edge_correction(dfy, 1. - dfx, bx + 2 == d, &egx1);
        def = def + ky0 * ecy0 + ky1 * ecy1 + kx0 * ecx0 + kx1 * ecx1;
        if (pgx) {
            double s = d / 2.;
            gx = gx + (ky0 * egy0.x + ky1 * egy1.x + kx0 * egx0.y - kx1 * egx1.y) * s;
            gy = gy + (ky0 * egy0.y - ky1 * egy1.y + kx0 * egx0.x + kx1 * egx1.x) * s;
        }
    }
    if (pgx) {
        *pgx = gx;
        *pgy = gy;
    }
    return def;
}

Vec3 CameraModel::pixel_to_ray(Vec2 v, Vec3* pgx, Vec3* pgy) const {
    double x = v.x + 0.5, y = v.y + 0.5;
    int plane = y >= H ? 1 : 0;
    int x0, w;
    int xz = locate_column(x, x0, w);
    x = (x - x0) / w;
    y = (y - (plane ? H : 0)) / H;
    const int* block = TILE_BLOCKS[plane * 5 + xz];
    Vec3 ray = spline(block[0], block[1], x, y, pgx, pgy);
    if (pgx) {
        *pgx = *pgx * (1. / w);
        *pgy = *pgy * (1. / H);
    }
    return ray;
}

Vec2 CameraModel::box_ray_to_pixel(Vec3 ray, int side) {
    if (side < 0)
        side = ray.z > 0 ? 0 : 1;
    else if (side == 0)
        ray.z = std::max(0.0, ray.z);
    else
        ray.z = std::min(0.0, ray.z);
    if (side == 1) {
        std::swap(ray.x, ray.y);
        ray.z *= -1;
    }
    if (ray.z > std::fabs(ray.x) && ray.z > std::fabs(ray.y)) {
        ray.x /= ray.z;
        ray.y /= ray.z;
        return Vec2(2 * Wz + H * (ray.x + 1) / 2. - 0.5, (side ? H : 0) + H * (1 - ray.y) / 2. - 0.5);
    }
    if (std::fabs(ray.x) > std::fabs(ray.y)) {
        double fy = clampd((1 - ray.y / std::fabs(ray.x)) / 2., 1e-8, 1 - 1e-8);
        double fz = clampd(ray.z / std::fabs(ray.x), 1e-8, 1 - 1e-8);
        Vec2 px(Wz + Wz * fz - 0.5, (side ? H : 0) + H * fy - 0.5);
        if (ray.x > 0)
            px.x = W - 1 - px.x;
        return px;
    }
    double fx = clampd((1 - ray.x / std::fabs(ray.y)) / 2., 1e-8, 1 - 1e-8);
    double fz = clampd(1. - ray.z / std::fabs(ray.y), 1e-8, 1 - 1e-8);
    Vec2 px(Wz * fz - 0.5, (side ? 0 : H) + H * fx - 0.5);
    if (ray.y > 0)
        px.x = W - 1 - px.x;
    return px;
}

Vec3 CameraModel::box_pixel_to_ray(Vec2 v) {
    double x = v.x + 0.5, y = v.y + 0.5;
    int plane = y >= H ? 1 : 0;
    int x0, w;
    int xz = locate_column(x, x0, w);
    x = (x - x0) / w;
    y = (y - (plane ? H : 0)) / H;
    Vec3 ray;
    if (xz == 2)
        ray = Vec3(-1 + x * 2, 1 - y * 2, 1);
    else if (xz == 1)
        ray = Vec3(-1, 1 - y * 2, x);
    else if (xz == 3)
        ray = Vec3(1, 1 - y * 2, 1 - x);
    else if (xz == 0)
        ray = Vec3(-1, 1 - y * 2, -1 + x);
    else
        ray = Vec3(1, 1 - y * 2, -x);
    if (plane) {
        std::swap(ray.x, ray.y);
        ray.z *= -1;
    }
    return ray;
}

Vec2 CameraModel::ray_to_pixel(Vec3 ray, int side) const {
    if (side < 0) {
        Vec2 px = ray_to_pixel_on(ray, 0);
        return px.x >= 0 ? px : ray_to_pixel_on(ray, 1);
    }
    return ray_to_pixel_on(ray, side);
}

Vec2 CameraModel::ray_to_pixel_on(Vec3 ray, int side) const {
    ray = ray.normalized();
    if ((side == 0 && ray.z < -0.20) || (side == 1 && ray.z > 0.20))
        return Vec2(-1e6, -1e6);
    Vec3 intl = ray;
    intl.z = side == 0 ? std::max(intl.z, 0.0) : std::min(intl.z, 0.0);
    Vec2 px = box_ray_to_pixel(intl, side);
    for (int pass = 0; pass <= 100; pass++) {
        Vec3 gx, gy;
        Vec3 def = pixel_to_ray(px, &gx, &gy);
        double len = def.norm();
        gx = gx * (1. / len);
        gy = gy * (1. / len);
        def = def * (1. / len);
        gx = gx - def * gx.dot(def);
        gy = gy - def * gy.dot(def);
        Vec3 err = ray - def;
        if (err.norm() < 1e-8)
            return px;
        Vec2 rhs(err.dot(gx), err.dot(gy));
        double xx = gx.dot(gx), xy = gx.dot(gy), yy = gy.dot(gy);
        double det = xx * yy - xy * xy;
        Vec2 step((rhs.x * yy - rhs.y * xy) / det, -(rhs.x * xy - rhs.y * xx) / det);
        Vec2 next = px + step;
        if (next.x < -0.5 || next.y < -0.5 || next.x > 4095.5 || next.y > 2687.5 || zone_index(px) != zone_index(next)) {
            // Newton would leave this tile, where the model has a kink: step on the box instead.
            intl = box_pixel_to_ray(px) + err;
            intl.z = side == 0 ? std::max(intl.z, 0.0) : std::min(intl.z, 0.0);
            intl = intl.normalized();
            px = box_ray_to_pixel(intl, side);
        } else {
            px = next;
        }
    }
    return Vec2(-1e5, -1e5);
}

// ---- StitchWarp -------------------------------------------------------------
Vec2 StitchWarp::recorded_to_model(Vec2 px, const float* disp) {
    if (px.x >= 688 * 2 && px.x < 3408 - 688)
        return px;
    int cam = lens_index(px);
    return from_stripped(warp_stripped(to_stripped(px), disp, cam), cam);
}

Vec2 StitchWarp::model_to_recorded(Vec2 px, const float* disp) {
    if (px.x >= 688 * 2 - 0.5 && px.x < 3408 - 688 - 0.5)
        return px;
    int cam = lens_index(px);
    Vec2 target = to_stripped(px);
    Vec2 guess = target, last_error, last_move, last_guess = guess;
    int unroll = 0;
    for (int pass = 0; pass <= 100; pass++) {
        Vec2 error = target - warp_stripped(guess, disp, cam);
        if (last_error.norm2() != 0 && error.norm2() > last_error.norm2()) {
            guess = last_guess;
            unroll++;
            continue;
        } else if (error.norm2() < last_error.norm2()) {
            unroll = 0;
        }
        if (error.norm2() < 1e-8)
            break;
        Vec2 move = error;
        if (pass > 0 && last_error.norm2() < 1000 && error.norm2() < 1000) {
            double x_scale = 0.0, y_scale = 0.0;
            if (std::fabs(last_move.x) > 1e-6)
                x_scale = (error.x - last_error.x) / (-last_move.x);
            if (std::fabs(last_move.y) > 1e-6)
                y_scale = (error.y - last_error.y) / (-last_move.y);
            if (std::fabs(x_scale - y_scale) < 0.2) {
                Vec2 r1 = last_move, r2 = error - last_error;
                double n1 = std::sqrt(r1.norm2()), n2 = std::sqrt(r2.norm2());
                // Either can be zero, where the warp is locally flat: no rotation to correct.
                if (n1 > 0 && n2 > 0) {
                    r1 = Vec2(r1.x / n1, r1.y / n1);
                    r2 = Vec2(r2.x / n2, r2.y / n2);
                    double phi = clampd(r2.x * r1.y - r2.y * r1.x, -0.25, 0.25);
                    double cp = std::sqrt(1. - phi * phi);
                    move = Vec2(move.x * cp - move.y * phi, move.x * phi + move.y * cp);
                }
            }
            x_scale = clampd(x_scale, 0.5, 6.0);
            y_scale = clampd(y_scale, 0.5, 6.0);
            move.x /= x_scale * (unroll ? 2 : 1);
            move.y /= y_scale * (unroll ? 2 : 1);
        }
        last_guess = guess;
        guess = guess + move;
        last_move = move;
        last_error = error;
        if (guess.y < 0)
            guess.y += 1344 * 4;
        else if (guess.y >= 1344 * 4)
            guess.y -= 1344 * 4;
        if (guess.x < 0)
            guess.x = 0;
    }
    return from_stripped(guess, cam);
}

}  // namespace maxcam
