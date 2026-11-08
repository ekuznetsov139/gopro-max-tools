#include "maxcam.h"

#include <algorithm>
#include <cstdio>

namespace maxcam {

namespace {

const int H = TILE_H, W = FRAME_W, Wz = SIDE_W;
const double PI = 3.14159265358979323846;
const double HW = 656.;     // px from a C edge to the plane between the lenses

// z of a side tile's outer edge on the cube, 32 px past that plane.
double cube_outer_z() { return std::tan((1 - Wz / HW) * PI / 4); }

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
// All four quadrants of a lens's outer tiles unrolled into one strip: x is the distance
// from the stitch, y runs along the whole stitch (4 * 1344 long, wrapping).  x is signed: a
// warped position can lie past the tile's outer edge, still on the same lens.
Vec2 to_stripped(Vec2 px, int cam) {
    double sign = cam ? 1 : -1;
    if (px.x >= 2048.) {
        if (px.y >= 1343.5)
            return Vec2((px.x - 3407.5) * -sign, 3 * 1344 + 2687.5 - px.y);
        return Vec2((px.x - 3407.5) * sign, 0.5 + px.y);
    }
    if (px.y >= 1343.5)
        return Vec2((px.x - 687.5) * sign, 1344 + px.y - 1343.5);
    return Vec2((px.x - 687.5) * -sign, 2 * 1344 + 1343.5 - px.y);
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

// The warp of docs/FORMAT.md.
const double WARP_WIDTH = 344., WARP_NEUTRAL = 66.77, WARP_GAIN = 88.4, INTERP = 0.333;
// Z(s), the warp's falloff with distance s from the stitch: a parabola reaching zero at Z_END,
// continued below Z_KNEE by its tangent line.  About 1 at the stitch.
const double Z_END = 332., Z_KNEE = 247., Z_CURV = 2.03e-5;

double falloff(double s) {
    if (s >= Z_END)
        return 0;
    if (s < 0)
        s = 0;      // past the outer edge, the falloff stays at its edge value
    if (s >= Z_KNEE)
        return Z_CURV * (Z_END - s) * (Z_END - s);
    return Z_CURV * (Z_END - Z_KNEE) * ((Z_END - Z_KNEE) + 2 * (Z_KNEE - s));
}

Vec2 warp_stripped(Vec2 s, const float* disp, int cam, double rho) {
    if (s.x >= WARP_WIDTH)
        return s;
    double z = falloff(s.x);

    double phi = s.y / (4 * 1344) * DISP_LEN + 80.;
    if (phi >= DISP_LEN)
        phi -= DISP_LEN;
    int i = std::min(std::max(int(phi), 0), DISP_LEN - 1);
    auto v = [&](int j) {
        return INTERP * disp[(j + DISP_LEN) % DISP_LEN]
            + (1. - INTERP) * (disp[(j - 1 + DISP_LEN) % DISP_LEN] + disp[(j + 1) % DISP_LEN]) / 2;
    };
    double d = v(i) * (1. - (phi - i)) + v(i + 1) * (phi - i);
    double a = (WARP_NEUTRAL - d) / 80. * z;
    double psi = 2 * PI * i / DISP_LEN;

    Vec2 out(s.x - WARP_GAIN * a / std::max(std::fabs(std::cos(psi)), std::fabs(std::sin(psi))),
             s.y + (cam ? -1 : 1) * WARP_GAIN * rho * a * std::cos(psi));
    if (out.y < 0)
        out.y += 1344 * 4;
    if (out.y >= 1344 * 4)
        out.y -= 1344 * 4;
    return out;
}

// Which tile a pixel falls in on a given lens.  A pixel across a stitch from its own lens's tile
// stays on that lens: the tile's fractions then run outside 0..1.
struct TilePos {
    int plane, col, x0, w;
    double fx, fy;
};

TilePos locate(Vec2 v, int lens) {
    double x = v.x + 0.5, y = v.y + 0.5;
    TilePos t;
    t.plane = y >= H ? 1 : 0;
    t.col = locate_column(x, t.x0, t.w);
    if (t.col != 2 && TILE_BLOCKS[t.plane * 5 + t.col][0] != lens) {
        static const int across[5] = {1, 0, 2, 4, 3};
        t.col = across[t.col];
        static const int left[5] = {0, Wz, 2 * Wz, W - 2 * Wz, W - Wz};
        t.x0 = left[t.col];
    }
    t.fx = (x - t.x0) / t.w;
    t.fy = (y - (t.plane ? H : 0)) / H;
    return t;
}

// Tile fractions -> the fractions at which a bilinear patch between the cube's corners gives the
// cube's ray, with their derivatives.
void cube_fractions(const TilePos& t, double& gx, double& gy, double& dgx, double& dgy) {
    double angle = (1 - 2 * t.fy) * PI / 4, tn = std::tan(angle);
    gy = (1 - tn) / 2;
    dgy = (PI / 4) * (1 + tn * tn);
    if (t.col == 2) {
        tn = std::tan((2 * t.fx - 1) * PI / 4);
        gx = (1 + tn) / 2;
        dgx = (PI / 4) * (1 + tn * tn);
        return;
    }
    bool inner_left = t.col == 0 || t.col == 3;     // the C edge is on the tile's left
    double n = Wz * (inner_left ? t.fx : 1 - t.fx);
    double q = std::tan((1 - n / HW) * PI / 4), qe = cube_outer_z();
    gx = inner_left ? (1 - q) / (1 - qe) : (q - qe) / (1 - qe);
    dgx = (PI / (4 * HW)) * (1 + q * q) * Wz / (1 - qe);
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

Vec3 CameraModel::spline(int side, int tile, double cx, double cy, Vec3* pgx, Vec3* pgy) const {
    // (cx, cy) drive the bilinear part, which carries on along the cube face past the tile's edge;
    // the bumps have no support there, so they stop at the edge.
    double fx = clampd(cx, 0, 1), fy = clampd(cy, 0, 1);
    const std::vector<int>& t0 = tiles_[side][0][tile];
    auto at = [&](const std::vector<int>& t, int i) -> const Vec3& { return data_[t[i]]; };

    // Level 0: bilinear between the tile's corners, and quadratic bumps from its edge
    // midpoints and centre.
    const Vec3 &c0 = at(t0, 0), &c1 = at(t0, 2), &c2 = at(t0, 6), &c3 = at(t0, 8);
    const Vec3 &cx0 = at(t0, 3), &cx1 = at(t0, 5), &cy0 = at(t0, 1), &cy1 = at(t0, 7), &cc = at(t0, 4);
    Vec3 def = (c0 * (1 - cx) + c1 * cx) * (1 - cy) + (c2 * (1 - cx) + c3 * cx) * cy;
    Vec3 gx, gy;
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
        // A clamped coordinate no longer moves the bumps.
        *pgx = (gx * (fx == cx ? 1.0 : 0.0)) + (c1 - c0) * (1 - cy) + (c3 - c2) * cy;
        *pgy = (gy * (fy == cy ? 1.0 : 0.0)) + (c2 - c0) * (1 - cx) + (c3 - c1) * cx;
    }
    return def;
}

Vec3 CameraModel::pixel_to_ray(Vec2 v, Vec3* pgx, Vec3* pgy) const {
    return pixel_to_ray(v, lens_index(v), pgx, pgy);
}

Vec3 CameraModel::pixel_to_ray(Vec2 v, int lens, Vec3* pgx, Vec3* pgy) const {
    TilePos t = locate(v, lens);
    double gx, gy, dgx, dgy;
    cube_fractions(t, gx, gy, dgx, dgy);
    const int* block = TILE_BLOCKS[t.plane * 5 + t.col];
    Vec3 ray = spline(block[0], block[1], gx, gy, pgx, pgy);
    if (pgx) {
        *pgx = *pgx * (dgx / t.w);
        *pgy = *pgy * (dgy / H);
    }
    return ray;
}

Vec2 CameraModel::cube_ray_to_pixel(Vec3 ray, int side) {
    if (side < 0)
        side = ray.z > 0 ? 0 : 1;
    if (side == 1) {
        std::swap(ray.x, ray.y);
        ray.z *= -1;
    }
    const double to_u = 4 / PI;
    if (ray.z > std::fabs(ray.x) && ray.z > std::fabs(ray.y)) {
        double fx = (to_u * std::atan(ray.x / ray.z) + 1) / 2, fy = (1 - to_u * std::atan(ray.y / ray.z)) / 2;
        return Vec2(2 * Wz + H * fx - 0.5, (side ? H : 0) + H * fy - 0.5);
    }
    if (std::fabs(ray.x) > std::fabs(ray.y)) {
        double ax = std::fabs(ray.x), fy = (1 - to_u * std::atan(ray.y / ax)) / 2;
        double n = HW * (1 - to_u * std::atan(ray.z / ax));
        return Vec2(ray.x > 0 ? W - 2 * Wz + n - 0.5 : 2 * Wz - n - 0.5, (side ? H : 0) + H * fy - 0.5);
    }
    double ay = std::fabs(ray.y), fy = (1 - to_u * std::atan(ray.x / ay)) / 2;
    double n = HW * (1 - to_u * std::atan(ray.z / ay));
    return Vec2(ray.y > 0 ? W - n - 0.5 : n - 0.5, (side ? 0 : H) + H * fy - 0.5);
}

Vec3 CameraModel::cube_pixel_to_ray(Vec2 v, int lens) {
    TilePos t = locate(v, lens);
    double b = std::tan((1 - 2 * t.fy) * PI / 4);
    Vec3 ray;
    if (t.col == 2) {
        ray = Vec3(std::tan((2 * t.fx - 1) * PI / 4), b, 1);
    } else {
        double n = Wz * ((t.col == 0 || t.col == 3) ? t.fx : 1 - t.fx);
        double q = std::tan((1 - n / HW) * PI / 4);
        ray = Vec3(t.col < 2 ? -1 : 1, b, (t.col == 1 || t.col == 3) ? q : -q);
    }
    if (t.plane) {
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
    Vec2 px = cube_ray_to_pixel(ray, side);
    for (int pass = 0; pass <= 100; pass++) {
        Vec3 gx, gy;
        Vec3 def = pixel_to_ray(px, side, &gx, &gy);
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
        TilePos a = locate(px, side), b = locate(next, side);
        if (next.x < -0.5 || next.y < -0.5 || next.x > 4095.5 || next.y > 2687.5 || a.plane != b.plane || a.col != b.col) {
            // Newton would leave this tile, where the model has a kink: step on the cube instead.
            px = cube_ray_to_pixel((cube_pixel_to_ray(px, side).normalized() + err).normalized(), side);
        } else {
            px = next;
        }
    }
    return Vec2(-1e5, -1e5);
}

Vec3 CameraModel::recorded_to_ray(Vec2 px, const float* disp) const {
    return pixel_to_ray(StitchWarp::recorded_to_model(px, disp, warp_ratio()), lens_index(px)).normalized();
}

Vec2 CameraModel::ray_to_recorded(Vec3 ray, int lens, const float* disp) const {
    Vec2 m = ray_to_pixel(ray, lens);
    return m.x < -10 ? m : StitchWarp::model_to_recorded(m, lens, disp, warp_ratio());
}

// ---- StitchWarp -------------------------------------------------------------
Vec2 StitchWarp::recorded_to_model(Vec2 px, const float* disp, double rho) {
    int cam = lens_index(px);
    if (locate(px, cam).col == 2)
        return px;
    return from_stripped(warp_stripped(to_stripped(px, cam), disp, cam, rho), cam);
}

Vec2 StitchWarp::model_to_recorded(Vec2 px, int lens, const float* disp, double rho) {
    if (locate(px, lens).col == 2)
        return px;
    Vec2 target = to_stripped(px, lens), guess = target;
    double gain = 1, last = 1e30;
    for (int pass = 0; pass < 200; pass++) {
        Vec2 error = target - warp_stripped(guess, disp, lens, rho);
        if (error.y > 2 * 1344)
            error.y -= 4 * 1344;
        else if (error.y < -2 * 1344)
            error.y += 4 * 1344;
        if (error.norm2() < 1e-16)
            break;
        if (error.norm2() > last)
            gain = std::max(gain / 2, 0.05);    // steep DISP steps can make the plain iteration overshoot
        last = error.norm2();
        guess = guess + Vec2(error.x * gain, error.y * gain);
        if (guess.y < 0)
            guess.y += 1344 * 4;
        else if (guess.y >= 1344 * 4)
            guess.y -= 1344 * 4;
        if (guess.x < 0)
            guess.x = 0;    // recorded pixels stop at the tile's edge
    }
    return from_stripped(guess, lens);
}

}  // namespace maxcam
