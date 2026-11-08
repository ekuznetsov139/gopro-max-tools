#include "maxcam.h"

#include <algorithm>
#include <cstdio>
#include <cstdint>
#include <cstring>

namespace maxcam {

namespace {

const int H = TILE_H, W = FRAME_W, Wz = SIDE_W;
const double PI = 3.14159265358979323846;
const double HW = 656.;     // px from a C edge to the plane between the lenses
// The baseline when no .360 file is given: lens 0 relative to lens 1, cm, as the camera tested
// reports it (the negative of its header's SHFX/Y/Z, in mm), with the 0.02 mm of y dropped.
const Vec3 DEFAULT_BASELINE(-2.336, 0., 3.580);

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

// The cube's ray at tile fractions (fx, fy) of a frame tile, as cube_pixel_to_ray gives it.
Vec3 cube_ray_at(int plane, int col, double fx, double fy) {
    double b = std::tan((1 - 2 * fy) * PI / 4);
    Vec3 ray;
    if (col == 2) {
        ray = Vec3(std::tan((2 * fx - 1) * PI / 4), b, 1);
    } else {
        double n = Wz * ((col == 0 || col == 3) ? fx : 1 - fx);
        double q = std::tan((1 - n / HW) * PI / 4);
        ray = Vec3(col < 2 ? -1 : 1, b, (col == 1 || col == 3) ? q : -q);
    }
    if (plane) {
        std::swap(ray.x, ray.y);
        ray.z *= -1;
    }
    return ray;
}

// ---- The .360 header: moov/udta/GPMF, a KLV tree; the DEVC whose DVID is GEOC holds SHFX/Y/Z --
uint32_t be32(const unsigned char* p) { return (uint32_t(p[0]) << 24) | (uint32_t(p[1]) << 16) | (uint32_t(p[2]) << 8) | p[3]; }
uint64_t be64(const unsigned char* p) { return (uint64_t(be32(p)) << 32) | be32(p + 4); }
double be_double(const unsigned char* p) {
    uint64_t u = be64(p);
    double d;
    std::memcpy(&d, &u, 8);
    return d;
}

// Finds a box by type among the children in [start, end) of a file; returns its payload range.
bool find_box(FILE* f, uint64_t start, uint64_t end, const char* type, uint64_t& pstart, uint64_t& pend) {
    unsigned char hdr[16];
    for (uint64_t off = start; off + 8 <= end;) {
        if (std::fseek(f, (long)off, SEEK_SET) != 0 || std::fread(hdr, 1, 8, f) != 8)
            return false;
        uint64_t size = be32(hdr), hs = 8;
        if (size == 1) {
            if (std::fread(hdr + 8, 1, 8, f) != 8)
                return false;
            size = be64(hdr + 8);
            hs = 16;
        } else if (size == 0) {
            size = end - off;
        }
        if (size < hs)
            return false;
        if (std::memcmp(hdr + 4, type, 4) == 0) {
            pstart = off + hs;
            pend = off + size;
            return true;
        }
        off += size;
    }
    return false;
}

// Walks a KLV buffer for the GEOC device's SHFX/Y/Z.  Returns how many of the three were found.
int find_shf(const unsigned char* b, size_t n, double shf[3], bool in_geoc) {
    int found = 0;
    for (size_t off = 0; off + 8 <= n;) {
        const unsigned char* k = b + off;
        unsigned type = k[4], size = k[5], rep = (k[6] << 8) | k[7];
        size_t len = size_t(size) * rep, pad = (len + 3) & ~size_t(3);
        if (off + 8 + len > n)
            break;
        const unsigned char* pay = k + 8;
        if (type == 0) {
            // A container: a DEVC whose DVID is GEOC is the one; look inside either way.
            bool geoc = in_geoc;
            for (size_t o = 0; o + 8 <= len;) {
                const unsigned char* kk = pay + o;
                size_t ll = size_t(kk[5]) * ((kk[6] << 8) | kk[7]);
                if (std::memcmp(kk, "DVID", 4) == 0 && kk[4] == 'F' && ll == 4 && std::memcmp(kk + 8, "GEOC", 4) == 0)
                    geoc = true;
                o += 8 + ((ll + 3) & ~size_t(3));
            }
            found += find_shf(pay, len, shf, geoc);
        } else if (in_geoc && type == 'd' && len == 8) {
            if (std::memcmp(k, "SHFX", 4) == 0) { shf[0] = be_double(pay); found++; }
            if (std::memcmp(k, "SHFY", 4) == 0) { shf[1] = be_double(pay); found++; }
            if (std::memcmp(k, "SHFZ", 4) == 0) { shf[2] = be_double(pay); found++; }
        }
        off += 8 + pad;
    }
    return found;
}

// The camera's baseline from a .360 file, cm: the negative of the header's SHFX/Y/Z in mm.
bool read_baseline(const std::string& path, Vec3& baseline, std::string& error) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) {
        error = "cannot open " + path;
        return false;
    }
    std::fseek(f, 0, SEEK_END);
    uint64_t end = (uint64_t)std::ftell(f), a, b, c, d, e, g;
    bool ok = find_box(f, 0, end, "moov", a, b) && find_box(f, a, b, "udta", c, d) && find_box(f, c, d, "GPMF", e, g);
    std::vector<unsigned char> buf;
    if (ok) {
        buf.resize(g - e);
        ok = std::fseek(f, (long)e, SEEK_SET) == 0 && std::fread(buf.data(), 1, buf.size(), f) == buf.size();
    }
    std::fclose(f);
    double shf[3] = {0, 0, 0};
    if (!ok || find_shf(buf.data(), buf.size(), shf, false) != 3) {
        error = path + ": no GEOC calibration (SHFX/Y/Z) in moov/udta/GPMF";
        return false;
    }
    baseline = Vec3(-shf[0] / 10, -shf[1] / 10, -shf[2] / 10);
    return true;
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
    base_.assign(vars, Vec3());
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
    // Level 0's corners carry the cube: spline() reads a tile's corners at its remapped fractions
    // (0,0), (1,0), (0,1), (1,1), which are the tile's own fractions (fx, fy) in {0,1}.
    for (int plane = 0; plane < 2; plane++)
        for (int col = 0; col < 5; col++) {
            const int* block = TILE_BLOCKS[plane * 5 + col];
            const std::vector<int>& t0 = tiles_[block[0]][0][block[1]];
            const int corner[4] = {0, 2, 6, 8};
            for (int q = 0; q < 4; q++)
                base_[t0[corner[q]]] = cube_ray_at(plane, col, q & 1, q >> 1);
        }
    reset();
}

void CameraModel::reset() {
    data_ = base_;
    lattice_zero_ = true;
    baseline_ = DEFAULT_BASELINE;
    set_rotation(Vec3());
    error_.clear();
}

void CameraModel::set_rotation(Vec3 deg) {
    rot1_deg_ = deg;
    Vec3 r = deg * (PI / 180);
    double th = r.norm();
    double R[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    if (th > 0) {
        Vec3 k = r * (1 / th);
        double cs = std::cos(th), sn = std::sin(th), v = 1 - cs;
        double K[9] = {0, -k.z, k.y, k.z, 0, -k.x, -k.y, k.x, 0};
        double kk[3] = {k.x, k.y, k.z};
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                R[3 * i + j] = (i == j ? cs : 0) + sn * K[3 * i + j] + v * kk[i] * kk[j];
    }
    std::memcpy(rot1_, R, sizeof R);
}

Vec3 CameraModel::rotate1(Vec3 v, bool inverse) const {
    const double* R = rot1_;
    if (inverse)
        return Vec3(R[0] * v.x + R[3] * v.y + R[6] * v.z, R[1] * v.x + R[4] * v.y + R[7] * v.z, R[2] * v.x + R[5] * v.y + R[8] * v.z);
    return Vec3(R[0] * v.x + R[1] * v.y + R[2] * v.z, R[3] * v.x + R[4] * v.y + R[5] * v.z, R[6] * v.x + R[7] * v.y + R[8] * v.z);
}

bool CameraModel::load(const std::string& clip_path, const std::string& corrections_path) {
    reset();
    if (!clip_path.empty() && !read_baseline(clip_path, baseline_, error_)) {
        reset();
        return false;
    }
    if (corrections_path.empty())
        return true;
    FILE* f = std::fopen(corrections_path.c_str(), "r");
    if (!f) {
        error_ = "cannot open " + corrections_path;
        reset();
        return false;
    }
    bool ok = true;
    std::vector<Vec3> delta(data_.size());
    for (Vec3& v : delta)
        ok = ok && std::fscanf(f, "%lf %lf %lf\n", &v.x, &v.y, &v.z) == 3;
    Vec3 db, dr;
    ok = ok && std::fscanf(f, "%lf %lf %lf\n", &db.x, &db.y, &db.z) == 3;
    ok = ok && std::fscanf(f, "%lf %lf %lf\n", &dr.x, &dr.y, &dr.z) == 3;
    double extra;
    bool more = std::fscanf(f, "%lf", &extra) == 1;
    std::fclose(f);
    if (!ok || more) {
        error_ = corrections_path + ": not a corrections file of " + std::to_string(data_.size() + 2) + " rows (calibration/README.md)";
        reset();
        return false;
    }
    Vec3 saved_baseline = baseline_;
    for (size_t i = 0; i < data_.size(); i++) {
        data_[i] = base_[i] + delta[i];
        if (delta[i].x != 0 || delta[i].y != 0 || delta[i].z != 0)
            lattice_zero_ = false;
    }
    baseline_ = saved_baseline + db;
    set_rotation(dr);
    return true;
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
    for (int k = 1; k <= 2 && !lattice_zero_; k++) {
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
    Vec3 ray;
    if (lattice_zero_) {
        // The plain cube.  Its derivative is that of the bilinear patch between the tile's corners.
        ray = cube_ray_at(t.plane, t.col, t.fx, t.fy);
        if (pgx) {
            const std::vector<int>& t0 = tiles_[block[0]][0][block[1]];
            const Vec3 &c0 = base_[t0[0]], &c1 = base_[t0[2]], &c2 = base_[t0[6]], &c3 = base_[t0[8]];
            *pgx = (c1 - c0) * (1 - gy) + (c3 - c2) * gy;
            *pgy = (c2 - c0) * (1 - gx) + (c3 - c1) * gx;
        }
    } else {
        ray = spline(block[0], block[1], gx, gy, pgx, pgy);
    }
    if (pgx) {
        *pgx = *pgx * (dgx / t.w);
        *pgy = *pgy * (dgy / H);
    }
    if (lens == 1) {
        ray = rotate1(ray, false);
        if (pgx) {
            *pgx = rotate1(*pgx, false);
            *pgy = rotate1(*pgy, false);
        }
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
    if (side == 1)
        ray = rotate1(ray, true);
    if ((side == 0 && ray.z < -0.20) || (side == 1 && ray.z > 0.20))
        return Vec2(-1e6, -1e6);
    Vec2 px = cube_ray_to_pixel(ray, side);
    if (lattice_zero_)
        return px;      // the cube in closed form
    for (int pass = 0; pass <= 100; pass++) {
        Vec3 gx, gy;
        Vec3 def = pixel_to_ray(px, side, &gx, &gy);
        if (side == 1) {        // back to the unrotated frame, where `ray` is
            def = rotate1(def, true);
            gx = rotate1(gx, true);
            gy = rotate1(gy, true);
        }
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
