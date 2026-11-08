// The camera model and warp: self-consistency, and agreement with docs/FORMAT.md.  No footage needed;
// a .360 file as the second argument also checks reading the baseline from its header.
//   make -C cpp test
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "../cpp/maxcam.h"
using namespace maxcam;

static int fails = 0;
static void check(const char* name, bool ok, double value) {
    std::printf("%s %s  (%g)\n", ok ? "PASS" : "FAIL", name, value);
    if (!ok) fails++;
}
static double frand(double hi) { return (std::rand() % 100000) / 100000.0 * hi; }
static const double PI_ = 3.14159265358979323846;

// The ray table of docs/FORMAT.md, as written there.
static Vec3 format_ray(double x, double y) {
    double X = x + 0.5, top = y + 0.5 >= 1344 ? 1344 : 0;
    double w = std::tan((1 - 2 * (y + 0.5 - top) / 1344) * PI_ / 4);
    auto q = [](double n) { return std::tan((1 - n / 656) * PI_ / 4); };
    Vec3 d;
    if (X < 688) d = Vec3(-1, w, -q(X));
    else if (X < 1376) d = Vec3(-1, w, q(1376 - X));
    else if (X < 2720) d = Vec3(std::tan((2 * (X - 1376) / 1344 - 1) * PI_ / 4), w, 1);
    else if (X < 3408) d = Vec3(1, w, q(X - 2720));
    else d = Vec3(1, w, -q(4096 - X));
    if (top) d = Vec3(d.y, d.x, -d.z);
    return d.normalized();
}

int main(int argc, char** argv) {
    const char* corrections = argc > 1 ? argv[1] : "../calibration/model4.txt";
    CameraModel m;
    check("loads the corrections file", m.load("", corrections), 0);
    check("a missing corrections file is an error", !CameraModel().load("", "/nonexistent"), 0);
    check("a missing clip is an error", !CameraModel().load("/nonexistent.360", ""), 0);
    check("the shipped file has no lattice corrections", !m.has_lattice(), 0);
    Vec3 b = m.baseline();
    check("default baseline (-2.336, 0, 3.58) cm", std::fabs(b.x + 2.336) < 1e-9 && b.y == 0 && std::fabs(b.z - 3.58) < 1e-9, b.norm());
    {   // the default-constructed model is the same as the shipped file
        CameraModel d;
        double worst = 0;
        std::srand(7);
        for (int i = 0; i < 5000; i++) {
            Vec2 px(frand(4095), frand(2687));
            worst = std::max(worst, (d.pixel_to_ray(px) - m.pixel_to_ray(px)).norm());
        }
        check("the default model equals the shipped file", worst == 0, worst);
    }
    if (argc > 2) {
        CameraModel c;
        bool ok = c.load(argv[2], corrections);
        check("reads the baseline from a .360 header", ok, 0);
        if (!ok) std::printf("  %s\n", c.error().c_str());
        Vec3 cb = c.baseline();
        check("that baseline is a few cm, mostly along z", cb.norm() > 2 && cb.norm() < 8 && std::fabs(cb.z) > std::fabs(cb.x), cb.norm());
        std::printf("  header baseline (%.4f, %.4f, %.4f) cm\n", cb.x, cb.y, cb.z);
    }

    check("zone of a C0 pixel", zone_index(Vec2(2048, 672)) == 2, zone_index(Vec2(2048, 672)));
    check("C0 is lens 0", lens_index(Vec2(2048, 672)) == 0, 0);
    check("C1 is lens 1", lens_index(Vec2(2048, 2016)) == 1, 0);

    // The middle of each centre tile looks along its lens's axis.
    Vec3 c0 = m.pixel_to_ray(Vec2(2047.5, 671.5)).normalized(), c1 = m.pixel_to_ray(Vec2(2047.5, 2015.5)).normalized();
    check("C0 looks along +z", c0.z > 0.999999, c0.z);
    check("C1 looks along -z", c1.z < -0.999999, c1.z);

    // The supplied file is the cube itself: all ten tiles, against the document's table.
    std::srand(1);
    double worst = 0;
    for (int i = 0; i < 40000; i++) {
        Vec2 px(frand(4095), frand(2687));
        worst = std::max(worst, (m.pixel_to_ray(px).normalized() - format_ray(px.x, px.y)).norm());
    }
    check("with no corrections the model is FORMAT.md's ray table", worst < 1e-12, worst);
    Vec3 mid = m.pixel_to_ray(Vec2(2720 + 656 - 0.5, 671.5)).normalized();
    check("R0's mid-plane, mid-height is (1, 0, 0)", (mid - Vec3(1, 0, 0)).norm() < 1e-8, (mid - Vec3(1, 0, 0)).norm());

    // pixel -> ray -> pixel, including past the tile's outer edge
    worst = 0;
    for (int i = 0; i < 20000; i++) {
        Vec2 px(frand(4095), frand(2687));
        int lens = lens_index(px);
        Vec2 back = m.ray_to_pixel(m.pixel_to_ray(px).normalized(), lens);
        worst = std::max(worst, std::sqrt((back - px).norm2()));
    }
    check("pixel -> ray -> pixel within 0.001 px", worst < 1e-3, worst);
    Vec2 past(3407.5 + 40, 500);      // on R1's pixels, 40 px past R0's outer edge, still lens 0
    Vec3 past_ray = m.pixel_to_ray(past, 0).normalized();
    check("past the outer edge, the cube carries on", past_ray.z < 0 && past_ray.x > 0.9, past_ray.z);
    check("...and comes back to the same pixel", std::sqrt((m.ray_to_pixel(past_ray, 0) - past).norm2()) < 1e-4, std::sqrt((m.ray_to_pixel(past_ray, 0) - past).norm2()));

    // the derivative against finite differences, on the cube and with corrections
    double wg = 0;
    for (int i = 0; i < 2000; i++) {
        Vec2 px(20 + frand(4050), 20 + frand(2640));
        Vec3 gx, gy;
        Vec3 r = m.pixel_to_ray(px, &gx, &gy);
        (void)r;
        double e = 1e-3;
        if (zone_index(px + Vec2(e, 0)) != zone_index(px - Vec2(e, 0))) continue;
        Vec3 fd = (m.pixel_to_ray(px + Vec2(e, 0)) - m.pixel_to_ray(px - Vec2(e, 0))) * (1 / (2 * e));
        wg = std::max(wg, (fd - gx).norm());
    }
    check("d(ray)/dx matches finite differences", wg < 1e-4, wg);

    float disp[DISP_LEN];
    for (int i = 0; i < DISP_LEN; i++) disp[i] = 50 + 15 * std::sin(i * 0.05);
    float neutral[DISP_LEN], near_[DISP_LEN];
    for (int i = 0; i < DISP_LEN; i++) { neutral[i] = 66.77f; near_[i] = 15; }
    double rho = m.warp_ratio();
    check("the warp's ratio comes from the baseline", std::fabs(rho - 0.6527) < 0.001, rho);

    Vec2 centre(2048, 672);
    check("the centre tiles are not warped", std::sqrt((StitchWarp::recorded_to_model(centre, near_, rho) - centre).norm2()) < 1e-9, 0);
    Vec2 away(687.5 - 400, 672), inside(687.5 - 100, 672);
    check("beyond 344 px from the stitch, DISP has no effect", std::sqrt((StitchWarp::recorded_to_model(away, near_, rho) - away).norm2()) < 1e-9, 0);
    check("at the neutral value, nothing moves", std::sqrt((StitchWarp::recorded_to_model(inside, neutral, rho) - inside).norm2()) < 1e-4, 0);
    double moved = std::sqrt((StitchWarp::recorded_to_model(inside, near_, rho) - inside).norm2());
    check("within it, DISP moves the pixel", moved > 1, moved);
    // FORMAT.md: for DISP = 15, about 57 px across at the middle of a face, right at the tile's edge
    Vec2 edge(3407, 671.5), w = StitchWarp::recorded_to_model(edge, near_, rho);
    check("closest objects move about 57 px across, at mid-face", std::fabs(std::fabs(w.x - edge.x) - 57) < 2, std::fabs(w.x - edge.x));
    check("that pixel goes past the edge, on the same lens", w.x > 3407.5 && lens_index(edge) == 0, w.x);

    {   // Corrections.  A tiny lattice value switches to the spline path, which must still give the cube.
        FILE* f = std::fopen("test_corrections.txt", "w");
        for (int i = 0; i < 902; i++) std::fprintf(f, "%s\n", i == 500 ? "1e-13 0 0" : "0 0 0");
        std::fprintf(f, "0.1 0 -0.2\n0 0.5 0\n");      // baseline +(0.1, 0, -0.2) cm; lens 1 rotated 0.5 deg about y
        std::fclose(f);
        CameraModel c;
        check("loads a corrections file", c.load("", "test_corrections.txt"), 0);
        check("it has lattice corrections", c.has_lattice(), 0);
        Vec3 cb = c.baseline();
        check("baseline correction applied", std::fabs(cb.x + 2.236) < 1e-9 && std::fabs(cb.z - 3.38) < 1e-9, cb.norm());
        double worst = 0, worst1 = 0, worstb = 0;
        std::srand(3);
        for (int i = 0; i < 20000; i++) {
            Vec2 px(frand(4095), frand(2687));
            int lens = lens_index(px);
            Vec3 r = c.pixel_to_ray(px).normalized(), r0 = m.pixel_to_ray(px).normalized();
            if (lens == 0) worst = std::max(worst, (r - r0).norm());
            else {      // r0 rotated by 0.5 deg about +y, right-handed
                double th = 0.5 * PI_ / 180, cs = std::cos(th), sn = std::sin(th);
                Vec3 e(r0.x * cs + r0.z * sn, r0.y, -r0.x * sn + r0.z * cs);
                worst1 = std::max(worst1, (r - e).norm());
            }
            worstb = std::max(worstb, std::sqrt((c.ray_to_pixel(r, lens) - px).norm2()));
        }
        check("spline path with ~zero corrections is still the cube (lens 0)", worst < 1e-9, worst);
        check("lens 1 rays are rotated by 0.5 deg about +y", worst1 < 1e-12, worst1);
        check("pixel -> ray -> pixel with corrections within 0.001 px", worstb < 1e-3, worstb);
        Vec3 c1 = c.pixel_to_ray(Vec2(2047.5, 2015.5)).normalized();
        check("rotation about +y tilts C1's axis towards -x", c1.x < -0.008 && c1.x > -0.009, c1.x);
        std::remove("test_corrections.txt");
    }

    // recorded -> ray -> recorded.  A few recorded pixels right at a stitch corner do not come back.
    int n = 0, bad = 0;
    for (int i = 0; i < 40000; i++) {
        Vec2 r(frand(4095), frand(2687));
        if (r.x >= 688 * 2 && r.x < 3408 - 688) continue;
        Vec2 back = m.ray_to_recorded(m.recorded_to_ray(r, disp), lens_index(r), disp);
        n++;
        if (std::sqrt((back - r).norm2()) > 0.1) bad++;
    }
    check("recorded -> ray -> recorded within 0.1 px on 99 % of side-tile pixels", bad < n / 100, double(bad) / n);
    std::printf(fails ? "%d FAILED\n" : "all passed\n", fails);
    return fails ? 1 : 0;
}
