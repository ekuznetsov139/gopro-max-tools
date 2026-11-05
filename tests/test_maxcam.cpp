// The extracted camera model and warp: self-consistency, on the example calibration.  No footage needed.
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

int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: %s model4.txt\n", argv[0]); return 2; }
    CameraModel m;
    check("loads", m.load(argv[1]), 0);
    check("a missing file is an error", !CameraModel().load("/nonexistent"), 0);
    Vec3 b = m.baseline();
    check("baseline about (-2.3, 0, 3.0) cm", std::fabs(b.x + 2.33) < 0.1 && std::fabs(b.y) < 0.1 && std::fabs(b.z - 3.04) < 0.1, b.norm());

    check("zone of a C0 pixel", zone_index(Vec2(2048, 672)) == 2, zone_index(Vec2(2048, 672)));
    check("C0 is lens 0", lens_index(Vec2(2048, 672)) == 0, 0);
    check("C1 is lens 1", lens_index(Vec2(2048, 2016)) == 1, 0);

    // The middle of each centre tile looks along its lens's axis.
    Vec3 c0 = m.pixel_to_ray(Vec2(2047.5, 671.5)).normalized(), c1 = m.pixel_to_ray(Vec2(2047.5, 2015.5)).normalized();
    check("C0 looks along +z", c0.z > 0.99, c0.z);
    check("C1 looks along -z", c1.z < -0.99, c1.z);

    // pixel -> ray -> pixel
    std::srand(1);
    double worst = 0;
    for (int i = 0; i < 20000; i++) {
        Vec2 px(frand(4095), frand(2687));
        Vec2 back = m.ray_to_pixel(m.pixel_to_ray(px).normalized(), lens_index(px));
        worst = std::max(worst, std::sqrt((back - px).norm2()));
    }
    check("pixel -> ray -> pixel within 0.001 px", worst < 1e-3, worst);

    // the derivative against finite differences
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

    check("the centre tiles are not warped", std::fabs(StitchWarp::recorded_to_model(Vec2(2048, 672), disp).x - 2048) < 1e-9, 0);
    float flat[DISP_LEN];
    for (int i = 0; i < DISP_LEN; i++) flat[i] = 80;
    Vec2 s = StitchWarp::recorded_to_model(Vec2(687.5 - 300, 672), flat), s2 = StitchWarp::recorded_to_model(Vec2(687.5 - 300, 672), disp);
    check("beyond 344 px from the stitch, DISP has no effect", (s - StitchWarp::recorded_to_model(Vec2(687.5 - 400, 672), flat)).norm2() > 1 && std::fabs(StitchWarp::recorded_to_model(Vec2(687.5 - 400, 672), flat).x - StitchWarp::recorded_to_model(Vec2(687.5 - 400, 672), disp).x) < 1e-6, 0);
    check("within it, DISP moves the pixel", (s - s2).norm2() > 1, std::sqrt((s - s2).norm2()));

    // recorded -> model -> recorded.  A few recorded pixels right at a stitch corner do not come back.
    int n = 0, bad = 0;
    for (int i = 0; i < 60000; i++) {
        Vec2 r(frand(4095), frand(2687));
        if (r.x >= 688 * 2 && r.x < 3408 - 688) continue;
        Vec2 back = StitchWarp::model_to_recorded(StitchWarp::recorded_to_model(r, disp), disp);
        n++;
        if (std::sqrt((back - r).norm2()) > 0.1) bad++;
    }
    check("recorded -> model -> recorded within 0.1 px on 99 % of side-tile pixels", bad < n / 100, double(bad) / n);
    std::printf(fails ? "%d FAILED\n" : "all passed\n", fails);
    return fails ? 1 : 0;
}
