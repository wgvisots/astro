/**********************************************************************
 * conjunction_bench.cpp
 *
 * Copyright 2024 W. Greg Visotski
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 *         limitations under the License.
 **********************************************************************/

// Conjunction throughput (docs/performance.md): closest-approach searches over a window, and
// collision probabilities.
//
//   astro_conjunction_bench [window seconds]   (default 86400)

#include <chrono>
#include <cstdio>
#include <cstdlib>

#include "astro/core/angle.h"
#include "astro/orbit/close_approach.h"

using namespace astro;

namespace {

double sink = 0.0; // keeps results observable so the work is not optimized away

template <class F>
double secondsPerRun(F&& run, double minSeconds = 0.5) {
    run(); // warm up
    using Clock  = std::chrono::steady_clock;
    int  runs    = 0;
    auto t0      = Clock::now();
    auto elapsed = 0.0;
    do {
        run();
        ++runs;
        elapsed = std::chrono::duration<double>(Clock::now() - t0).count();
    } while (elapsed < minSeconds);
    return elapsed / runs;
}

void report(const char* name, double seconds) {
    std::printf("%-58s %10.2f us\n", name, seconds * 1e6);
}

} // namespace

int main(int argc, char** argv) {
    const double window = argc > 1 ? std::atof(argv[1]) : 86400.0;
    if (!(window > 0.0)) {
        std::fprintf(stderr, "usage: %s [window seconds]\n", argv[0]);
        return 2;
    }
    const Time epoch = Time::fromIso8601("2026-03-01T00:00:00Z");
    const Time end   = epoch + window;

    // Two LEO objects, and a LEO object against a GTO.
    const KeplerianElements leo1{7000e3, 0.001, 0.9, 0.3, 0.0, 0.0};
    const KeplerianElements leo2{7100e3, 0.002, 1.7, 2.0, 1.0, 2.0};
    const KeplerianElements gto{24400e3, 0.73, deg2rad(27.0), 1.0, 3.0, 0.5};

    std::printf("Window %.0f s\n\n", window);
    std::printf("%-58s %13s\n", "call", "time");
    const auto search = [&](const char* name, const Trajectory& a, const Trajectory& b) {
        report(name, secondsPerRun([&] { sink += closestApproach(a, b, epoch, end).distance; }));
    };
    search("closestApproach, LEO-LEO, two-body", {epoch, leo1}, {epoch, leo2});
    search("closestApproach, LEO-LEO, secular J2/J4",
           {epoch, leo1, TrajectoryModel::SecularJ2J4},
           {epoch, leo2, TrajectoryModel::SecularJ2J4});
    search("closestApproach, LEO-GTO, two-body", {epoch, leo1}, {epoch, gto});
    report("closeApproaches (all minima), LEO-LEO, two-body", secondsPerRun([&] {
               sink += static_cast<double>(closeApproaches({epoch, leo1}, {epoch, leo2}, epoch, end).size());
           }));

    // probabilityInCircle(mx, my, xx, xy, yy, R) on the geometries of docs/performance.md.
    struct Geometry {
        const char* name;
        double      mx, my, xx, xy, yy, R;
    };
    // clang-format off
    const Geometry geometries[] = {
        {"probabilityInCircle, typical (miss ~ sigma, R << sigma)", 100.0,  50.0, 400.0 * 400.0, 0.3 * 400.0 * 100.0, 100.0 * 100.0, 20.0},
        {"probabilityInCircle, aspect ratio 1000:1",                 30.0,  -5.0,           1e6,                 0.0,           1.0, 10.0},
        {"probabilityInCircle, R >> sigma",                           0.0,  40.0,          25.0,                -7.5,           9.0, 50.0},
        {"probabilityInCircle, 123 km x 19 m, R = 30 km",         19168.0,   0.0,      1.518e10,                 0.0,         347.6,  3e4},
    };
    // clang-format on
    std::printf("\n");
    for (const Geometry& g : geometries) {
        report(g.name, secondsPerRun([&] { sink += probabilityInCircle(g.mx, g.my, g.xx, g.xy, g.yy, g.R); }));
    }

    // The whole collisionProbability call, from 3x3 covariances: a 200 m miss, crossing at
    // ~10 km/s.
    CloseApproach ca;
    ca.first = elementsToState(leo1);
    const Vec3 v2{-ca.first.v.x, ca.first.v.z, ca.first.v.y};
    ca.second = {ca.first.r + 200.0 * normalized(cross(v2 - ca.first.v, ca.first.r)), v2};
    const Mat3 C{1e4, 2e3, 0.0, 2e3, 4e5, 0.0, 0.0, 0.0, 1e4};
    report("collisionProbability (projection + probabilityInCircle)",
           secondsPerRun([&] { sink += collisionProbability(ca, C, C, 10.0).probability; }));
    std::printf("  (Pc %.2e)\n", collisionProbability(ca, C, C, 10.0).probability);

    std::printf("\n(sink %g)\n", sink);
    return 0;
}
