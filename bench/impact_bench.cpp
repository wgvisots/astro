/**********************************************************************
 * impact_bench.cpp
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

// Impact-prediction timings (docs/performance.md): one core.

#include <chrono>
#include <cstdio>

#include "astro/core/angle.h"
#include "astro/geodesy/local.h"
#include "astro/orbit/impact.h"

using namespace astro;

namespace {

double sink = 0.0;

template <class F>
double microseconds(F&& run, double minSeconds = 0.5) {
    run();
    using Clock = std::chrono::steady_clock;
    int    n    = 0;
    auto   t0   = Clock::now();
    double el   = 0.0;
    do {
        run();
        ++n;
        el = std::chrono::duration<double>(Clock::now() - t0).count();
    } while (el < minSeconds);
    return el / n * 1e6;
}

} // namespace

int main() {
    const Time        epoch = Time::fromIso8601("2026-03-01T12:00:00Z");
    const ConstantEop eop;
    const Geodetic    site{deg2rad(30.0), deg2rad(-100.0), 0.0};
    const Mat3        enu = ecefToEnuMatrix(site.lat, site.lon);
    Mat6              P{};
    for (std::size_t i = 0; i < 3; ++i) {
        P[i][i]         = 1e4;
        P[i + 3][i + 3] = 1e-2;
    }

    for (const auto model : {PrecessionNutationModel::IAU2006_2000B, PrecessionNutationModel::IAU2006_2000A}) {
        const double us = microseconds([&] {
            const ItrsRotationSeries s(epoch, 7200.0, eop, {model, true});
            sink += s.duration();
        });
        std::printf("ItrsRotationSeries, 2 h window (3 nodes), %s: %8.1f us\n",
                    model == PrecessionNutationModel::IAU2006_2000B ? "IAU 2006/2000B" : "IAU 2006/2000A",
                    us);
    }

    const ItrsRotationSeries series(epoch, 7200.0, eop, ImpactOptions{}.frameOptions);
    ImpactOptions            shared;
    shared.rotations = &series;
    std::printf("\n%-24s %10s %14s %14s %16s\n",
                "shot",
                "flight (s)",
                "per call (us)",
                "shared (us)",
                "+ covariance (us)");
    for (const double V : {3000.0, 5000.0, 7000.0}) {
        const double      g = deg2rad(35.0), a = deg2rad(45.0);
        const Vec3        v{V * std::cos(g) * std::sin(a), V * std::cos(g) * std::cos(a), V * std::sin(g)};
        const StateVector s0  = itrsToGcrs(StateVector{geodeticToEcef(site), enu.transpose() * v}, epoch, eop);
        const double      own = microseconds([&] { sink += predictImpact(epoch, s0, Frame::GCRS, eop)->flightTime; });
        const double sh = microseconds([&] { sink += predictImpact(epoch, s0, Frame::GCRS, eop, shared)->flightTime; });
        const double cov = microseconds(
            [&] { sink += predictImpact(epoch, s0, P, Frame::GCRS, eop, shared)->uncertainty->ellipse.semiMajor; });
        char label[48];
        std::snprintf(label, sizeof label, "%.0f km/s, 35 deg", V / 1000.0);
        std::printf("%-24s %10.1f %14.1f %14.1f %16.1f\n",
                    label,
                    predictImpact(epoch, s0, Frame::GCRS, eop)->flightTime,
                    own,
                    sh,
                    cov);
    }
    return sink == 0.123 ? 1 : 0;
}
