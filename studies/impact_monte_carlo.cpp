/**********************************************************************
 * impact_monte_carlo.cpp
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

// Impact study (docs/accuracy.md):
//   1. impact points against a brute-force reference (tight tolerance, the full IAU 2006/2000A
//      transformation at every evaluation, Brent's method on the propagated height);
//   2. the linear 1-sigma impact ellipse against a Monte Carlo of the initial covariance.
//
//   astro_study_impact [samples]     (default 2000; the published table used 4000)

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <random>
#include <vector>

#include "astro/core/angle.h"
#include "astro/core/rootfind.h"
#include "astro/geodesy/local.h"
#include "astro/orbit/impact.h"

using namespace astro;

namespace {

const Time kEpoch = Time::fromIso8601("2026-03-01T12:00:00Z");

struct Shot {
    const char* name;
    double      speed; // m/s, relative to the ground
    double      gamma; // deg above the horizon
    double      azimuth;
};

StateVector launch(const Shot& s, const EopProvider& eop) {
    const Geodetic site{deg2rad(30.0), deg2rad(-100.0), 0.0};
    const Mat3     enu = ecefToEnuMatrix(site.lat, site.lon);
    const double   g = deg2rad(s.gamma), a = deg2rad(s.azimuth);
    const Vec3     v{s.speed * std::cos(g) * std::sin(a), s.speed * std::cos(g) * std::cos(a), s.speed * std::sin(g)};
    return itrsToGcrs(StateVector{geodeticToEcef(site), enu.transpose() * v}, kEpoch, eop);
}

double referenceError(const StateVector& s0, const Impact& imp, const EopProvider& eop) {
    const ImpactOptions o;
    AdaptiveOptions     tight;
    tight.relTol = 1e-13;
    tight.absTol = 1e-9;
    const NumericalPropagator prop(
        kEpoch,
        s0,
        {std::make_shared<GravityField>(o.J, o.mu, o.re, poleOfDateMatrix(Frame::GCRS, kEpoch))},
        tight);
    const auto height = [&](double dt) { return ecefToGeodetic(gcrsToItrs(prop.propagate(dt).r, kEpoch + dt, eop)).h; };
    const double dt   = brentRoot(height, imp.flightTime - 2.0, imp.flightTime + 2.0, 1e-9);
    return norm(gcrsToItrs(prop.propagate(dt).r, kEpoch + dt, eop) - imp.itrs.r);
}

} // namespace

int main(int argc, char** argv) {
    const int samples = argc > 1 ? std::atoi(argv[1]) : 2000;
    EopValues ev;
    ev.dut1 = 0.3;
    ev.xp   = 0.10 * constants::arcsecToRad;
    ev.yp   = 0.35 * constants::arcsecToRad;
    const ConstantEop eop(ev);

    const Shot shots[] = {
        {"3 km/s, 35 deg", 3000.0, 35.0, 45.0},
        {"5 km/s, 35 deg", 5000.0, 35.0, 45.0},
        {"7 km/s, 35 deg", 7000.0, 35.0, 45.0},
        {"7 km/s, 15 deg", 7000.0, 15.0, 45.0},
    };

    std::printf("Impact point against the brute-force reference\n");
    std::printf("%-16s %10s %12s %10s\n", "shot", "flight (s)", "range (km)", "error (m)");
    for (const Shot& s : shots) {
        const StateVector s0    = launch(s, eop);
        const auto        imp   = predictImpact(kEpoch, s0, Frame::GCRS, eop);
        const double      range = constants::wgs84A
                           * std::acos(dot(normalized(imp->itrs.r),
                                           normalized(geodeticToEcef({deg2rad(30.0), deg2rad(-100.0), 0.0}))));
        std::printf("%-16s %10.1f %12.1f %10.4f\n",
                    s.name,
                    imp->flightTime,
                    range / 1000.0,
                    referenceError(s0, *imp, eop));
    }

    std::printf("\nLinear 1-sigma impact ellipse against a Monte Carlo of %d samples\n", samples);
    std::printf("%-16s %-14s %21s %21s %15s %11s\n",
                "shot",
                "1-sigma r, v",
                "semi-major lin / MC",
                "semi-minor lin / MC",
                "rotation lin/MC",
                "in 95% ell.");
    const double sigmas[][2] = {
        {10.0,   0.01},
        {100.0,  0.1 },
        {1000.0, 1.0 }
    };
    std::mt19937_64                  rng(11);
    std::normal_distribution<double> n01;
    for (const Shot& s : shots) {
        const StateVector s0 = launch(s, eop);
        for (const auto& sg : sigmas) {
            Mat6 P{};
            for (std::size_t i = 0; i < 3; ++i) {
                P[i][i]         = sg[0] * sg[0];
                P[i + 3][i + 3] = sg[1] * sg[1];
            }
            const auto               nominal = predictImpact(kEpoch, s0, P, Frame::GCRS, eop);
            const ImpactUncertainty& u       = *nominal->uncertainty;
            const Mat3               enu     = ecefToEnuMatrix(nominal->point.lat, nominal->point.lon);
            // Inverse of the linear covariance, for the fraction inside the linear 95 % ellipse.
            const double      det = u.eastVariance * u.northVariance - u.eastNorthCovariance * u.eastNorthCovariance;
            const double      k95 = ellipseScale(0.95);
            double            me = 0.0, mn = 0.0, cee = 0.0, cnn = 0.0, cen = 0.0;
            int               inside = 0;
            std::vector<Vec3> d(static_cast<std::size_t>(samples));
            for (int m = 0; m < samples; ++m) {
                StateVector x                   = s0;
                x.r                            += Vec3{n01(rng), n01(rng), n01(rng)} * sg[0];
                x.v                            += Vec3{n01(rng), n01(rng), n01(rng)} * sg[1];
                const auto imp                  = predictImpact(kEpoch, x, Frame::GCRS, eop);
                d[static_cast<std::size_t>(m)]  = enu * (imp->itrs.r - nominal->itrs.r);
                const Vec3& e                   = d[static_cast<std::size_t>(m)];
                me                             += e.x;
                mn                             += e.y;
                const double q =
                    (u.northVariance * e.x * e.x - 2.0 * u.eastNorthCovariance * e.x * e.y + u.eastVariance * e.y * e.y)
                    / det;
                inside += q <= k95 * k95 ? 1 : 0;
            }
            me /= samples;
            mn /= samples;
            for (const Vec3& e : d) {
                cee += (e.x - me) * (e.x - me);
                cnn += (e.y - mn) * (e.y - mn);
                cen += (e.x - me) * (e.y - mn);
            }
            const ErrorEllipse mc = errorEllipse(cee / (samples - 1), cnn / (samples - 1), cen / (samples - 1));
            char               label[64];
            std::snprintf(label, sizeof label, "%g m, %g m/s", sg[0], sg[1]);
            std::printf("%-16s %-14s %10.1f / %8.1f %10.1f / %8.1f %7.1f / %5.1f %10.1f%%\n",
                        s.name,
                        label,
                        u.ellipse.semiMajor,
                        mc.semiMajor,
                        u.ellipse.semiMinor,
                        mc.semiMinor,
                        rad2deg(u.ellipse.rotation),
                        rad2deg(mc.rotation),
                        100.0 * inside / samples);
        }
    }
}
