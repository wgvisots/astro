/**********************************************************************
 * conjunctions.cpp
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

// Conjunction study (docs/accuracy.md, docs/guide/propagation.md):
//  1. closest-approach search: minima missed against a brute-force scan, by sampling step;
//  2. collision probability against NASA CARA's reference values and Alfano's (2009) test cases,
//     with a brute-force grid where the references disagree.
//
//   astro_study_conjunctions [number of orbit pairs]   (default 300)

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <utility>
#include <vector>

#include "astro/core/angle.h"
#include "astro/orbit/close_approach.h"
#include "orbit/conjunction_reference.h"

using namespace astro;
using namespace astro::test;

namespace {

const Time kEpoch = Time::fromIso8601("2026-03-01T00:00:00Z");

/// Local minima of the distance found by scanning it every `dt` seconds over [0, T].
int bruteForceMinima(const Trajectory& a, const Trajectory& b, double T, double dt) {
    const auto fa = a.secondsAfter(kEpoch), fb = b.secondsAfter(kEpoch);
    const auto d2 = [&](double s) { return norm2(fb(s).r - fa(s).r); };
    int        n  = 0;
    double     f0 = d2(0.0), f1 = d2(dt);
    for (double s = 2.0 * dt; s <= T; s += dt) {
        const double f2 = d2(s);
        if (f1 < f0 && f1 <= f2) {
            ++n;
        }
        f0 = f1;
        f1 = f2;
    }
    return n;
}

void searchStepStudy(int nPairs) {
    std::printf("1. Closest-approach search: %d random two-body orbit pairs (perigee >= 6600 km, a up to\n"
                "   66 700 km, e up to 0.9), 1-day window, against a 0.5 s brute-force scan of the distance\n\n",
                nPairs);
    std::mt19937                           rng(7);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    const auto                             orbit = [&] {
        const double a = (6700.0 + u(rng) * 60000.0) * 1e3;
        const double e = u(rng) * std::min(0.9, 1.0 - 6600e3 / a);
        return KeplerianElements{a,
                                 e,
                                 u(rng) * constants::pi,
                                 u(rng) * constants::twoPi,
                                 u(rng) * constants::twoPi,
                                 u(rng) * constants::twoPi};
    };
    const double                                   T = 86400.0;
    std::vector<std::pair<Trajectory, Trajectory>> pairs;
    std::vector<int>                               brute;
    int                                            total = 0;
    for (int k = 0; k < nPairs; ++k) {
        pairs.emplace_back(Trajectory{kEpoch, orbit()}, Trajectory{kEpoch, orbit()});
        brute.push_back(bruteForceMinima(pairs.back().first, pairs.back().second, T, 0.5));
        total += brute.back();
    }

    std::printf("| stepFraction | Minima missed | Extra minima | Time per search |\n");
    std::printf("|---|---|---|---|\n");
    for (double fraction : {0.05, 0.1, 0.2, 0.4}) {
        CloseApproachOptions options;
        options.stepFraction = fraction;
        int        missed = 0, extra = 0;
        const auto t0 = std::chrono::steady_clock::now();
        for (std::size_t k = 0; k < pairs.size(); ++k) {
            const int n = static_cast<int>(
                closeApproaches(pairs[k].first, pairs[k].second, kEpoch, kEpoch + T, 1e300, options).size());
            missed += std::max(0, brute[k] - n);
            extra  += std::max(0, n - brute[k]);
        }
        const double us = std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count()
                        / static_cast<double>(pairs.size());
        std::printf("| %.2f%s | %d of %d | %d | %.0f us |\n",
                    fraction,
                    fraction == 0.05 ? " (default)" : "",
                    missed,
                    total,
                    extra,
                    us);
    }
    std::printf("\n");
}

void row(const char* what, double hbr, double pc, double reference) {
    std::printf("| %s | %g | %.15e | %.15e | %.1e |\n", what, hbr, pc, reference, pc / reference - 1.0);
}

void collisionProbabilityStudy() {
    std::printf("2. Collision probability against the references (relative difference = ours / reference - 1;\n"
                "   CARA's miss-distance convention, see caraPc)\n\n");
    std::printf("| Case | HBR | collisionProbability | Reference | Relative difference |\n");
    std::printf("|---|---|---|---|---|\n");
    row("CARA Omitron 1 (km)", kOmitronHbr, caraPc(kCara[0], kOmitronHbr), kOmitronPc[0]);
    row("CARA Omitron 2 (km)", kOmitronHbr, caraPc(kCara[1], kOmitronHbr), kOmitronPc[1]);
    for (int k = 0; k < 6; ++k) {
        row("CARA multiple radii (m)", kMultiRadiusHbr[k], caraPc(kCara[2], kMultiRadiusHbr[k]), kMultiRadiusPc[k]);
    }
    for (int k = 0; k < 3; ++k) {
        char what[64];
        std::snprintf(what, sizeof what, "Alfano (2009) case %d, linear Pc", 5 + k);
        row(what, kAlfanoHbr, caraPc(kCara[3 + k], kAlfanoHbr), kAlfanoLinearPc[k]);
    }
    for (int k = 0; k < 3; ++k) {
        char what[64];
        std::snprintf(what, sizeof what, "Alfano (2009) case %d, 1e8-sample Monte Carlo", 5 + k);
        row(what, kAlfanoHbr, caraPc(kCara[3 + k], kAlfanoHbr), kAlfanoMonteCarloPc[k]);
    }

    // CARA's multiple-radius geometry, where its 64-point Gauss-Chebyshev values disagree: a
    // brute-force polar grid at two resolutions (their difference shows the grid's own error).
    std::printf("\nCARA multiple radii, against a brute-force polar grid (Simpson in r x midpoint in theta)\n\n");
    std::printf(
        "| HBR (m) | collisionProbability | Grid 2000 x 40000 | Grid 4000 x 80000 | Ours vs grid | CARA vs grid |\n");
    std::printf("|---|---|---|---|---|---|\n");
    const Vec3    dr = kCara[2].r2 - kCara[2].r1, z = normalized(kCara[2].v2 - kCara[2].v1);
    CloseApproach c;
    c.first  = {Vec3{}, kCara[2].v1};
    c.second = {norm(dr) * normalized(dr - dot(dr, z) * z), kCara[2].v2};
    for (int k = 2; k < 6; ++k) {
        const double hbr    = kMultiRadiusHbr[k];
        const auto   q      = collisionProbability(c, kCara[2].c1, kCara[2].c2, hbr);
        const double coarse = polarReference(q.missDistance, 0.0, q.xx, q.xy, q.yy, hbr, 2000, 40000);
        const double fine   = polarReference(q.missDistance, 0.0, q.xx, q.xy, q.yy, hbr, 4000, 80000);
        std::printf("| %g | %.12e | %.12e | %.12e | %.1e | %.1e |\n",
                    hbr,
                    q.probability,
                    coarse,
                    fine,
                    q.probability / fine - 1.0,
                    kMultiRadiusPc[k] / fine - 1.0);
    }
    std::printf("\nEncounter-plane sigmas there: %.0f m x %.1f m, miss %.0f m\n",
                collisionProbability(c, kCara[2].c1, kCara[2].c2, 10.0).sigmaMajor,
                collisionProbability(c, kCara[2].c1, kCara[2].c2, 10.0).sigmaMinor,
                collisionProbability(c, kCara[2].c1, kCara[2].c2, 10.0).missDistance);
}

} // namespace

int main(int argc, char** argv) {
    const int nPairs = argc > 1 ? std::atoi(argv[1]) : 300;
    searchStepStudy(nPairs);
    collisionProbabilityStudy();
    return 0;
}
