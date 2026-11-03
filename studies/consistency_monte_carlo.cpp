/**********************************************************************
 * consistency_monte_carlo.cpp
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

// Consistency-check study (docs/accuracy.md): with priors sampled from their covariance, how
// often does each ConsistencyTest raise a false alarm, and how often does it catch a burn?
//
//   astro_study_consistency [samples]     (default 1500; the published tables used 1500 and 400)

#include <cstdio>
#include <cstdlib>
#include <random>

#include "astro/core/angle.h"
#include "astro/orbit/elements.h"
#include "astro/satellite/estimate_tracker.h"
#include "astro/satellite/orbital_frames.h"

using namespace astro;

namespace {

const Time kEpoch = Time::fromIso8601("2026-03-01T00:00:00Z");

struct Orbit {
    const char*       name;
    KeplerianElements el;
};

const Orbit kOrbits[] = {
    {"LEO 400 km i=51.6",     {6778e3, 0.0005, deg2rad(51.6), 0.3, 0.5, 0.7}},
    {"SSO 800 km e=0.01",     {7178e3, 0.01, deg2rad(98.6), 1.3, 2.5, 0.2}  },
    {"HEO 1000 x 20000 km",   {16878e3, 0.56, deg2rad(63.4), 0.1, 4.7, 0.0} },
    {"near-equatorial i=0.5", {7000e3, 0.002, deg2rad(0.5), 0.0, 0.0, 1.0}  },
    {"polar i=90",            {7000e3, 0.001, deg2rad(90.0), 2.0, 0.3, 0.1} },
    {"retrograde i=140",      {7000e3, 0.001, deg2rad(140.0), 0.5, 1.0, 2.0}},
};

/// Reference ("truth") propagator: J2 + J4 with a tight tolerance.
TrackerOptions tightOptions() {
    TrackerOptions o;
    o.integrator.relTol  = 1e-14;
    o.integrator.absTol  = 1e-11;
    o.integrator.maxStep = 20.0;
    return o;
}

[[maybe_unused]] Mat6 diagonal(double positionSigma, double velocitySigma) {
    Mat6 P{};
    for (std::size_t i = 0; i < 3; ++i) {
        P[i][i]         = positionSigma * positionSigma;
        P[i + 3][i + 3] = velocitySigma * velocitySigma;
    }
    return P;
}

} // namespace

namespace {

StateVector sample(const StateVector& s0, double sr, double sv, std::mt19937_64& rng) {
    std::normal_distribution<double> n(0.0, 1.0);
    StateVector                      s  = s0;
    s.r                                += Vec3{n(rng), n(rng), n(rng)} * sr;
    s.v                                += Vec3{n(rng), n(rng), n(rng)} * sv;
    return s;
}

void falseAlarms(int samples) {
    std::printf("False-alarm rate at a nominal 1e-3 (mean d^2 in parentheses; ideal 6 / 3 / 6)\n\n");
    std::printf("%-24s %-16s %5s   %-19s %-19s %-19s\n",
                "orbit",
                "prior sigma",
                "gap",
                "PositionVelocity",
                "Position",
                "Elements");
    const Mat6 tiny = diagonal(1e-3, 1e-6);
    for (const auto& o : kOrbits) {
        const StateVector s0 = elementsToState(o.el);
        for (auto [sr, sv] : {
                 std::pair{20.0,  0.02},
                 std::pair{100.0, 0.1 }
        }) {
            for (double dt : {1800.0, 7200.0, 14400.0}) {
                const StateEstimate truth{kEpoch + dt,
                                          makePropagator({kEpoch, s0, Mat6{}}, tightOptions()).propagate(dt),
                                          tiny};
                std::mt19937_64     rng(11);
                double              m6 = 0, m3 = 0, me = 0;
                int                 a6 = 0, a3 = 0, ae = 0;
                for (int k = 0; k < samples; ++k) {
                    const auto c =
                        checkConsistency({kEpoch, sample(s0, sr, sv, rng), diagonal(sr, sv)}, truth, tightOptions());
                    m6 += c.mahalanobis2;
                    m3 += c.positionMahalanobis2;
                    me += c.elementsMahalanobis2;
                    a6 += c.pValue < 1e-3;
                    a3 += c.positionPValue < 1e-3;
                    ae += c.elementsPValue < 1e-3;
                }
                char prior[32];
                std::snprintf(prior, sizeof prior, "%.0f m, %.2f m/s", sr, sv);
                std::printf("%-24s %-16s %4.1fh   %6.4f (%8.2f)   %6.4f (%6.2f)   %6.4f (%6.2f)\n",
                            o.name,
                            prior,
                            dt / 3600,
                            double(a6) / samples,
                            m6 / samples,
                            double(a3) / samples,
                            m3 / samples,
                            double(ae) / samples,
                            me / samples);
            }
        }
    }
}

void detection(int samples) {
    std::printf("\nDetection rate at alarm probability 1e-3: LEO 400 km, prior 20 m / 2 cm/s (sampled)\n\n");
    std::printf("%-6s %-10s %-26s %10s %10s %10s\n",
                "gap",
                "burn",
                "when / direction",
                "PosVel",
                "Position",
                "Elements");
    const StateVector s0   = elementsToState(kOrbits[0].el);
    const Mat6        tiny = diagonal(1e-3, 1e-6);
    for (double gap : {1800.0, 14400.0}) {
        for (double dv : {0.02, 0.1, 0.5}) {
            for (double before : {20.0, 900.0}) {
                for (int dir = 0; dir < 3; ++dir) {
                    StateVector b  = makePropagator({kEpoch, s0, Mat6{}}, tightOptions()).propagate(gap - before);
                    const Mat3  R  = ricMatrix(b);
                    b.v           += Vec3{R(dir, 0), R(dir, 1), R(dir, 2)} * dv;
                    const StateEstimate next{
                        kEpoch + gap,
                        makePropagator({kEpoch + gap - before, b, Mat6{}}, tightOptions()).propagate(before),
                        tiny};
                    std::mt19937_64 rng(5);
                    int             d6 = 0, d3 = 0, de = 0;
                    for (int k = 0; k < samples; ++k) {
                        const auto c  = checkConsistency({kEpoch, sample(s0, 20.0, 0.02, rng), diagonal(20.0, 0.02)},
                                                        next,
                                                        tightOptions());
                        d6           += c.pValue < 1e-3;
                        d3           += c.positionPValue < 1e-3;
                        de           += c.elementsPValue < 1e-3;
                    }
                    const char* dn[] = {"radial", "in-track", "cross-track"};
                    char        when[40], burn[16];
                    std::snprintf(when, sizeof when, "%4.0f s before, %s", before, dn[dir]);
                    std::snprintf(burn, sizeof burn, "%.2f m/s", dv);
                    std::printf("%4.1fh  %-10s %-26s %9.0f%% %9.0f%% %9.0f%%\n",
                                gap / 3600,
                                burn,
                                when,
                                100.0 * d6 / samples,
                                100.0 * d3 / samples,
                                100.0 * de / samples);
                }
            }
        }
    }
}

} // namespace

int main(int argc, char** argv) {
    const int samples = argc > 1 ? std::atoi(argv[1]) : 1500;
    falseAlarms(samples);
    detection(std::max(100, samples / 4));
    return 0;
}
