/**********************************************************************
 * propagation_bench.cpp
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

// Propagation throughput: one state propagated from its epoch out X seconds at Y-second steps,
// repeated over many states. Compares the per-time calls with the grid APIs. Also times Lambert's
// problem and single calls of propagateKepler and SGP4.
//
//   astro_propagation_bench [X seconds] [Y seconds] [number of states]   (defaults 86400 60 100)

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "astro/core/angle.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/numerical.h"
#include "astro/orbit/secular.h"
#include "astro/orbit/sgp4.h"
#include "astro/orbit/tle.h"
#include "astro/orbit/two_body.h"

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

} // namespace

int main(int argc, char** argv) {
    const double X       = argc > 1 ? std::atof(argv[1]) : 86400.0;
    const double Y       = argc > 2 ? std::atof(argv[2]) : 60.0;
    const int    nStates = argc > 3 ? std::atoi(argv[3]) : 100;
    if (!(X > 0.0) || !(Y > 0.0) || nStates < 1) {
        std::fprintf(stderr, "usage: %s [X seconds] [Y seconds] [number of states]\n", argv[0]);
        return 2;
    }

    // LEO states spread over node, perigee and anomaly.
    std::vector<StateVector> states;
    for (int k = 0; k < nStates; ++k) {
        states.push_back(elementsToState({7000e3, 0.001, deg2rad(51.6), 0.1 * k, 0.3 * k, 0.05 * k}));
    }
    const std::size_t   count = static_cast<std::size_t>(X / Y) + 1;
    std::vector<double> dts(count);
    for (std::size_t k = 0; k < count; ++k) {
        dts[k] = Y * static_cast<double>(k);
    }
    std::vector<StateVector>   out(count);
    std::vector<StateVector12> out12(count);

    std::printf("X = %.0f s, Y = %g s: %zu outputs per state, %d states\n\n", X, Y, count, nStates);
    std::printf("%-46s %12s %12s\n", "method", "us/state", "ns/output");
    auto report = [&](const char* name, double seconds) {
        const double perState = seconds / nStates;
        std::printf("%-46s %12.1f %12.1f\n", name, perState * 1e6, perState * 1e9 / static_cast<double>(count));
    };

    report("secular: propagateSecular(state, dt) per time", secondsPerRun([&] {
               for (const auto& s : states) {
                   for (double dt : dts) {
                       sink += propagateSecular(s, dt).r.x;
                   }
               }
           }));
    report("secular: SecularPropagator::propagateGrid", secondsPerRun([&] {
               for (const auto& s : states) {
                   SecularPropagator(s).propagateGrid(0.0, Y, out);
                   sink += out.back().r.x;
               }
           }));
    report("secular: SecularPropagatorCorrected::propagateGrid", secondsPerRun([&] {
               for (const auto& s : states) {
                   SecularPropagatorCorrected(s).propagateGrid(0.0, Y, out);
                   sink += out.back().r.x;
               }
           }));
    report("secular: SecularPropagator::propagateGrid12", secondsPerRun([&] {
               for (const auto& s : states) {
                   SecularPropagator(s).propagateGrid12(0.0, Y, out12);
                   sink += out12.back().r.x;
               }
           }));
    report("numerical: j2j4 propagate(dts), stops at each", secondsPerRun([&] {
               for (const auto& s : states) {
                   sink += NumericalPropagator::j2j4(Time(), s).propagate(dts).back().r.x;
               }
           }));
    report("numerical: j2j4 propagateGrid", secondsPerRun([&] {
               for (const auto& s : states) {
                   NumericalPropagator::j2j4(Time(), s).propagateGrid(0.0, Y, out);
                   sink += out.back().r.x;
               }
           }));
    report("numerical: j2j4 propagateGrid12", secondsPerRun([&] {
               for (const auto& s : states) {
                   NumericalPropagator::j2j4(Time(), s).propagateGrid12(0.0, Y, out12);
                   sink += out12.back().r.x;
               }
           }));
    AdaptiveOptions loose; // the 4 h / few-metre setting from docs/performance.md
    loose.relTol = 1e-9;
    loose.absTol = 1e-6;
    report("numerical: propagateGrid, relTol 1e-9", secondsPerRun([&] {
               for (const auto& s : states) {
                   NumericalPropagator::j2j4(Time(), s, Frame::TOD, loose).propagateGrid(0.0, Y, out);
                   sink += out.back().r.x;
               }
           }));
    report("numerical: propagateGrid12, relTol 1e-9", secondsPerRun([&] {
               for (const auto& s : states) {
                   NumericalPropagator::j2j4(Time(), s, Frame::TOD, loose).propagateGrid12(0.0, Y, out12);
                   sink += out12.back().r.x;
               }
           }));
    // Lambert's problem between pairs of the states above, half an orbit to two days apart.
    std::printf("\n%-46s %12s\n", "Lambert (per problem)", "us");
    for (int revolutions : {0, 3}) {
        const double seconds = secondsPerRun([&] {
            for (std::size_t k = 0; k + 1 < states.size(); ++k) {
                const double tof = 3000.0 + 1700.0 * static_cast<double>(k);
                for (const LambertSolution& l :
                     lambertIzzo(states[k].r, 2.0 * states[k + 1].r, tof, constants::gmEarth, true, revolutions)) {
                    sink += l.v1.x;
                }
            }
        });
        std::printf("%-46s %12.2f\n",
                    revolutions == 0 ? "lambertIzzo, no complete revolution" : "lambertIzzo, up to 3 revolutions",
                    seconds / static_cast<double>(states.size() - 1) * 1e6);
    }
    // Single calls of the closed-form propagators: one state or element set, many times.
    std::printf("\n%-46s %12s\n", "Analytic propagators (per call)", "ns");
    {
        const double seconds = secondsPerRun([&] {
            for (const StateVector& s : states) {
                for (int k = 1; k <= 100; ++k) {
                    sink += propagateKepler(s, 60.0 * k).r.x;
                }
            }
        });
        std::printf("%-46s %12.0f\n", "propagateKepler", seconds / (100.0 * static_cast<double>(states.size())) * 1e9);
    }
    // Element sets of the SGP4 verification set: the ISS (near-Earth) and a 12 h Molniya orbit
    // (deep-space, resonant; its cost grows with the time from the epoch, here up to 1 day).
    const struct {
        const char* label;
        const char* line1;
        const char* line2;
    } sets[] = {
        {"Sgp4::propagate, near-Earth",
         "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927", "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537"},
        {"Sgp4::propagate, deep-space (within 1 day)",
         "1 08195U 75081A   06176.33215444  .00000099  00000-0  11873-3 0   813", "2 08195  64.1586 279.0717 6877146 264.7651  20.2257  2.00491383225656"},
    };
    for (const auto& set : sets) {
        const Sgp4   sgp4(Tle::parse(set.line1, set.line2));
        const double seconds = secondsPerRun([&] {
            for (int k = 0; k < 1440; ++k) {
                sink += sgp4.propagate(static_cast<double>(k)).r.x;
            }
        });
        std::printf("%-46s %12.0f\n", set.label, seconds / 1440.0 * 1e9);
    }
    {
        // The resonance terms are integrated in 12 h steps from the epoch on every call, so a
        // resonant orbit costs more the further the time is from its epoch.
        const Sgp4   sgp4(Tle::parse(sets[1].line1, sets[1].line2));
        const double seconds = secondsPerRun([&] {
            for (int k = 0; k < 100; ++k) {
                sink += sgp4.propagate(525600.0 + static_cast<double>(k)).r.x;
            }
        });
        std::printf("%-46s %12.0f\n", "Sgp4::propagate, deep-space, 1 year out", seconds / 100.0 * 1e9);
    }
    std::printf("\n(checksum %g)\n", sink);
    return 0;
}
