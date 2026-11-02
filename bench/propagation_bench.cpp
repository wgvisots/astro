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
// repeated over many states. Compares the per-time calls with the grid APIs.
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
    std::printf("\n(checksum %g)\n", sink);
    return 0;
}
