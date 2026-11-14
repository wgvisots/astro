/**********************************************************************
 * covariance_bench.cpp
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

// Covariance, consistency-check and frame timings (docs/performance.md): LEO, one core.

#include <chrono>
#include <cstdio>
#include <vector>

#include "astro/core/angle.h"
#include "astro/frames/eop.h"
#include "astro/frames/transform.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/numerical.h"
#include "astro/orbit/secular.h"
#include "astro/orbit/two_body.h"
#include "astro/satellite/estimate_tracker.h"

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
    const Time        t0 = Time::fromIso8601("2026-03-01T00:00:00Z");
    const StateVector s0 = elementsToState({6778e3, 0.0005, deg2rad(51.6), 0.3, 0.5, 0.7});
    Mat6              P0{};
    for (std::size_t i = 0; i < 3; ++i) {
        P0[i][i]         = 100.0;
        P0[i + 3][i + 3] = 1e-4;
    }

    std::printf("4 h at Y = 10 s (1441 outputs)\n\n%-8s %16s %16s %24s\n",
                "relTol",
                "states only",
                "+ covariance",
                "+ covariance + noise");
    std::vector<StateVector>     st(1441);
    std::vector<StateCovariance> pc(1441);
    for (double rt : {1e-12, 1e-9}) {
        AdaptiveOptions o;
        o.relTol          = rt;
        o.absTol          = rt * 1e3;
        const auto   prop = NumericalPropagator::j2j4(t0, s0, Frame::GCRS, o);
        const double a    = microseconds([&] {
            prop.propagateGrid(0.0, 10.0, st);
            sink += st.back().r.x;
        });
        const double b    = microseconds([&] {
            prop.propagateGridWithCovariance(0.0, 10.0, P0, pc);
            sink += pc.back().state.r.x;
        });
        const double c    = microseconds([&] {
            prop.propagateGridWithCovariance(0.0, 10.0, P0, pc, 1e-9);
            sink += pc.back().state.r.x;
        });
        std::printf("%-8.0e %13.0f us %13.0f us %21.0f us\n", rt, a, b, c);
    }

    std::printf("\nAnalytic propagators with covariance (same orbit and P0)\n\n");
    int        n       = 0;
    const auto secular = SecularPropagatorCorrected(s0);
    std::printf("  propagateKeplerWithCovariance                     %7.2f us\n",
                microseconds([&] { sink += propagateKeplerWithCovariance(s0, 100.0 + n++, P0).covariance[0][0]; }));
    std::printf("  SecularPropagatorCorrected::propagateWithCovariance %5.2f us\n",
                microseconds([&] { sink += secular.propagateWithCovariance(100.0 + n++, P0).covariance[0][0]; }));
    std::printf("  SecularPropagatorCorrected::propagateGridWithCovariance, 1441 outputs %7.0f us\n", microseconds([&] {
                    secular.propagateGridWithCovariance(0.0, 10.0, P0, pc);
                    sink += pc.back().covariance[0][0];
                }));

    std::printf("\nConsistency check (propagate the previous estimate with covariance, compare)\n\n");
    for (double gap : {600.0, 1800.0, 14400.0}) {
        const StateEstimate prior{t0, s0, P0};
        const StateEstimate next{t0 + gap, makePropagator(prior).propagate(gap), P0};
        std::printf("  gap %6.0f s: %7.1f us\n", gap, microseconds([&] {
                        sink += checkConsistency(prior, next).mahalanobis2;
                    }));
    }

    std::printf("\nFrame matrices\n\n");
    FrameOptions b;
    b.model = PrecessionNutationModel::IAU2006_2000B;
    int k   = 0;
    std::printf("  FrameTransformer, IAU 2006/2000A   %7.1f us\n",
                microseconds([&] { sink += FrameTransformer(t0 + 7.0 * k++, ConstantEop{}).gast(); }));
    std::printf("  FrameTransformer, IAU 2006/2000B   %7.1f us\n",
                microseconds([&] { sink += FrameTransformer(t0 + 7.0 * k++, ConstantEop{}, b).gast(); }));
    std::printf("  poleOfDateMatrix(GCRS)             %7.1f us\n",
                microseconds([&] { sink += poleOfDateMatrix(Frame::GCRS, t0 + 7.0 * k++)(2, 2); }));
    std::printf("\n(checksum %g)\n", sink);
    return 0;
}
