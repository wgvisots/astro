/**********************************************************************
 * tracker.cpp
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

// Following an upstream filter: latest-wins tracking with a consistency check on each delivery.
// Deliveries are simulated from a "truth" orbit that makes a 0.3 m/s burn between deliveries.

#include <iomanip>
#include <iostream>

#include "astro/core/angle.h"
#include "astro/orbit/elements.h"
#include "astro/satellite/estimate_tracker.h"

using namespace astro;

namespace {

void print(const ConsistencyCheck& c) {
    std::cout << std::fixed << std::setprecision(1) << "  residual R/I/C " << c.positionResidualRic << " m  vs 1-sigma "
              << c.positionSigmaRic << " m\n  semi-major axis residual " << c.semiMajorAxisResidual << " m (sigma "
              << c.semiMajorAxisSigma << " m), element p-value " << std::setprecision(3) << c.elementsPValue << " -> "
              << (c.consistent ? "consistent" : "FLAGGED") << "\n";
}

} // namespace

int main() {
    const Time t0 = Time::fromIso8601("2026-03-01T00:00:00Z");
    Mat6       P{};
    for (std::size_t i = 0; i < 3; ++i) {
        P[i][i]         = 25.0; // 5 m
        P[i + 3][i + 3] = 1e-5; // 3 mm/s
    }
    // Truth for the simulation: the tracker's own model.
    const StateVector s0    = elementsToState({6778e3, 0.0005, degToRad(51.6), 0.3, 0.5, 0.7});
    auto              truth = [](const Time& t, const StateVector& s, double dt) {
        return trackerPropagator({t, s, Mat6{}}).propagate(dt);
    };

    //! [tracker]
    TrackerOptions options;
    options.frame = Frame::GCRS;      // frame of the deliveries (also TOD, TEME, J2000, ...)
    EstimateTracker tracker(options); // J2 + J4, relTol 1e-9, element consistency test

    (void)tracker.update({t0, s0, P}); // the first delivery has nothing to be checked against
    if (auto check = tracker.update({t0 + 1800.0, truth(t0, s0, 1800.0), P})) {
        std::cout << "delivery 2 (no burn):\n";
        print(*check);
    }
    //! [tracker]

    // A 0.3 m/s in-track burn 5 minutes before the third delivery.
    StateVector burn  = truth(t0, s0, 3300.0);
    burn.v           += normalized(burn.v) * 0.3;
    //! [flagged]
    if (auto check = tracker.update({t0 + 3600.0, truth(t0 + 3300.0, burn, 300.0), P}); check && !check->consistent) {
        std::cout << "delivery 3 (after a burn):\n";
        print(*check);
    }
    // Latest wins: predictions now start from delivery 3.
    const auto ahead = tracker.predictGrid(0.0, 60.0, 61); // states + covariances for the next hour
    std::cout << "predicted 1 h ahead: " << ahead.back().state.r << " m\n";
    //! [flagged]
    return 0;
}
