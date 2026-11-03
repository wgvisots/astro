/**********************************************************************
 * covariance.cpp
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

// Covariance propagation along an output grid, and the uncertainty in radial / in-track /
// cross-track coordinates.

#include <cmath>
#include <iomanip>
#include <iostream>

#include "astro/core/angle.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/numerical.h"
#include "astro/satellite/orbital_frames.h"

using namespace astro;

int main() {
    const Time        epoch = Time::fromIso8601("2026-03-01T00:00:00Z");
    const StateVector s0    = elementsToState({6778e3, 0.0005, deg2rad(51.6), 0.3, 0.5, 0.7});

    //! [covariance]
    // 10 m and 1 cm/s (1-sigma) per axis at the epoch, in the state's frame (GCRS).
    Mat6 P0{};
    for (std::size_t i = 0; i < 3; ++i) {
        P0[i][i]         = 10.0 * 10.0;
        P0[i + 3][i + 3] = 0.01 * 0.01;
    }
    AdaptiveOptions opt;
    opt.relTol      = 1e-9;
    opt.absTol      = 1e-6;
    const auto prop = NumericalPropagator::j2j4(epoch, s0, Frame::GCRS, opt);
    // States and covariances every 10 min for 4 h, with a little process noise (m^2/s^3).
    const std::vector<StateCovariance> out = prop.propagateGridWithCovariance(0.0, 600.0, 25, P0, 1e-10);
    //! [covariance]

    //! [ric]
    std::cout << " t (h)   sigma R (m)   sigma I (m)   sigma C (m)\n" << std::fixed << std::setprecision(1);
    for (std::size_t k = 0; k < out.size(); k += 6) {
        const Mat6 ric = rotateCovariance(out[k].covariance, ricMatrix(out[k].state));
        std::cout << std::setw(5) << k / 6.0 << std::setw(14) << std::sqrt(ric[0][0]) << std::setw(14)
                  << std::sqrt(ric[1][1]) << std::setw(14) << std::sqrt(ric[2][2]) << '\n';
    }
    //! [ric]
    return 0;
}
