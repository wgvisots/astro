/**********************************************************************
 * secular.h
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

#pragma once

#include "astro/core/constants.h"
#include "astro/orbit/elements.h"

/// Analytic secular (mean-element) drift due to the Earth's zonal harmonics,
/// after Brouwer (1959): first order in J2, plus J2^2 and J4.
namespace astro {

struct ZonalParameters {
    double mu{constants::gmEarth};
    double re{constants::earthRadiusEgm2008};
    double J2{constants::J2};
    double J4{constants::J4};
};

enum class SecularOrder {
    J2,         ///< first-order J2 only
    J2SquaredJ4 ///< J2 + J2^2 + J4 (Brouwer)
};

/// Secular rates (rad/s) of the mean elements.
struct SecularRates {
    double raanDot{0.0};
    double argpDot{0.0};
    double meanAnomalyDot{0.0};  ///< total rate including the mean motion
};

/// Rates for mean semi-major axis a, eccentricity e and inclination i.
[[nodiscard]] SecularRates secularRates(double a, double e, double i, const ZonalParameters& zp = {},
                                        SecularOrder order = SecularOrder::J2SquaredJ4) noexcept;

/// Propagate mean elements by dt seconds with the secular rates (a, e, i constant).
/// The true anomaly is advanced through the mean anomaly.
[[nodiscard]] KeplerianElements propagateSecular(const KeplerianElements& mean, double dt, const ZonalParameters& zp = {},
                                                 SecularOrder order = SecularOrder::J2SquaredJ4);

/// Inclination (rad) of a sun-synchronous orbit with mean a and e (first-order J2;
/// nodal rate 360 deg per tropical year). Throws astro::Error if no solution exists.
[[nodiscard]] double sunSynchronousInclination(double a, double e = 0.0, const ZonalParameters& zp = {});

/// Mean semi-major axis of a repeat-ground-track orbit that completes `revs`
/// revolutions in `days` nodal days, for given e and i (J2, iterative).
[[nodiscard]] double repeatGroundTrackSemiMajorAxis(int revs, int days, double e, double i,
                                                    const ZonalParameters& zp = {});

/// Critical inclination (rad) where the J2 apsidal rotation vanishes: acos(1/sqrt(5)).
inline constexpr double kCriticalInclination = 1.1071487177940904;

} // namespace astro
