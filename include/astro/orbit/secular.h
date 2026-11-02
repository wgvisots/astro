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

#include <iosfwd>

#include <boost/describe/enum.hpp>

#include "astro/core/constants.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/state12.h"

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

BOOST_DESCRIBE_ENUM(SecularOrder, J2, J2SquaredJ4)

/// Writes "SecularOrder::<enumerator>", or "SecularOrder::INVALID_VALUE" for any other value.
std::ostream& operator<<(std::ostream& os, SecularOrder value);

/// Secular rates (rad/s) of the mean elements.
struct SecularRates {
    double raanDot{0.0};
    double argpDot{0.0};
    double meanAnomalyDot{0.0}; ///< total rate including the mean motion
};

/// Rates for mean semi-major axis a, eccentricity e and inclination i.
[[nodiscard]] SecularRates secularRates(double                 a,
                                        double                 e,
                                        double                 i,
                                        const ZonalParameters& zp    = {},
                                        SecularOrder           order = SecularOrder::J2SquaredJ4) noexcept;

/// Propagate mean elements by dt seconds with the secular rates (a, e, i constant).
/// The true anomaly is advanced through the mean anomaly.
[[nodiscard]] KeplerianElements propagateSecular(const KeplerianElements& mean,
                                                 double                   dt,
                                                 const ZonalParameters&   zp    = {},
                                                 SecularOrder             order = SecularOrder::J2SquaredJ4);

/// propagateSecular for an inertial state whose z axis is the Earth's pole. The osculating
/// elements of s0 are treated as mean elements (no short-period correction), so expect errors at the
/// short-period level (kilometres in LEO); use NumericalPropagator::j2j4 when that matters.
[[nodiscard]] StateVector propagateSecular(const StateVector&     s0,
                                           double                 dt,
                                           const ZonalParameters& zp    = {},
                                           SecularOrder           order = SecularOrder::J2SquaredJ4);

/// propagateSecular for a 12-component state (position, velocity, acceleration, jerk), with the
/// same mean-element caveat. The force model for acceleration, jerk and residuals is the point mass
/// plus J2 (SecularOrder::J2) or plus J2 and J4 (SecularOrder::J2SquaredJ4).
///  - ModelOnly: secular propagation from position and velocity.
///  - UnmodeledForce: the unmodeled force da0 + dj0 t is superposed on the secular trajectory,
///    dr = da0 dt^2/2 + dj0 dt^3/6, dv = da0 dt + dj0 dt^2/2, ignoring how the displacement feeds
///    back into gravity. That is fine while dr stays small compared with the orbit; otherwise use
///    NumericalPropagator::j2j4 with HigherOrderMode::UnmodeledForce.
[[nodiscard]] Propagation12 propagateSecular(const StateVector12&   s0,
                                             double                 dt,
                                             const ZonalParameters& zp    = {},
                                             SecularOrder           order = SecularOrder::J2SquaredJ4,
                                             HigherOrderMode        mode  = HigherOrderMode::ModelOnly);

/// Inclination (rad) of a sun-synchronous orbit with mean a and e (first-order J2;
/// nodal rate 360 deg per tropical year). Throws astro::Error if no solution exists.
[[nodiscard]] double sunSynchronousInclination(double a, double e = 0.0, const ZonalParameters& zp = {});

/// Mean semi-major axis of a repeat-ground-track orbit that completes `revs`
/// revolutions in `days` nodal days, for given e and i (J2, iterative).
[[nodiscard]] double
repeatGroundTrackSemiMajorAxis(int revs, int days, double e, double i, const ZonalParameters& zp = {});

/// Critical inclination (rad) where the J2 apsidal rotation vanishes: acos(1/sqrt(5)).
inline constexpr double kCriticalInclination = 1.1071487177940904;

} // namespace astro
