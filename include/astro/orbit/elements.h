/**********************************************************************
 * elements.h
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
#include "astro/core/mat6.h"
#include "astro/core/vec.h"

/// Classical and equinoctial orbital elements. SI units, radians.
namespace astro {

/// Classical (Keplerian) osculating elements.
///
/// Special cases (e and i below 1e-11) follow Vallado: for circular orbits the
/// argument of perigee is 0 and `nu` holds the argument of latitude; for
/// equatorial orbits the RAAN is 0 and `argp` holds the longitude of
/// periapsis; for circular equatorial orbits both are 0 and `nu` holds the
/// true longitude. Hyperbolic orbits have a < 0. Exactly parabolic orbits
/// (a infinite) are not representable here; use states.
struct KeplerianElements {
    double a{0.0};    ///< semi-major axis (m)
    double e{0.0};    ///< eccentricity
    double i{0.0};    ///< inclination
    double raan{0.0}; ///< right ascension of the ascending node
    double argp{0.0}; ///< argument of periapsis
    double nu{0.0};   ///< true anomaly

    [[nodiscard]] double semiLatusRectum() const noexcept { return a * (1.0 - e * e); }
    [[nodiscard]] double meanAnomaly() const noexcept;
    [[nodiscard]] double eccentricAnomaly() const noexcept;
};

[[nodiscard]] KeplerianElements stateToElements(const StateVector& s, double mu = constants::gmEarth);
[[nodiscard]] StateVector       elementsToState(const KeplerianElements& el, double mu = constants::gmEarth);
/// Elements from mean anomaly instead of true anomaly.
[[nodiscard]] KeplerianElements
elementsFromMeanAnomaly(double a, double e, double i, double raan, double argp, double M);

/// Modified equinoctial-style (prograde) equinoctial elements:
///   h = e sin(raan + argp), k = e cos(raan + argp),
///   p = tan(i/2) sin(raan), q = tan(i/2) cos(raan), lambda = mean longitude.
/// Non-singular for circular and equatorial orbits (singular at i = 180 deg).
struct EquinoctialElements {
    double a{0.0};
    double h{0.0};
    double k{0.0};
    double p{0.0};
    double q{0.0};
    double lambda{0.0};
};

[[nodiscard]] EquinoctialElements keplerianToEquinoctial(const KeplerianElements& el) noexcept;
[[nodiscard]] KeplerianElements   equinoctialToKeplerian(const EquinoctialElements& eq);
[[nodiscard]] EquinoctialElements stateToEquinoctial(const StateVector& s, double mu = constants::gmEarth);
[[nodiscard]] StateVector         equinoctialToState(const EquinoctialElements& eq, double mu = constants::gmEarth);

/// Partial derivatives of stateToEquinoctial: row i is element i (a, h, k, p, q, lambda), column j
/// is state component j (x, y, z, vx, vy, vz). Central differences with steps of 1e-7 |r| and
/// 1e-7 |v| (relative accuracy ~1e-8); lambda differences are wrapped. Use it to map a position /
/// velocity covariance P into elements as J P J^T. Same singularity as stateToEquinoctial.
[[nodiscard]] Mat6 stateToEquinoctialJacobian(const StateVector& s, double mu = constants::gmEarth);

/// Orbital period (s) for a > 0.
[[nodiscard]] double orbitalPeriod(double a, double mu = constants::gmEarth) noexcept;
/// Mean motion (rad/s).
[[nodiscard]] double meanMotion(double a, double mu = constants::gmEarth) noexcept;
/// Semi-major axis from mean motion (rad/s).
[[nodiscard]] double semiMajorAxisFromMeanMotion(double n, double mu = constants::gmEarth) noexcept;
/// Specific orbital energy (J/kg).
[[nodiscard]] double specificEnergy(const StateVector& s, double mu = constants::gmEarth) noexcept;
/// Flight-path angle (rad): angle of the velocity above the local horizontal.
[[nodiscard]] double        flightPathAngle(const StateVector& s) noexcept;
[[nodiscard]] inline double periapsisRadius(double a, double e) noexcept {
    return a * (1.0 - e);
}
[[nodiscard]] inline double apoapsisRadius(double a, double e) noexcept {
    return a * (1.0 + e);
}

} // namespace astro
