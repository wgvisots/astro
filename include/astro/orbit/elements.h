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

/// \file
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

    /// Semi-latus rectum p = a (1 - e^2) (m); positive for elliptic and hyperbolic orbits.
    [[nodiscard]] double semiLatusRectum() const noexcept { return a * (1.0 - e * e); }
    /// Mean anomaly (rad) from `nu`: M = E - e sin E (elliptic) or M = e sinh F - F (hyperbolic).
    [[nodiscard]] double meanAnomaly() const noexcept;
    /// Eccentric anomaly E (rad) for e < 1; the hyperbolic anomaly F for e > 1.
    [[nodiscard]] double eccentricAnomaly() const noexcept;
};

/// Osculating elements of an inertial state (m, m/s) about a body with gravitational parameter
/// `mu` (m^3/s^2). Inclination and RAAN refer to the z axis and x-y plane of the state's frame.
/// RAAN, argument of periapsis and elliptic true anomalies are in [0, 2 pi); hyperbolic true
/// anomalies are in [-pi, pi). The angles are accurate to rounding everywhere, including
/// inclinations and anomalies close to 0 and pi, so elementsToState gives the state back to
/// rounding (an inclination or eccentricity that is not zero but below the 1e-11 of the special
/// cases is the exception: the angle that is set to 0 is lost). Throws astro::Error for a
/// rectilinear state (r = 0 or r parallel to v) or a parabolic orbit (|e - 1| < 1e-14).
[[nodiscard]] KeplerianElements stateToElements(const StateVector& s, double mu = constants::gmEarth);
/// Inertial position (m) and velocity (m/s) from elements; the inverse of stateToElements.
/// Throws astro::Error if the semi-latus rectum is not positive (inconsistent a and e).
[[nodiscard]] StateVector elementsToState(const KeplerianElements& el, double mu = constants::gmEarth);
/// Elements from mean anomaly instead of true anomaly. Throws astro::Error as meanToTrueAnomaly
/// does: for an e that is negative, infinite or NaN.
[[nodiscard]] KeplerianElements
elementsFromMeanAnomaly(double a, double e, double i, double raan, double argp, double M);

/// Modified equinoctial-style (prograde) equinoctial elements:
///   h = e sin(raan + argp), k = e cos(raan + argp),
///   p = tan(i/2) sin(raan), q = tan(i/2) cos(raan), lambda = mean longitude.
/// Non-singular for circular and equatorial orbits (singular at i = 180 deg).
/// Defined for elliptic orbits only (a > 0, e < 1), as in Broucke & Cefola (1972): the mean
/// longitude is an angle only when the mean anomaly is one, and a hyperbolic mean anomaly is not.
struct EquinoctialElements {
    double a{0.0};      ///< semi-major axis (m)
    double h{0.0};      ///< e sin(raan + argp)
    double k{0.0};      ///< e cos(raan + argp)
    double p{0.0};      ///< tan(i/2) sin(raan)
    double q{0.0};      ///< tan(i/2) cos(raan)
    double lambda{0.0}; ///< mean longitude raan + argp + M (rad, [0, 2pi))
};

/// Equinoctial elements from Keplerian ones (lambda from the mean anomaly). If the orbit is not
/// elliptic (a <= 0 or e >= 1, or either is NaN) every element returned is NaN.
[[nodiscard]] EquinoctialElements keplerianToEquinoctial(const KeplerianElements& el) noexcept;
/// Keplerian elements from equinoctial ones (RAAN and argp set to 0 where undefined). Throws
/// astro::Error unless the elements are those of an ellipse (a > 0, h^2 + k^2 < 1).
[[nodiscard]] KeplerianElements equinoctialToKeplerian(const EquinoctialElements& eq);
/// Equinoctial elements of an inertial state (m, m/s), computed directly from the vectors (no
/// classical singularities). Throws astro::Error for a retrograde equatorial orbit (i = 180 deg)
/// and for a state that is not on an elliptic orbit (parabolic, hyperbolic or rectilinear).
[[nodiscard]] EquinoctialElements stateToEquinoctial(const StateVector& s, double mu = constants::gmEarth);
/// Inertial position (m) and velocity (m/s) from equinoctial elements, via equinoctialToKeplerian
/// and elementsToState; throws astro::Error as they do.
[[nodiscard]] StateVector equinoctialToState(const EquinoctialElements& eq, double mu = constants::gmEarth);

/// Partial derivatives of stateToEquinoctial: row i is element i (a, h, k, p, q, lambda), column j
/// is state component j (x, y, z, vx, vy, vz). Central differences with steps of 1e-7 |r| and
/// 1e-7 |v| (relative accuracy ~1e-8); lambda differences are wrapped. Use it to map a position /
/// velocity covariance P into elements as J P J^T. Same singularity as stateToEquinoctial, and
/// throws astro::Error where it does, for the state or any of the twelve perturbed states.
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
[[nodiscard]] double flightPathAngle(const StateVector& s) noexcept;
/// Periapsis radius a (1 - e) (m).
[[nodiscard]] inline double periapsisRadius(double a, double e) noexcept {
    return a * (1.0 - e);
}
/// Apoapsis radius a (1 + e) (m); meaningful for elliptic orbits only.
[[nodiscard]] inline double apoapsisRadius(double a, double e) noexcept {
    return a * (1.0 + e);
}

} // namespace astro
