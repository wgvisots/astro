/**********************************************************************
 * kepler.h
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

/// \file
/// Anomaly conversions and Kepler's equation for elliptic (e < 1), parabolic
/// (e = 1) and hyperbolic (e > 1) orbits. Angles in radians.
namespace astro {

/// Elliptic: true -> eccentric anomaly.
[[nodiscard]] double trueToEccentricAnomaly(double nu, double e) noexcept;
/// Elliptic: eccentric -> true anomaly.
[[nodiscard]] double eccentricToTrueAnomaly(double E, double e) noexcept;
/// Kepler's equation M = E - e sin E.
[[nodiscard]] double eccentricToMeanAnomaly(double E, double e) noexcept;
/// Solve Kepler's equation M = E - e sin E for E, for e in [0, 1) and any M (Newton-Raphson from a
/// starting value that guarantees convergence). E is on the same revolution as M (|E - M| <= e)
/// and satisfies the equation to rounding for every e < 1, including e within 1e-12 of 1 near
/// perigee. A NaN or infinite M gives NaN. Throws astro::Error unless 0 <= e < 1, or if the
/// iteration fails to converge (a guard: it is not known to happen).
[[nodiscard]] double meanToEccentricAnomaly(double M, double e);

/// Hyperbolic: true -> hyperbolic anomaly F (tanh(F/2) = sqrt((e-1)/(e+1)) tan(nu/2)).
[[nodiscard]] double trueToHyperbolicAnomaly(double nu, double e) noexcept;
/// Hyperbolic: hyperbolic anomaly F -> true anomaly.
[[nodiscard]] double hyperbolicToTrueAnomaly(double F, double e) noexcept;
/// M = e sinh F - F.
[[nodiscard]] double hyperbolicToMeanAnomaly(double F, double e) noexcept;
/// Solve M = e sinh F - F for F, for e > 1 and any M (Newton-Raphson from a starting value that
/// guarantees convergence). F has the sign of M and satisfies the equation to rounding for every
/// e > 1, including e within 1e-12 of 1. An infinite M gives an infinite F, a NaN M gives NaN.
/// Throws astro::Error unless e is finite and greater than 1, or if the iteration fails to converge
/// (a guard: it is not known to happen).
[[nodiscard]] double meanToHyperbolicAnomaly(double M, double e);

/// Parabolic: D = tan(nu/2); Barker's equation M = D + D^3/3.
[[nodiscard]] double trueToParabolicAnomaly(double nu) noexcept;
/// Parabolic: D -> true anomaly, nu = 2 atan(D).
[[nodiscard]] double parabolicToTrueAnomaly(double D) noexcept;
/// Barker's equation M = D + D^3/3. Here D = tan(nu / 2) and M is the parabolic mean anomaly
/// sqrt(mu / (2 q^3)) (t - T), with q the perifocal distance and T the time of perifocal passage.
[[nodiscard]] double parabolicToMeanAnomaly(double D) noexcept;
/// Solve Barker's equation M = D + D^3/3 for D (closed form), to rounding for every M of either
/// sign. D has the sign of M and is finite for every finite M; an infinite M gives an infinite D
/// and a NaN M gives NaN.
[[nodiscard]] double meanToParabolicAnomaly(double M) noexcept;

/// Mean anomaly <-> true anomaly for any eccentricity. Throws astro::Error as
/// meanToEccentricAnomaly (e < 1) and meanToHyperbolicAnomaly (e > 1) do: for an e that is negative,
/// infinite or NaN.
[[nodiscard]] double meanToTrueAnomaly(double M, double e);
/// True -> mean anomaly for any eccentricity (elliptic, parabolic or hyperbolic M).
[[nodiscard]] double trueToMeanAnomaly(double nu, double e) noexcept;

} // namespace astro
