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
[[nodiscard]] double eccentricToTrueAnomaly(double E, double e) noexcept;
/// Kepler's equation M = E - e sin E.
[[nodiscard]] double eccentricToMeanAnomaly(double E, double e) noexcept;
/// Solve Kepler's equation for E (Newton-Raphson with a robust starter), e in [0, 1).
[[nodiscard]] double meanToEccentricAnomaly(double M, double e);

/// Hyperbolic: true -> hyperbolic anomaly F (tanh(F/2) = sqrt((e-1)/(e+1)) tan(nu/2)).
[[nodiscard]] double trueToHyperbolicAnomaly(double nu, double e) noexcept;
[[nodiscard]] double hyperbolicToTrueAnomaly(double F, double e) noexcept;
/// M = e sinh F - F.
[[nodiscard]] double hyperbolicToMeanAnomaly(double F, double e) noexcept;
[[nodiscard]] double meanToHyperbolicAnomaly(double M, double e);

/// Parabolic: D = tan(nu/2); Barker's equation M = D + D^3/3.
[[nodiscard]] double trueToParabolicAnomaly(double nu) noexcept;
[[nodiscard]] double parabolicToTrueAnomaly(double D) noexcept;
[[nodiscard]] double parabolicToMeanAnomaly(double D) noexcept;
[[nodiscard]] double meanToParabolicAnomaly(double M) noexcept;

/// Mean anomaly <-> true anomaly for any eccentricity.
[[nodiscard]] double meanToTrueAnomaly(double M, double e);
[[nodiscard]] double trueToMeanAnomaly(double nu, double e) noexcept;

} // namespace astro
