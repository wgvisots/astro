/**********************************************************************
 * refraction.h
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
/// Atmospheric refraction. Elevations in radians; refraction is the amount by
/// which the observed (apparent) elevation exceeds the true (geometric) one.
namespace astro {

/// Meteorological conditions at the observer.
struct Atmosphere {
    double pressureHpa{1013.25};    ///< air pressure at the observer, hPa (0 = no refraction)
    double temperatureC{15.0};      ///< air temperature, degrees Celsius
    double relativeHumidity{0.5};   ///< 0..1
    double wavelengthMicron{0.574}; ///< > 100 means radio
};

/// Refraction constants for the model  dZ = A tan Z + B tan^3 Z  (rad), computed as in SOFA/ERFA refco
/// (Green 1987 / Hohenkerk & Sinclair). Valid for zenith distances up to ~75-80 deg.
struct RefractionConstants {
    double A{0.0}; ///< coefficient of tan Z (rad)
    double B{0.0}; ///< coefficient of tan^3 Z (rad)
};
/// Refraction constants for the given conditions (inputs are clamped to safe
/// ranges, as in SOFA refco).
[[nodiscard]] RefractionConstants refractionConstants(const Atmosphere& atm) noexcept;

/// Observed elevation from true elevation using the A/B model, with the SOFA
/// treatment near the horizon (refraction is held constant below ~3 deg).
[[nodiscard]] double refractTrueToObserved(double trueElevation, const RefractionConstants& k) noexcept;
/// True elevation from observed elevation (A/B model).
[[nodiscard]] double refractObservedToTrue(double observedElevation, const RefractionConstants& k) noexcept;

/// Bennett (1982): refraction (rad) for an *observed* elevation; scaled for
/// pressure and temperature. Good to ~0.1 arcmin down to the horizon.
[[nodiscard]] double
refractionBennett(double observedElevation, double pressureHpa = 1010.0, double temperatureC = 10.0) noexcept;
/// Saemundsson (1986): refraction (rad) for a *true* elevation.
[[nodiscard]] double
refractionSaemundsson(double trueElevation, double pressureHpa = 1010.0, double temperatureC = 10.0) noexcept;

} // namespace astro
