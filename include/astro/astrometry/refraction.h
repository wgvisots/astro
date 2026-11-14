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

/// Refraction constants for the model  dZ = A tan Z + B tan^3 Z  (rad), computed as in ERFA's eraRefco
/// (Green 1987 / Hohenkerk & Sinclair). Valid for zenith distances up to ~75-80 deg. Nearer the
/// horizon it is not the conventional refraction (about 11' at a true elevation of 0, where the
/// usual figure is 34'), and refractTrueToObserved and refractObservedToTrue are no longer exact
/// inverses of each other (27" apart at 0-2 deg): use refractionBennett or refractionSaemundsson
/// there.
struct RefractionConstants {
    double A{0.0}; ///< coefficient of tan Z (rad)
    double B{0.0}; ///< coefficient of tan^3 Z (rad)
};
/// Refraction constants for the given conditions (inputs are clamped to safe
/// ranges, as in ERFA's eraRefco).
[[nodiscard]] RefractionConstants refractionConstants(const Atmosphere& atm) noexcept;

/// Observed elevation from true elevation using the A/B model, with ERFA's (eraAtioq)
/// treatment near the horizon (refraction is held constant below ~3 deg).
[[nodiscard]] double refractTrueToObserved(double trueElevation, const RefractionConstants& k) noexcept;
/// True elevation from observed elevation (A/B model).
[[nodiscard]] double refractObservedToTrue(double observedElevation, const RefractionConstants& k) noexcept;

/// Bennett (1982): refraction (rad) for an *observed* elevation; scaled for
/// pressure and temperature. Good to ~0.1 arcmin down to the horizon.
///
/// The formula is  R = cot(h + 7.31 / (h + 4.4)) arcmin  for h in degrees, less its value at
/// h = 90 deg (-0.08 arcsec), so that the refraction is exactly zero at the zenith (the correction
/// given by Meeus 1998, ch. 16). It is meaningful for observed elevations from the horizon (0) to the
/// zenith. Below the horizon it is an extrapolation, and it is largest at -1.70 deg (56.9 arcmin
/// at 1010 hPa and 10 C); for every elevation below that, where the formula as printed turns back
/// down and then divides by zero, that largest value is returned. The result is therefore finite,
/// continuous and never negative for any elevation; a NaN elevation gives NaN.
[[nodiscard]] double
refractionBennett(double observedElevation, double pressureHpa = 1010.0, double temperatureC = 10.0) noexcept;
/// Saemundsson (1986): refraction (rad) for a *true* elevation; scaled for pressure and
/// temperature as refractionBennett(), with which it agrees to ~0.1 arcmin.
///
/// The formula is  R = 1.02 cot(h + 10.3 / (h + 5.11)) arcmin  for h in degrees, less its value at
/// h = 90 deg (-0.12 arcsec), so that the refraction is exactly zero at the zenith (Meeus 1998,
/// ch. 16). It is meaningful for true elevations from about -0.6 deg (a body seen on the horizon)
/// to the zenith. Below that it is an extrapolation, and it is largest at -1.90 deg (44.6 arcmin at
/// 1010 hPa and 10 C); for every elevation below that the largest value is returned. The result is
/// finite, continuous and never negative for any elevation; a NaN elevation gives NaN.
[[nodiscard]] double
refractionSaemundsson(double trueElevation, double pressureHpa = 1010.0, double temperatureC = 10.0) noexcept;

} // namespace astro
