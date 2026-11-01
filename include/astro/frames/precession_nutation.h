/**********************************************************************
 * precession_nutation.h
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

#include <array>

#include "astro/core/vec.h"

/// IAU 2006 precession and IAU 2000A/2000B nutation (IERS Conventions 2010, Ch. 5).
///
/// All functions take `t`, Julian centuries of TT since J2000.0
/// (astro::Time::to(TimeScale::TT).centuriesSinceJ2000()), and return radians.
namespace astro {

/// Precession-nutation model selection.
enum class PrecessionNutationModel {
    IAU2006_2000A, ///< full model, ~microarcsecond (default)
    IAU2006_2000B, ///< truncated nutation, ~1 mas, ~20x faster
};

/// The 14 fundamental arguments of nutation theory (IERS Conventions 2003/2010, eq. 5.43),
/// radians: l, l', F, D, Om, L_Me, L_Ve, L_E, L_Ma, L_J, L_Sa, L_U, L_Ne, p_A.
[[nodiscard]] std::array<double, 14> fundamentalArguments(double t) noexcept;

/// Nutation in longitude and obliquity (rad).
struct Nutation {
    double dpsi{0.0};
    double deps{0.0};
};

/// IAU 2000A nutation with the IAU 2006 (R06) adjustments: 1358 + 1056 terms.
[[nodiscard]] Nutation nutationIau2000A(double t) noexcept;
/// IAU 2000B nutation: 77 luni-solar terms plus fixed planetary offsets.
[[nodiscard]] Nutation nutationIau2000B(double t) noexcept;
[[nodiscard]] Nutation nutation(double t, PrecessionNutationModel model) noexcept;

/// IAU 2006 mean obliquity of the ecliptic.
[[nodiscard]] double meanObliquity06(double t) noexcept;

/// IAU 2006 Fukushima-Williams bias-precession angles (rad).
struct FukushimaWilliams {
    double gammaBar{0.0};
    double phiBar{0.0};
    double psiBar{0.0};
    double epsA{0.0};  ///< mean obliquity of date
};
[[nodiscard]] FukushimaWilliams fukushimaWilliams06(double t) noexcept;
/// Rotation matrix from Fukushima-Williams angles: R1(-eps) R3(-psi) R1(phi) R3(gamma).
[[nodiscard]] Mat3 fwToMatrix(double gamma, double phi, double psi, double eps) noexcept;

/// GCRS -> mean J2000.0 (EME2000) frame bias matrix (IAU 2006).
[[nodiscard]] Mat3 frameBiasMatrix06() noexcept;
/// GCRS -> mean equator and equinox of date (bias + precession).
[[nodiscard]] Mat3 biasPrecessionMatrix06(double t) noexcept;
/// Mean J2000.0 -> mean of date (precession only).
[[nodiscard]] Mat3 precessionMatrix06(double t) noexcept;
/// Nutation matrix: mean of date -> true of date.
[[nodiscard]] Mat3 nutationMatrix(double epsA, const Nutation& n) noexcept;
/// GCRS -> true equator and equinox of date (bias-precession-nutation), with
/// optional extra nutation corrections (rad) added to dpsi/deps.
[[nodiscard]] Mat3 npbMatrix06(double t, PrecessionNutationModel model = PrecessionNutationModel::IAU2006_2000A,
                               double ddpsi = 0.0, double ddeps = 0.0) noexcept;

/// Coordinates of the Celestial Intermediate Pole in the GCRS (rad).
struct CipXY {
    double x{0.0};
    double y{0.0};
};
/// CIP X, Y from the IAU 2006/2000A series (IERS Tables 5.2a/b).
[[nodiscard]] CipXY cipXY06(double t) noexcept;
/// CIP X, Y extracted from a bias-precession-nutation matrix (used for 2000B).
[[nodiscard]] CipXY cipXYFromMatrix(const Mat3& npb) noexcept;
[[nodiscard]] CipXY cipXY(double t, PrecessionNutationModel model) noexcept;
/// CIO locator s, given the CIP X, Y (IERS Table 5.2d).
[[nodiscard]] double cioLocatorS06(double t, const CipXY& xy) noexcept;
/// TIO locator s' (IERS 2010 eq. 5.13).
[[nodiscard]] double tioLocatorSp(double t) noexcept;
/// GCRS -> CIRS matrix from X, Y, s.
[[nodiscard]] Mat3 celestialToIntermediateMatrix(const CipXY& xy, double s) noexcept;

/// Convert celestial pole offsets dX, dY to nutation corrections ddpsi, ddeps
/// (first-order, IERS 2010 eq. 5.27 inverse).
[[nodiscard]] Nutation poleOffsetsToNutation(double t, double dX, double dY) noexcept;

} // namespace astro
