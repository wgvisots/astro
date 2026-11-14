/**********************************************************************
 * precession_nutation.cpp
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

// nutationIau2000B, fukushimaWilliams06, fwToMatrix and celestialToIntermediateMatrix follow the
// ERFA routines eraNut00b, eraPfw06, eraFw2m and eraC2ixys. See THIRD_PARTY_NOTICES.md.

#include "astro/frames/precession_nutation.h"

#include <cmath>

#include "../core/enum_io.h"
#include "../data/iau2006_series.h"
#include "astro/core/constants.h"
#include "astro/core/rotation.h"

namespace astro {

using constants::twoPi;
using constants::uasToRad;

namespace {

constexpr double kTurnArcsec = 1296000.0;

double delaunay(double c0, double c1, double c2, double c3, double c4, double t) {
    return std::fmod(c0 + t * (c1 + t * (c2 + t * (c3 + t * c4))), kTurnArcsec) * constants::arcsecToRad;
}

// Evaluate a Poisson series: returns sum over powers of (sum of terms) * t^p, in the table's units.
double evalSeries(const detail::SeriesTerm* terms, std::size_t count, const std::array<double, 14>& fa, double t) {
    double sums[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    for (std::size_t k = count; k-- > 0;) { // smallest terms first
        const auto& term = terms[k];
        double      arg  = 0.0;
        for (int i = 0; i < 14; ++i) {
            if (term.n[i] != 0) {
                arg += term.n[i] * fa[static_cast<std::size_t>(i)];
            }
        }
        sums[term.power] += term.sinCoef * std::sin(arg) + term.cosCoef * std::cos(arg);
    }
    return sums[0] + t * (sums[1] + t * (sums[2] + t * (sums[3] + t * sums[4])));
}

double evalPoly(const double* c, double t) {
    return c[0] + t * (c[1] + t * (c[2] + t * (c[3] + t * (c[4] + t * c[5]))));
}

} // namespace

std::array<double, 14> fundamentalArguments(double t) noexcept {
    std::array<double, 14> f{};
    f[0]  = delaunay(485868.249036, 1717915923.2178, 31.8792, 0.051635, -0.00024470, t);  // l
    f[1]  = delaunay(1287104.793048, 129596581.0481, -0.5532, 0.000136, -0.00001149, t);  // l'
    f[2]  = delaunay(335779.526232, 1739527262.8478, -12.7512, -0.001037, 0.00000417, t); // F
    f[3]  = delaunay(1072260.703692, 1602961601.2090, -6.3706, 0.006593, -0.00003169, t); // D
    f[4]  = delaunay(450160.398036, -6962890.5431, 7.4722, 0.007702, -0.00005939, t);     // Om
    f[5]  = std::fmod(4.402608842 + 2608.7903141574 * t, twoPi);                          // Mercury
    f[6]  = std::fmod(3.176146697 + 1021.3285546211 * t, twoPi);                          // Venus
    f[7]  = std::fmod(1.753470314 + 628.3075849991 * t, twoPi);                           // Earth
    f[8]  = std::fmod(6.203480913 + 334.0612426700 * t, twoPi);                           // Mars
    f[9]  = std::fmod(0.599546497 + 52.9690962641 * t, twoPi);                            // Jupiter
    f[10] = std::fmod(0.874016757 + 21.3299104960 * t, twoPi);                            // Saturn
    f[11] = std::fmod(5.481293872 + 7.4781598567 * t, twoPi);                             // Uranus
    f[12] = std::fmod(5.311886287 + 3.8133035638 * t, twoPi);                             // Neptune
    f[13] = (0.02438175 + 0.00000538691 * t) * t; // general precession in longitude
    return f;
}

Nutation nutationIau2000A(double t) noexcept {
    const auto fa = fundamentalArguments(t);
    return {evalSeries(detail::kNutLongitude, detail::kNutLongitudeCount, fa, t) * uasToRad,
            evalSeries(detail::kNutObliquity, detail::kNutObliquityCount, fa, t) * uasToRad};
}

Nutation nutationIau2000B(double t) noexcept {
    // Simon et al. (1994) linear Delaunay arguments, as used by IAU 2000B.
    const double el  = std::fmod(485868.249036 + 1717915923.2178 * t, kTurnArcsec) * constants::arcsecToRad;
    const double elp = std::fmod(1287104.79305 + 129596581.0481 * t, kTurnArcsec) * constants::arcsecToRad;
    const double f   = std::fmod(335779.526232 + 1739527262.8478 * t, kTurnArcsec) * constants::arcsecToRad;
    const double d   = std::fmod(1072260.70369 + 1602961601.2090 * t, kTurnArcsec) * constants::arcsecToRad;
    const double om  = std::fmod(450160.398036 - 6962890.5431 * t, kTurnArcsec) * constants::arcsecToRad;
    double       dp = 0.0, de = 0.0;
    for (int i = 76; i >= 0; --i) {
        const auto&  x   = detail::kNut2000B[i];
        const double arg = std::fmod(x.n[0] * el + x.n[1] * elp + x.n[2] * f + x.n[3] * d + x.n[4] * om, twoPi);
        const double s = std::sin(arg), c = std::cos(arg);
        dp += (x.ps + x.pst * t) * s + x.pc * c;
        de += (x.ec + x.ect * t) * c + x.es * s;
    }
    constexpr double u2r    = constants::arcsecToRad / 1e7;
    constexpr double dpplan = -0.135 * constants::masToRad;
    constexpr double deplan = 0.388 * constants::masToRad;
    return {dp * u2r + dpplan, de * u2r + deplan};
}

Nutation nutation(double t, PrecessionNutationModel model) noexcept {
    return model == PrecessionNutationModel::IAU2006_2000A ? nutationIau2000A(t) : nutationIau2000B(t);
}

double meanObliquity06(double t) noexcept {
    return (84381.406 + t * (-46.836769 + t * (-0.0001831 + t * (0.00200340 + t * (-0.000000576 + t * -0.0000000434)))))
         * constants::arcsecToRad;
}

FukushimaWilliams fukushimaWilliams06(double t) noexcept {
    FukushimaWilliams fw;
    fw.gammaBar =
        (-0.052928 + t * (10.556378 + t * (0.4932044 + t * (-0.00031238 + t * (-0.000002788 + t * 0.0000000260)))))
        * constants::arcsecToRad;
    fw.phiBar =
        (84381.412819 + t * (-46.811016 + t * (0.0511268 + t * (0.00053289 + t * (-0.000000440 + t * -0.0000000176)))))
        * constants::arcsecToRad;
    fw.psiBar =
        (-0.041775 + t * (5038.481484 + t * (1.5584175 + t * (-0.00018522 + t * (-0.000026452 + t * -0.0000000148)))))
        * constants::arcsecToRad;
    fw.epsA = meanObliquity06(t);
    return fw;
}

Mat3 fwToMatrix(double gamma, double phi, double psi, double eps) noexcept {
    return R1(-eps) * R3(-psi) * R1(phi) * R3(gamma);
}

Mat3 frameBiasMatrix06() noexcept {
    const FukushimaWilliams fw = fukushimaWilliams06(0.0);
    return fwToMatrix(fw.gammaBar, fw.phiBar, fw.psiBar, fw.epsA);
}

Mat3 biasPrecessionMatrix06(double t) noexcept {
    const FukushimaWilliams fw = fukushimaWilliams06(t);
    return fwToMatrix(fw.gammaBar, fw.phiBar, fw.psiBar, fw.epsA);
}

Mat3 precessionMatrix06(double t) noexcept {
    return biasPrecessionMatrix06(t) * frameBiasMatrix06().transpose();
}

Mat3 nutationMatrix(double epsA, const Nutation& n) noexcept {
    return R1(-(epsA + n.deps)) * R3(-n.dpsi) * R1(epsA);
}

Mat3 npbMatrix06(double t, PrecessionNutationModel model, double ddpsi, double ddeps) noexcept {
    const FukushimaWilliams fw = fukushimaWilliams06(t);
    const Nutation          n  = nutation(t, model);
    return fwToMatrix(fw.gammaBar, fw.phiBar, fw.psiBar + n.dpsi + ddpsi, fw.epsA + n.deps + ddeps);
}

CipXY cipXY06(double t) noexcept {
    const auto fa = fundamentalArguments(t);
    return {(evalPoly(detail::kCipXPoly, t) + evalSeries(detail::kCipX, detail::kCipXCount, fa, t)) * uasToRad,
            (evalPoly(detail::kCipYPoly, t) + evalSeries(detail::kCipY, detail::kCipYCount, fa, t)) * uasToRad};
}

CipXY matrixToCipXY(const Mat3& npb) noexcept {
    return {npb(2, 0), npb(2, 1)};
}

CipXY cipXY(double t, PrecessionNutationModel model) noexcept {
    return model == PrecessionNutationModel::IAU2006_2000A ? cipXY06(t) : matrixToCipXY(npbMatrix06(t, model));
}

double cioLocatorS06(double t, const CipXY& xy) noexcept {
    const auto   fa = fundamentalArguments(t);
    const double sxy2 =
        (evalPoly(detail::kSXY2Poly, t) + evalSeries(detail::kSXY2, detail::kSXY2Count, fa, t)) * uasToRad;
    return sxy2 - 0.5 * xy.x * xy.y;
}

double tioLocatorSp(double t) noexcept {
    return -47e-6 * t * constants::arcsecToRad;
}

Mat3 celestialToIntermediateMatrix(const CipXY& xy, double s) noexcept {
    const double r2 = xy.x * xy.x + xy.y * xy.y;
    const double e  = (r2 > 0.0) ? std::atan2(xy.y, xy.x) : 0.0;
    const double d  = std::atan(std::sqrt(r2 / (1.0 - r2)));
    return R3(-(e + s)) * R2(d) * R3(e);
}

Nutation poleOffsetsToNutation(double t, double dX, double dY) noexcept {
    // dX = dpsi sin(epsA) + c deps,  dY = deps - c dpsi sin(epsA),  c = psiA cos(eps0) - chiA
    const double psiA       = (5038.481507 * t - 1.0790069 * t * t) * constants::arcsecToRad;
    const double chiA       = (10.556403 * t - 2.3814292 * t * t) * constants::arcsecToRad;
    const double eps0       = 84381.406 * constants::arcsecToRad;
    const double c          = psiA * std::cos(eps0) - chiA;
    const double k          = 1.0 + c * c;
    const double dpsiSinEps = (dX - c * dY) / k;
    const double deps       = (dY + c * dX) / k;
    return {dpsiSinEps / std::sin(meanObliquity06(t)), deps};
}

std::ostream& operator<<(std::ostream& os, PrecessionNutationModel value) {
    return detail::writeEnum(os, "PrecessionNutationModel", value);
}

} // namespace astro
