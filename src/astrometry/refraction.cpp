/**********************************************************************
 * refraction.cpp
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

// refractionConstants follows ERFA's eraRefco; refractTrueToObserved and refractObservedToTrue follow
// the refraction steps of eraAtioq and eraAtoiq. See THIRD_PARTY_NOTICES.md.

#include "astro/astrometry/refraction.h"

#include <algorithm>
#include <cmath>

#include "astro/core/constants.h"

namespace astro {

RefractionConstants refractionConstants(const Atmosphere& atm) noexcept {
    const bool   optic = atm.wavelengthMicron <= 100.0;
    const double t     = std::clamp(atm.temperatureC, -150.0, 200.0);
    const double p     = std::clamp(atm.pressureHpa, 0.0, 10000.0);
    const double r     = std::clamp(atm.relativeHumidity, 0.0, 1.0);
    const double w     = std::clamp(atm.wavelengthMicron, 0.1, 1e6);
    double       pw    = 0.0;
    if (p > 0.0) {
        const double ps =
            std::pow(10.0, (0.7859 + 0.03477 * t) / (1.0 + 0.00412 * t)) * (1.0 + p * (4.5e-6 + 6e-10 * t * t));
        pw = r * ps / (1.0 - (1.0 - r) * ps / p);
    }
    const double tk = t + 273.15;
    double       gamma;
    if (optic) {
        const double wlsq = w * w;
        gamma             = ((77.53484e-6 + (4.39108e-7 + 3.666e-9 / wlsq) / wlsq) * p - 11.2684e-6 * pw) / tk;
    } else {
        gamma = (77.6890e-6 * p - (6.3938e-6 - 0.375463 / tk) * pw) / tk;
    }
    double beta = 4.4474e-6 * tk;
    if (!optic) {
        beta -= 0.0074 * pw * beta;
    }
    return {gamma * (1.0 - beta), -gamma * (beta - gamma / 2.0)};
}

double refractTrueToObserved(double trueElevation, const RefractionConstants& k) noexcept {
    // ERFA eraAtioq: work with the topocentric unit vector (r = cos el, z = sin el).
    constexpr double celmin = 1e-6, selmin = 0.05;
    const double     rr  = std::max(std::cos(trueElevation), celmin);
    const double     z   = std::max(std::sin(trueElevation), selmin);
    const double     tz  = rr / z;
    const double     w   = k.B * tz * tz;
    const double     del = (k.A + w) * tz / (1.0 + (k.A + 3.0 * w) / (z * z));
    // As in ERFA, the clamped r and z are used in the update.
    const double cosdel = 1.0 - del * del / 2.0;
    const double f      = cosdel - del * z / rr;
    const double horiz  = std::cos(trueElevation) * f;
    const double up     = cosdel * std::sin(trueElevation) + del * rr;
    return std::atan2(up, horiz);
}

double refractObservedToTrue(double observedElevation, const RefractionConstants& k) noexcept {
    // ERFA eraAtoiq: dZ = (A + B tan^2 Z) tan Z applied to the observed zenith distance.
    constexpr double selmin = 0.05;
    const double     c      = std::cos(observedElevation);
    const double     s      = std::max(std::sin(observedElevation), selmin);
    const double     tz     = c / s;
    const double     dref   = (k.A + k.B * tz * tz) * tz;
    return observedElevation - dref;
}

namespace {

// Refraction in arcminutes, at 1010 hPa and 10 C, from a formula of the form
//     R = c cot(h + p / (h + q)),   h = elevation in degrees:
// Bennett's (c = 1, p = 7.31, q = 4.4; h observed) and Saemundsson's (c = 1.02, p = 10.3, q = 5.11;
// h true). As printed, the formula is not usable over the whole sky, so two things are done here:
//  - At the zenith it is slightly negative (-0.08 and -0.12 arcsec). Meeus (1998), Astronomical
//    Algorithms, 2nd ed., ch. 16, adds 0.0013515 and 0.0019279 arcmin so that R is zero at h = 90
//    deg: those are minus the formula's value at 90 deg, which is what is subtracted here.
//  - The cotangent's argument is smallest, and the refraction largest, at h = sqrt(p) - q (-1.70 and
//    -1.90 deg). Below that the formula turns back down, and at h = -q it divides by zero. The
//    elevation is held at the turning point, so the largest value is returned for everything below.
// The result is finite, continuous, never negative and never rising with elevation; NaN gives NaN.
double cotangentFormulaArcmin(double elevation, double c, double p, double q) noexcept {
    if (elevation >= constants::halfPi) {
        return 0.0; // the zenith, and anything past it
    }
    const auto   formula      = [&](double h) { return c / std::tan((h + p / (h + q)) * constants::degToRad); };
    const double turningPoint = std::sqrt(p) - q;
    const double degrees      = elevation * constants::radToDeg;
    const double held         = degrees < turningPoint ? turningPoint : degrees; // NaN passes through
    const double r            = formula(held) - formula(90.0);
    return r < 0.0 ? 0.0 : r; // negative only by rounding, next to the zenith
}

} // namespace

double refractionBennett(double observedElevation, double pressureHpa, double temperatureC) noexcept {
    const double rArcmin = cotangentFormulaArcmin(observedElevation, 1.0, 7.31, 4.4);
    return rArcmin / 60.0 * constants::degToRad * (pressureHpa / 1010.0) * (283.0 / (273.0 + temperatureC));
}

double refractionSaemundsson(double trueElevation, double pressureHpa, double temperatureC) noexcept {
    const double rArcmin = cotangentFormulaArcmin(trueElevation, 1.02, 10.3, 5.11);
    return rArcmin / 60.0 * constants::degToRad * (pressureHpa / 1010.0) * (283.0 / (273.0 + temperatureC));
}

} // namespace astro
