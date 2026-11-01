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

#include "astro/astrometry/refraction.h"

#include <algorithm>
#include <cmath>

#include "astro/core/constants.h"

namespace astro {

RefractionConstants refractionConstants(const Atmosphere& atm) noexcept {
    const bool optic = atm.wavelengthMicron <= 100.0;
    const double t = std::clamp(atm.temperatureC, -150.0, 200.0);
    const double p = std::clamp(atm.pressureHpa, 0.0, 10000.0);
    const double r = std::clamp(atm.relativeHumidity, 0.0, 1.0);
    const double w = std::clamp(atm.wavelengthMicron, 0.1, 1e6);
    double pw = 0.0;
    if (p > 0.0) {
        const double ps = std::pow(10.0, (0.7859 + 0.03477 * t) / (1.0 + 0.00412 * t)) * (1.0 + p * (4.5e-6 + 6e-10 * t * t));
        pw = r * ps / (1.0 - (1.0 - r) * ps / p);
    }
    const double tk = t + 273.15;
    double gamma;
    if (optic) {
        const double wlsq = w * w;
        gamma = ((77.53484e-6 + (4.39108e-7 + 3.666e-9 / wlsq) / wlsq) * p - 11.2684e-6 * pw) / tk;
    } else {
        gamma = (77.6890e-6 * p - (6.3938e-6 - 0.375463 / tk) * pw) / tk;
    }
    double beta = 4.4474e-6 * tk;
    if (!optic) beta -= 0.0074 * pw * beta;
    return {gamma * (1.0 - beta), -gamma * (beta - gamma / 2.0)};
}

double refractTrueToObserved(double trueElevation, const RefractionConstants& k) noexcept {
    // SOFA atioq: work with the topocentric unit vector (r = cos el, z = sin el).
    constexpr double celmin = 1e-6, selmin = 0.05;
    const double rr = std::max(std::cos(trueElevation), celmin);
    const double z = std::max(std::sin(trueElevation), selmin);
    const double tz = rr / z;
    const double w = k.B * tz * tz;
    const double del = (k.A + w) * tz / (1.0 + (k.A + 3.0 * w) / (z * z));
    // As in SOFA, the clamped r and z are used in the update.
    const double cosdel = 1.0 - del * del / 2.0;
    const double f = cosdel - del * z / rr;
    const double horiz = std::cos(trueElevation) * f;
    const double up = cosdel * std::sin(trueElevation) + del * rr;
    return std::atan2(up, horiz);
}

double refractObservedToTrue(double observedElevation, const RefractionConstants& k) noexcept {
    // SOFA atoiq: dZ = (A + B tan^2 Z) tan Z applied to the observed zenith distance.
    constexpr double selmin = 0.05;
    const double c = std::cos(observedElevation);
    const double s = std::max(std::sin(observedElevation), selmin);
    const double tz = c / s;
    const double dref = (k.A + k.B * tz * tz) * tz;
    return observedElevation - dref;
}

double refractionBennett(double observedElevation, double pressureHpa, double temperatureC) noexcept {
    const double h = observedElevation * constants::radToDeg;
    const double rArcmin = 1.0 / std::tan((h + 7.31 / (h + 4.4)) * constants::degToRad);
    return rArcmin / 60.0 * constants::degToRad * (pressureHpa / 1010.0) * (283.0 / (273.0 + temperatureC));
}

double refractionSaemundsson(double trueElevation, double pressureHpa, double temperatureC) noexcept {
    const double h = trueElevation * constants::radToDeg;
    const double rArcmin = 1.02 / std::tan((h + 10.3 / (h + 5.11)) * constants::degToRad);
    return rArcmin / 60.0 * constants::degToRad * (pressureHpa / 1010.0) * (283.0 / (273.0 + temperatureC));
}

} // namespace astro
