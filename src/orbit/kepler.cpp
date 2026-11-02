/**********************************************************************
 * kepler.cpp
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

#include "astro/orbit/kepler.h"

#include <cmath>

#include "astro/core/angle.h"
#include "astro/core/error.h"

namespace astro {

double trueToEccentricAnomaly(double nu, double e) noexcept {
    return std::atan2(std::sqrt(1.0 - e * e) * std::sin(nu), e + std::cos(nu));
}

double eccentricToTrueAnomaly(double E, double e) noexcept {
    return std::atan2(std::sqrt(1.0 - e * e) * std::sin(E), std::cos(E) - e);
}

double eccentricToMeanAnomaly(double E, double e) noexcept {
    return E - e * std::sin(E);
}

double meanToEccentricAnomaly(double M, double e) {
    if (e < 0.0 || e >= 1.0) {
        throw Error("meanToEccentricAnomaly: requires 0 <= e < 1");
    }
    // Solve on M in [-pi, pi), then restore the full revolutions.
    const double Mw   = wrapPi(M);
    const double revs = M - Mw;
    // Starter (Danby): E0 = M + 0.85 e sign(sin M).
    double E = Mw + 0.85 * e * (std::sin(Mw) >= 0.0 ? 1.0 : -1.0);
    for (int i = 0; i < 50; ++i) {
        const double s = std::sin(E), c = std::cos(E);
        const double f    = E - e * s - Mw;
        const double fp   = 1.0 - e * c;
        const double fpp  = e * s;
        const double fppp = e * c;
        // Third-order (Householder-like) correction.
        double d  = -f / fp;
        d         = -f / (fp + 0.5 * d * fpp);
        d         = -f / (fp + 0.5 * d * fpp + d * d * fppp / 6.0);
        E        += d;
        if (std::abs(d) < 1e-15) {
            break;
        }
    }
    return E + revs;
}

double trueToHyperbolicAnomaly(double nu, double e) noexcept {
    return 2.0 * std::atanh(std::sqrt((e - 1.0) / (e + 1.0)) * std::tan(0.5 * nu));
}

double hyperbolicToTrueAnomaly(double F, double e) noexcept {
    return 2.0 * std::atan(std::sqrt((e + 1.0) / (e - 1.0)) * std::tanh(0.5 * F));
}

double hyperbolicToMeanAnomaly(double F, double e) noexcept {
    return e * std::sinh(F) - F;
}

double meanToHyperbolicAnomaly(double M, double e) {
    if (e <= 1.0) {
        throw Error("meanToHyperbolicAnomaly: requires e > 1");
    }
    double F = std::asinh(M / e); // good starter for large |M|
    if (std::abs(M) < 1.0) {
        F = M / (e - 1.0);
    }
    for (int i = 0; i < 100; ++i) {
        const double f   = e * std::sinh(F) - F - M;
        const double fp  = e * std::cosh(F) - 1.0;
        const double d   = f / fp;
        F               -= d;
        if (std::abs(d) < 1e-15 * (1.0 + std::abs(F))) {
            return F;
        }
    }
    throw Error("meanToHyperbolicAnomaly: did not converge");
}

double trueToParabolicAnomaly(double nu) noexcept {
    return std::tan(0.5 * nu);
}
double parabolicToTrueAnomaly(double D) noexcept {
    return 2.0 * std::atan(D);
}
double parabolicToMeanAnomaly(double D) noexcept {
    return D + D * D * D / 3.0;
}

double meanToParabolicAnomaly(double M) noexcept {
    // Closed-form solution of Barker's equation.
    const double w = std::cbrt(1.5 * M + std::sqrt(1.0 + 2.25 * M * M));
    return w - 1.0 / w;
}

double meanToTrueAnomaly(double M, double e) {
    if (e < 1.0) {
        return eccentricToTrueAnomaly(meanToEccentricAnomaly(M, e), e);
    }
    if (e == 1.0) {
        return parabolicToTrueAnomaly(meanToParabolicAnomaly(M));
    }
    return hyperbolicToTrueAnomaly(meanToHyperbolicAnomaly(M, e), e);
}

double trueToMeanAnomaly(double nu, double e) noexcept {
    if (e < 1.0) {
        return eccentricToMeanAnomaly(trueToEccentricAnomaly(nu, e), e);
    }
    if (e == 1.0) {
        return parabolicToMeanAnomaly(trueToParabolicAnomaly(nu));
    }
    return hyperbolicToMeanAnomaly(trueToHyperbolicAnomaly(nu, e), e);
}

} // namespace astro
