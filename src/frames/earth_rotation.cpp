/**********************************************************************
 * earth_rotation.cpp
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

#include "astro/frames/earth_rotation.h"

#include <cmath>

#include "astro/core/angle.h"
#include "astro/core/rotation.h"

namespace astro {

double earthRotationAngle(const JulianDate& ut1) noexcept {
    const double t = ut1.daysSinceJ2000();
    // Fractional part of the date, taken from both parts to keep precision.
    const double f = std::fmod(ut1.jd1, 1.0) + std::fmod(ut1.jd2, 1.0);
    return wrapTwoPi(constants::twoPi * (f + 0.7790572732640 + 0.00273781191135448 * t));
}

double gmst06(const JulianDate& ut1, const JulianDate& tt) noexcept {
    const double t = tt.centuriesSinceJ2000();
    return wrapTwoPi(earthRotationAngle(ut1) +
                     (0.014506 + t * (4612.156534 + t * (1.3915817 + t * (-0.00000044 + t * (-0.000029956 + t * -0.0000000368))))) *
                         constants::arcsecToRad);
}

double equationOfOrigins(const Mat3& r, double s) noexcept {
    const double x = r(2, 0);
    const double ax = x / (1.0 + r(2, 2));
    const double xs = 1.0 - ax * x;
    const double ys = -ax * r(2, 1);
    const double zs = -x;
    const double p = r(0, 0) * xs + r(0, 1) * ys + r(0, 2) * zs;
    const double q = r(1, 0) * xs + r(1, 1) * ys + r(1, 2) * zs;
    return (p != 0.0 || q != 0.0) ? s - std::atan2(q, p) : s;
}

double gast06(const JulianDate& ut1, const JulianDate& tt, PrecessionNutationModel model) noexcept {
    const double t = tt.centuriesSinceJ2000();
    const Mat3 npb = npbMatrix06(t, model);
    const double s = cioLocatorS06(t, cipXYFromMatrix(npb));
    return wrapTwoPi(earthRotationAngle(ut1) - equationOfOrigins(npb, s));
}

double gmst82(const JulianDate& ut1) noexcept {
    constexpr double A = 24110.54841 - constants::secondsPerDay / 2.0;
    constexpr double B = 8640184.812866, C = 0.093104, D = -6.2e-6;
    const double t = ut1.centuriesSinceJ2000();
    const double f = constants::secondsPerDay * (std::fmod(ut1.jd1, 1.0) + std::fmod(ut1.jd2, 1.0));
    return wrapTwoPi(constants::twoPi / constants::secondsPerDay * ((A + (B + (C + D * t) * t) * t) + f));
}

Mat3 polarMotionMatrix(double xp, double yp, double sp) noexcept {
    return R1(-yp) * R2(-xp) * R3(sp);
}

} // namespace astro
