/**********************************************************************
 * spherical.cpp
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

#include "astro/astrometry/spherical.h"

#include <cmath>

#include "astro/core/angle.h"
#include "astro/core/rotation.h"

namespace astro {

RaDec raDecFromVector(const Vec3& v) noexcept {
    const double rho = std::hypot(v.x, v.y);
    return {rho == 0.0 ? 0.0 : wrapTwoPi(std::atan2(v.y, v.x)), std::atan2(v.z, rho), norm(v)};
}

Vec3 vectorFromRaDec(const RaDec& rd) noexcept {
    const double cd = std::cos(rd.dec);
    return {rd.distance * cd * std::cos(rd.ra), rd.distance * cd * std::sin(rd.ra), rd.distance * std::sin(rd.dec)};
}

double angularSeparation(double ra1, double dec1, double ra2, double dec2) noexcept {
    return angleBetween(unitVectorFromRaDec(ra1, dec1), unitVectorFromRaDec(ra2, dec2));
}

double positionAngle(double ra1, double dec1, double ra2, double dec2) noexcept {
    const double dra = ra2 - ra1;
    const double y   = std::sin(dra) * std::cos(dec2);
    const double x   = std::cos(dec1) * std::sin(dec2) - std::sin(dec1) * std::cos(dec2) * std::cos(dra);
    return wrapTwoPi(std::atan2(y, x));
}

Mat3 equatorialToEclipticMatrix(double eps) noexcept {
    return R1(eps);
}

RaDec equatorialToEcliptic(const Vec3& eq, double eps) noexcept {
    return raDecFromVector(R1(eps) * eq);
}

Vec3 eclipticToEquatorial(const RaDec& ecl, double eps) noexcept {
    return R1(-eps) * vectorFromRaDec(ecl);
}

const Mat3& icrsToGalacticMatrix() noexcept {
    // clang-format off
    static const Mat3 m{
        -0.054875560416215368492398900454, -0.873437090234885048760383168409, -0.483835015548713226831774175116,
        +0.494109427875583673525222371358, -0.444829629960011178146614061616, +0.746982244497218890527388004556,
        -0.867666149019004701181616534570, -0.198076373431201528180486091412, +0.455983776175066922272100478348};
    // clang-format on
    return m;
}

RaDec icrsToGalactic(double ra, double dec) noexcept {
    return raDecFromVector(icrsToGalacticMatrix() * unitVectorFromRaDec(ra, dec));
}

RaDec galacticToIcrs(double l, double b) noexcept {
    return raDecFromVector(icrsToGalacticMatrix().transpose() * unitVectorFromRaDec(l, b));
}

} // namespace astro
