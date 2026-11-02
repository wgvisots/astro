/**********************************************************************
 * local.cpp
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

#include "astro/geodesy/local.h"

#include <cmath>

#include "astro/core/angle.h"

namespace astro {

Mat3 ecefToEnuMatrix(double lat, double lon) noexcept {
    const double sp = std::sin(lat), cp = std::cos(lat), sl = std::sin(lon), cl = std::cos(lon);
    // clang-format off
    return {     -sl,       cl, 0.0,
            -sp * cl, -sp * sl,  cp,
             cp * cl,  cp * sl,  sp};
    // clang-format on
}

Mat3 ecefToNedMatrix(double lat, double lon) noexcept {
    const double sp = std::sin(lat), cp = std::cos(lat), sl = std::sin(lon), cl = std::cos(lon);
    // clang-format off
    return {-sp * cl, -sp * sl,  cp,
                 -sl,       cl, 0.0,
            -cp * cl, -cp * sl, -sp};
    // clang-format on
}

Vec3 ecefToEnu(const Vec3& ecef, const Geodetic& o, const Ellipsoid& ell) noexcept {
    return ecefToEnuMatrix(o.lat, o.lon) * (ecef - geodeticToEcef(o, ell));
}

Vec3 enuToEcef(const Vec3& enu, const Geodetic& o, const Ellipsoid& ell) noexcept {
    return geodeticToEcef(o, ell) + ecefToEnuMatrix(o.lat, o.lon).transpose() * enu;
}

Vec3 ecefToNed(const Vec3& ecef, const Geodetic& o, const Ellipsoid& ell) noexcept {
    return ecefToNedMatrix(o.lat, o.lon) * (ecef - geodeticToEcef(o, ell));
}

Vec3 nedToEcef(const Vec3& ned, const Geodetic& o, const Ellipsoid& ell) noexcept {
    return geodeticToEcef(o, ell) + ecefToNedMatrix(o.lat, o.lon).transpose() * ned;
}

Aer enuToAer(const Vec3& enu) noexcept {
    const double horiz = std::hypot(enu.x, enu.y);
    return {wrapTwoPi(std::atan2(enu.x, enu.y)), std::atan2(enu.z, horiz), std::hypot(horiz, enu.z)};
}

Vec3 aerToEnu(const Aer& a) noexcept {
    const double ce = std::cos(a.el);
    return {a.range * ce * std::sin(a.az), a.range * ce * std::cos(a.az), a.range * std::sin(a.el)};
}

Aer ecefToAer(const Vec3& target, const Geodetic& observer, const Ellipsoid& ell) noexcept {
    return enuToAer(ecefToEnu(target, observer, ell));
}

Vec3 aerToEcef(const Aer& aer, const Geodetic& observer, const Ellipsoid& ell) noexcept {
    return enuToEcef(aerToEnu(aer), observer, ell);
}

} // namespace astro
