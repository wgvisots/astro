/**********************************************************************
 * low_precision.cpp
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

#include "astro/ephemeris/low_precision.h"

#include <cmath>
#include <string>

#include "../core/enum_io.h"
#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/core/rotation.h"
#include "astro/frames/precession_nutation.h"

namespace astro {

std::string_view toString(Body b) noexcept {
    return detail::enumName(b, "INVALID_VALUE");
}

Vec3 sunPositionLowPrecision(const Time& time) {
    const double T  = time.to(TimeScale::TT).centuriesSinceJ2000();
    const double L0 = (280.46646 + T * (36000.76983 + T * 0.0003032)) * constants::degToRad;
    const double M  = (357.52911 + T * (35999.05029 - T * 0.0001537)) * constants::degToRad;
    const double e  = 0.016708634 - T * (0.000042037 + T * 0.0000001267);
    const double Cc = ((1.914602 - T * (0.004817 + T * 0.000014)) * std::sin(M)
                       + (0.019993 - T * 0.000101) * std::sin(2.0 * M) + 0.000289 * std::sin(3.0 * M))
                    * constants::degToRad;
    const double lon = L0 + Cc; // geometric longitude, mean equinox of date
    const double nu  = M + Cc;
    const double R   = 1.000001018 * (1.0 - e * e) / (1.0 + e * std::cos(nu)) * constants::au;
    const double eps = meanObliquity06(T);
    const Vec3   mod{R * std::cos(lon), R * std::sin(lon) * std::cos(eps), R * std::sin(lon) * std::sin(eps)};
    // The series describes the Earth-Moon barycenter's orbit; shift to the Earth's center.
    const Vec3 emb = moonPositionLowPrecision(time) / (1.0 + constants::earthMoonMassRatio);
    return biasPrecessionMatrix06(T).transpose() * mod + emb;
}

Vec3 moonPositionLowPrecision(const Time& time) {
    const double T = time.to(TimeScale::TT).centuriesSinceJ2000();
    // Mean elements (deg); L0 includes -1.3972 T to refer longitudes to the J2000 equinox.
    const double L0 = (218.31617 + 481267.88088 * T - 1.3972 * T) * constants::degToRad;
    const double l  = (134.96292 + 477198.86753 * T) * constants::degToRad;
    const double lp = (357.52543 + 35999.04944 * T) * constants::degToRad;
    const double F  = (93.27283 + 483202.01873 * T) * constants::degToRad;
    const double D  = (297.85027 + 445267.11135 * T) * constants::degToRad;
    const double as = constants::arcsecToRad;
    const double lambda =
        L0
        + as
              * (22640.0 * std::sin(l) + 769.0 * std::sin(2.0 * l) - 4586.0 * std::sin(l - 2.0 * D)
                 + 2370.0 * std::sin(2.0 * D) - 668.0 * std::sin(lp) - 412.0 * std::sin(2.0 * F)
                 - 212.0 * std::sin(2.0 * l - 2.0 * D) - 206.0 * std::sin(l + lp - 2.0 * D)
                 + 192.0 * std::sin(l + 2.0 * D) - 165.0 * std::sin(lp - 2.0 * D) + 148.0 * std::sin(l - lp)
                 - 125.0 * std::sin(D) - 110.0 * std::sin(l + lp) - 55.0 * std::sin(2.0 * F - 2.0 * D));
    const double beta =
        as
        * (18520.0 * std::sin(F + lambda - L0 + as * (412.0 * std::sin(2.0 * F) + 541.0 * std::sin(lp)))
           - 526.0 * std::sin(F - 2.0 * D) + 44.0 * std::sin(l + F - 2.0 * D) - 31.0 * std::sin(-l + F - 2.0 * D)
           - 25.0 * std::sin(-2.0 * l + F) - 23.0 * std::sin(lp + F - 2.0 * D) + 21.0 * std::sin(-l + F)
           + 11.0 * std::sin(-lp + F - 2.0 * D));
    const double r = (385000.0 - 20905.0 * std::cos(l) - 3699.0 * std::cos(2.0 * D - l) - 2956.0 * std::cos(2.0 * D)
                      - 570.0 * std::cos(2.0 * l) + 246.0 * std::cos(2.0 * l - 2.0 * D) - 205.0 * std::cos(lp - 2.0 * D)
                      - 171.0 * std::cos(l + 2.0 * D) - 152.0 * std::cos(l + lp - 2.0 * D))
                   * 1000.0;
    // Ecliptic of J2000 -> equator of J2000 (~GCRS; the 23 mas frame bias is negligible here).
    const Vec3 ecl{r * std::cos(beta) * std::cos(lambda), r * std::cos(beta) * std::sin(lambda), r * std::sin(beta)};
    return R1(-meanObliquity06(0.0)) * ecl;
}

StateVector LowPrecisionEphemeris::state(Body target, Body center, const Time& t) const {
    auto geocentric = [](Body b, const Time& tt) -> Vec3 {
        switch (b) {
        case Body::Earth: return {};
        case Body::Sun:   return sunPositionLowPrecision(tt);
        case Body::Moon:  return moonPositionLowPrecision(tt);
        default:
            throw Error("LowPrecisionEphemeris::state: supports only Sun, Moon and Earth (got Body::"
                        + std::string(toString(b)) + ")");
        }
    };
    auto             pos = [&](const Time& tt) { return geocentric(target, tt) - geocentric(center, tt); };
    constexpr double h   = 60.0;
    StateVector      s;
    s.r = pos(t);
    s.v = (pos(t + h) - pos(t - h)) / (2.0 * h);
    return s;
}

std::ostream& operator<<(std::ostream& os, Body value) {
    return detail::writeEnum(os, "Body", value);
}

} // namespace astro
