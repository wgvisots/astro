/**********************************************************************
 * main.cpp
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

// Uses an installed astro: headers, library and CMake package files. Exits with 0 if a time
// conversion, an orbit propagation and an element set all work.

#include <cmath>
#include <cstdio>

#include "astro/orbit/sgp4.h"
#include "astro/orbit/tle.h"
#include "astro/orbit/two_body.h"
#include "astro/time/time.h"

int main() {
    using namespace astro;
    // J2000.0 is 2000-01-01T11:58:55.816 UTC: the TT Julian date must then be 2451545.0.
    const Time   j2000 = Time::fromIso8601("2000-01-01T11:58:55.816Z");
    const double jdTt  = j2000.to(TimeScale::TT).jd().jd1 + j2000.to(TimeScale::TT).jd().jd2;
    // A circular orbit comes back to its start after one period.
    const StateVector s0{
        {7000e3, 0.0,                                    0.0},
        {0.0,    std::sqrt(constants::gmEarth / 7000e3), 0.0}
    };
    const double period = 2.0 * constants::pi * std::sqrt(7000e3 * 7000e3 * 7000e3 / constants::gmEarth);
    const double miss   = norm(propagateKepler(s0, period).r - s0.r);
    // SGP4 on the ISS element set gives a state about 6700 km from the centre.
    const Sgp4   sgp4(Tle::parse("1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
                               "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537"));
    const double radius = norm(sgp4.propagate(0.0).r);

    std::printf("TT Julian date of J2000.0: %.9f\nmiss after one period: %.3g m\nISS radius: %.1f km\n",
                jdTt,
                miss,
                radius / 1e3);
    const bool ok = std::abs(jdTt - 2451545.0) < 1e-9 && miss < 1e-6 && radius > 6.6e6 && radius < 6.8e6;
    std::printf("%s\n", ok ? "installed astro works" : "installed astro does NOT work");
    return ok ? 0 : 1;
}
