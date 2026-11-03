/**********************************************************************
 * time_scales.cpp
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

// Time scales: UTC with leap seconds, TAI, TT, TDB, UT1, GPS; Unix time; GNSS week and
// time of week; ISO-8601.

#include <iomanip>
#include <iostream>

#include "astro/frames/eop.h"
#include "astro/time/time.h"

using namespace astro;

int main() {
    //! [convert]
    // A UTC instant inside the leap second at the end of 2016.
    const Time utc = Time::fromIso8601("2016-12-31T23:59:60.5Z");
    std::cout << utc.iso8601(3) << '\n'                     // 2016-12-31T23:59:60.500Z
              << utc.to(TimeScale::TAI).iso8601(3) << '\n'  // 2017-01-01T00:00:36.500 TAI
              << utc.to(TimeScale::TT).iso8601(3) << '\n'   // TAI + 32.184 s
              << utc.to(TimeScale::TDB).iso8601(6) << '\n'  // TT + periodic terms (< 2 ms)
              << utc.to(TimeScale::GPS).iso8601(3) << '\n'; // TAI - 19 s
    //! [convert]

    //! [arithmetic]
    // Arithmetic is in SI seconds, so it counts the leap second.
    const Time before = Time::fromIso8601("2016-12-31T23:59:59Z");
    const Time after  = Time::fromIso8601("2017-01-01T00:00:00Z");
    std::cout << "elapsed: " << (after - before) << " s\n"; // 2 s
    std::cout << (before + 1.0).iso8601(0) << '\n';         // 2016-12-31T23:59:60Z
    //! [arithmetic]

    //! [unix_gnss]
    const Time        t   = Time::fromUnix(1700000000.0); // Unix time ignores leap seconds
    const GnssWeekTow gps = t.gnssWeekTow(GnssSystem::GPS);
    std::cout << t.iso8601(0) << "  GPS week " << gps.week << ", TOW " << gps.tow << " s\n";
    std::cout << "back to Unix: " << std::setprecision(12)
              << Time::fromGnssWeekTow(GnssSystem::GPS, gps.week, gps.tow).unixSeconds() << '\n';
    //! [unix_gnss]

    //! [ut1]
    // UT1 needs Earth orientation data (UT1 - UTC); here a constant -0.1 s instead of an IERS file.
    EopValues v;
    v.dut1 = -0.1;
    const ConstantEop eop(v);
    std::cout << t.to(TimeScale::UT1, eop).iso8601(3) << '\n';
    //! [ut1]
    return 0;
}
