/**********************************************************************
 * constants_test.cpp
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

#include "astro/core/constants.h"

#include <gtest/gtest.h>

#include <cmath>
#include <string>

#include "astro/version.h"

using namespace astro;

// Defining constants against the documents that define them, written out again here so that a mistyped
// digit in constants.h fails. All are exact decimal values, so the comparison is for equality of the doubles.
TEST(Constants, DefiningValuesAsPublished) {
    EXPECT_EQ(constants::c, 299792458.0);      // SI definition of the metre
    EXPECT_EQ(constants::au, 149597870700.0);  // IAU 2012 Resolution B2
    EXPECT_EQ(constants::ttMinusTai, 32.184);  // IAU 1991 Resolution A4
    EXPECT_EQ(constants::taiMinusGps, 19.0);   // IS-GPS-200
    EXPECT_EQ(constants::taiMinusBdt, 33.0);   // BeiDou ICD: BDT = UTC on 2006-01-01, when TAI - UTC was 33 s
    EXPECT_EQ(constants::LG, 6.969290134e-10); // IAU 2000 Resolution B1.9
    EXPECT_EQ(constants::LB, 1.550519768e-8);  // IAU 2006 Resolution B3
    EXPECT_EQ(constants::TDB0, -6.55e-5);      // IAU 2006 Resolution B3
    EXPECT_EQ(constants::sunRadius, 6.957e8);  // IAU 2015 Resolution B3
    EXPECT_EQ(constants::wgs84A, 6378137.0);   // NIMA TR8350.2
    EXPECT_EQ(constants::wgs84InvF, 298.257223563);
    EXPECT_EQ(constants::grs80A, 6378137.0); // Moritz, Geodetic Reference System 1980
    EXPECT_EQ(constants::grs80InvF, 298.257222101);
    EXPECT_EQ(constants::iers2010A, 6378136.6); // IERS Conventions (2010), Table 1.1
    EXPECT_EQ(constants::iers2010InvF, 298.25642);
    EXPECT_EQ(constants::gmEarth, 3.986004418e14);
    EXPECT_EQ(constants::earthRotationRate, 7.292115e-5);
}

// Epochs, against their published Modified Julian Dates: J2000.0 is MJD 51544.5, the Unix epoch 1970-01-01 is
// MJD 40587, the GPS epoch 1980-01-06 is MJD 44244; MJD = JD - 2400000.5.
TEST(Constants, EpochsAgreeWithTheirModifiedJulianDates) {
    EXPECT_EQ(constants::mjdOffset, 2400000.5);
    EXPECT_EQ(constants::jdJ2000, 2400000.5 + 51544.5);
    EXPECT_EQ(constants::mjdJ2000, 51544.5);
    EXPECT_EQ(constants::jdUnixEpoch, 2400000.5 + 40587.0);
    EXPECT_EQ(constants::jdGpsEpoch, 2400000.5 + 44244.0);
    // T0 is 1977-01-01T00:00:00 TAI, which is JD 2443144.5 in TAI; in TT it is 32.184 s later. constants.h says
    // one double holds it to 14 microseconds: the nearest double is 1.37e-5 s off (worked out in exact rational
    // arithmetic); the spacing of doubles there is 4.7e-10 day, 40 microseconds.
    EXPECT_NEAR((constants::jdT0 - 2443144.5) * 86400.0, 32.184, 1.4e-5);
    EXPECT_EQ(constants::daysPerJulianCentury, 100.0 * 365.25);
    EXPECT_EQ(constants::daysPerJulianMillennium, 1000.0 * 365.25);
}

// The angle constants against pi in long double: one rounding error each (2.3e-16 relative at most).
TEST(Constants, AngleUnitsAgainstPi) {
    const long double piL = std::acos(-1.0L);
    const auto rel = [](long double got, long double want) { return static_cast<double>(std::abs(got - want) / want); };
    EXPECT_LE(rel(constants::pi, piL), 1.2e-16);
    EXPECT_LE(rel(constants::twoPi, 2 * piL), 1.2e-16);
    EXPECT_LE(rel(constants::halfPi, piL / 2), 1.2e-16);
    EXPECT_LE(rel(constants::degToRad, piL / 180), 2.3e-16);
    EXPECT_LE(rel(constants::radToDeg, 180 / piL), 2.3e-16);
    EXPECT_LE(rel(constants::arcsecToRad, piL / 648000), 2.3e-16);
    EXPECT_LE(rel(constants::masToRad, piL / 648000000), 3.4e-16);
    EXPECT_LE(rel(constants::uasToRad, piL / 648000000000), 3.4e-16);
}

TEST(Version, MacrosAndStringAgree) {
    // The pieces, the combined number and the text all come from CMake's project version.
    EXPECT_EQ(ASTRO_VERSION, ASTRO_VERSION_MAJOR * 10000 + ASTRO_VERSION_MINOR * 100 + ASTRO_VERSION_PATCH);
    const std::string expected = std::to_string(ASTRO_VERSION_MAJOR) + "." + std::to_string(ASTRO_VERSION_MINOR) + "."
                               + std::to_string(ASTRO_VERSION_PATCH);
    EXPECT_EQ(std::string(ASTRO_VERSION_STRING), expected);
    EXPECT_EQ(std::string(astro::versionString), expected);
}
