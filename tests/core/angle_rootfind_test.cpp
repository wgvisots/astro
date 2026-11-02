/**********************************************************************
 * angle_rootfind_test.cpp
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

#include <gtest/gtest.h>

#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/core/rootfind.h"

using namespace astro;
using constants::pi;

TEST(Angle, Wrap) {
    EXPECT_NEAR(wrapTwoPi(-0.5), 2 * pi - 0.5, 1e-15);
    EXPECT_NEAR(wrapTwoPi(7.0), 7.0 - 2 * pi, 1e-15);
    EXPECT_GE(wrapTwoPi(-1e-300), 0.0);
    EXPECT_LT(wrapTwoPi(-1e-300), 2 * pi);
    EXPECT_NEAR(wrapPi(pi + 0.1), -pi + 0.1, 1e-15);
    EXPECT_NEAR(wrapPi(-pi - 0.1), pi - 0.1, 1e-15);
    EXPECT_DOUBLE_EQ(wrap360(-90.0), 270.0);
    EXPECT_DOUBLE_EQ(wrap180(270.0), -90.0);
}

TEST(Angle, Sexagesimal) {
    EXPECT_NEAR(parseHms("12:00:00"), pi, 1e-15);
    EXPECT_NEAR(parseHms("06h30m"), deg2rad(97.5), 1e-15);
    EXPECT_NEAR(parseDms("-00 30 00"), deg2rad(-0.5), 1e-15);
    EXPECT_NEAR(parseDms("+45d15'36\""), deg2rad(45.26), 1e-15);
    EXPECT_NEAR(parseDms("45\xC2\xB0 15' 36\""), deg2rad(45.26), 1e-15);
    EXPECT_NEAR(parseDms("12.5"), deg2rad(12.5), 1e-15);
    EXPECT_THROW((void)parseDms("12:61:00"), Error);
    EXPECT_THROW((void)parseDms("abc"), Error);
    EXPECT_THROW((void)parseDms("1:2:3:4"), Error);

    EXPECT_EQ(formatHms(parseHms("01:02:03.456")), "01:02:03.456");
    EXPECT_EQ(formatDms(deg2rad(-0.5)), "-00:30:00.00");
    EXPECT_EQ(formatDms(deg2rad(10.0)), "+10:00:00.00");
    // Rounding carries: 59.9999 s rounds up to the next minute.
    EXPECT_EQ(formatDms(dmsToRad(10, 59, 59.9999), 2), "+11:00:00.00");
    EXPECT_EQ(formatHms(hmsToRad(23, 59, 59.9996), 3), "24:00:00.000");
}

TEST(RootFind, BrentAndBisect) {
    auto         f    = [](double x) { return std::cos(x) - x; };
    const double root = 0.7390851332151607;
    EXPECT_NEAR(brentRoot(f, 0.0, 1.0, 1e-14), root, 1e-13);
    EXPECT_NEAR(bisect(f, 0.0, 1.0, 1e-14), root, 1e-13);
    EXPECT_THROW((void)brentRoot(f, 1.0, 2.0), Error);
    auto fdf = [](double x) { return std::pair{std::cos(x) - x, -std::sin(x) - 1.0}; };
    EXPECT_NEAR(newtonRoot(fdf, 1.0), root, 1e-14);
}

TEST(RootFind, BrentMinimize) {
    auto r = brentMinimize([](double x) { return (x - 1.234) * (x - 1.234) + 3.0; }, -10.0, 10.0, 1e-12);
    EXPECT_NEAR(r.x, 1.234, 5e-8); // minimizers are limited to ~sqrt(eps) in x
    EXPECT_NEAR(r.fx, 3.0, 1e-14);
}
