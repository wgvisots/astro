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

#include <cmath>
#include <limits>
#include <random>

#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/core/rootfind.h"

using namespace astro;
using constants::pi;
using constants::twoPi;

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInf = std::numeric_limits<double>::infinity();

/// True if a and b are the same double, bit for bit (so +0 and -0 differ).
bool sameBits(double a, double b) {
    return a == b && std::signbit(a) == std::signbit(b);
}

} // namespace

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

// Wrapping is the identity on its own range, so a value already in range must come back bit for bit.
TEST(Angle, WrapLeavesValuesAlreadyInRangeUnchanged) {
    const double belowPi    = std::nextafter(pi, 0.0);
    const double belowTwoPi = std::nextafter(twoPi, 0.0);
    const double below180   = std::nextafter(180.0, 0.0);
    const double below360   = std::nextafter(360.0, 0.0);
    for (double x : {0.0, -0.0, 0.1, -0.1, 1e-20, -1e-17, 5e-324, -5e-324, 3.0, -3.0, belowPi, -belowPi, -pi}) {
        EXPECT_TRUE(sameBits(wrapPi(x), x)) << "wrapPi(" << x << ") = " << wrapPi(x);
    }
    for (double x : {0.0, -0.0, 0.1, -0.1, 1e-20, -1e-17, 5e-324, -5e-324, 179.0, below180, -below180, -180.0}) {
        EXPECT_TRUE(sameBits(wrap180(x), x)) << "wrap180(" << x << ") = " << wrap180(x);
    }
    for (double x : {0.0, -0.0, 0.1, 1e-20, 5e-324, 3.0, pi, 6.0, belowTwoPi}) {
        EXPECT_TRUE(sameBits(wrapTwoPi(x), x)) << "wrapTwoPi(" << x << ") = " << wrapTwoPi(x);
    }
    for (double x : {0.0, -0.0, 0.1, 1e-20, 5e-324, 180.0, 359.0, below360}) {
        EXPECT_TRUE(sameBits(wrap360(x), x)) << "wrap360(" << x << ") = " << wrap360(x);
    }

    // Random values, each scaled by the largest double inside the range so that it is in range itself.
    std::mt19937_64                        rng(5);
    std::uniform_real_distribution<double> half(-1.0, 1.0);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    int                                    changed = 0;
    for (int i = 0; i < 100000; ++i) {
        const double h = half(rng), u = unit(rng);
        changed += wrapPi(h * belowPi) != h * belowPi;
        changed += wrap180(h * below180) != h * below180;
        changed += wrapTwoPi(u * belowTwoPi) != u * belowTwoPi;
        changed += wrap360(u * below360) != u * below360;
    }
    EXPECT_EQ(changed, 0);
}

// Out of range, the result is the input less a whole number of turns. Checked against the same sum in
// long double, where it is exact: one turn (the double twoPi, or 360) is subtracted from a double within a
// factor of two of it, so no rounding is involved and the comparison is for equality.
TEST(Angle, WrapPiAnd180SubtractWholeTurnsExactly) {
    const long double turn = twoPi;
    for (double x : {pi, pi + 0.1, 4.0, 6.0, std::nextafter(pi, 4.0)}) {
        EXPECT_EQ(wrapPi(x), static_cast<double>(x - turn)) << x;
    }
    for (double x : {-pi - 0.1, -4.0, -6.0, std::nextafter(-pi, -4.0)}) {
        EXPECT_EQ(wrapPi(x), static_cast<double>(x + turn)) << x;
    }
    for (double x : {180.0, 180.1, 270.0, 359.9, std::nextafter(180.0, 400.0)}) {
        EXPECT_EQ(wrap180(x), static_cast<double>(x - 360.0L)) << x;
    }
    for (double x : {-180.1, -270.0, -359.9, std::nextafter(-180.0, -400.0)}) {
        EXPECT_EQ(wrap180(x), static_cast<double>(x + 360.0L)) << x;
    }
    // The upper end of each range is excluded: it maps to the lower end.
    EXPECT_EQ(wrapPi(pi), -pi);
    EXPECT_EQ(wrap180(180.0), -180.0);
    EXPECT_EQ(wrapTwoPi(twoPi), 0.0);
    EXPECT_EQ(wrap360(360.0), 0.0);
    // Not-a-number and infinity have no wrapped value.
    for (double x : {kNaN, kInf, -kInf}) {
        EXPECT_TRUE(std::isnan(wrapPi(x)));
        EXPECT_TRUE(std::isnan(wrap180(x)));
        EXPECT_TRUE(std::isnan(wrapTwoPi(x)));
        EXPECT_TRUE(std::isnan(wrap360(x)));
    }
}

// A negative first argument makes the whole angle negative, as "-5 deg 30 min" does on paper. Truth: the
// same angle written as a decimal, -(5 + 30/60) = -5.5 degrees; both sides round the same way, so the
// comparison is for equality.
TEST(Angle, DmsAndHmsToRadTakeTheSignOfTheFirstArgument) {
    EXPECT_EQ(dmsToRad(-5, 30, 0.0), deg2rad(-5.5));
    EXPECT_EQ(hmsToRad(-1, 30, 0.0), deg2rad(-22.5));
    // `negative` is still how to write an angle between 0 and -1: an int has no negative zero.
    EXPECT_EQ(dmsToRad(0, 30, 0.0, true), deg2rad(-0.5));
    EXPECT_EQ(hmsToRad(0, 30, 0.0, true), deg2rad(-7.5));
    // The two ways of giving the sign agree; they do not cancel.
    EXPECT_EQ(dmsToRad(-5, 30, 0.0, true), deg2rad(-5.5));
    EXPECT_EQ(hmsToRad(-1, 30, 0.0, true), deg2rad(-22.5));
    // Positive angles, and the flag alone, are as before.
    EXPECT_EQ(dmsToRad(5, 30, 0.0), deg2rad(5.5));
    EXPECT_EQ(dmsToRad(5, 30, 0.0, true), deg2rad(-5.5));
    // Minutes and seconds are magnitudes: their signs do not change the sign of the angle.
    EXPECT_EQ(dmsToRad(5, -30, -0.0), deg2rad(5.5));
    EXPECT_EQ(dmsToRad(-5, -30, -0.0), deg2rad(-5.5));
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
