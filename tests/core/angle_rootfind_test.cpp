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

// Every form the grammar in angle.h allows. Truth: the same angle written as a decimal number of hours or
// degrees. The parser adds the fields in the same order, so the two agree to a rounding error or two
// (EXPECT_DOUBLE_EQ allows four).
TEST(Angle, ParseAcceptsTheDocumentedSexagesimalForms) {
    const double hms = deg2rad(15.0 * (12.0 + 34.0 / 60.0 + 56.7 / 3600.0));
    for (const char* text : {"12:34:56.7",
                             "12 34 56.7",
                             "12h34m56.7s",
                             "12h 34m 56.7s",
                             "12 h 34 m 56.7 s",
                             "12H34M56.7S",
                             "12,34,56.7",
                             "12: 34: 56.7",
                             "12h34:56.7",
                             "  12:34:56.7\t",
                             "+12:34:56.7",
                             "+ 12 34 56.7",
                             "012:034:056.700"}) {
        EXPECT_DOUBLE_EQ(parseHms(text), hms) << text;
    }
    const double dms = deg2rad(-(5.0 + 6.0 / 60.0 + 7.25 / 3600.0));
    for (const char* text : {"-5:06:07.25",
                             "-05 06 07.25",
                             "- 5  6  7.25",
                             "-5d6m7.25s",
                             "-5D 6M 7.25S",
                             "-5d6'7.25\"",
                             "-5\xC2\xB0"
                             "6'7.25\"", // degree sign in UTF-8
                             "-5\xB0"
                             "6'7.25\"", // degree sign in ISO 8859-1 (Latin-1)
                             "-5\xC2\xB0 6' 7.25\" "}) {
        EXPECT_DOUBLE_EQ(parseDms(text), dms) << text;
    }
    // One or two fields; only the last may have a fraction, and it may end or begin with the point.
    EXPECT_DOUBLE_EQ(parseDms("12.5"), deg2rad(12.5));
    EXPECT_DOUBLE_EQ(parseDms("12.5d"), deg2rad(12.5));
    EXPECT_DOUBLE_EQ(parseDms("12."), deg2rad(12.0));
    EXPECT_DOUBLE_EQ(parseDms(".5"), deg2rad(0.5));
    EXPECT_DOUBLE_EQ(parseDms("12:30.5"), deg2rad(12.0 + 30.5 / 60.0));
    EXPECT_DOUBLE_EQ(parseDms("12d30'"), deg2rad(12.5));
    EXPECT_DOUBLE_EQ(parseHms("12.5h"), deg2rad(187.5));
    EXPECT_DOUBLE_EQ(parseHms("12h30"), deg2rad(187.5));
    // The first field has no upper limit; minutes and seconds just below 60 are in range.
    EXPECT_DOUBLE_EQ(parseDms("360:00:00"), deg2rad(360.0));
    EXPECT_DOUBLE_EQ(parseHms("24:00:00"), twoPi);
    EXPECT_DOUBLE_EQ(parseDms("0:59:59.999"), deg2rad(59.0 / 60.0 + 59.999 / 3600.0));
    // The sign applies to the whole angle, including an angle between 0 and -1.
    EXPECT_DOUBLE_EQ(parseDms("-0:30"), deg2rad(-0.5));
    EXPECT_TRUE(sameBits(parseDms("-0"), -0.0));
}

// Each string breaks one rule of the grammar in angle.h, and must throw as other malformed text does.
// Before the parser was made strict, most of these returned a number (or NaN, or infinity).
TEST(Angle, ParseRejectsMalformedSexagesimalText) {
    for (const char* text : {
             "nan",         "NaN",        "inf",     "-inf",        "infinity", // not numbers
             "1e2",         "1:2e1:3",    "1:2:3e1", "1E2",         "0x10",     // exponents, hexadecimal
             "5::30",       ":5:30",      "5:30:",   "5: :30",      "5,,30",    // an empty field
             "12.5:30",     "12:30.5:30", "12.5 30", "12.5.5",                  // a fraction before the last field
             "5:-30",       "--5",        "+-5",     "5 +30",       "5-",       // a sign that is not in front
             "12:34:56.7x", "12x",        "x12",     "12:34:56.7\n"             // stray characters
         }) {
        EXPECT_THROW((void)parseDms(text), Error) << "parseDms(\"" << text << "\")";
        EXPECT_THROW((void)parseHms(text), Error) << "parseHms(\"" << text << "\")";
    }
    for (const char* text : {"", " ", "+", "-", ".", "12:.", "1 2 3 4", "1:2:3:4", "12:60:00", "12:00:60"}) {
        EXPECT_THROW((void)parseDms(text), Error) << "parseDms(\"" << text << "\")";
        EXPECT_THROW((void)parseHms(text), Error) << "parseHms(\"" << text << "\")";
    }
    // A unit mark must fit the quantity: no degrees or arc marks in hours, no hours in degrees.
    for (const char* text : {"12d", "12D", "12\xC2\xB0", "12\xB0", "12h30'", "12h30m15\"", "12d30m"}) {
        EXPECT_THROW((void)parseHms(text), Error) << "parseHms(\"" << text << "\")";
    }
    for (const char* text : {"12h", "12H", "12h30m", "12h30'"}) {
        EXPECT_THROW((void)parseDms(text), Error) << "parseDms(\"" << text << "\")";
    }
    // A unit mark must fit its field: the first number is hours or degrees, the second minutes, ...
    for (const char* text : {"30m", "15s", "12h30h", "12h30s", "12m30s", "12:30:15m", "12h:30", "12hh", "h"}) {
        EXPECT_THROW((void)parseHms(text), Error) << "parseHms(\"" << text << "\")";
    }
    for (const char* text : {"30m", "30'", "15\"", "12d30d", "12d30\"", "12'30\"", "12:30:15'", "12\xC2", "d"}) {
        EXPECT_THROW((void)parseDms(text), Error) << "parseDms(\"" << text << "\")";
    }
}

// An angle that is not finite has no sexagesimal digits. The result is the one stated in angle.h; before,
// the hours or degrees were converted to int unchecked, which is undefined behaviour (it printed
// "+-2147483648:-2147483648:nan").
TEST(Angle, FormatAndDecomposeOfNonFiniteAnglesAreDefined) {
    EXPECT_EQ(formatHms(kNaN), "nan");
    EXPECT_EQ(formatDms(kNaN), "nan");
    EXPECT_EQ(formatHms(kInf), "inf");
    EXPECT_EQ(formatHms(-kInf), "-inf");
    EXPECT_EQ(formatDms(kInf), "+inf");
    EXPECT_EQ(formatDms(-kInf), "-inf");

    for (const Sexagesimal& s : {radToHms(kNaN), radToDms(kNaN)}) {
        EXPECT_FALSE(s.negative);
        EXPECT_EQ(s.whole, 0);
        EXPECT_EQ(s.minutes, 0);
        EXPECT_TRUE(std::isnan(s.seconds));
    }
    for (double sign : {1.0, -1.0}) {
        for (const Sexagesimal& s : {radToHms(sign * kInf), radToDms(sign * kInf)}) {
            EXPECT_EQ(s.negative, sign < 0.0);
            EXPECT_EQ(s.whole, 0);
            EXPECT_EQ(s.minutes, 0);
            EXPECT_EQ(s.seconds, kInf);
        }
    }
}

// The whole part is an int. The largest one still decomposes and prints; one degree more does not fit
// and is reported like infinity, as angle.h says.
TEST(Angle, FormatAndDecomposeAtTheLargestWholePart) {
    const double      largest = std::numeric_limits<int>::max(); // 2147483647, exact in a double
    const Sexagesimal fits    = radToDms(deg2rad(-largest));
    EXPECT_TRUE(fits.negative);
    EXPECT_EQ(fits.whole, std::numeric_limits<int>::max());
    EXPECT_EQ(fits.minutes, 0);
    EXPECT_EQ(fits.seconds, 0.0);
    EXPECT_EQ(formatDms(deg2rad(largest)), "+2147483647:00:00.00");

    const Sexagesimal tooLarge = radToDms(deg2rad(-(largest + 1.0)));
    EXPECT_TRUE(tooLarge.negative);
    EXPECT_EQ(tooLarge.whole, 0);
    EXPECT_EQ(tooLarge.minutes, 0);
    EXPECT_EQ(tooLarge.seconds, kInf);
    EXPECT_EQ(formatDms(deg2rad(largest + 1.0)), "+inf");
    EXPECT_EQ(formatHms(-1e300), "-inf");
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
