/**********************************************************************
 * tle_test.cpp
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

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "astro/core/angle.h"
#include "astro/core/constants.h"
#include "astro/core/error.h"
#include "astro/orbit/tle.h"
#include "test_util.h"

using namespace astro;

namespace {

// The element set used as the example in CelesTrak's description of the format.
const std::string kLine1 = "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927";
const std::string kLine2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";

// `line` cut or padded to 68 columns, with its checksum added in column 69: the sum of the digits, each minus sign
// counting 1, modulo 10. Written out here so that the tests do not lean on astro::tleChecksum.
std::string withChecksum(std::string line) {
    line.resize(68, ' ');
    int sum = 0;
    for (const char c : line) {
        if (c >= '0' && c <= '9') {
            sum += c - '0';
        } else if (c == '-') {
            sum += 1;
        }
    }
    line += static_cast<char>('0' + sum % 10);
    return line;
}

// `line` with `text` written from column `column` (counted from 1, as in the format description), and the
// checksum recomputed.
std::string withField(std::string line, std::size_t column, const std::string& text) {
    line.replace(column - 1, text.size(), text);
    return withChecksum(line);
}

// Columns `first` to `last` of `line` (counted from 1, inclusive).
std::string columns(const std::string& line, std::size_t first, std::size_t last) {
    return line.substr(first - 1, last - first + 1);
}

// Expects `f` to throw an astro::Error whose message contains `name`: the field, line or function at fault.
template <class F>
void expectErrorNaming(F&& f, const std::string& name) {
    try {
        f();
        ADD_FAILURE() << "no astro::Error was thrown; expected one naming '" << name << "'";
    } catch (const Error& e) {
        EXPECT_NE(std::string(e.what()).find(name), std::string::npos)
            << "expected the message to name '" << name << "'; it was: " << e.what();
    }
}

const double kNan = std::numeric_limits<double>::quiet_NaN();
const double kInf = std::numeric_limits<double>::infinity();

} // namespace

// --- The epoch ------------------------------------------------------------------------------------------------------

TEST(Tle, EpochDayFractionIsOfAn86400SecondDayAlsoOnALeapSecondDay) {
    // The epoch field holds the day of the year and a fraction of an 86400 s day counted from 0h UTC. That is how
    // element sets are produced and how the SGP4 reference code reads them (twoline2rv in Vallado et al. 2006
    // turns the fraction into hours, minutes and seconds of a 24 h day). 2016-12-31 (day 366 of 2016) and
    // 2015-06-30 (day 181 of 2015) end with a leap second and are 86401 s long.
    struct Case {
        const char* field;
        int         year, month, day, hour, minute;
        double      second;
    };
    // clang-format off
    const Case cases[] = {
        {"16366.50000000", 2016, 12, 31, 12,  0,  0.0     }, // 0.5        x 86400 s = 43200 s
        {"16366.75000000", 2016, 12, 31, 18,  0,  0.0     }, // 0.75       x 86400 s = 64800 s
        {"16366.99998843", 2016, 12, 31, 23, 59, 59.000352}, // 0.99998843 x 86400 s = 86399.000352 s
        {"15181.50000000", 2015,  6, 30, 12,  0,  0.0     },
        {"16365.75000000", 2016, 12, 30, 18,  0,  0.0     }, // ordinary days
        {"17001.00000000", 2017,  1,  1,  0,  0,  0.0     },
        {"08264.51782528", 2008,  9, 20, 12, 25, 40.104192}, // 0.51782528 x 86400 s = 44740.104192 s
    };
    // clang-format on
    for (const Case& c : cases) {
        const std::string line1 = withField(kLine1, 19, c.field);
        const Tle         tle   = Tle::parse(line1, kLine2);
        const Time expected = Time::fromCalendar(TimeScale::UTC, c.year, c.month, c.day, c.hour, c.minute, c.second);
        // 0.1 microsecond: the expected seconds above are exact products, and a two-part Julian date carries the
        // time of day to about 1e-11 s. The defect was up to 1 s.
        EXPECT_NEAR(tle.epoch - expected, 0.0, 1e-7) << c.field << " was read as " << tle.epoch.iso8601(6);
        // Writing is the inverse of reading.
        EXPECT_EQ(tle.line1(), line1);
    }
}

TEST(Tle, Line1WritesTheEpochAsAFractionOfAn86400SecondDay) {
    Tle tle = Tle::parse(kLine1, kLine2);
    // 18:00:00 is 64800 s after 0h, 0.75 of 86400 s, whether or not the day ends with a leap second.
    tle.epoch = Time::fromCalendar(TimeScale::UTC, 2016, 12, 31, 18, 0, 0.0); // leap-second day
    EXPECT_EQ(columns(tle.line1(), 19, 32), "16366.75000000");
    tle.epoch = Time::fromCalendar(TimeScale::UTC, 2016, 12, 30, 18, 0, 0.0); // ordinary day
    EXPECT_EQ(columns(tle.line1(), 19, 32), "16365.75000000");
    // An epoch held in another time scale is the same instant: 2017-01-01 00:00:37 TAI is 2016-12-31 23:59:60
    // UTC (TAI - UTC was 36 s until then), so 18:00:36 TAI on the 31st is 18:00:00 UTC.
    tle.epoch = Time::fromCalendar(TimeScale::TAI, 2016, 12, 31, 18, 0, 36.0);
    EXPECT_EQ(columns(tle.line1(), 19, 32), "16366.75000000");
}

TEST(Tle, ParseRejectsAnEpochDayThatIsNotInTheYear) {
    // Days are numbered from 1; 2008 is a leap year (366 days), 2007 is not.
    expectErrorNaming([] { (void)Tle::parse(withField(kLine1, 19, "08000.50000000"), kLine2); }, "epoch day");
    expectErrorNaming([] { (void)Tle::parse(withField(kLine1, 19, "08367.00000000"), kLine2); }, "epoch day");
    expectErrorNaming([] { (void)Tle::parse(withField(kLine1, 19, "07366.00000000"), kLine2); }, "epoch day");
    expectErrorNaming([] { (void)Tle::parse(withField(kLine1, 19, "08964.51782528"), kLine2); }, "epoch day");
    const Tle last = Tle::parse(withField(kLine1, 19, "08366.99999999"), kLine2);
    EXPECT_EQ(last.epoch.iso8601(3), "2008-12-31T23:59:59.999Z"); // 0.99999999 x 86400 s = 86399.999136 s
    const Tle first = Tle::parse(withField(kLine1, 19, "07001.00000000"), kLine2);
    EXPECT_EQ(first.epoch.iso8601(0), "2007-01-01T00:00:00Z");
}

// --- Writing B* and the second derivative of the mean motion --------------------------------------------------------

TEST(Tle, Line1WritesExponentFieldsAtTheNearestValueTheFieldCanHold) {
    // The field is a sign, five mantissa digits with "0." implied in front, and a signed one-digit exponent:
    // +-0.NNNNN x 10^(+-N). Below 0.1e-9 the exponent stays at -9 and the mantissa gets leading zeros; the
    // smallest step is 0.00001e-9 = 1e-14, and anything nearer to zero than half of that is written as zero.
    struct Case {
        double      value;
        const char* field;
        double      readBack;
    };
    // clang-format off
    const Case cases[] = {
        { 0.0,         " 00000-0",  0.0       },
        { 1e-5,        " 10000-4",  1e-5      },
        {-0.11606e-4,  "-11606-4", -0.11606e-4},
        { 0.999996e-3, " 10000-2",  1e-3      }, // the rounding carries into the exponent
        { 0.5,         " 50000+0",  0.5       },
        { 12.0,        " 12000+2",  12.0      },
        { 0.99999e9,   " 99999+9",  0.99999e9 }, // the largest value
        { 1e-10,       " 10000-9",  1e-10     }, // the smallest value with a full mantissa
        { 1e-11,       " 01000-9",  1e-11     },
        { 3.2e-12,     " 00320-9",  3.2e-12   },
        {-4.5e-11,     "-04500-9", -4.5e-11   },
        { 1e-13,       " 00010-9",  1e-13     },
        { 1.26e-13,    " 00013-9",  1.3e-13   }, // rounded to the last digit that is left
        { 1e-14,       " 00001-9",  1e-14     }, // the smallest value
        { 4e-15,       " 00000-0",  0.0       },
        {-4e-15,       " 00000-0",  0.0       },
        { 1e-320,      " 00000-0",  0.0       },
    };
    // clang-format on
    for (const Case& c : cases) {
        Tle tle                 = Tle::parse(kLine1, kLine2);
        tle.bstar               = c.value;
        tle.meanMotionDdot      = c.value;
        const std::string line1 = tle.line1();
        EXPECT_EQ(columns(line1, 54, 61), c.field) << "B* = " << c.value;
        EXPECT_EQ(columns(line1, 45, 52), c.field) << "n-ddot = " << c.value;
        // Reading the line back gives the value of the field, to the rounding of the two multiplications that
        // form it (a relative 1e-15).
        const Tle back = Tle::parse(line1, kLine2);
        EXPECT_NEAR(back.bstar, c.readBack, 1e-14 * std::abs(c.readBack)) << "B* = " << c.value;
        EXPECT_NEAR(back.meanMotionDdot, c.readBack, 1e-14 * std::abs(c.readBack)) << "n-ddot = " << c.value;
    }
}

TEST(Tle, ExponentFieldsReadBackWithinHalfAStepOverTheirWholeRange) {
    // "Nearest value" as a property: whatever is written, reading it back is within half a step of the last
    // mantissa digit. For a value d.dddd x 10^p that step is 10^(p-4), and never less than 1e-14.
    for (int power = -16; power <= 7; ++power) {
        for (const double digits : {1.0, 1.00004, 1.00005001, 1.23456789, 5.55555, 9.99994, 9.99996}) {
            for (const double sign : {1.0, -1.0}) {
                const double value = sign * digits * std::pow(10.0, power);
                const double step  = std::max(1e-14, std::pow(10.0, power - 4));
                Tle          tle   = Tle::parse(kLine1, kLine2);
                tle.bstar          = value;
                tle.meanMotionDdot = value;
                const Tle back     = Tle::parse(tle.line1(), kLine2);
                // The 1e-9 allows for the rounding of the powers of ten here and in the reader.
                EXPECT_LE(std::abs(back.bstar - value), 0.5 * step * (1.0 + 1e-9)) << "B* = " << value;
                EXPECT_LE(std::abs(back.meanMotionDdot - value), 0.5 * step * (1.0 + 1e-9)) << "n-ddot = " << value;
            }
        }
    }
}

TEST(Tle, Line1RefusesExponentFieldValuesTooLargeForTheField) {
    // The largest value the field holds is 0.99999e9; from 0.999995e9 a value rounds beyond it.
    for (const double value : {1e9, -1e9, 0.999996e9, 2.5e12, kInf, -kInf, kNan}) {
        Tle tle   = Tle::parse(kLine1, kLine2);
        tle.bstar = value;
        expectErrorNaming([&] { (void)tle.line1(); }, "B*");
        tle                = Tle::parse(kLine1, kLine2);
        tle.meanMotionDdot = value;
        expectErrorNaming([&] { (void)tle.line1(); }, "n-ddot");
    }
}

// --- Reading number fields ------------------------------------------------------------------------------------------

TEST(Tle, BlanksInImpliedDecimalFieldsCountAsZeros) {
    // The decimal point of these fields is implied before the first column, so a blank column is a zero digit
    // and not padding: Vallado's reference reader twoline2rv fills them with zeros before it reads the number.
    struct Ecc {
        const char* field;
        double      value;
    };
    for (const Ecc& c : {
             Ecc{"0006703", 0.0006703},
             Ecc{" 006703", 0.0006703},
             Ecc{"   6703", 0.0006703},
             Ecc{"6703   ", 0.6703   },
             Ecc{"      1", 0.0000001},
             Ecc{"       ", 0.0      }
    }) {
        const Tle tle = Tle::parse(kLine1, withField(kLine2, 27, c.field));
        // 1e-16: the value is one correctly rounded decimal-to-binary conversion.
        EXPECT_NEAR(tle.eccentricity, c.value, 1e-16) << "eccentricity field '" << c.field << "'";
    }
    // B* (columns 54-61) and the second derivative of the mean motion (columns 45-52): 0.NNNNN x 10^N.
    struct Exp {
        const char* field;
        double      value;
    };
    for (const Exp& c : {
             Exp{" 11606-4", 0.11606e-4 },
             Exp{"  1606-4", 0.01606e-4 },
             Exp{"-  606-4", -0.00606e-4},
             Exp{" 116  -4", 0.11600e-4 },
             Exp{" 11606- ", 0.11606    },
             Exp{"        ", 0.0        }
    }) {
        const Tle tle = Tle::parse(withField(withField(kLine1, 54, c.field), 45, c.field), kLine2);
        // A relative 1e-15: the mantissa times a power of ten, each rounded once.
        EXPECT_NEAR(tle.bstar, c.value, 1e-14 * std::abs(c.value)) << "B* field '" << c.field << "'";
        EXPECT_NEAR(tle.meanMotionDdot, c.value, 1e-14 * std::abs(c.value)) << "n-ddot field '" << c.field << "'";
    }
}

TEST(Tle, ParseRejectsNumberFieldsThatAreNotPlainDecimalNumbers) {
    // A field holds digits, at most one sign and one decimal point, and blanks. Each of these was accepted
    // because the C library's number parsing also takes "nan", "inf", hexadecimal and exponents.
    struct Case {
        int         line;
        std::size_t column;
        const char* text;
        const char* name;
    };
    // clang-format off
    const Case cases[] = {
        {2,  9, "     nan",     "inclination"        },
        {2,  9, "     inf",     "inclination"        },
        {2,  9, "  0x1p+5",     "inclination"        },
        {2,  9, "   5e+01",     "inclination"        },
        {2, 18, "    -nan",     "RAAN"               },
        {2, 27, "1e-3   ",      "eccentricity"       },
        {2, 35, "infinity",     "argument of perigee"},
        {2, 44, "    1e2 ",     "mean anomaly"       },
        {2, 53, "        inf",  "mean motion"        },
        {2, 53, " 0x10.0p0  ",  "mean motion"        },
        {2, 64, "  -84",        "revolution number"  },
        {1, 19, "-8",           "epoch year"         },
        {1, 21, "       1e2  ", "epoch day"          },
        {1, 21, "         nan", "epoch day"          },
        {1, 34, "  1e-5    ",   "n-dot"              },
        {1, 45, " 11e2-5 ",     "n-ddot"             },
        {1, 54, " 11e2-4 ",     "B*"                 },
        {1, 65, "  -3",         "element set number" },
        {1, 65, "  +3",         "element set number" },
    };
    // clang-format on
    for (const Case& c : cases) {
        const std::string line1 = c.line == 1 ? withField(kLine1, c.column, c.text) : kLine1;
        const std::string line2 = c.line == 2 ? withField(kLine2, c.column, c.text) : kLine2;
        expectErrorNaming([&] { (void)Tle::parse(line1, line2); }, c.name);
    }
}

TEST(Tle, ParseAcceptsSignsDecimalPointsAndBlankPadding) {
    // What a fixed-format field may hold: blank padding, one sign, one decimal point. A blank field is zero,
    // as in a Fortran formatted read (the STR#3 test element sets leave the ephemeris type blank).
    const std::string line1 = withField(withField(withField(kLine1, 34, "+.00002182"), 63, " "), 65, "    ");
    const std::string line2 = withField(withField(kLine2, 9, "  51.64 "), 44, "   5.   ");
    const Tle         tle   = Tle::parse(line1, line2);
    EXPECT_EQ(tle.meanMotionDot, 0.00002182);
    EXPECT_EQ(tle.ephemerisType, 0);
    EXPECT_EQ(tle.elementSetNumber, 0);
    EXPECT_EQ(tle.inclination, deg2rad(51.64));
    EXPECT_EQ(tle.meanAnomaly, deg2rad(5.0));
    // A blank-padded day of the year, as some writers produce.
    EXPECT_EQ(Tle::parse(withField(kLine1, 19, "08  4.51782528"), kLine2).epoch.iso8601(3), "2008-01-04T12:25:40.104Z");
}

// --- Writing values at and beyond the edges of their fields ---------------------------------------------------------

TEST(Tle, Line2WritesAnglesThatRoundTo360AsZero) {
    // 2 pi - 1e-9 rad is 359.99999994 deg, which rounds to 360.0000 at the four decimals of the field. The
    // range of these angles is [0, 360), so that is 0.0000.
    const double almostTwoPi = constants::twoPi - 1e-9;
    Tle          tle         = Tle::parse(kLine1, kLine2);
    tle.raan                 = almostTwoPi;
    tle.argPerigee           = almostTwoPi;
    tle.meanAnomaly          = almostTwoPi;
    const std::string line2  = tle.line2();
    EXPECT_EQ(columns(line2, 18, 25), "  0.0000");
    EXPECT_EQ(columns(line2, 35, 42), "  0.0000");
    EXPECT_EQ(columns(line2, 44, 51), "  0.0000");
    // The largest angle that does not round up, and angles outside [0, 2 pi) that wrap into it.
    tle.raan        = deg2rad(359.99994);
    tle.argPerigee  = deg2rad(-90.0);
    tle.meanAnomaly = deg2rad(725.5);
    EXPECT_EQ(columns(tle.line2(), 18, 25), "359.9999");
    EXPECT_EQ(columns(tle.line2(), 35, 42), "270.0000");
    EXPECT_EQ(columns(tle.line2(), 44, 51), "  5.5000");
}

TEST(Tle, Line1CarriesAnEpochThatRoundsUpToMidnightIntoTheNextDay) {
    // The field resolves 1e-8 day (0.864 ms). 0.1 ms before midnight is 0.9999999988 of the day and rounds to
    // the next day's 0h: the day number, and at the end of December the year, must follow.
    Tle tle   = Tle::parse(kLine1, kLine2);
    tle.epoch = Time::fromCalendar(TimeScale::UTC, 2020, 12, 31, 23, 59, 59.9999); // day 366 of a leap year
    EXPECT_EQ(columns(tle.line1(), 19, 32), "21001.00000000");
    tle.epoch = Time::fromCalendar(TimeScale::UTC, 2021, 3, 4, 23, 59, 59.9999); // day 63
    EXPECT_EQ(columns(tle.line1(), 19, 32), "21064.00000000");
    tle.epoch = Time::fromCalendar(TimeScale::UTC, 1999, 12, 31, 23, 59, 59.9999);
    EXPECT_EQ(columns(tle.line1(), 19, 32), "00001.00000000");
    // 1 ms before midnight is 0.99999998843 of the day and stays in its day.
    tle.epoch = Time::fromCalendar(TimeScale::UTC, 2020, 12, 31, 23, 59, 59.999);
    EXPECT_EQ(columns(tle.line1(), 19, 32), "20366.99999999");
}

TEST(Tle, Line1DoesNotCarryAnEpochAcrossALeapSecond) {
    // On 2016-12-31 the next day's 0h is a whole second after 23:59:60, so 0.1 ms before 23:59:60 must not be
    // written as 17001.00000000, which would be 1.0001 s late. The nearest epoch the field can hold is
    // 16366.99999999 (23:59:59.999136, 0.76 ms early).
    Tle tle   = Tle::parse(kLine1, kLine2);
    tle.epoch = Time::fromCalendar(TimeScale::UTC, 2016, 12, 31, 23, 59, 59.9999);
    EXPECT_EQ(columns(tle.line1(), 19, 32), "16366.99999999");
    // Inside the leap second there is no fraction of an 86400 s day that names the instant.
    tle.epoch = Time::fromCalendar(TimeScale::UTC, 2016, 12, 31, 23, 59, 60.5);
    expectErrorNaming([&] { (void)tle.line1(); }, "epoch");
    // 0.1 ms before the end of the leap second rounds to the next day's 0h, as on any other day.
    tle.epoch = Time::fromCalendar(TimeScale::UTC, 2016, 12, 31, 23, 59, 60.9999);
    EXPECT_EQ(columns(tle.line1(), 19, 32), "17001.00000000");
}

TEST(Tle, Line1RefusesAnEpochYearTheTwoDigitFieldCannotHold) {
    // Two-digit years 57-99 are 1957-1999 and 00-56 are 2000-2056, so 2057 would be read back as 1957.
    Tle tle   = Tle::parse(kLine1, kLine2);
    tle.epoch = Time::fromCalendar(TimeScale::UTC, 2057, 3, 1);
    expectErrorNaming([&] { (void)tle.line1(); }, "epoch");
    tle.epoch = Time::fromCalendar(TimeScale::UTC, 2056, 12, 31, 12);
    EXPECT_EQ(columns(tle.line1(), 19, 32), "56366.50000000");
    // A non-finite epoch has no date at all.
    tle.epoch = Time(JulianDate{kNan, 0.0}, TimeScale::UTC);
    expectErrorNaming([&] { (void)tle.line1(); }, "epoch");
}

TEST(Tle, EncodeAlpha5RefusesNumbersOutsideItsRange) {
    // Five characters hold 0 to 99999 as digits and 100000 to 339999 with a leading letter.
    expectErrorNaming([] { (void)encodeAlpha5(-5); }, "catalog number");
    expectErrorNaming([] { (void)encodeAlpha5(-1); }, "catalog number");
    expectErrorNaming([] { (void)encodeAlpha5(340000); }, "catalog number");
    EXPECT_EQ(encodeAlpha5(0), "00000");
    EXPECT_EQ(encodeAlpha5(99999), "99999");
    EXPECT_EQ(encodeAlpha5(100000), "A0000");
    EXPECT_EQ(encodeAlpha5(339999), "Z9999");
}

TEST(Tle, Line1RefusesValuesItsFieldsCannotHold) {
    const Tle good = Tle::parse(kLine1, kLine2);
    struct Case {
        const char* name;
        void (*set)(Tle&);
    };
    // clang-format off
    const Case cases[] = {
        {"catalog number",           [](Tle& t) { t.satnum = -7; }                           },
        {"catalog number",           [](Tle& t) { t.satnum = 340000; }                       },
        {"classification",           [](Tle& t) { t.classification = '\0'; }                 },
        {"classification",           [](Tle& t) { t.classification = '\n'; }                 },
        {"international designator", [](Tle& t) { t.intlDesignator = "1998-067ABCDEF"; }     }, // 14 characters for 8 columns
        {"international designator", [](Tle& t) { t.intlDesignator = "98067\tA"; }           },
        {"n-dot",                    [](Tle& t) { t.meanMotionDot = -1.5; }                  }, // the field has no digit before the point
        {"n-dot",                    [](Tle& t) { t.meanMotionDot = 1.0; }                   },
        {"n-dot",                    [](Tle& t) { t.meanMotionDot = 0.999999996; }           }, // rounds to 1.00000000
        {"n-dot",                    [](Tle& t) { t.meanMotionDot = kNan; }                  },
        {"n-dot",                    [](Tle& t) { t.meanMotionDot = kInf; }                  },
        {"ephemeris type",           [](Tle& t) { t.ephemerisType = 12; }                    },
        {"ephemeris type",           [](Tle& t) { t.ephemerisType = -1; }                    },
        {"element set number",       [](Tle& t) { t.elementSetNumber = -3; }                 },
    };
    // clang-format on
    for (const Case& c : cases) {
        Tle tle = good;
        c.set(tle);
        expectErrorNaming([&] { (void)tle.line1(); }, c.name);
    }
}

TEST(Tle, Line2RefusesValuesItsFieldsCannotHold) {
    const Tle good = Tle::parse(kLine1, kLine2);
    struct Case {
        const char* name;
        void (*set)(Tle&);
    };
    // clang-format off
    const Case cases[] = {
        {"catalog number",      [](Tle& t) { t.satnum = -7; }                  },
        {"inclination",         [](Tle& t) { t.inclination = kNan; }           },
        {"inclination",         [](Tle& t) { t.inclination = deg2rad(1000.0); }}, // four digits before the point
        {"RAAN",                [](Tle& t) { t.raan = kNan; }                  },
        {"RAAN",                [](Tle& t) { t.raan = kInf; }                  },
        {"eccentricity",        [](Tle& t) { t.eccentricity = 0.99999996; }    }, // rounds to 1.0000000
        {"eccentricity",        [](Tle& t) { t.eccentricity = 1.2; }           },
        {"eccentricity",        [](Tle& t) { t.eccentricity = -0.001; }        },
        {"eccentricity",        [](Tle& t) { t.eccentricity = kNan; }          },
        {"argument of perigee", [](Tle& t) { t.argPerigee = -kInf; }           },
        {"mean anomaly",        [](Tle& t) { t.meanAnomaly = kNan; }           },
        {"mean motion",         [](Tle& t) { t.meanMotion = 100.5; }           }, // three digits before the point
        {"mean motion",         [](Tle& t) { t.meanMotion = 99.999999996; }    }, // rounds to 100.00000000
        {"mean motion",         [](Tle& t) { t.meanMotion = kNan; }            },
        {"revolution number",   [](Tle& t) { t.revolutionNumber = -8; }        },
    };
    // clang-format on
    for (const Case& c : cases) {
        Tle tle = good;
        c.set(tle);
        expectErrorNaming([&] { (void)tle.line2(); }, c.name);
    }
}

TEST(Tle, ValuesAtTheEdgesOfTheirFieldsAreWrittenAndReadBack) {
    // The largest (or otherwise extreme) value each field can hold: the lines must be well formed, 69 columns
    // with a valid checksum, and give the values back.
    Tle tle;
    tle.satnum              = 339999;
    tle.classification      = 'S';
    tle.intlDesignator      = "99999ZZZ";
    tle.epoch               = Time::fromCalendar(TimeScale::UTC, 1960, 1, 1);
    tle.meanMotionDot       = -0.99999999;
    tle.meanMotionDdot      = -0.99999e9;
    tle.bstar               = 0.99999e9;
    tle.ephemerisType       = 9;
    tle.elementSetNumber    = 9999;
    tle.inclination         = deg2rad(180.0);
    tle.raan                = deg2rad(359.9999);
    tle.eccentricity        = 0.9999999;
    tle.argPerigee          = deg2rad(0.0001);
    tle.meanAnomaly         = deg2rad(359.9999);
    tle.meanMotion          = 99.99999999;
    tle.revolutionNumber    = 99999;
    const std::string line1 = tle.line1(), line2 = tle.line2();
    EXPECT_EQ(line1, withChecksum("1 Z9999S 99999ZZZ 60001.00000000 -.99999999 -99999+9  99999+9 9 9999"));
    EXPECT_EQ(line2, withChecksum("2 Z9999 180.0000 359.9999 9999999   0.0001 359.9999 99.9999999999999"));
    const Tle back = Tle::parse(line1, line2);
    EXPECT_EQ(back.satnum, tle.satnum);
    EXPECT_EQ(back.classification, tle.classification);
    EXPECT_EQ(back.intlDesignator, tle.intlDesignator);
    EXPECT_TRUE(back.epoch == tle.epoch) << back.epoch.iso8601(9);
    EXPECT_EQ(back.meanMotionDot, tle.meanMotionDot);
    EXPECT_EQ(back.ephemerisType, tle.ephemerisType);
    EXPECT_EQ(back.elementSetNumber, tle.elementSetNumber);
    EXPECT_EQ(back.eccentricity, tle.eccentricity);
    EXPECT_EQ(back.meanMotion, tle.meanMotion);
    EXPECT_EQ(back.revolutionNumber, tle.revolutionNumber);
    // 1e-15: a degree-radian-degree round trip of a four-decimal angle, and a mantissa times a power of ten.
    EXPECT_NEAR(back.meanMotionDdot, tle.meanMotionDdot, 1e-15 * 1e9);
    EXPECT_NEAR(back.bstar, tle.bstar, 1e-15 * 1e9);
    EXPECT_NEAR(back.inclination, tle.inclination, 1e-15);
    EXPECT_NEAR(back.raan, tle.raan, 1e-14);
    EXPECT_NEAR(back.argPerigee, tle.argPerigee, 1e-15);
    EXPECT_NEAR(back.meanAnomaly, tle.meanAnomaly, 1e-14);
    // The element set and revolution numbers wrap, as documented, when they outgrow their fields.
    tle.elementSetNumber = 12345;
    tle.revolutionNumber = 123456;
    EXPECT_EQ(columns(tle.line1(), 65, 68), "2345");
    EXPECT_EQ(columns(tle.line2(), 64, 68), "23456");
}

// --- Checksums, catalog numbers and files of element sets -----------------------------------------------------------

TEST(Tle, ChecksumVerificationNeedsAll69ColumnsAndADigit) {
    // Without column 69 there is nothing to verify, which is not the same as a checksum that matches.
    expectErrorNaming([] { (void)Tle::parse(kLine1.substr(0, 68), kLine2); }, "line 1");
    expectErrorNaming([] { (void)Tle::parse(kLine1, kLine2.substr(0, 68)); }, "line 2");
    std::string notADigit = kLine1;
    notADigit[68]         = 'X';
    expectErrorNaming([&] { (void)Tle::parse(notADigit, kLine2); }, "checksum");
    notADigit     = kLine2;
    notADigit[68] = '-';
    expectErrorNaming([&] { (void)Tle::parse(kLine1, notADigit); }, "checksum");
    // Without verification the lines may stop short of column 69, as before; a field that is missing is zero.
    const Tle tle = Tle::parse(kLine1.substr(0, 68), kLine2.substr(0, 63), "", false);
    EXPECT_EQ(tle.satnum, 25544);
    EXPECT_EQ(tle.elementSetNumber, 292);
    EXPECT_EQ(tle.meanMotion, 15.72125391);
    EXPECT_EQ(tle.revolutionNumber, 0);
}

TEST(Tle, EpochBeforeTheStartOfUtcIsATleError) {
    // Two-digit years 57-59 are 1957-1959, before UTC and the leap-second table begin (1960-01-01).
    for (const char* year : {"57", "58", "59"}) {
        try {
            (void)Tle::parse(withField(kLine1, 19, year), kLine2);
            ADD_FAILURE() << "epoch year " << year << " was accepted";
        } catch (const Error& e) {
            const std::string what = e.what();
            EXPECT_EQ(what.rfind("TLE: epoch", 0), 0u) << what;
            EXPECT_NE(what.find("19" + std::string(year)), std::string::npos) << what;
        }
    }
    // The pivot of the two-digit year: 57-99 are 1957-1999 and 00-56 are 2000-2056.
    EXPECT_EQ(Tle::parse(withField(kLine1, 19, "60"), kLine2).epoch.iso8601(0), "1960-09-20T12:25:40Z");
    EXPECT_EQ(Tle::parse(withField(kLine1, 19, "99"), kLine2).epoch.iso8601(0), "1999-09-21T12:25:40Z");
    EXPECT_EQ(Tle::parse(withField(kLine1, 19, "00"), kLine2).epoch.iso8601(0), "2000-09-20T12:25:40Z");
    EXPECT_EQ(Tle::parse(withField(kLine1, 19, "56"), kLine2).epoch.iso8601(0), "2056-09-20T12:25:40Z");
}

TEST(Tle, DecodeAlpha5AcceptsOnlyTheAlpha5Alphabet) {
    // Five characters: all digits, or a capital letter other than I and O followed by four digits.
    for (const char* bad : {"-1234",
                            "+1234",
                            "a0001",
                            "1234567",
                            "A12345",
                            "A-123",
                            "A 123",
                            "A12",
                            "1234",
                            "    5",
                            "12 34",
                            "I0001",
                            "O0001",
                            "0000A",
                            "A000B",
                            "     ",
                            ""}) {
        expectErrorNaming([&] { (void)decodeAlpha5(bad); }, "catalog number");
    }
    // The letters count on from 10, skipping I and O: A = 10 ... H = 17, J = 18 ... N = 22, P = 23 ... Z = 33.
    EXPECT_EQ(decodeAlpha5("00000"), 0);
    EXPECT_EQ(decodeAlpha5("00005"), 5);
    EXPECT_EQ(decodeAlpha5("99999"), 99999);
    EXPECT_EQ(decodeAlpha5("A0000"), 100000);
    EXPECT_EQ(decodeAlpha5("H9999"), 179999);
    EXPECT_EQ(decodeAlpha5("J0000"), 180000);
    EXPECT_EQ(decodeAlpha5("N9999"), 229999);
    EXPECT_EQ(decodeAlpha5("P0000"), 230000);
    EXPECT_EQ(decodeAlpha5("Z9999"), 339999);
}

TEST(Tle, ParseReadsAlpha5AndBlankPaddedCatalogNumbers) {
    // Alpha-5 through parse and back: B = 11, so "B1234" is 111234.
    const std::string line1 = withField(kLine1, 3, "B1234"), line2 = withField(kLine2, 3, "B1234");
    const Tle         tle = Tle::parse(line1, line2);
    EXPECT_EQ(tle.satnum, 111234);
    EXPECT_EQ(tle.line1(), line1);
    EXPECT_EQ(tle.line2(), line2);
    // The number is right-justified in its five columns; some writers pad it with blanks instead of zeros.
    EXPECT_EQ(Tle::parse(withField(kLine1, 3, "    5"), withField(kLine2, 3, "    5")).satnum, 5);
    // Anything else in the field is refused, and so are two lines for different objects.
    expectErrorNaming([] { (void)Tle::parse(withField(kLine1, 3, "b1234"), withField(kLine2, 3, "b1234")); },
                      "catalog number");
    expectErrorNaming([] { (void)Tle::parse(withField(kLine1, 3, "-5544"), withField(kLine2, 3, "-5544")); },
                      "catalog number");
    EXPECT_THROW((void)Tle::parse(kLine1, withField(kLine2, 3, "25545")), Error);
}

TEST(Tle, ParseManyReportsUnpairedLinesWithTheirLineNumber) {
    // A line 1 at the end of the text.
    expectErrorNaming([] { (void)Tle::parseMany(kLine1 + "\n" + kLine2 + "\n" + kLine1 + "\n"); }, "text line 3");
    // Two line 1s in a row: the first one (text line 2, after the name) has no line 2.
    expectErrorNaming([] { (void)Tle::parseMany("ISS (ZARYA)\n" + kLine1 + "\n" + kLine1 + "\n" + kLine2 + "\n"); },
                      "text line 2");
    // A line 2 on its own, here after a complete set, a blank line and a comment.
    expectErrorNaming([] { (void)Tle::parseMany(kLine1 + "\n" + kLine2 + "\n\n# note\n" + kLine2 + "\n"); },
                      "text line 5");
    // Something else between the two lines.
    expectErrorNaming([] { (void)Tle::parseMany(kLine1 + "\nISS (ZARYA)\n" + kLine2 + "\n"); }, "text line 1");
}

TEST(Tle, ParseManyAcceptsNamesBlankLinesCommentsAndBothLineEndings) {
    const std::vector<Tle> sets =
        Tle::parseMany("# two-line and three-line sets, mixed\n"
                       "\n"
                       + kLine1 + "\n" + kLine2 + "\n" + "ISS (ZARYA)\n" + kLine1 + "\r\n" + kLine2 + "\r\n" + "   \n"
                       + "0 ISS (ZARYA)\n" + kLine1 + "\n" + "# a comment between the lines of a set\n" + kLine2 + "\n"
                       + "1 KUNS-PF\n" + kLine1 + "\n" + kLine2 + "\n" + "2 SOPS\n" + kLine1 + "\n" + kLine2);
    ASSERT_EQ(sets.size(), 5u);
    EXPECT_EQ(sets[0].name, "");
    EXPECT_EQ(sets[1].name, "ISS (ZARYA)");
    EXPECT_EQ(sets[2].name, "ISS (ZARYA)"); // the "0 " of the three-line format is not part of the name
    EXPECT_EQ(sets[3].name, "1 KUNS-PF");   // a short line is a name even if it starts like a data line
    EXPECT_EQ(sets[4].name, "2 SOPS");
    for (const Tle& tle : sets) {
        EXPECT_EQ(tle.satnum, 25544);
        EXPECT_EQ(tle.line1(), kLine1);
        EXPECT_EQ(tle.line2(), kLine2);
    }
    EXPECT_TRUE(Tle::parseMany("").empty());
    EXPECT_TRUE(Tle::parseMany("# nothing here\n\n").empty());
}
