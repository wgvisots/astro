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

#include "astro/core/angle.h"
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
