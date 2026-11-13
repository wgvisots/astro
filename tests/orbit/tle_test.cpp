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

#include <cmath>
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

} // namespace

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
