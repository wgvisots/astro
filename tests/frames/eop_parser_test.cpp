/**********************************************************************
 * eop_parser_test.cpp
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

// The Earth orientation file readers on damaged input. A damaged field must be reported, with its
// line: a value that is silently wrong, or a line that is silently dropped, turns into metres of
// error with nothing to show for it. (Found by fuzzing the parsers; see fuzz/eop.cpp.)

#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/frames/eop.h"

using namespace astro;

namespace {

/// One finals2000A line with the fields this library reads, each right-aligned in its columns
/// (counted from 0): MJD 7-14, x 18-26, y 37-45, UT1-UTC 58-67, LOD 79-85, dX 97-105, dY 116-124.
std::string finalsLine(const std::string& mjd,
                       const std::string& x,
                       const std::string& y,
                       const std::string& dut1,
                       const std::string& lod = "0.3000",
                       const std::string& dX  = "0.100",
                       const std::string& dY  = "-0.050") {
    std::string line(130, ' ');
    const auto  put = [&](std::size_t begin, std::size_t end, const std::string& text) {
        line.replace(end - text.size(), text.size(), text);
        (void)begin;
    };
    put(7, 15, mjd);
    put(18, 27, x);
    put(37, 46, y);
    put(58, 68, dut1);
    put(79, 86, lod);
    put(97, 106, dX);
    put(116, 125, dY);
    return line + "\n";
}

std::string threeGoodFinalsLines() {
    return finalsLine("60310.00", "0.130000", "0.200000", "-0.0100000")
         + finalsLine("60311.00", "0.131000", "0.198000", "-0.0103000")
         + finalsLine("60312.00", "0.132000", "0.196000", "-0.0106000");
}

/// The message of the astro::Error that parsing `text` throws; empty if it does not throw.
template <class Parse>
std::string errorOf(Parse&& parse, const std::string& text) {
    try {
        (void)parse(text);
    } catch (const Error& e) {
        return e.what();
    }
    return {};
}

const auto parseFinals = [](const std::string& text) { return EopTable::parseFinals2000A(text); };
const auto parseC04    = [](const std::string& text) { return EopTable::parseC04(text); };

} // namespace

TEST(EopParser, FinalsFieldsAreReadFromTheirColumns) {
    // The values written above, in the units of EopValues: arcseconds and milliarcseconds to
    // radians, milliseconds of length of day to seconds.
    const EopTable  table = EopTable::parseFinals2000A(threeGoodFinalsLines());
    const EopValues v     = table.at(60311.0);
    EXPECT_EQ(table.firstMjd(), 60310.0);
    EXPECT_EQ(table.lastMjd(), 60312.0);
    EXPECT_NEAR(v.xp, arcsecToRad(0.131), 1e-18);
    EXPECT_NEAR(v.yp, arcsecToRad(0.198), 1e-18);
    EXPECT_NEAR(v.dut1, -0.0103, 1e-15);
    EXPECT_NEAR(v.lod, 0.3e-3, 1e-18);
    EXPECT_NEAR(v.dX, arcsecToRad(0.1e-3), 1e-21);
    EXPECT_NEAR(v.dY, arcsecToRad(-0.05e-3), 1e-21);
}

TEST(EopParser, FinalsRowsWithoutDataAreSkippedNotRejected) {
    // The published file runs on past its last prediction with rows that hold a date and nothing
    // else. Those are not errors.
    const std::string text = threeGoodFinalsLines() + finalsLine("60313.00", "", "", "", "", "", "");
    EXPECT_EQ(EopTable::parseFinals2000A(text).lastMjd(), 60312.0);
}

TEST(EopParser, FinalsFieldThatIsNotAPlainNumberIsReportedWithItsLine) {
    // Each of these was accepted: "nan" and "inf" as values (every interpolation near them then
    // returned NaN), a number followed by other text as the number, and a value with an exponent.
    const struct {
        const char* x;
        const char* what;
    } cases[] = {
        {"nan",      "not a number"         },
        {"inf",      "infinity"             },
        {"0.2N10",   "text after the number"},
        {"1.3e-1",   "an exponent"          },
        {"0x1p-3",   "hexadecimal"          },
        {"0.13.000", "two decimal points"   },
    };
    for (const auto& c : cases) {
        const std::string text = finalsLine("60310.00", "0.130000", "0.200000", "-0.0100000")
                               + finalsLine("60311.00", c.x, "0.198000", "-0.0103000")
                               + finalsLine("60312.00", "0.132000", "0.196000", "-0.0106000");
        const std::string message = errorOf(parseFinals, text);
        EXPECT_NE(message.find("line 2"), std::string::npos) << c.what << ": '" << message << "'";
        EXPECT_NE(message.find("polar motion x"), std::string::npos) << c.what << ": '" << message << "'";
    }
    // The same for a field the row could do without: a damaged value is not a missing one.
    const std::string lod = finalsLine("60310.00", "0.130000", "0.200000", "-0.0100000", "nan")
                          + finalsLine("60311.00", "0.131000", "0.198000", "-0.0103000");
    EXPECT_NE(errorOf(parseFinals, lod).find("LOD"), std::string::npos);
}

TEST(EopParser, C04DataLineWithADamagedNumberIsReportedWithItsLine) {
    // Header lines (no 4-digit year first, or fewer than ten columns) are skipped as before. A data
    // line with a damaged number used to be dropped without a word, leaving a gap in the table.
    const std::string header =
        "# EOP C04\n  Date      MJD      x          y        UT1-UTC       LOD         dX        dY\n";
    const auto line = [](int day, const std::string& x) {
        return "2024   1  " + std::to_string(day) + "  " + std::to_string(60309 + day) + "   " + x
             + "   0.200000  -0.0100000   0.0003000   0.000100  -0.000050\n";
    };
    const EopTable good = EopTable::parseC04(header + line(1, "0.130000") + line(2, "0.131000") + line(3, "0.132000"));
    EXPECT_EQ(good.firstMjd(), 60310.0);
    EXPECT_NEAR(good.at(60311.0).xp, arcsecToRad(0.131), 1e-18);
    for (const char* damaged : {"0.13l000", "nan", "inf", "1.3e-1", "--0.13"}) {
        const std::string message =
            errorOf(parseC04, header + line(1, "0.130000") + line(2, damaged) + line(3, "0.132000"));
        EXPECT_NE(message.find("line 4"), std::string::npos) << damaged << ": '" << message << "'";
        EXPECT_NE(message.find("polar motion x"), std::string::npos) << damaged << ": '" << message << "'";
    }
}

TEST(EopTable, ValuesThatAreNotFiniteAreRejected) {
    for (double bad : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()}) {
        EopRecord a, b, c;
        a.mjd = 60310.0;
        b.mjd = 60311.0;
        c.mjd = 60312.0;
        for (int which = 0; which < 6; ++which) {
            EopRecord damaged  = b;
            double*   fields[] = {&damaged.values.xp,
                                  &damaged.values.yp,
                                  &damaged.values.dut1,
                                  &damaged.values.lod,
                                  &damaged.values.dX,
                                  &damaged.values.dY};
            *fields[which]     = bad;
            EXPECT_THROW(EopTable({a, damaged, c}), Error) << "value " << which << " = " << bad;
        }
    }
}

TEST(EopTable, InterpolationThatOverflowsThrowsInsteadOfReturningInfinity) {
    // Finite records can still interpolate to infinity: values near the top of the range of a double
    // with the cubic's weights, which are larger than 1 between unevenly spaced records.
    std::vector<EopRecord> records(4);
    const double           mjd[] = {60310.0, 60310.001, 60312.0, 60315.0};
    for (std::size_t k = 0; k < 4; ++k) {
        records[k].mjd       = mjd[k];
        records[k].values.xp = (k % 2 == 0 ? 1.0 : -1.0) * 1.5e308;
    }
    const EopTable table(records, EopInterpolation::Lagrange4);
    EXPECT_THROW((void)table.at(60311.0), Error);
    // Ordinary values on the same epochs interpolate as usual.
    for (EopRecord& r : records) {
        r.values.xp = 1e-6;
    }
    EXPECT_NEAR(EopTable(records, EopInterpolation::Lagrange4).at(60311.0).xp, 1e-6, 1e-18);
}
