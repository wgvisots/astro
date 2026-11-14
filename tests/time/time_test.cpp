/**********************************************************************
 * time_test.cpp
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

#include "astro/time/time.h"

#include <chrono>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "astro/core/error.h"
#include "astro/frames/eop.h"
#include "astro/time/leap_seconds.h"
#include "test_util.h"

using namespace astro;
using astro::test::dataPath;
using astro::test::readCsv;

namespace {
// Difference in seconds between two two-part JDs.
double diffSec(const JulianDate& a, double b1, double b2) {
    return ((a.jd1 - b1) + (a.jd2 - b2)) * 86400.0;
}

// `call` must throw astro::Error, and the message must contain `what`.
template <typename Call>
void expectErrorSaying(Call&& call, const std::string& what, const char* label) {
    try {
        call();
        ADD_FAILURE() << label << ": no error";
    } catch (const Error& e) {
        EXPECT_NE(std::string(e.what()).find(what), std::string::npos)
            << label << ": expected \"" << what << "\" in the message, which was: " << e.what();
    }
}

const double kNaN = std::numeric_limits<double>::quiet_NaN();
const double kInf = std::numeric_limits<double>::infinity();

// Puts the built-in leap-second table back when the test ends, also when an assertion or an exception ends it
// early, so that no later test runs with a table a test installed.
struct RestoreBuiltinLeapSecondTable {
    RestoreBuiltinLeapSecondTable()                                                = default;
    RestoreBuiltinLeapSecondTable(const RestoreBuiltinLeapSecondTable&)            = delete;
    RestoreBuiltinLeapSecondTable& operator=(const RestoreBuiltinLeapSecondTable&) = delete;
    ~RestoreBuiltinLeapSecondTable() { setLeapSecondTable(LeapSecondTable::builtin()); }
};

// TAI - UTC from 1960 to 1971, when UTC ran at a rate of its own: offset + (MJD_UTC - mjdRef) * rate seconds
// from 0h UTC on the first of the month given. Copied from USNO tai-utc.dat (from 1961) and from the SOFA/ERFA
// function dat (the 1960 line, which tai-utc.dat does not have); written out here so that the test does not
// read the library's own table.
struct DriftEraStep {
    int    year, month;
    double offset, mjdRef, rate;
};
// clang-format off
constexpr DriftEraStep kDriftEra[] = {
    {1960,  1, 1.4178180, 37300.0, 0.001296 },
    {1961,  1, 1.4228180, 37300.0, 0.001296 },
    {1961,  8, 1.3728180, 37300.0, 0.001296 },
    {1962,  1, 1.8458580, 37665.0, 0.0011232},
    {1963, 11, 1.9458580, 37665.0, 0.0011232},
    {1964,  1, 3.2401300, 38761.0, 0.001296 },
    {1964,  4, 3.3401300, 38761.0, 0.001296 },
    {1964,  9, 3.4401300, 38761.0, 0.001296 },
    {1965,  1, 3.5401300, 38761.0, 0.001296 },
    {1965,  3, 3.6401300, 38761.0, 0.001296 },
    {1965,  7, 3.7401300, 38761.0, 0.001296 },
    {1965,  9, 3.8401300, 38761.0, 0.001296 },
    {1966,  1, 4.3131700, 39126.0, 0.002592 },
    {1968,  2, 4.2131700, 39126.0, 0.002592 },
};
// clang-format on

// TAI - UTC in seconds at the UTC instant MJD `mjdUtc`, in the month (year, month) of 1960-1971.
double driftEraTaiMinusUtc(int year, int month, double mjdUtc) {
    const DriftEraStep* step = nullptr;
    for (const DriftEraStep& s : kDriftEra) {
        if (s.year * 12 + s.month <= year * 12 + month) {
            step = &s;
        }
    }
    if (step == nullptr) {
        throw std::logic_error("no drift-era step before the date");
    }
    return step->offset + (mjdUtc - step->mjdRef) * step->rate;
}
} // namespace

TEST(JulianDate, CalendarRoundTrip) {
    EXPECT_DOUBLE_EQ(JulianDate::fromCalendar(2000, 1, 1, 12).jd(), 2451545.0);
    EXPECT_DOUBLE_EQ(JulianDate::fromCalendar(1858, 11, 17).mjd(), 0.0);
    EXPECT_DOUBLE_EQ(JulianDate::fromCalendar(-4713, 11, 24, 12).jd(), 0.0); // proleptic Gregorian JD 0
    const CalendarDate c = JulianDate::fromCalendar(2024, 2, 29, 23, 59, 59.123456).calendar(6);
    EXPECT_EQ(c.year, 2024);
    EXPECT_EQ(c.month, 2);
    EXPECT_EQ(c.day, 29);
    EXPECT_EQ(c.hour, 23);
    EXPECT_EQ(c.minute, 59);
    EXPECT_NEAR(c.second, 59.123456, 1e-9);
    // Rounding carries into the next day.
    const CalendarDate n = JulianDate::fromCalendar(2023, 12, 31, 23, 59, 59.9999999).calendar(3);
    EXPECT_EQ(n.year, 2024);
    EXPECT_EQ(n.month, 1);
    EXPECT_EQ(n.day, 1);
    EXPECT_EQ(n.hour, 0);
}

TEST(JulianDate, Helpers) {
    EXPECT_EQ(dayOfYear(2024, 12, 31), 366);
    EXPECT_EQ(dayOfYear(2023, 3, 1), 60);
    EXPECT_EQ(dayOfYearToDate(2024, 60), (YearMonthDay{2024, 2, 29}));
    EXPECT_EQ(dayOfWeek(1980, 1, 6), 0); // GPS epoch was a Sunday
    EXPECT_EQ(dayOfWeek(2000, 1, 1), 6); // Saturday
    EXPECT_TRUE(isLeapYear(2000));
    EXPECT_FALSE(isLeapYear(1900));
    // Precision of large second offsets.
    const JulianDate j = JulianDate{constants::jdGpsEpoch, 0.0}.addSeconds(1.4e9 + 0.000001);
    EXPECT_NEAR(j.secondsSince(JulianDate{constants::jdGpsEpoch, 0.0}), 1.4e9 + 0.000001, 1e-9);
}

TEST(LeapSeconds, BuiltinValues) {
    const auto t = leapSecondTable();
    EXPECT_DOUBLE_EQ(t->taiMinusUtc(2017, 1, 1), 37.0);
    EXPECT_DOUBLE_EQ(t->taiMinusUtc(2016, 12, 31, 0.999), 36.0);
    EXPECT_DOUBLE_EQ(t->taiMinusUtc(1972, 1, 1), 10.0);
    EXPECT_NEAR(t->taiMinusUtc(1965, 3, 1, 0.5), 3.717242, 1e-6); // drift era (matches ERFA)
    EXPECT_THROW((void)t->taiMinusUtc(1959, 12, 31), Error);
}

TEST(LeapSeconds, ParseFileFormats) {
    const char*           leapSecondDat = "#  File expires on 28 June 2026\n"
                                          "#    MJD        Date        TAI-UTC (s)\n"
                                          "    41317.0    1  1 1972       10\n"
                                          "    41499.0    1  7 1972       11\n"
                                          "    57754.0    1  1 2017       37\n"
                                          "    61041.0    1  1 2026       38\n";
    const LeapSecondTable a             = LeapSecondTable::parse(leapSecondDat);
    EXPECT_DOUBLE_EQ(a.taiMinusUtc(2026, 1, 1), 38.0);
    EXPECT_DOUBLE_EQ(a.taiMinusUtc(2025, 12, 31), 37.0);
    EXPECT_NEAR(a.taiMinusUtc(1965, 3, 1, 0.5), 3.717242, 1e-6); // pre-1972 retained
    ASSERT_TRUE(a.expiresMjd().has_value());
    EXPECT_DOUBLE_EQ(*a.expiresMjd(), JulianDate::fromCalendar(2026, 6, 28).mjd());

    const char*           list = "#@\t3960057600\n"
                                 "2272060800\t10\t# 1 Jan 1972\n"
                                 "2287785600\t11\t# 1 Jul 1972\n"
                                 "3692217600\t37\t# 1 Jan 2017\n";
    const LeapSecondTable b    = LeapSecondTable::parse(list);
    EXPECT_DOUBLE_EQ(b.taiMinusUtc(1972, 6, 30), 10.0);
    EXPECT_DOUBLE_EQ(b.taiMinusUtc(1972, 7, 1), 11.0);
    EXPECT_DOUBLE_EQ(b.taiMinusUtc(2017, 1, 1), 37.0);
    ASSERT_TRUE(b.expiresMjd().has_value());
    EXPECT_NEAR(*b.expiresMjd(), JulianDate::fromCalendar(2025, 6, 28).mjd(), 1e-9);

    const char*           taiUtc = " 1961 JAN  1 =JD 2437300.5  TAI-UTC=   1.4228180 S + (MJD - 37300.) X 0.001296 S\n"
                                   " 1972 JAN  1 =JD 2441317.5  TAI-UTC=  10.0       S + (MJD - 41317.) X 0.0      S\n";
    const LeapSecondTable c      = LeapSecondTable::parse(taiUtc);
    EXPECT_NEAR(c.taiMinusUtc(1961, 1, 1), 1.4228180, 1e-12);
    EXPECT_DOUBLE_EQ(c.taiMinusUtc(1980, 1, 1), 10.0);
}

TEST(LeapSeconds, MalformedFilesThrow) {
    EXPECT_THROW((void)LeapSecondTable::parse("    41317.0    1 13 1972       10\n"), Error); // month 13
    EXPECT_THROW((void)LeapSecondTable::parse("    41317.0    1  1 99999  10\n"), Error);     // year > 9999
    EXPECT_THROW((void)LeapSecondTable::parse(" 1972 JAN  1 =JD 2441317.5  TAI-UTC= abc S + (MJD - 41317.) X 0.0 S\n"),
                 Error);
    EXPECT_THROW((void)LeapSecondTable::parse("#@ 3960057600\n1e300\t10\n"), Error);
    EXPECT_THROW((void)LeapSecondTable::builtin().taiMinusUtc(std::nan("")), Error);
}

TEST(TimeScales, MatchesErfa) {
    const auto csv = readCsv(dataPath("time_scales.csv"));
    ASSERT_GT(csv.rows.size(), 300u);
    const ConstantEop eop(EopValues{.dut1 = 0.1234567});
    auto num = [&](const std::vector<std::string>& r, const char* name) { return std::stod(r[csv.col(name)]); };
    // Tolerances. Scales: 1e-10 s, against a worst difference of 1.9e-11 s (TDB, TCG and TCB; TAI and TT agree
    // to the bit), which is a few units in the last place of a day fraction (4.8e-12 s). TDB - TT: the series is
    // ERFA's, summed in the same order, and agrees to the bit; 1e-15 s (of 1.7e-3 s) allows another compiler's
    // rounding. These were 1e-9 s and 1e-12 s.
    for (const auto& r : csv.rows) {
        const Time        utc   = Time::fromCalendar(TimeScale::UTC,
                                            std::stoi(r[0]),
                                            std::stoi(r[1]),
                                            std::stoi(r[2]),
                                            std::stoi(r[3]),
                                            std::stoi(r[4]),
                                            num(r, "second"));
        const std::string label = utc.iso8601(6);
        EXPECT_NEAR(diffSec(utc.jd(), num(r, "utc1"), num(r, "utc2")), 0.0, 1e-10) << label;
        EXPECT_NEAR(diffSec(utc.to(TimeScale::TAI).jd(), num(r, "tai1"), num(r, "tai2")), 0.0, 1e-10) << label;
        EXPECT_NEAR(diffSec(utc.to(TimeScale::TT).jd(), num(r, "tt1"), num(r, "tt2")), 0.0, 1e-10) << label;
        EXPECT_NEAR(diffSec(utc.to(TimeScale::TDB).jd(), num(r, "tdb1"), num(r, "tdb2")), 0.0, 1e-10) << label;
        EXPECT_NEAR(diffSec(utc.to(TimeScale::TCG).jd(), num(r, "tcg1"), num(r, "tcg2")), 0.0, 1e-10) << label;
        EXPECT_NEAR(diffSec(utc.to(TimeScale::TCB).jd(), num(r, "tcb1"), num(r, "tcb2")), 0.0, 1e-10) << label;
        const int  year         = std::stoi(r[0]);
        const bool inLeapSecond = num(r, "second") >= 60.0;
        // ERFA's eraUtcut1 evaluates TAI - UTC at 0h of the day, which before 1972, while UTC drifted, is off by
        // up to 1.3 ms within a day; so its UT1 column is used from 1972 only. The tests below cover the years
        // before.
        if (year >= 1972) {
            EXPECT_NEAR(diffSec(utc.to(TimeScale::UT1, eop).jd(), num(r, "ut11"), num(r, "ut12")), 0.0, 1e-10) << label;
        }
        EXPECT_NEAR(tdbMinusTt(utc.to(TimeScale::TT).jd()), num(r, "tdb_minus_tt"), 1e-15) << label;

        // Round trips back to UTC from every scale. (UT1 is skipped inside a leap second: with a constant
        // dut1, which real EOP never has across a leap second, UT1 is discontinuous there.)
        for (TimeScale s : {TimeScale::TAI,
                            TimeScale::TT,
                            TimeScale::TDB,
                            TimeScale::TCG,
                            TimeScale::TCB,
                            TimeScale::GPS,
                            TimeScale::GST,
                            TimeScale::BDT,
                            TimeScale::UT1}) {
            if (s == TimeScale::UT1 && inLeapSecond) {
                continue;
            }
            const Time back = utc.to(s, eop).to(TimeScale::UTC, eop);
            EXPECT_NEAR(back - utc, 0.0, 1e-10) << label << " via " << toString(s);
            EXPECT_EQ(back.iso8601(6), label) << " via " << toString(s);
        }
    }
}

namespace {

// UT1 before 1972 against UT1 = TAI - (TAI - UTC) + (UT1 - UTC), with TAI from ERFA (the tai1, tai2 columns of
// time_scales.csv, which agree with the library to 1e-11 s) and TAI - UTC at the instant from the published
// coefficients (kDriftEra) at the plain UTC time of day. ERFA's utc1 is the midnight that starts the UTC day,
// so utc1 - 2400000.5 is its MJD. `onShortenedDay` picks the rows on 1961-07-31 and 1968-01-31, the two UTC
// days that were shorter than 86400 s, or the other rows.
void checkDriftEraUt1(bool onShortenedDay) {
    const auto   csv = readCsv(dataPath("time_scales.csv"));
    auto         num = [&](const std::vector<std::string>& r, const char* name) { return std::stod(r[csv.col(name)]); };
    const double dut1 = 0.1234567;
    const ConstantEop eop(EopValues{.dut1 = dut1});
    int               rows = 0;
    for (const auto& r : csv.rows) {
        const int  year = std::stoi(r[0]), month = std::stoi(r[1]), day = std::stoi(r[2]);
        const bool shortened = (year == 1961 && month == 7 && day == 31) || (year == 1968 && month == 1 && day == 31);
        if (year >= 1972 || shortened != onShortenedDay) {
            continue;
        }
        ++rows;
        const Time utc =
            Time::fromCalendar(TimeScale::UTC, year, month, day, std::stoi(r[3]), std::stoi(r[4]), num(r, "second"));
        const double fractionOfDay = (std::stoi(r[3]) * 3600.0 + std::stoi(r[4]) * 60.0 + num(r, "second")) / 86400.0;
        const double dat           = driftEraTaiMinusUtc(year, month, (num(r, "utc1") - 2400000.5) + fractionOfDay);
        // Tolerance 1e-10 s. ERFA's TAI agrees with the published TAI - UTC to 2.3e-11 s (worked out in exact
        // rational arithmetic), and on the ordinary days the worst difference seen here was 4.8e-12 s.
        EXPECT_NEAR(diffSec(utc.to(TimeScale::UT1, eop).jd(), num(r, "tai1"), num(r, "tai2")), dut1 - dat, 1e-10)
            << utc.iso8601(6);
    }
    // 48 ordinary rows and 6 on the shortened days (three on each: 0h, 12h and half a second before the end).
    EXPECT_EQ(rows, onShortenedDay ? 6 : 48);
}

} // namespace

TEST(TimeScales, Ut1Before1972IsUtcPlusDut1) {
    checkDriftEraUt1(false);
}

// On the two shortened UTC days (1961-07-31 and 1968-01-31, 0.05 s and 0.1 s short) TAI - UTC must be
// evaluated at the clock time of day, (seconds of day) / 86400, not at the quasi-JD fraction, which is
// scaled to the short day. The conversion to UT1 used the second, and was off by the drift rate times the
// shift in the fraction: 3.0e-9 s at 1968-01-31T23:59:59.5, 1.5e-9 s at 12h, zero at 0h.
TEST(TimeScales, Ut1OnTheShortenedUtcDaysIsUtcPlusDut1) {
    checkDriftEraUt1(true);
}

TEST(TimeScales, TopocentricTdb) {
    // ERFA eraDtdb with the observer terms, 50 random epochs and observers (tools/gen_reference.py). An empty or
    // cut-short file must fail rather than pass with nothing compared. Tolerance 1e-15 s: the arithmetic is ERFA's
    // and agrees to the bit here (this was 1e-13 s).
    const auto csv = readCsv(dataPath("dtdb_topo.csv"));
    ASSERT_EQ(csv.rows.size(), 50u);
    for (const auto& r : csv.rows) {
        EXPECT_NEAR(tdbMinusTt(JulianDate{std::stod(r[0]), 0.0},
                               std::stod(r[1]),
                               std::stod(r[2]),
                               std::stod(r[3]),
                               std::stod(r[4])),
                    std::stod(r[5]),
                    1e-15);
    }
}

TEST(Time, LeapSecondHandling) {
    const Time t = Time::fromCalendar(TimeScale::UTC, 2016, 12, 31, 23, 59, 60.5);
    EXPECT_EQ(t.iso8601(3), "2016-12-31T23:59:60.500Z");
    EXPECT_EQ(t.to(TimeScale::TAI).iso8601(3), "2017-01-01T00:00:36.500 TAI");
    // Elapsed time across the leap second counts it.
    const Time a = Time::fromCalendar(TimeScale::UTC, 2016, 12, 31, 23, 59, 59.0);
    const Time b = Time::fromCalendar(TimeScale::UTC, 2017, 1, 1, 0, 0, 0.0);
    EXPECT_NEAR(b - a, 2.0, 1e-9);
    EXPECT_EQ((a + 1.0).iso8601(3), "2016-12-31T23:59:60.000Z");
    EXPECT_EQ((a + 2.0).iso8601(3), "2017-01-01T00:00:00.000Z");
    // No leap second on an ordinary day.
    EXPECT_THROW((void)Time::fromCalendar(TimeScale::UTC, 2016, 12, 30, 23, 59, 60.0), Error);
    EXPECT_THROW((void)Time::fromCalendar(TimeScale::TT, 2016, 12, 31, 23, 59, 60.0), Error);
    EXPECT_TRUE(a < b);
}

TEST(Time, LastMinuteOfAShortenedUtcDayIsShorter) {
    // TAI-UTC stepped down by 0.05 s at 0h on 1961-08-01 and by 0.1 s on 1968-02-01 (USNO tai-utc.dat), so the UTC
    // days before those dates were 86399.95 s and 86399.9 s long and their last minutes end at 23:59:59.95 and
    // 23:59:59.9. A later reading of the clock did not exist, just as 23:59:60 does not exist on an ordinary day.
    //
    // The readings 0.01 s before those ends are 0.01 s before the next midnight. Tolerance 1e-8 s: a UTC second
    // was up to 3e-8 longer than an SI second then (3e-10 s here), and the conversion, like ERFA's, treats the
    // step and the rate separately, which leaves step * rate = 0.1 s * 3e-8 = 3e-9 s at such a midnight.
    const Time aug1961 = Time::fromCalendar(TimeScale::UTC, 1961, 8, 1);
    EXPECT_NEAR(aug1961 - Time::fromCalendar(TimeScale::UTC, 1961, 7, 31, 23, 59, 59.94), 0.01, 1e-8);
    EXPECT_THROW((void)Time::fromCalendar(TimeScale::UTC, 1961, 7, 31, 23, 59, 59.96), Error);
    EXPECT_THROW((void)Time::fromIso8601("1961-07-31T23:59:59.999Z"), Error);
    const Time feb1968 = Time::fromCalendar(TimeScale::UTC, 1968, 2, 1);
    EXPECT_NEAR(feb1968 - Time::fromCalendar(TimeScale::UTC, 1968, 1, 31, 23, 59, 59.89), 0.01, 1e-8);
    EXPECT_THROW((void)Time::fromCalendar(TimeScale::UTC, 1968, 1, 31, 23, 59, 59.91), Error);
    // Only the last minute of such a day is shorter, and the days either side are ordinary.
    EXPECT_NO_THROW((void)Time::fromCalendar(TimeScale::UTC, 1961, 7, 31, 23, 58, 59.99));
    EXPECT_NO_THROW((void)Time::fromCalendar(TimeScale::UTC, 1961, 7, 30, 23, 59, 59.99));
    EXPECT_NO_THROW((void)Time::fromCalendar(TimeScale::UTC, 1961, 8, 1, 23, 59, 59.99));
}

TEST(Time, NegativeLeapSecondRemovesTheLastSecondOfTheDay) {
    // None has happened yet. If TAI-UTC went from 37 s to 36 s at the start of 2030, the UTC day 2029-12-31 would
    // be 86399 s long and end at 23:59:59 (Recommendation ITU-R TF.460-6): 23:59:58.5 is in its last second.
    const RestoreBuiltinLeapSecondTable restore;
    std::vector<LeapSecondEntry>        entries = LeapSecondTable::builtin().entries();
    entries.push_back({2030, 1, 36.0});
    setLeapSecondTable(LeapSecondTable(entries));

    const Time last = Time::fromCalendar(TimeScale::UTC, 2029, 12, 31, 23, 59, 58.5);
    // 23:59:58.5 UTC + 37 s is 00:00:35.5 TAI on 2030-01-01, half a second before 00:00:00 UTC (00:00:36 TAI).
    EXPECT_NEAR(last.to(TimeScale::TAI).jd().secondsSince(JulianDate::fromCalendar(2030, 1, 1)), 35.5, 1e-9);
    EXPECT_NEAR(Time::fromCalendar(TimeScale::UTC, 2030, 1, 1) - last, 0.5, 1e-9);
    EXPECT_EQ(last.iso8601(1), "2029-12-31T23:59:58.5Z");
    EXPECT_EQ((last + 0.5).iso8601(1), "2030-01-01T00:00:00.0Z");
    // 23:59:59 is not a time on that day.
    EXPECT_THROW((void)Time::fromCalendar(TimeScale::UTC, 2029, 12, 31, 23, 59, 59.5), Error);
    EXPECT_THROW((void)Time::fromIso8601("2029-12-31T23:59:59Z"), Error);
    EXPECT_THROW((void)Time::fromIso8601("2030-01-01T00:59:59+01:00"), Error);
}

TEST(Time, Unix) {
    EXPECT_EQ(Time::fromUnix(0.0).iso8601(3), "1970-01-01T00:00:00.000Z");
    const Time t = Time::fromUnix(1483228800.25); // 2017-01-01T00:00:00.25Z
    EXPECT_EQ(t.iso8601(3), "2017-01-01T00:00:00.250Z");
    EXPECT_NEAR(t.unixSeconds(), 1483228800.25, 1e-6);
    EXPECT_NEAR(Time::fromUnix(1483228799.5).unixSeconds(), 1483228799.5, 1e-6);
    EXPECT_EQ(Time::fromUnix(1483228799.5).iso8601(1), "2016-12-31T23:59:59.5Z");
    // Leap second maps onto the next Unix second.
    EXPECT_NEAR(Time::fromCalendar(TimeScale::UTC, 2016, 12, 31, 23, 59, 60.5).unixSeconds(), 1483228800.5, 1e-6);
    // Full precision path.
    const Time p = Time::fromUnix(1700000000LL, 0.123456789);
    EXPECT_EQ(p.iso8601(9), "2023-11-14T22:13:20.123456789Z");
    const auto tp = p.sysTime();
    EXPECT_EQ(Time::fromSysTime(tp).iso8601(6), "2023-11-14T22:13:20.123457Z");
}

TEST(Time, Gnss) {
    // GPS epoch.
    const Time e = Time::fromCalendar(TimeScale::UTC, 1980, 1, 6);
    EXPECT_NEAR(e.gpsSeconds(), 0.0, 1e-9);
    // 2017-01-01T00:00:00Z: GPS = TAI - 19 s = UTC + 18 s. It is 13510 days = 1930 weeks after the GPS epoch
    // 1980-01-06 (JD 2457754.5 - 2444244.5), so GPS week 1930, time of week 18 s.
    const Time        t = Time::fromCalendar(TimeScale::UTC, 2017, 1, 1);
    const GnssWeekTow g = t.gnssWeekTow(GnssSystem::GPS);
    EXPECT_EQ(g.week, 1930);
    EXPECT_NEAR(g.tow, 18.0, 1e-9);
    const GnssWeekTow gal = t.gnssWeekTow(GnssSystem::Galileo);
    EXPECT_EQ(gal.week, 1930 - 1024);
    EXPECT_NEAR(gal.tow, 18.0, 1e-9);
    // BeiDou: BDT = TAI - 33 s = UTC + 4 s in 2017. 2017-01-01 is 4018 days = 574 weeks after the BDT epoch
    // 2006-01-01 (JD 2457754.5 - 2453736.5), so BDT week 574 starts at 0h BDT that day.
    const GnssWeekTow bd = t.gnssWeekTow(GnssSystem::BeiDou);
    EXPECT_EQ(bd.week, 574);
    EXPECT_NEAR(bd.tow, 4.0, 1e-9);
    EXPECT_NEAR(Time::fromGnssWeekTow(GnssSystem::GPS, 1930, 18.0) - t, 0.0, 1e-9);
    EXPECT_NEAR(Time::fromGnssWeekTow(GnssSystem::Galileo, 906, 18.0) - t, 0.0, 1e-9);
    EXPECT_NEAR(Time::fromGnssWeekTow(GnssSystem::BeiDou, 574, 4.0) - t, 0.0, 1e-9);
    // Rollover: 10-bit week 1930 % 1024 = 906 near 2017.
    EXPECT_EQ(resolveWeekRollover(906, 10, GnssSystem::GPS, t), 1930);
    EXPECT_EQ(resolveWeekRollover(2, 10, GnssSystem::GPS, Time::fromCalendar(TimeScale::UTC, 2019, 4, 6)), 2050);
    EXPECT_NEAR(Time::fromGpsSeconds(t.gpsSeconds()) - t, 0.0, 1e-9);
}

TEST(Time, GnssTimeOfWeekStaysBelowOneWeek) {
    // Week 0 of each system starts at 0h on these Julian dates, in the system's own time scale: GPS 1980-01-06
    // (IS-GPS-200), Galileo 1999-08-22, which is GPS week 1024 (Galileo OS SIS ICD), BeiDou 2006-01-01 (BeiDou ICD).
    struct System {
        GnssSystem system;
        TimeScale  scale;
        double     epochJd;
    };
    const System systems[] = {
        {GnssSystem::GPS,     TimeScale::GPS, 2444244.5},
        {GnssSystem::Galileo, TimeScale::GST, 2451412.5},
        {GnssSystem::BeiDou,  TimeScale::BDT, 2453736.5}
    };
    // Fractions of a day: the last one is the largest below 1, a tenth of a nanosecond before midnight.
    const double fractions[] = {0.0, 1.0e-17, 0.25, 0.75, std::nextafter(1.0, 0.0)};
    for (const System& s : systems) {
        for (const int week : {0, 1, 1023, 2000, 2399}) {
            for (int day = 0; day < 7; ++day) {
                for (const double fraction : fractions) {
                    const Time        t = Time::fromJd(s.epochJd + 7.0 * week + day, fraction, s.scale);
                    const GnssWeekTow g = t.gnssWeekTow(s.system);
                    EXPECT_GE(g.tow, 0.0) << week << " " << day << " " << fraction;
                    EXPECT_LT(g.tow, 604800.0) << week << " " << day << " " << fraction;
                    // Week and time of week together are the seconds since the epoch, here in extended precision.
                    // 1.2e-10 s is the spacing of doubles just below 604800: the time of week cannot be closer.
                    const long double want = (7.0L * week + day + static_cast<long double>(fraction)) * 86400.0L;
                    const long double got  = 604800.0L * g.week + static_cast<long double>(g.tow);
                    EXPECT_NEAR(static_cast<double>(got - want), 0.0, 1.2e-10)
                        << week << " " << day << " " << fraction << ": week " << g.week << ", tow " << g.tow;
                }
            }
        }
    }
    // The last instant of GPS week 2000 that a Time can hold is week 2001, time of week 0, to the nearest double.
    const GnssWeekTow last =
        Time::fromJd(2444244.5 + 7.0 * 2000 + 6, std::nextafter(1.0, 0.0), TimeScale::GPS).gnssWeekTow(GnssSystem::GPS);
    EXPECT_EQ(last.week, 2001);
    EXPECT_EQ(last.tow, 0.0);
}

TEST(Time, Iso8601) {
    EXPECT_EQ(Time::fromIso8601("2024-03-01T12:34:56.789Z").iso8601(3), "2024-03-01T12:34:56.789Z");
    EXPECT_EQ(Time::fromIso8601("2024-03-01 12:34").iso8601(0), "2024-03-01T12:34:00Z");
    EXPECT_EQ(Time::fromIso8601("2024-03-01").iso8601(0), "2024-03-01T00:00:00Z");
    EXPECT_EQ(Time::fromIso8601("2024-061T06:00:00").iso8601(0), "2024-03-01T06:00:00Z");
    EXPECT_EQ(Time::fromIso8601("2024-03-01T12:00:00+05:30").iso8601(0), "2024-03-01T06:30:00Z");
    EXPECT_EQ(Time::fromIso8601("2024-03-01T00:00:00-01:00").iso8601(0), "2024-03-01T01:00:00Z");
    const Time tt = Time::fromIso8601("2000-01-01T12:00:00 TT");
    EXPECT_EQ(tt.scale(), TimeScale::TT);
    EXPECT_DOUBLE_EQ(tt.jd().jd(), constants::jdJ2000);
    EXPECT_EQ(Time::fromIso8601("2000-01-01T12:00:00", TimeScale::TDB).scale(), TimeScale::TDB);
    EXPECT_EQ(Time::fromIso8601("2016-12-31T23:59:60.5Z").iso8601(1), "2016-12-31T23:59:60.5Z");
    EXPECT_THROW((void)Time::fromIso8601("2024-3-1"), Error);
    EXPECT_THROW((void)Time::fromIso8601("2024-03-01T12:00:00 XYZ"), Error);
    EXPECT_THROW((void)Time::fromIso8601("2024-03-01T12:00:00Zjunk"), Error);
}

TEST(Time, Iso8601OffsetIsCivilClockArithmetic) {
    // A clock keeping "UTC+01:00" shows the same seconds as a UTC clock and is one hour ahead in its hours and
    // minutes, also on a day with a leap second. The expected UTC calendar fields below are worked out by hand.
    struct Case {
        const char* text;
        int         year, month, day, hour, minute;
        double      second;
    };
    // clang-format off
    const Case cases[] = {
        // The leap-second day 2016-12-31 (it ends with 23:59:60) and the days either side of it.
        {"2016-12-30T23:30:00+01:00",    2016, 12, 30, 22, 30,  0.0 },
        {"2016-12-31T00:30:00+01:00",    2016, 12, 30, 23, 30,  0.0 }, // back into the day before
        {"2016-12-30T23:30:00-01:00",    2016, 12, 31,  0, 30,  0.0 }, // forward into the leap-second day
        {"2016-12-31T12:00:00+05:30",    2016, 12, 31,  6, 30,  0.0 }, // inside the leap-second day
        {"2017-01-01T00:30:00+01:00",    2016, 12, 31, 23, 30,  0.0 }, // back across the leap second
        {"2016-12-31T18:30:00-06:00",    2017,  1,  1,  0, 30,  0.0 }, // forward across the leap second
        {"2017-01-01T00:59:59.25+01:00", 2016, 12, 31, 23, 59, 59.25}, // the second before the leap second
        {"2017-01-01T00:59:60+01:00",    2016, 12, 31, 23, 59, 60.0 }, // the leap second itself
        {"2016-12-31T18:59:60.5-05:00",  2016, 12, 31, 23, 59, 60.5 },
        {"2017-01-01T05:29:60.75+05:30", 2016, 12, 31, 23, 59, 60.75}, // a zone with a half-hour offset
        {"2017-01-01T01:00:00+01:00",    2017,  1,  1,  0,  0,  0.0 }, // the second after the leap second
        {"2017-01-02T00:30:00+01:00",    2017,  1,  1, 23, 30,  0.0 },
        {"2017-01-01T23:30:00-01:00",    2017,  1,  2,  0, 30,  0.0 },
        // Month and year boundaries, with and without a leap day.
        {"2024-03-01T00:30:00+01:00",    2024,  2, 29, 23, 30,  0.0 },
        {"2024-02-29T23:30:00-01:00",    2024,  3,  1,  0, 30,  0.0 },
        {"2023-03-01T00:30:00+01:00",    2023,  2, 28, 23, 30,  0.0 },
        {"2024-01-01T00:30:00+14:00",    2023, 12, 31, 10, 30,  0.0 },
        {"2023-12-31T23:30:00-12:00",    2024,  1,  1, 11, 30,  0.0 },
        // The largest offsets.
        {"2024-03-01T00:00:00+23:59",    2024,  2, 29,  0,  1,  0.0 },
        {"2024-03-01T23:59:00-23:59",    2024,  3,  2, 23, 58,  0.0 },
        // Other spellings: no colon, hours only, an ordinal date, a zero offset.
        {"2024-03-01T12:00:00+0530",     2024,  3,  1,  6, 30,  0.0 },
        {"2024-03-01T12:00:00+05",       2024,  3,  1,  7,  0,  0.0 },
        {"2024-061T12:00:00-00:45",      2024,  3,  1, 12, 45,  0.0 },
        {"2016-12-31T23:59:60.5+00:00",  2016, 12, 31, 23, 59, 60.5 },
        // Before 1972 a UTC second was not an SI second; the clock arithmetic is the same.
        {"1965-01-01T05:00:00+06:00",    1964, 12, 31, 23,  0,  0.0 },
    };
    // clang-format on
    for (const Case& c : cases) {
        const Time got  = Time::fromIso8601(c.text);
        const Time want = Time::fromCalendar(TimeScale::UTC, c.year, c.month, c.day, c.hour, c.minute, c.second);
        EXPECT_EQ(got.scale(), TimeScale::UTC) << c.text;
        // The same calendar fields give the same two-part Julian date, bit for bit.
        EXPECT_EQ(got.jd().jd1, want.jd().jd1) << c.text << " gave " << got.iso8601(6);
        EXPECT_EQ(got.jd().jd2, want.jd().jd2) << c.text << " gave " << got.iso8601(6);
    }

    // One hour on the "+01:00" clock that contains the leap second lasts 3601 s. (1e-9 s: about a hundred times
    // the resolution of a two-part Julian date.)
    EXPECT_NEAR(Time::fromIso8601("2017-01-01T01:30:00+01:00") - Time::fromIso8601("2017-01-01T00:30:00+01:00"),
                3601.0,
                1e-9);
    // TAI - UTC was 36 s until the end of that leap second, so it began at 00:00:36 TAI.
    EXPECT_EQ(Time::fromIso8601("2017-01-01T00:59:60+01:00").to(TimeScale::TAI).iso8601(3),
              "2017-01-01T00:00:36.000 TAI");
}

TEST(Time, Iso8601OffsetRejectsInvalidFields) {
    // A leap second exists only at 23:59:60 UTC; with these offsets ":60" falls in another UTC minute.
    EXPECT_THROW((void)Time::fromIso8601("2017-01-01T00:59:60+00:30"), Error); // 00:29:60 UTC
    EXPECT_THROW((void)Time::fromIso8601("2016-12-31T23:59:60+01:00"), Error); // 22:59:60 UTC
    // The local fields must be a real date and time of day, whatever the offset.
    EXPECT_THROW((void)Time::fromIso8601("2024-02-30T12:00:00+01:00"), Error);
    EXPECT_THROW((void)Time::fromIso8601("2024-03-01T24:00:00+01:00"), Error);
    EXPECT_THROW((void)Time::fromIso8601("2024-03-01T12:60:00+01:00"), Error);
    // Offset hours are 00-23 and offset minutes 00-59 (RFC 3339, section 5.6).
    EXPECT_THROW((void)Time::fromIso8601("2024-03-01T12:00:00+99:99"), Error);
    EXPECT_THROW((void)Time::fromIso8601("2024-03-01T12:00:00+24:00"), Error);
    EXPECT_THROW((void)Time::fromIso8601("2024-03-01T12:00:00-24:00"), Error);
    EXPECT_THROW((void)Time::fromIso8601("2024-03-01T12:00:00+23:60"), Error);
    EXPECT_THROW((void)Time::fromIso8601("2024-03-01T12:00:00+0560"), Error);
}

TEST(Time, Iso8601RejectsExponentInSeconds) {
    // ISO-8601 seconds are digits with an optional decimal fraction: "05e1" is not 50 s.
    EXPECT_THROW((void)Time::fromIso8601("2024-03-01T12:00:05e1"), Error);
    EXPECT_THROW((void)Time::fromIso8601("2024-03-01T12:00:05e0"), Error);
    EXPECT_THROW((void)Time::fromIso8601("2024-03-01T12:00:05E1Z"), Error);
    EXPECT_THROW((void)Time::fromIso8601("2024-03-01T12:00:05.5e1 TT"), Error);
    // A decimal fraction is still read in full.
    EXPECT_EQ(Time::fromIso8601("2024-03-01T12:00:05.125Z").iso8601(3), "2024-03-01T12:00:05.125Z");
    EXPECT_EQ(Time::fromIso8601("2024-03-01T12:00:05.125 TT").iso8601(3), "2024-03-01T12:00:05.125 TT");
}

TEST(Time, Iso8601TextCoversYears0000To9999) {
    // Inside the four-digit years, what iso8601() writes is read back by fromIso8601() as the same instant.
    for (int year : {0, 1, 999, 1000, 9999}) {
        const Time t    = Time::fromCalendar(TimeScale::TT, year, 12, 31, 23, 59, 58.5);
        const Time back = Time::fromIso8601(t.iso8601(1));
        EXPECT_EQ(back.scale(), TimeScale::TT) << year;
        EXPECT_EQ(back.jd().jd1, t.jd().jd1) << year;
        EXPECT_EQ(back.jd().jd2, t.jd().jd2) << year;
    }
    EXPECT_EQ(Time::fromCalendar(TimeScale::TT, 0, 1, 1).iso8601(0), "0000-01-01T00:00:00 TT");
    EXPECT_EQ(Time::fromCalendar(TimeScale::TT, 9999, 12, 31, 23, 59, 59.0).iso8601(0), "9999-12-31T23:59:59 TT");

    // Outside them there is no four-digit year to write, and both directions refuse.
    EXPECT_THROW((void)Time::fromCalendar(TimeScale::TT, -5, 1, 1).iso8601(), Error);
    EXPECT_THROW((void)Time::fromCalendar(TimeScale::TT, 12000, 1, 1).iso8601(), Error);
    EXPECT_THROW((void)Time::fromCalendar(TimeScale::TT, 9999, 12, 31, 23, 59, 59.9996).iso8601(3), Error); // rounds up
    EXPECT_THROW((void)Time::fromIso8601("-005-01-01T00:00:00 TT"), Error);
    EXPECT_THROW((void)Time::fromIso8601("-0005-01-01T00:00:00 TT"), Error);
    EXPECT_THROW((void)Time::fromIso8601("12000-01-01T00:00:00 TT"), Error);
    EXPECT_THROW((void)Time::fromIso8601("9999-12-31T23:30:00-01:00"), Error); // 00:30 UTC in the year 10000

    // Only the text is limited: the calendar fields of such a time are still available.
    EXPECT_EQ(Time::fromCalendar(TimeScale::TT, 12000, 1, 1).calendar().year, 12000);
    EXPECT_EQ(Time::fromCalendar(TimeScale::TT, -5, 1, 1).calendar().year, -5);
}

TEST(Time, Ut1RequiresEop) {
    const Time t = Time::fromCalendar(TimeScale::UTC, 2020, 1, 1);
    EXPECT_THROW((void)t.to(TimeScale::UT1), Error);
    const ConstantEop eop(EopValues{.dut1 = -0.2});
    EXPECT_NEAR(t.to(TimeScale::UT1, eop).jd().secondsSince(t.jd()), -0.2, 1e-9);
}

TEST(Time, GlobalLeapSecondTableOverride) {
    // Install a table with a hypothetical leap second. The built-in table is put back however the test ends.
    const RestoreBuiltinLeapSecondTable restore;
    std::vector<LeapSecondEntry>        entries = LeapSecondTable::builtin().entries();
    entries.push_back({2030, 1, 38.0});
    setLeapSecondTable(LeapSecondTable(entries));
    const Time t = Time::fromCalendar(TimeScale::UTC, 2029, 12, 31, 23, 59, 60.5);
    EXPECT_NEAR(t.to(TimeScale::TAI).jd().secondsSince(JulianDate::fromCalendar(2030, 1, 1)), 37.5, 1e-9);
    setLeapSecondTable(LeapSecondTable::builtin());
    EXPECT_THROW((void)Time::fromCalendar(TimeScale::UTC, 2029, 12, 31, 23, 59, 60.5), Error);
}

TEST(Time, NonFiniteSecondsThrowAndSaySo) {
    // A NaN or an infinity is refused where it comes in, with a message that names it: it must never reach a
    // conversion to an integer (undefined behaviour) or surface later as an unrelated error.
    for (const TimeScale scale : {TimeScale::UTC, TimeScale::TT}) {
        for (const double bad : {kNaN, kInf, -kInf}) {
            expectErrorSaying([&] { (void)Time::fromCalendar(scale, 2024, 1, 1, 0, 0, bad); },
                              "not finite",
                              "fromCalendar");
            expectErrorSaying([&] { (void)Time::fromCalendar(scale, CalendarDate{2024, 1, 1, 0, 0, bad}); },
                              "not finite",
                              "fromCalendar(CalendarDate)");
        }
    }
    for (const double bad : {kNaN, kInf, -kInf}) {
        expectErrorSaying([&] { (void)Time::fromUnix(bad); }, "not finite", "fromUnix(double)");
        expectErrorSaying([&] { (void)Time::fromUnix(1700000000LL, bad); }, "not finite", "fromUnix(whole, fraction)");
    }
    // Finite, but too large to be turned into whole seconds and days.
    for (const double big : {1.0e30, -1.0e30}) {
        expectErrorSaying([&] { (void)Time::fromUnix(big); }, "out of range", "fromUnix(double)");
        expectErrorSaying([&] { (void)Time::fromUnix(0LL, big); }, "out of range", "fromUnix(whole, fraction)");
    }
    expectErrorSaying([] { (void)Time::fromUnix(std::numeric_limits<long long>::max(), 0.5); },
                      "out of range",
                      "fromUnix(whole, fraction)");
    expectErrorSaying([] { (void)Time::fromUnix(std::numeric_limits<long long>::min(), 0.5); },
                      "out of range",
                      "fromUnix(whole, fraction)");
    // A fraction outside [0, 1) is still carried into the whole seconds: 1483228799 + 1.25 s.
    EXPECT_EQ(Time::fromUnix(1483228799LL, 1.25).iso8601(2), "2017-01-01T00:00:00.25Z");
}

TEST(Time, CalendarOfANonFiniteTimeThrows) {
    // A Time made from a NaN Julian date cannot be written as a date in any scale.
    for (const TimeScale scale : {TimeScale::UTC, TimeScale::TT, TimeScale::GPS}) {
        for (const double bad : {kNaN, kInf, -kInf}) {
            const Time t = Time::fromJd(bad, 0.0, scale);
            expectErrorSaying([&] { (void)t.calendar(); }, "not finite", "calendar");
            expectErrorSaying([&] { (void)t.iso8601(); }, "not finite", "iso8601");
        }
    }
}

TEST(JulianDate, ToCalendarRejectsNonFiniteAndHugeDates) {
    for (const double bad : {kNaN, kInf, -kInf}) {
        expectErrorSaying([&] { (void)JulianDate{bad, 0.0}.calendar(); }, "not finite", "jd1");
        expectErrorSaying([&] { (void)JulianDate{2451545.0, bad}.calendar(); }, "not finite", "jd2");
    }
    // Finite, but beyond the years an int can hold.
    expectErrorSaying([] { (void)JulianDate{1.0e30, 0.0}.calendar(); }, "out of range", "1e30");
    expectErrorSaying([] { (void)JulianDate{-1.0e30, 0.0}.calendar(); }, "out of range", "-1e30");
    // Far dates inside the range are unaffected.
    const CalendarDate far = JulianDate::fromCalendar(99999, 12, 31, 12).calendar(3);
    EXPECT_EQ(far.year, 99999);
    EXPECT_EQ(far.month, 12);
    EXPECT_EQ(far.day, 31);
    EXPECT_EQ(far.hour, 12);
    const CalendarDate old = JulianDate::fromCalendar(-99999, 1, 1).calendar(3);
    EXPECT_EQ(old.year, -99999);
    EXPECT_EQ(old.month, 1);
    EXPECT_EQ(old.day, 1);
}

// Time::fromMjd is the Modified Julian Date in the scale given: MJD 51544.5 is JD 2451545.0 (J2000.0), and
// MJD 57754 is 2017-01-01 (the day of the last leap second, IERS Bulletin C).
TEST(Time, FromMjd) {
    const Time j2000 = Time::fromMjd(51544.5, TimeScale::TT);
    EXPECT_EQ(j2000.scale(), TimeScale::TT);
    EXPECT_EQ(j2000.jd().jd(), constants::jdJ2000);
    EXPECT_EQ(j2000.mjd(), 51544.5);
    EXPECT_EQ(Time::fromMjd(57754.0, TimeScale::UTC).iso8601(3), "2017-01-01T00:00:00.000Z");
    EXPECT_EQ(Time::fromMjd(57754.25, TimeScale::TAI).iso8601(3), "2017-01-01T06:00:00.000 TAI");
}

// calendar() gives the fields in the time's own scale, with the seconds rounded to `decimals` places and
// carried into the next day; in UTC the leap second is second 60.
TEST(Time, CalendarFieldsInTheTimesOwnScale) {
    const auto expectFields = [](const CalendarDate& c, int y, int mo, int d, int h, int mi, double s) {
        EXPECT_EQ(c.year, y);
        EXPECT_EQ(c.month, mo);
        EXPECT_EQ(c.day, d);
        EXPECT_EQ(c.hour, h);
        EXPECT_EQ(c.minute, mi);
        EXPECT_DOUBLE_EQ(c.second, s);
    };
    const Time leap = Time::fromCalendar(TimeScale::UTC, 2016, 12, 31, 23, 59, 60.26);
    expectFields(leap.calendar(), 2016, 12, 31, 23, 59, 60.26);
    expectFields(leap.calendar(1), 2016, 12, 31, 23, 59, 60.3);
    // The same instant in TAI. TAI - UTC was 36 s until the end of the leap second, which therefore began at
    // 00:00:36 TAI on 2017-01-01.
    expectFields(leap.to(TimeScale::TAI).calendar(6), 2017, 1, 1, 0, 0, 36.26);
    // Rounding to whole seconds carries into the next day, in UTC after the leap second.
    expectFields(Time::fromCalendar(TimeScale::UTC, 2016, 12, 31, 23, 59, 60.6).calendar(0), 2017, 1, 1, 0, 0, 0.0);
    expectFields(Time::fromCalendar(TimeScale::TT, 2023, 12, 31, 23, 59, 59.6).calendar(0), 2024, 1, 1, 0, 0, 0.0);
}

TEST(Time, Iso8601WithoutTheScale) {
    const Time tt = Time::fromCalendar(TimeScale::TT, 2000, 1, 1, 12);
    EXPECT_EQ(tt.iso8601(3, false), "2000-01-01T12:00:00.000");
    EXPECT_EQ(tt.iso8601(0, false), "2000-01-01T12:00:00");
    EXPECT_EQ(tt.iso8601(3, true), "2000-01-01T12:00:00.000 TT");
    const Time utc = Time::fromCalendar(TimeScale::UTC, 2016, 12, 31, 23, 59, 60.5);
    EXPECT_EQ(utc.iso8601(1, false), "2016-12-31T23:59:60.5");
}

// Time::now() is the system clock read as UTC: it lies between two readings of the clock taken either side of
// it. (sysTime() rounds to the nanosecond; 1 microsecond allows for the rounding of the day fraction.)
TEST(Time, NowIsTheSystemClockAsUtc) {
    using std::chrono::system_clock;
    const system_clock::time_point before = system_clock::now();
    const Time                     t      = Time::now();
    const system_clock::time_point after  = system_clock::now();
    EXPECT_EQ(t.scale(), TimeScale::UTC);
    EXPECT_GE(t.sysTime(), before - std::chrono::microseconds(1));
    EXPECT_LE(t.sysTime(), after + std::chrono::microseconds(1));
    EXPECT_GE(t.calendar().year, 2026); // this test was written in November 2026
}

// Every scale's short name, in any mix of case; "GPST" is GPS time. Anything else is an error that quotes the
// name.
TEST(Time, TimeScaleFromString) {
    // clang-format off
    const std::pair<const char*, TimeScale> names[] = {
        {"UTC", TimeScale::UTC}, {"tai", TimeScale::TAI}, {"Tt",  TimeScale::TT},  {"tDb",  TimeScale::TDB},
        {"TCG", TimeScale::TCG}, {"tcb", TimeScale::TCB}, {"UT1", TimeScale::UT1}, {"gps",  TimeScale::GPS},
        {"GST", TimeScale::GST}, {"bdt", TimeScale::BDT}, {"GPST", TimeScale::GPS}, {"gpst", TimeScale::GPS},
    };
    // clang-format on
    for (const auto& [name, scale] : names) {
        EXPECT_EQ(parseTimeScale(name), scale) << name;
    }
    for (const char* bad : {"", "XYZ", "UTC ", " TT", "UT", "GLONASS", "TT\n"}) {
        expectErrorSaying([&] { (void)parseTimeScale(bad); }, "'" + std::string(bad) + "'", bad);
    }
}

// += and -= add SI seconds as + and - do: in UTC through TAI, so that the leap second is counted.
TEST(Time, AddAndSubtractInPlace) {
    Time t  = Time::fromCalendar(TimeScale::UTC, 2016, 12, 31, 23, 59, 59.0);
    t      += 1.0;
    EXPECT_EQ(t.iso8601(3), "2016-12-31T23:59:60.000Z");
    t += 1.0;
    EXPECT_EQ(t.iso8601(3), "2017-01-01T00:00:00.000Z");
    t -= 2.0;
    EXPECT_EQ(t.iso8601(3), "2016-12-31T23:59:59.000Z");
    EXPECT_EQ((t - 86400.0).iso8601(3), "2016-12-30T23:59:59.000Z");
    // In TT a day is 86400 s, whatever UTC does.
    Time tt  = Time::fromCalendar(TimeScale::TT, 2016, 12, 31, 12);
    tt      += 86400.0;
    EXPECT_EQ(tt.iso8601(3), "2017-01-01T12:00:00.000 TT");
    tt -= 0.5;
    EXPECT_EQ(tt.iso8601(3), "2017-01-01T11:59:59.500 TT");
}

// Differences and ordering across scales, against the defined offsets: TT = TAI + 32.184 s (IAU 1991
// Resolution A4), GPS = TAI - 19 s (IS-GPS-200), BDT = TAI - 33 s (BeiDou ICD), and UTC = TAI - 37 s in 2020
// (IERS Bulletin C). 1e-9 s is about a hundred times the resolution of a two-part Julian date.
TEST(Time, DifferencesAndOrderingAcrossScales) {
    const Time tai = Time::fromCalendar(TimeScale::TAI, 2020, 6, 1, 12);
    const Time tt  = Time::fromCalendar(TimeScale::TT, 2020, 6, 1, 12, 0, 32.184);
    const Time gps = Time::fromCalendar(TimeScale::GPS, 2020, 6, 1, 11, 59, 41.0);
    const Time bdt = Time::fromCalendar(TimeScale::BDT, 2020, 6, 1, 11, 59, 27.0);
    const Time utc = Time::fromCalendar(TimeScale::UTC, 2020, 6, 1, 11, 59, 23.0);
    for (const Time* t : {&tt, &gps, &bdt, &utc}) {
        EXPECT_NEAR(*t - tai, 0.0, 1e-9) << t->iso8601(6);
        EXPECT_NEAR(tai - *t, 0.0, 1e-9) << t->iso8601(6);
        // A microsecond later in the other scale is later, whichever side it is on.
        EXPECT_TRUE(tai < *t + 1e-6) << t->iso8601(6);
        EXPECT_TRUE(*t + 1e-6 > tai) << t->iso8601(6);
        EXPECT_TRUE(*t - 1e-6 < tai) << t->iso8601(6);
    }
    EXPECT_NEAR(tt - utc, 0.0, 1e-9);
    EXPECT_NEAR(Time::fromCalendar(TimeScale::TT, 2020, 6, 1, 12) - tai, -32.184, 1e-9);
}

// UTC is undefined before 1960-01-01, where the leap-second table starts (time.h). 1960-01-01 is 3653 days
// (ten years, three of them leap years) before 1970-01-01, so Unix time -315619200.
TEST(Time, UtcBefore1960IsAnError) {
    EXPECT_THROW((void)Time::fromCalendar(TimeScale::UTC, 1959, 12, 31, 23, 59, 59.0), Error);
    EXPECT_THROW((void)Time::fromIso8601("1959-12-31T23:59:59Z"), Error);
    EXPECT_THROW((void)Time::fromUnix(-315619200.5), Error);
    EXPECT_EQ(Time::fromUnix(-315619200.0).iso8601(0), "1960-01-01T00:00:00Z");
    const Time first = Time::fromCalendar(TimeScale::UTC, 1960, 1, 1);
    EXPECT_THROW((void)(first - 1.0), Error);
    const Time tt1959 = Time::fromCalendar(TimeScale::TT, 1959, 12, 31);
    EXPECT_THROW((void)tt1959.to(TimeScale::UTC), Error);
    EXPECT_THROW((void)tt1959.unixSeconds(), Error);
    // Other scales have no such limit.
    EXPECT_NEAR(tt1959.to(TimeScale::TAI) - tt1959, 0.0, 1e-9);
}

// time.h: the differences, comparisons and Unix or GPS seconds of a UT1 time need Earth orientation data, so
// they throw unless both sides are UT1.
TEST(Time, Ut1WithoutEopThrowsWhereItMustConvert) {
    const Time ut1 = Time::fromCalendar(TimeScale::UT1, 2020, 1, 1);
    const Time tt  = Time::fromCalendar(TimeScale::TT, 2020, 1, 1);
    for (const auto& call : std::initializer_list<std::function<void()>>{
             [&] { (void)(ut1 - tt); },
             [&] { (void)(tt - ut1); },
             [&] { (void)(ut1 < tt); },
             [&] { (void)ut1.unixSeconds(); },
             [&] { (void)ut1.gpsSeconds(); },
             [&] { (void)ut1.gnssWeekTow(GnssSystem::GPS); },
         }) {
        expectErrorSaying(call, "Earth orientation", "UT1");
    }
    EXPECT_NEAR((ut1 + 1.0) - ut1, 1.0, 1e-9); // both UT1: no conversion
}

// Before 1970 the Unix time is negative, and the time point is before the system clock's epoch.
TEST(Time, SysTimeBefore1970) {
    using namespace std::chrono;
    const Time t = Time::fromIso8601("1969-12-31T23:59:59.5Z");
    EXPECT_EQ(duration_cast<nanoseconds>(t.sysTime().time_since_epoch()).count(), -500000000LL);
    EXPECT_NEAR(t.unixSeconds(), -0.5, 1e-9);
    EXPECT_EQ(Time::fromSysTime(t.sysTime()).iso8601(3), "1969-12-31T23:59:59.500Z");
    EXPECT_EQ(
        duration_cast<seconds>(Time::fromCalendar(TimeScale::UTC, 1960, 1, 1).sysTime().time_since_epoch()).count(),
        -315619200LL);
}

// time.h: resolveWeekRollover throws for `bits` outside 1-30. Within it, the truncated week is resolved to the
// full week nearest the reference. GPS week 2048 began on 2019-04-07 (JD 2444244.5 + 2048 x 7 = 2458580.5), the
// second rollover of the 10-bit week of the legacy navigation message, so the reference 2019-04-10 is in week
// 2048 and Galileo week 1024.
TEST(Time, ResolveWeekRollover) {
    const Time ref = Time::fromCalendar(TimeScale::UTC, 2019, 4, 10);
    EXPECT_EQ(resolveWeekRollover(2, 10, GnssSystem::GPS, ref), 2050);    // 2, 1026 or 2050: 2050 is nearest
    EXPECT_EQ(resolveWeekRollover(1023, 10, GnssSystem::GPS, ref), 2047); // the week before the rollover
    EXPECT_EQ(resolveWeekRollover(2050, 13, GnssSystem::GPS, ref), 2050); // 13 bits (CNAV): no rollover yet
    EXPECT_EQ(resolveWeekRollover(3, 2, GnssSystem::GPS, ref), 2047);     // 2 bits: 2047 or 2051; 2047 is nearer
    EXPECT_EQ(resolveWeekRollover(2, 10, GnssSystem::Galileo, ref), 1026);
    for (int bits : {0, -1, 31, 32}) {
        expectErrorSaying([&] { (void)resolveWeekRollover(2, bits, GnssSystem::GPS, ref); },
                          "bits",
                          ("bits " + std::to_string(bits)).c_str());
    }
}

// Each field out of range is refused; in scales other than UTC a minute has 60 seconds.
TEST(Time, FromCalendarRejectsInvalidFields) {
    // clang-format off
    struct Case { int month, day, hour, minute; double second; };
    const Case cases[] = {
        {13,  1,  0,  0,  0.0 }, {0,  1,  0,  0,  0.0 }, {2, 30,  0,  0,  0.0 }, {4, 31,  0,  0,  0.0 },
        { 1,  0,  0,  0,  0.0 }, {1,  1, 24,  0,  0.0 }, {1,  1, -1,  0,  0.0 }, {1,  1,  0, 60,  0.0 },
        { 1,  1,  0, -1,  0.0 }, {1,  1,  0,  0, -0.5 }, {1,  1,  0,  0, 60.0 }, {1,  1, 12, 30, 61.0 },
    };
    // clang-format on
    for (const TimeScale scale : {TimeScale::UTC, TimeScale::TT, TimeScale::GPS}) {
        for (const Case& c : cases) {
            EXPECT_THROW((void)Time::fromCalendar(scale, 2024, c.month, c.day, c.hour, c.minute, c.second), Error)
                << toString(scale) << " " << c.month << "-" << c.day << " " << c.hour << ":" << c.minute << ":"
                << c.second;
        }
    }
    // The largest valid fields.
    EXPECT_NO_THROW((void)Time::fromCalendar(TimeScale::TT, 2024, 2, 29, 23, 59, std::nextafter(60.0, 0.0)));
}
