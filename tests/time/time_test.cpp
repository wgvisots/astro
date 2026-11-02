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
} // namespace

TEST(JulianDate, CalendarRoundTrip) {
    EXPECT_DOUBLE_EQ(JulianDate::fromCalendar(2000, 1, 1, 12).jd(), 2451545.0);
    EXPECT_DOUBLE_EQ(JulianDate::fromCalendar(1858, 11, 17).mjd(), 0.0);
    EXPECT_DOUBLE_EQ(JulianDate::fromCalendar(-4713, 11, 24, 12).jd(), 0.0); // proleptic Gregorian JD 0
    const CalendarDate c = JulianDate::fromCalendar(2024, 2, 29, 23, 59, 59.123456).toCalendar(6);
    EXPECT_EQ(c.year, 2024);
    EXPECT_EQ(c.month, 2);
    EXPECT_EQ(c.day, 29);
    EXPECT_EQ(c.hour, 23);
    EXPECT_EQ(c.minute, 59);
    EXPECT_NEAR(c.second, 59.123456, 1e-9);
    // Rounding carries into the next day.
    const CalendarDate n = JulianDate::fromCalendar(2023, 12, 31, 23, 59, 59.9999999).toCalendar(3);
    EXPECT_EQ(n.year, 2024);
    EXPECT_EQ(n.month, 1);
    EXPECT_EQ(n.day, 1);
    EXPECT_EQ(n.hour, 0);
}

TEST(JulianDate, Helpers) {
    EXPECT_EQ(dayOfYear(2024, 12, 31), 366);
    EXPECT_EQ(dayOfYear(2023, 3, 1), 60);
    int m, d;
    monthDayFromDayOfYear(2024, 60, m, d);
    EXPECT_EQ(m, 2);
    EXPECT_EQ(d, 29);
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
                                          "    61000.0    1  1 2026       38\n";
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

TEST(TimeScales, MatchesErfa) {
    const auto csv = readCsv(dataPath("time_scales.csv"));
    ASSERT_GT(csv.rows.size(), 300u);
    const ConstantEop eop(EopValues{.dut1 = 0.1234567});
    auto num = [&](const std::vector<std::string>& r, const char* name) { return std::stod(r[csv.col(name)]); };
    for (const auto& r : csv.rows) {
        const Time        utc   = Time::fromCalendar(TimeScale::UTC,
                                            std::stoi(r[0]),
                                            std::stoi(r[1]),
                                            std::stoi(r[2]),
                                            std::stoi(r[3]),
                                            std::stoi(r[4]),
                                            num(r, "second"));
        const std::string label = utc.iso8601(6);
        EXPECT_NEAR(diffSec(utc.jd(), num(r, "utc1"), num(r, "utc2")), 0.0, 1e-9) << label;
        EXPECT_NEAR(diffSec(utc.to(TimeScale::TAI).jd(), num(r, "tai1"), num(r, "tai2")), 0.0, 1e-9) << label;
        EXPECT_NEAR(diffSec(utc.to(TimeScale::TT).jd(), num(r, "tt1"), num(r, "tt2")), 0.0, 1e-9) << label;
        EXPECT_NEAR(diffSec(utc.to(TimeScale::TDB).jd(), num(r, "tdb1"), num(r, "tdb2")), 0.0, 1e-9) << label;
        EXPECT_NEAR(diffSec(utc.to(TimeScale::TCG).jd(), num(r, "tcg1"), num(r, "tcg2")), 0.0, 1e-9) << label;
        EXPECT_NEAR(diffSec(utc.to(TimeScale::TCB).jd(), num(r, "tcb1"), num(r, "tcb2")), 0.0, 1e-9) << label;
        const int  year         = std::stoi(r[0]);
        const bool inLeapSecond = num(r, "second") >= 60.0;
        if (year >= 1972) {
            EXPECT_NEAR(diffSec(utc.to(TimeScale::UT1, eop).jd(), num(r, "ut11"), num(r, "ut12")), 0.0, 1e-9) << label;
        } else {
            // ERFA's utcut1 evaluates TAI-UTC at 0h of the day, which is off by up to ~1 ms/day during the
            // 1961-1971 drift era. We evaluate it at the instant, so UT1 - UTC == dut1 exactly.
            // (Compare against the plain time of day, since the quasi-JD fraction is scaled on step days.)
            const JulianDate plain = JulianDate::fromCalendar(year,
                                                              std::stoi(r[1]),
                                                              std::stoi(r[2]),
                                                              std::stoi(r[3]),
                                                              std::stoi(r[4]),
                                                              num(r, "second"));
            EXPECT_NEAR(utc.to(TimeScale::UT1, eop).jd().secondsSince(plain), 0.1234567, 1e-8) << label;
        }
        EXPECT_NEAR(tdbMinusTt(utc.to(TimeScale::TT).jd()), num(r, "tdb_minus_tt"), 1e-12) << label;

        // Round trips back to UTC from every scale. (UT1 is skipped inside a leap second: with a constant
        // dut1, which real EOP never has across a leap second, UT1 is discontinuous there.)
        for (TimeScale s : {TimeScale::TAI,
                            TimeScale::TT,
                            TimeScale::TDB,
                            TimeScale::TCG,
                            TimeScale::TCB,
                            TimeScale::GPS,
                            TimeScale::BDT,
                            TimeScale::UT1}) {
            if (s == TimeScale::UT1 && inLeapSecond) {
                continue;
            }
            const Time back = utc.to(s, eop).to(TimeScale::UTC, eop);
            EXPECT_NEAR(back - utc, 0.0, 1e-9) << label << " via " << toString(s);
            EXPECT_EQ(back.iso8601(6), label) << " via " << toString(s);
        }
    }
}

TEST(TimeScales, TopocentricTdb) {
    const auto csv = readCsv(dataPath("dtdb_topo.csv"));
    for (const auto& r : csv.rows) {
        EXPECT_NEAR(tdbMinusTt(JulianDate{std::stod(r[0]), 0.0},
                               std::stod(r[1]),
                               std::stod(r[2]),
                               std::stod(r[3]),
                               std::stod(r[4])),
                    std::stod(r[5]),
                    1e-13);
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
    const auto tp = p.toSysTime();
    EXPECT_EQ(Time::fromSysTime(tp).iso8601(6), "2023-11-14T22:13:20.123457Z");
}

TEST(Time, Gnss) {
    // GPS epoch.
    const Time e = Time::fromCalendar(TimeScale::UTC, 1980, 1, 6);
    EXPECT_NEAR(e.gpsSeconds(), 0.0, 1e-9);
    // 2017-01-01T00:00:00Z: GPS = UTC + 18 s; GPS week 1930, TOW = 0 + 18 (Sunday).
    const Time        t = Time::fromCalendar(TimeScale::UTC, 2017, 1, 1);
    const GnssWeekTow g = t.gnssWeekTow(GnssSystem::GPS);
    EXPECT_EQ(g.week, 1930);
    EXPECT_NEAR(g.tow, 18.0, 1e-9);
    const GnssWeekTow gal = t.gnssWeekTow(GnssSystem::Galileo);
    EXPECT_EQ(gal.week, 1930 - 1024);
    EXPECT_NEAR(gal.tow, 18.0, 1e-9);
    // BeiDou: BDT = UTC + 4 s in 2017; BDT week 573 starts 2016-12-25.
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

TEST(Time, Ut1RequiresEop) {
    const Time t = Time::fromCalendar(TimeScale::UTC, 2020, 1, 1);
    EXPECT_THROW((void)t.to(TimeScale::UT1), Error);
    const ConstantEop eop(EopValues{.dut1 = -0.2});
    EXPECT_NEAR(t.to(TimeScale::UT1, eop).jd().secondsSince(t.jd()), -0.2, 1e-9);
}

TEST(Time, GlobalLeapSecondTableOverride) {
    // Install a table with a hypothetical leap second, then restore.
    std::vector<LeapSecondEntry> entries = LeapSecondTable::builtin().entries();
    entries.push_back({2030, 1, 38.0});
    setLeapSecondTable(LeapSecondTable(entries));
    const Time t = Time::fromCalendar(TimeScale::UTC, 2029, 12, 31, 23, 59, 60.5);
    EXPECT_NEAR(t.to(TimeScale::TAI).jd().secondsSince(JulianDate::fromCalendar(2030, 1, 1)), 37.5, 1e-9);
    setLeapSecondTable(LeapSecondTable::builtin());
    EXPECT_THROW((void)Time::fromCalendar(TimeScale::UTC, 2029, 12, 31, 23, 59, 60.5), Error);
}
