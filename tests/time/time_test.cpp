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
    // Install a table with a hypothetical leap second, then restore.
    std::vector<LeapSecondEntry> entries = LeapSecondTable::builtin().entries();
    entries.push_back({2030, 1, 38.0});
    setLeapSecondTable(LeapSecondTable(entries));
    const Time t = Time::fromCalendar(TimeScale::UTC, 2029, 12, 31, 23, 59, 60.5);
    EXPECT_NEAR(t.to(TimeScale::TAI).jd().secondsSince(JulianDate::fromCalendar(2030, 1, 1)), 37.5, 1e-9);
    setLeapSecondTable(LeapSecondTable::builtin());
    EXPECT_THROW((void)Time::fromCalendar(TimeScale::UTC, 2029, 12, 31, 23, 59, 60.5), Error);
}
