/**********************************************************************
 * julian_date_test.cpp
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

#include "astro/time/julian_date.h"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <string>

#include "astro/core/error.h"

using namespace astro;

namespace {

// `call` must throw astro::Error, and the message must contain `what`.
template <typename Call>
void expectErrorSaying(Call&& call, const std::string& what, const std::string& label) {
    try {
        call();
        ADD_FAILURE() << label << ": no error";
    } catch (const Error& e) {
        EXPECT_NE(std::string(e.what()).find(what), std::string::npos)
            << label << ": expected \"" << what << "\" in the message, which was: " << e.what();
    }
}

// The Gregorian rules, written out here so that the tests below do not lean on the library's isLeapYear.
bool gregorianLeap(int year) {
    return (year % 4 == 0 && year % 100 != 0) || year % 400 == 0;
}

int gregorianMonthLength(int year, int month) {
    constexpr int lengths[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    return month == 2 && gregorianLeap(year) ? 29 : lengths[month - 1];
}

} // namespace

// Julian day numbers (the JD at noon of the date) of dates published as Julian dates. JD 2451545.0 is J2000.0,
// noon on 2000-01-01; MJD 0 is JD 2400000.5, the midnight that starts 1858-11-17; the Unix epoch 1970-01-01 is
// JD 2440587.5 and the GPS epoch 1980-01-06 is JD 2444244.5 at 0h; the Gregorian calendar began on 1582-10-15,
// JD 2299160.5 at 0h (Meeus, Astronomical Algorithms, chapter 7). JD 0 is noon on 24 November 4714 BC in the
// proleptic Gregorian calendar, year -4713 in astronomical numbering. 2024-02-29: 2024-01-01 is 24 years of 365
// days and 6 leap days after 2000-01-01, so JDN 2451545 + 8766 = 2460311, and 31 + 28 days later 2460370.
TEST(JulianDayNumber, PublishedDates) {
    // clang-format off
    struct Case { int year, month, day; long long jdn; };
    const Case cases[] = {
        { 2000,  1,  1, 2451545},
        { 1858, 11, 17, 2400001},
        { 1970,  1,  1, 2440588},
        { 1980,  1,  6, 2444245},
        { 1582, 10, 15, 2299161},
        {-4713, 11, 24,       0},
        { 2024,  2, 29, 2460370},
    };
    // clang-format on
    for (const Case& c : cases) {
        EXPECT_EQ(julianDayNumber(c.year, c.month, c.day), c.jdn) << c.year << "-" << c.month << "-" << c.day;
        int y = 0, m = 0, d = 0;
        dateFromJulianDayNumber(c.jdn, y, m, d);
        EXPECT_EQ(y, c.year) << c.jdn;
        EXPECT_EQ(m, c.month) << c.jdn;
        EXPECT_EQ(d, c.day) << c.jdn;
    }
}

// Brute force: walk from JD 0 (-4713-11-24) one day at a time to the year 2200, stepping the date with the
// Gregorian rules written out above. Each date must have the next day number, and the day number must give
// the date back. 2.5 million days cover every kind of year and the years before 0. (Takes about 35 ms.)
TEST(JulianDayNumber, EveryDayFromJdZeroToTheYear2200) {
    int       year = -4713, month = 11, day = 24;
    long long wrong = 0;
    for (long long jdn = 0; year < 2200; ++jdn) {
        int y = 0, m = 0, d = 0;
        dateFromJulianDayNumber(jdn, y, m, d);
        if (julianDayNumber(year, month, day) != jdn || y != year || m != month || d != day) {
            if (++wrong <= 5) {
                ADD_FAILURE() << "JDN " << jdn << " should be " << year << "-" << month << "-" << day << "; got "
                              << julianDayNumber(year, month, day) << " and " << y << "-" << m << "-" << d;
            }
        }
        if (++day > gregorianMonthLength(year, month)) {
            day = 1;
            if (++month > 12) {
                month = 1;
                ++year;
            }
        }
    }
    EXPECT_EQ(wrong, 0);
}

TEST(Calendar, DaysInMonth) {
    const int common[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    for (int month = 1; month <= 12; ++month) {
        EXPECT_EQ(daysInMonth(2023, month), common[month - 1]) << month;
        EXPECT_EQ(daysInMonth(2024, month), month == 2 ? 29 : common[month - 1]) << month;
    }
    // Century years are leap years only when divisible by 400.
    EXPECT_EQ(daysInMonth(1900, 2), 28);
    EXPECT_EQ(daysInMonth(2000, 2), 29);
    EXPECT_EQ(daysInMonth(2100, 2), 28);
    EXPECT_EQ(daysInMonth(-4, 2), 29); // astronomical year -4 is 5 BC, a leap year
    for (int month : {0, 13, -1}) {
        expectErrorSaying([&] { (void)daysInMonth(2024, month); }, "month", "month " + std::to_string(month));
    }
}

// Day of year both ways, for every day of a common and a leap year, against a count kept by hand.
TEST(Calendar, DayOfYearBothWays) {
    for (int year : {2023, 2024, 1900, 2000}) {
        int doy = 0;
        for (int month = 1; month <= 12; ++month) {
            for (int day = 1; day <= gregorianMonthLength(year, month); ++day) {
                ++doy;
                EXPECT_EQ(dayOfYear(year, month, day), doy) << year << "-" << month << "-" << day;
                int m = 0, d = 0;
                monthDayFromDayOfYear(year, doy, m, d);
                EXPECT_EQ(m, month) << year << " day " << doy;
                EXPECT_EQ(d, day) << year << " day " << doy;
            }
        }
        // julian_date.h: throws for a day of year outside 1-365 (1-366 in a leap year).
        const int last = gregorianLeap(year) ? 366 : 365;
        for (int bad : {0, -1, last + 1}) {
            int m = 0, d = 0;
            expectErrorSaying([&] { monthDayFromDayOfYear(year, bad, m, d); },
                              "day of year",
                              std::to_string(year) + " day " + std::to_string(bad));
        }
    }
}

// Days of the week that are known: JD 0 was a Monday; the Gregorian calendar began on Friday 15 October 1582;
// 1 January 1970 was a Thursday; 29 February 2024 was a Thursday.
TEST(Calendar, DayOfWeek) {
    EXPECT_EQ(dayOfWeek(-4713, 11, 24), 1);
    EXPECT_EQ(dayOfWeek(1582, 10, 15), 5);
    EXPECT_EQ(dayOfWeek(1970, 1, 1), 4);
    EXPECT_EQ(dayOfWeek(2024, 2, 29), 4);
    EXPECT_EQ(dayOfWeek(2024, 3, 3), 0); // the Sunday after
}

TEST(JulianDate, FromCalendarRejectsInvalidDates) {
    expectErrorSaying([] { (void)JulianDate::fromCalendar(2024, 0, 1); }, "month", "month 0");
    expectErrorSaying([] { (void)JulianDate::fromCalendar(2024, 13, 1); }, "month", "month 13");
    expectErrorSaying([] { (void)JulianDate::fromCalendar(2024, 4, 31); }, "day", "April 31");
    expectErrorSaying([] { (void)JulianDate::fromCalendar(2024, 1, 0); }, "day", "day 0");
    expectErrorSaying([] { (void)JulianDate::fromCalendar(2023, 2, 29); }, "day", "29 February 2023");
    expectErrorSaying([] { (void)JulianDate::fromCalendar(2100, 2, 29); }, "day", "29 February 2100");
    EXPECT_EQ(JulianDate::fromCalendar(2000, 2, 29).jd(), 2451603.5); // 2451544.5 + 31 + 28
    // julian_date.h: the time of day is not checked, so 25 h on 1 January is 1 h on 2 January (to a rounding
    // error of the fraction 25/24 of a day, 2.2e-16 day).
    EXPECT_NEAR(JulianDate::fromCalendar(2024, 1, 1, 25).daysSince(JulianDate::fromCalendar(2024, 1, 2, 1)),
                0.0,
                2.3e-16);
}

// The halves and quarters below are exact in binary, so the parts are compared for equality.
TEST(JulianDate, FromMjdAndTheDerivedDayCounts) {
    const JulianDate j2000 = JulianDate::fromMjd(51544.5); // MJD 51544.5 = JD 2451545.0
    EXPECT_EQ(j2000.jd1, 2451544.5);
    EXPECT_EQ(j2000.jd2, 0.5);
    EXPECT_EQ(j2000.jd(), 2451545.0);
    EXPECT_EQ(j2000.mjd(), 51544.5);
    EXPECT_EQ(j2000.daysSinceJ2000(), 0.0);
    const JulianDate before = JulianDate::fromMjd(-0.25); // a quarter of a day before MJD 0
    EXPECT_EQ(before.jd1, 2399999.5);
    EXPECT_EQ(before.jd2, 0.75);
    EXPECT_EQ(before.mjd(), -0.25);
    const JulianDate century = JulianDate{2451545.0 + 36525.0, 0.25};
    EXPECT_EQ(century.daysSinceJ2000(), 36525.25);
    EXPECT_EQ(JulianDate(2451545.0 + 36525.0).centuriesSinceJ2000(), 1.0);
    EXPECT_EQ(JulianDate(2451545.0 - 73050.0).centuriesSinceJ2000(), -2.0);
}

TEST(JulianDate, NormalizedPutsMidnightInTheFirstPart) {
    // clang-format off
    struct Case { double jd1, jd2, wantJd1, wantJd2; };
    const Case cases[] = {
        {2451545.0,   0.0,   2451544.5, 0.5  },
        {2451545.25,  0.875, 2451545.5, 0.625}, // 2451546.125 in all
        {2451545.5,  -0.25,  2451544.5, 0.75 },
        {2451545.5,   3.5,   2451548.5, 0.5  },
        {0.0,        -0.5,   -0.5,      0.0  }, // JD -0.5 is a midnight
    };
    // clang-format on
    for (const Case& c : cases) {
        const JulianDate n = JulianDate{c.jd1, c.jd2}.normalized();
        EXPECT_EQ(n.jd1, c.wantJd1) << c.jd1 << " + " << c.jd2;
        EXPECT_EQ(n.jd2, c.wantJd2) << c.jd1 << " + " << c.jd2;
    }
}

TEST(JulianDate, AddDaysAndDaysSince) {
    const JulianDate a = JulianDate{2451544.5, 0.75};
    const JulianDate b = a.addDays(1.5);
    // 2451544.5 + 0.75 + 1.5 = 2451546.75, which normalized is the midnight 2451546.5 and a quarter day.
    EXPECT_EQ(b.jd1, 2451546.5);
    EXPECT_EQ(b.jd2, 0.25);
    // 2451545.25 - 2.25 = 2451543.0, the midnight 2451542.5 and half a day.
    const JulianDate c = a.addDays(-2.25);
    EXPECT_EQ(c.jd1, 2451542.5);
    EXPECT_EQ(c.jd2, 0.5);
    EXPECT_EQ(b.daysSince(a), 1.5);
    EXPECT_EQ(c.daysSince(a), -2.25);
    // The two parts are subtracted separately, so a nanosecond difference between dates 2.4 million days from
    // JD 0 survives: one day plus 1e-9 day, against long double. (A single JD would keep only 4.7e-10 day.)
    const JulianDate p{2451545.0, 1e-9}, q{2451544.0, 0.0};
    EXPECT_NEAR(p.daysSince(q), static_cast<double>(1.0L + 1e-9L), 1.2e-16);
    EXPECT_NEAR(p.secondsSince(q), static_cast<double>((1.0L + 1e-9L) * 86400.0L), 1.2e-16 * 86400.0);
}

// Comparison is by elapsed days, whatever the split between the two parts.
TEST(JulianDate, ComparisonIgnoresTheSplit) {
    const JulianDate a{2451545.0, 0.25}, sameAsA{2451545.25, 0.0}, later{2451545.0, 0.25 + 1e-12};
    EXPECT_TRUE(a == sameAsA);
    EXPECT_TRUE(a < later);
    EXPECT_TRUE(later > sameAsA);
    EXPECT_FALSE(later == a);
    // A date that is not a number is unordered with every date.
    const JulianDate nan{std::numeric_limits<double>::quiet_NaN(), 0.0};
    EXPECT_FALSE(nan < a);
    EXPECT_FALSE(a < nan);
    EXPECT_FALSE(nan == a);
    EXPECT_EQ(nan <=> a, std::partial_ordering::unordered);
}
