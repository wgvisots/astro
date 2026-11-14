/**********************************************************************
 * julian_date.h
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

#pragma once

#include <compare>

#include "astro/core/constants.h"

/// \file
/// Two-part Julian dates and calendar conversions.
namespace astro {

/// Gregorian calendar date and time of day. `second` may reach 60.x on a UTC
/// leap-second day.
struct CalendarDate {
    int    year{2000};  ///< Gregorian year
    int    month{1};    ///< month, 1-12
    int    day{1};      ///< day of month, 1-31
    int    hour{0};     ///< hour, 0-23
    int    minute{0};   ///< minute, 0-59
    double second{0.0}; ///< seconds, [0, 60), or [60, 61) during a UTC leap second

    friend constexpr bool operator==(const CalendarDate&, const CalendarDate&) noexcept = default;
};

/// Two-part Julian Date for full double precision over long spans.
///
/// The value is jd1 + jd2 days. After normalize(), jd1 is a midnight
/// (integer + 0.5) and jd2 is the fraction of the day in [0, 1), which gives
/// ~10 picosecond resolution. The scale (UTC, TT, ...) is not stored here; see
/// astro::Time.
struct JulianDate {
    double jd1{constants::jdJ2000}; ///< first part, days (a midnight, integer + 0.5, after normalized())
    double jd2{0.0};                ///< second part, days (fraction of day in [0, 1) after normalized())

    /// J2000.0 (JD 2451545.0).
    constexpr JulianDate() noexcept = default;
    /// From the two parts as given (not normalized), in days.
    constexpr JulianDate(double jd1_, double jd2_ = 0.0) noexcept
        : jd1(jd1_)
        , jd2(jd2_) {}

    /// From Gregorian calendar date + time of day, assuming 86400-second days
    /// (i.e. not leap-second aware; use astro::Time for UTC). Throws
    /// astro::Error if `month` is outside 1-12 or `day` is not a day of that
    /// month. The time of day is not checked.
    [[nodiscard]] static JulianDate
    fromCalendar(int year, int month, int day, int hour = 0, int minute = 0, double second = 0.0);
    /// From a Modified Julian Date (JD - 2400000.5), split into whole day and fraction.
    [[nodiscard]] static JulianDate fromMjd(double mjd) noexcept;

    /// Single-double JD (loses precision; ~20 microseconds today).
    [[nodiscard]] constexpr double jd() const noexcept { return jd1 + jd2; }
    /// Modified Julian Date (JD - 2400000.5) as a single double.
    [[nodiscard]] constexpr double mjd() const noexcept { return (jd1 - constants::mjdOffset) + jd2; }
    /// Days since J2000.0 (2451545.0) in this scale.
    [[nodiscard]] constexpr double daysSinceJ2000() const noexcept { return (jd1 - constants::jdJ2000) + jd2; }
    /// Julian centuries since J2000.0 in this scale.
    [[nodiscard]] constexpr double centuriesSinceJ2000() const noexcept {
        return daysSinceJ2000() / constants::daysPerJulianCentury;
    }

    /// jd1 = midnight (integer + 0.5), jd2 in [0, 1).
    [[nodiscard]] JulianDate normalized() const noexcept;

    /// Shifted by a number of days / seconds (86400-second days); the result
    /// is normalized.
    /// @{
    [[nodiscard]] JulianDate addDays(double days) const noexcept;
    [[nodiscard]] JulianDate addSeconds(double seconds) const noexcept;
    /// @}
    /// (this - other) in days / seconds, computed without cancellation.
    /// @{
    [[nodiscard]] constexpr double daysSince(const JulianDate& o) const noexcept {
        return (jd1 - o.jd1) + (jd2 - o.jd2);
    }
    [[nodiscard]] constexpr double secondsSince(const JulianDate& o) const noexcept {
        return daysSince(o) * constants::secondsPerDay;
    }
    /// @}

    /// Calendar date assuming 86400-second days. Seconds are rounded to
    /// `decimals` places with carries into minutes/hours/days. Throws
    /// astro::Error if the date is not finite (NaN or infinity), or if either
    /// part is larger than 1e11 days in size.
    [[nodiscard]] CalendarDate toCalendar(int decimals = 9) const;

    friend constexpr std::partial_ordering operator<=>(const JulianDate& a, const JulianDate& b) noexcept {
        return a.daysSince(b) <=> 0.0;
    }
    friend constexpr bool operator==(const JulianDate& a, const JulianDate& b) noexcept {
        return a.daysSince(b) == 0.0;
    }
};

/// True for a Gregorian leap year.
[[nodiscard]] constexpr bool isLeapYear(int year) noexcept {
    return (year % 4 == 0 && year % 100 != 0) || year % 400 == 0;
}
/// Number of days in a month (1-12) of a Gregorian year. Throws astro::Error
/// if `month` is out of range.
[[nodiscard]] int daysInMonth(int year, int month);
/// Day of year, 1-based (Jan 1 = 1).
[[nodiscard]] int dayOfYear(int year, int month, int day);
/// Month and day from year and 1-based day of year. Throws astro::Error if
/// `doy` is outside 1-365 (1-366 in a leap year).
void dayOfYearToMonthDay(int year, int doy, int& month, int& day);
/// Day of week, 0 = Sunday ... 6 = Saturday.
[[nodiscard]] int dayOfWeek(int year, int month, int day);

/// Julian Day Number (integer, noon-based) for a Gregorian date.
[[nodiscard]] long long julianDayNumber(int year, int month, int day) noexcept;
/// Gregorian date from a Julian Day Number.
void julianDayNumberToDate(long long jdn, int& year, int& month, int& day) noexcept;

} // namespace astro
