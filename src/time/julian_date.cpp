/**********************************************************************
 * julian_date.cpp
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

#include <cmath>
#include <string>

#include "astro/core/error.h"

namespace astro {

namespace {
// Floor division for possibly-negative integers.
constexpr long long floorDiv(long long a, long long b) noexcept {
    long long q = a / b;
    if ((a % b != 0) && ((a < 0) != (b < 0))) --q;
    return q;
}
} // namespace

long long julianDayNumber(int year, int month, int day) noexcept {
    // Proleptic Gregorian calendar, valid for all years (floor division).
    const long long a = floorDiv(14 - month, 12);
    const long long y = static_cast<long long>(year) + 4800 - a;
    const long long m = month + 12 * a - 3;
    return day + floorDiv(153 * m + 2, 5) + 365 * y + floorDiv(y, 4) - floorDiv(y, 100) + floorDiv(y, 400) - 32045;
}

void dateFromJulianDayNumber(long long jdn, int& year, int& month, int& day) noexcept {
    const long long a = jdn + 32044;
    const long long b = floorDiv(4 * a + 3, 146097);
    const long long c = a - floorDiv(146097 * b, 4);
    const long long d = floorDiv(4 * c + 3, 1461);
    const long long e = c - floorDiv(1461 * d, 4);
    const long long m = floorDiv(5 * e + 2, 153);
    day = static_cast<int>(e - floorDiv(153 * m + 2, 5) + 1);
    month = static_cast<int>(m + 3 - 12 * floorDiv(m, 10));
    year = static_cast<int>(100 * b + d - 4800 + floorDiv(m, 10));
}

int daysInMonth(int year, int month) {
    static constexpr int dim[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    if (month < 1 || month > 12) throw Error("daysInMonth: month out of range: " + std::to_string(month));
    return (month == 2 && isLeapYear(year)) ? 29 : dim[month - 1];
}

int dayOfYear(int year, int month, int day) {
    return static_cast<int>(julianDayNumber(year, month, day) - julianDayNumber(year, 1, 1)) + 1;
}

void monthDayFromDayOfYear(int year, int doy, int& month, int& day) {
    const int ndays = isLeapYear(year) ? 366 : 365;
    if (doy < 1 || doy > ndays) throw Error("day of year out of range: " + std::to_string(doy));
    int y;
    dateFromJulianDayNumber(julianDayNumber(year, 1, 1) + doy - 1, y, month, day);
}

int dayOfWeek(int year, int month, int day) {
    // JDN 0 was a Monday; JDN + 1 mod 7 gives 0 = Sunday.
    return static_cast<int>(((julianDayNumber(year, month, day) + 1) % 7 + 7) % 7);
}

JulianDate JulianDate::fromCalendar(int year, int month, int day, int hour, int minute, double second) {
    if (month < 1 || month > 12) throw Error("JulianDate: month out of range");
    if (day < 1 || day > daysInMonth(year, month)) throw Error("JulianDate: day out of range");
    const double jd1 = static_cast<double>(julianDayNumber(year, month, day)) - 0.5;
    const double jd2 = (hour * 3600.0 + minute * 60.0 + second) / constants::secondsPerDay;
    return JulianDate{jd1, jd2}.normalized();
}

JulianDate JulianDate::fromMjd(double mjd) noexcept {
    const double whole = std::floor(mjd);
    return JulianDate{whole + constants::mjdOffset, mjd - whole};
}

JulianDate JulianDate::normalized() const noexcept {
    // Split both parts into integer-midnight and fraction, then recombine.
    const double a1 = std::floor(jd1 - 0.5);
    const double f1 = (jd1 - 0.5) - a1;
    const double a2 = std::floor(jd2);
    const double f2 = jd2 - a2;
    double whole = a1 + a2;
    double frac = f1 + f2;
    const double carry = std::floor(frac);
    whole += carry;
    frac -= carry;
    if (frac >= 1.0) {  // guard against rounding to exactly 1
        whole += 1.0;
        frac = 0.0;
    }
    return JulianDate{whole + 0.5, frac};
}

JulianDate JulianDate::addDays(double days) const noexcept {
    const double whole = std::floor(days);
    return JulianDate{jd1 + whole, jd2 + (days - whole)}.normalized();
}

JulianDate JulianDate::addSeconds(double seconds) const noexcept {
    // Split off whole days first so large offsets keep sub-microsecond precision.
    const double wholeDays = std::floor(seconds / constants::secondsPerDay);
    const double rem = seconds - wholeDays * constants::secondsPerDay;
    return JulianDate{jd1 + wholeDays, jd2 + rem / constants::secondsPerDay}.normalized();
}

CalendarDate JulianDate::toCalendar(int decimals) const {
    const JulianDate n = normalized();
    const double scale = std::pow(10.0, decimals);
    const double ticksPerDay = constants::secondsPerDay * scale;
    double ticks = std::round(n.jd2 * ticksPerDay);
    long long jdn = static_cast<long long>(n.jd1 + 0.5);
    if (ticks >= ticksPerDay) {
        ticks -= ticksPerDay;
        ++jdn;
    }
    CalendarDate c;
    dateFromJulianDayNumber(jdn, c.year, c.month, c.day);
    const double secs = ticks / scale;
    c.hour = static_cast<int>(secs / 3600.0);
    c.minute = static_cast<int>((secs - c.hour * 3600.0) / 60.0);
    c.second = (ticks - (c.hour * 3600.0 + c.minute * 60.0) * scale) / scale;
    return c;
}

} // namespace astro
