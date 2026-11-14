/**********************************************************************
 * time.h
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

#include <chrono>
#include <compare>
#include <iosfwd>
#include <string>
#include <string_view>

#include <boost/describe/enum.hpp>

#include "astro/time/julian_date.h"

/// \file
/// astro::Time, an instant in a stated time scale, and conversions between UTC, TAI, TT, TDB, TCG, TCB, UT1, GPS, Galileo and BeiDou time, Unix time and GNSS weeks.
namespace astro {

class EopProvider;

/// Supported time scales.
enum class TimeScale {
    UTC, ///< Coordinated Universal Time (with leap seconds)
    TAI, ///< International Atomic Time
    TT,  ///< Terrestrial Time = TAI + 32.184 s
    TDB, ///< Barycentric Dynamical Time
    TCG, ///< Geocentric Coordinate Time
    TCB, ///< Barycentric Coordinate Time
    UT1, ///< Universal Time (Earth rotation); requires EOP (UT1-UTC)
    GPS, ///< GPS Time = TAI - 19 s
    GST, ///< Galileo System Time (same as GPS time)
    BDT, ///< BeiDou Time = TAI - 33 s
};

BOOST_DESCRIBE_ENUM(TimeScale, UTC, TAI, TT, TDB, TCG, TCB, UT1, GPS, GST, BDT)

/// Writes "TimeScale::<enumerator>", or "TimeScale::INVALID_VALUE" for any other value.
std::ostream& operator<<(std::ostream& os, TimeScale value);

/// Short name of a time scale ("UTC", "TT", ...), or "?" for an invalid value.
[[nodiscard]] std::string_view toString(TimeScale scale) noexcept;
/// Parse "UTC", "TAI", "TT", ... (case-insensitive). Throws astro::Error for an unknown name.
[[nodiscard]] TimeScale parseTimeScale(std::string_view name);

/// Global navigation satellite systems with their own time scale and week count.
enum class GnssSystem {
    GPS,     ///< GPS (GPS time, weeks from 1980-01-06)
    Galileo, ///< Galileo (Galileo System Time, GST)
    BeiDou   ///< BeiDou (BeiDou Time, BDT)
};

BOOST_DESCRIBE_ENUM(GnssSystem, GPS, Galileo, BeiDou)

/// Writes "GnssSystem::<enumerator>", or "GnssSystem::INVALID_VALUE" for any other value.
std::ostream& operator<<(std::ostream& os, GnssSystem value);

/// GNSS week number and time of week, in the system's own time scale.
struct GnssWeekTow {
    int    week{0};  ///< full (non-rolled-over) week number
    double tow{0.0}; ///< seconds of week, [0, 604800)
};

/// An instant tagged with its time scale.
///
/// Internally a normalized two-part Julian Date in the given scale. For UTC
/// the fraction follows the "quasi-JD" convention of SOFA and ERFA: on a day containing a
/// leap second the fraction is (seconds of day) / 86401, so 23:59:60.5 is
/// representable and round-trips.
///
/// Conversions to and from UT1 need Earth orientation data (UT1-UTC) and throw
/// astro::Error if no EopProvider is given.
///
/// UTC is undefined before 1960-01-01, where the leap-second table starts. A
/// member that has to convert to or from UTC before then throws astro::Error;
/// each such member says so below.
class Time {
public:
    /// J2000.0 TT.
    Time() = default;
    /// From a two-part Julian Date in `scale` (normalized on construction).
    Time(const JulianDate& jd, TimeScale scale);

    /// From a two-part Julian Date (jd1 + jd2 days) in `scale`.
    [[nodiscard]] static Time fromJd(double jd1, double jd2, TimeScale scale);
    /// From a Modified Julian Date (JD - 2400000.5) in `scale`.
    [[nodiscard]] static Time fromMjd(double mjd, TimeScale scale);
    /// From calendar fields in the given scale. For UTC, second may be in
    /// [60, 61) on a day that ends in a leap second; on a day that ends early
    /// (a negative leap second, or the shortened days 1961-07-31 and 1968-01-31)
    /// the last minute is that much shorter. Throws astro::Error for an invalid
    /// date or time of day, which includes a second that is not finite (NaN or
    /// infinity) and a second at or after the end of the UTC day. Also throws
    /// for a UTC date before 1960.
    /// @{
    [[nodiscard]] static Time
    fromCalendar(TimeScale scale, int year, int month, int day, int hour = 0, int minute = 0, double second = 0.0);
    [[nodiscard]] static Time fromCalendar(TimeScale scale, const CalendarDate& cal);
    /// @}
    /// From POSIX/Unix time (seconds since 1970-01-01T00:00:00 UTC, leap seconds not counted). Result is UTC.
    /// Throws astro::Error if `unixSeconds` is not finite (NaN or infinity) or is more than 1e9 days from 1970.
    /// Also throws for a time before 1960-01-01 (`unixSeconds` below -315619200), where UTC is undefined.
    [[nodiscard]] static Time fromUnix(double unixSeconds);
    /// From whole + fractional Unix seconds (full precision). Result is UTC. `fraction` need not be in
    /// [0, 1): its whole part is carried into the seconds. Throws astro::Error if `fraction` is not finite,
    /// or if either argument is more than 1e9 days from 1970. Also throws for a time before 1960-01-01.
    [[nodiscard]] static Time fromUnix(long long wholeSeconds, double fraction);
    /// From a std::chrono::system_clock time point, read as Unix time
    /// (nanosecond resolution). Result is UTC. Throws astro::Error for a time
    /// point before 1960-01-01.
    [[nodiscard]] static Time fromSysTime(std::chrono::system_clock::time_point tp);
    /// The current system clock time, as UTC. Throws astro::Error only if the
    /// system clock is set to a date before 1960 (see fromSysTime).
    [[nodiscard]] static Time now();
    /// From seconds since the GPS epoch (1980-01-06T00:00:00 UTC). Result is GPS.
    [[nodiscard]] static Time fromGpsSeconds(double seconds);
    /// From a GNSS week number (full, not rolled over) and time of week. Result
    /// is in the system's own scale (GPS, GST or BDT).
    [[nodiscard]] static Time fromGnssWeekTow(GnssSystem system, int week, double tow);
    /// Parse an ISO-8601 date/time:
    ///   "2024-03-01", "2024-03-01T12:34:56.789", "2024-061T12:00:00" (ordinal),
    ///   with optional "Z" (UTC) or "+hh:mm"/"-hh:mm" (UTC offset), or a trailing
    ///   scale name such as " TT" or " TAI". Space is accepted in place of 'T'.
    ///   Without a zone/scale the `defaultScale` is used.
    ///
    /// A UTC offset (hours 00-23, minutes 00-59) is taken off the calendar fields,
    /// as on a civil clock, and not as elapsed seconds: "2017-01-01T00:30:00+01:00"
    /// is 2016-12-31T23:30:00 UTC although a leap second lies in between, and
    /// "2017-01-01T00:59:60+01:00" is the leap second 2016-12-31T23:59:60 UTC.
    /// Seconds are digits with an optional decimal fraction; an exponent ("05e1")
    /// is not accepted. The year has four digits, 0000 to 9999 (also after an
    /// offset has been taken off), the same range iso8601() writes.
    /// Throws astro::Error for malformed text, for an invalid date, time of day
    /// or offset, for an unknown scale name, and for a UTC time before 1960.
    [[nodiscard]] static Time fromIso8601(std::string_view text, TimeScale defaultScale = TimeScale::UTC);

    /// Normalized two-part Julian Date in this time's own scale.
    [[nodiscard]] const JulianDate& jd() const noexcept { return jd_; }
    /// The time scale this instant is expressed in.
    [[nodiscard]] TimeScale scale() const noexcept { return scale_; }
    /// Modified Julian Date in this time's own scale.
    [[nodiscard]] double mjd() const noexcept { return jd_.mjd(); }
    /// Julian centuries since J2000.0, in this time's own scale.
    [[nodiscard]] double centuriesSinceJ2000() const noexcept { return jd_.centuriesSinceJ2000(); }

    /// Convert to another scale. `eop` is required for UT1 (to or from); without
    /// it a UT1 conversion throws astro::Error. So does a UT1 conversion for
    /// which `eop` has no data at this time, and a conversion to or from UTC or
    /// UT1 (which goes through UTC) before 1960.
    /// @{
    [[nodiscard]] Time to(TimeScale target, const EopProvider* eop = nullptr) const;
    [[nodiscard]] Time to(TimeScale target, const EopProvider& eop) const { return to(target, &eop); }
    /// @}

    /// Calendar date in this scale (UTC: leap-second aware, second may be 60.x).
    /// Throws astro::Error if this time is not finite (it was made from a NaN or
    /// an infinity), if it is a UTC time before 1960, or, in the other scales,
    /// if its Julian date is larger than 1e11 days in size.
    [[nodiscard]] CalendarDate calendar(int decimals = 9) const;
    /// ISO-8601 text. UTC ends with 'Z'; other scales get a " TT"-style suffix
    /// unless `withScale` is false. fromIso8601() reads the text back.
    /// Throws astro::Error if the year (after rounding the seconds to `decimals`
    /// places) is outside 0000 to 9999: ISO-8601 text has a four-digit year.
    /// calendar() has no such limit. Also throws as calendar() does.
    [[nodiscard]] std::string iso8601(int decimals = 3, bool withScale = true) const;

    /// Unix seconds (converts to UTC). During a leap second, returns the value
    /// of the following second (POSIX time repeats). Throws astro::Error as
    /// to(TimeScale::UTC) does: for a time before 1960, and for a UT1 time
    /// (convert that with to() and an EopProvider first).
    [[nodiscard]] double unixSeconds() const;
    /// std::chrono::system_clock time point from the Unix time (converts to
    /// UTC), rounded to the nearest nanosecond. Throws astro::Error as
    /// unixSeconds() does.
    [[nodiscard]] std::chrono::system_clock::time_point toSysTime() const;
    /// Seconds since the GPS epoch (converts to GPS time). Throws astro::Error
    /// as to(TimeScale::GPS) does: for a UTC time before 1960, and for a UT1
    /// time (convert that with to() and an EopProvider first).
    [[nodiscard]] double gpsSeconds() const;
    /// Week number and time of week in the given GNSS system's scale. The time
    /// of week is always in [0, 604800): an instant so close to the end of a
    /// week that its time of week rounds to 604800 s is returned as time of
    /// week 0 of the next week. Throws astro::Error as gpsSeconds() does.
    [[nodiscard]] GnssWeekTow gnssWeekTow(GnssSystem system) const;

    /// Add / subtract SI seconds. UTC is advanced through TAI so leap seconds
    /// are counted; other scales advance their Julian Date directly. Throws
    /// astro::Error for a UTC time if it, or the result, is before 1960.
    /// @{
    [[nodiscard]] Time operator+(double seconds) const;
    [[nodiscard]] Time operator-(double seconds) const { return *this + (-seconds); }
    Time&              operator+=(double seconds) { return *this = *this + seconds; }
    Time&              operator-=(double seconds) { return *this = *this + (-seconds); }
    /// @}

    /// Elapsed seconds a - b, measured in a's scale (b is converted first; UTC
    /// differences are taken in TAI so leap seconds are counted). Scales that tick at
    /// different rates (TDB, TCG, TCB against the others) therefore give a slightly
    /// different number for b - a than for -(a - b): up to 1.5e-8 of the interval.
    /// Throws astro::Error if exactly one of the two is UT1 (the conversion needs
    /// Earth orientation data: use to() with an EopProvider first), or if a
    /// conversion to or from UTC falls before 1960.
    friend double operator-(const Time& a, const Time& b);
    /// Ordering after converting b to a's scale; throws as operator- does. Equality is
    /// exact: two times that are the same instant in different scales can compare
    /// unequal by a rounding error of the conversion (about 1 in 30 do). Compare
    /// std::abs(a - b) with a tolerance instead.
    friend std::partial_ordering operator<=>(const Time& a, const Time& b);
    friend bool                  operator==(const Time& a, const Time& b) { return (a <=> b) == 0; }

private:
    JulianDate jd_{constants::jdJ2000, 0.0};
    TimeScale  scale_{TimeScale::TT};
};

/// TDB - TT in seconds: the model of Fairhead & Bretagnon (1990) in its full form (787 terms)
/// with the topocentric terms of Moyer (1981) and Murray (1983), as distributed in ERFA's
/// eraDtdb, which this function follows; ~3 ns over 1950-2050. `tt` is TT (or TDB; the difference
/// is negligible). For the geocentric value leave the observer terms at zero.
/// @param ut1Fraction  UT1 fraction of day (only needed for topocentric terms)
/// @param elong        observer east longitude, rad
/// @param uKm          observer distance from Earth's spin axis, km
/// @param vKm          observer distance north of the equatorial plane, km
[[nodiscard]] double tdbMinusTt(const JulianDate& tt,
                                double            ut1Fraction = 0.0,
                                double            elong       = 0.0,
                                double            uKm         = 0.0,
                                double            vKm         = 0.0) noexcept;

/// Resolve a truncated (rolled-over) GNSS week number, e.g. 10-bit GPS LNAV
/// weeks, to the full week closest to `reference`. Throws astro::Error if
/// `bits` is outside 1-30, or as Time::gnssWeekTow does for `reference`.
[[nodiscard]] int resolveWeekRollover(int truncatedWeek, int bits, GnssSystem system, const Time& reference);

} // namespace astro
