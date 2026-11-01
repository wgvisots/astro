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
#include <string>
#include <string_view>

#include "astro/time/julian_date.h"

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

[[nodiscard]] std::string_view toString(TimeScale scale) noexcept;
/// Parse "UTC", "TAI", "TT", ... (case-insensitive). Throws astro::Error.
[[nodiscard]] TimeScale timeScaleFromString(std::string_view name);

enum class GnssSystem { GPS, Galileo, BeiDou };

struct GnssWeekTow {
    int week{0};      ///< full (non-rolled-over) week number
    double tow{0.0};  ///< seconds of week, [0, 604800)
};

/// An instant tagged with its time scale.
///
/// Internally a normalized two-part Julian Date in the given scale. For UTC
/// the fraction follows the SOFA "quasi-JD" convention: on a day containing a
/// leap second the fraction is (seconds of day) / 86401, so 23:59:60.5 is
/// representable and round-trips.
///
/// Conversions to and from UT1 need Earth orientation data (UT1-UTC) and throw
/// astro::Error if no EopProvider is given.
class Time {
public:
    /// J2000.0 TT.
    Time() = default;
    Time(const JulianDate& jd, TimeScale scale);

    [[nodiscard]] static Time fromJd(double jd1, double jd2, TimeScale scale);
    [[nodiscard]] static Time fromMjd(double mjd, TimeScale scale);
    /// From calendar fields in the given scale. For UTC, second may be in
    /// [60, 61) on a day that ends in a leap second.
    [[nodiscard]] static Time fromCalendar(TimeScale scale, int year, int month, int day, int hour = 0,
                                           int minute = 0, double second = 0.0);
    [[nodiscard]] static Time fromCalendar(TimeScale scale, const CalendarDate& cal);
    /// From POSIX/Unix time (seconds since 1970-01-01T00:00:00 UTC, leap seconds not counted). Result is UTC.
    [[nodiscard]] static Time fromUnix(double unixSeconds);
    /// From whole + fractional Unix seconds (full precision). Result is UTC.
    [[nodiscard]] static Time fromUnix(long long wholeSeconds, double fraction);
    [[nodiscard]] static Time fromSysTime(std::chrono::system_clock::time_point tp);
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
    [[nodiscard]] static Time fromIso8601(std::string_view text, TimeScale defaultScale = TimeScale::UTC);

    [[nodiscard]] const JulianDate& jd() const noexcept { return jd_; }
    [[nodiscard]] TimeScale scale() const noexcept { return scale_; }
    [[nodiscard]] double mjd() const noexcept { return jd_.mjd(); }
    /// Julian centuries since J2000.0, in this time's own scale.
    [[nodiscard]] double centuriesSinceJ2000() const noexcept { return jd_.centuriesSinceJ2000(); }

    /// Convert to another scale. `eop` is required for UT1.
    [[nodiscard]] Time to(TimeScale target, const EopProvider* eop = nullptr) const;
    [[nodiscard]] Time to(TimeScale target, const EopProvider& eop) const { return to(target, &eop); }

    /// Calendar date in this scale (UTC: leap-second aware, second may be 60.x).
    [[nodiscard]] CalendarDate calendar(int decimals = 9) const;
    /// ISO-8601 text. UTC ends with 'Z'; other scales get a " TT"-style suffix
    /// unless `withScale` is false.
    [[nodiscard]] std::string iso8601(int decimals = 3, bool withScale = true) const;

    /// Unix seconds (converts to UTC). During a leap second, returns the value
    /// of the following second (POSIX time repeats).
    [[nodiscard]] double unixSeconds() const;
    [[nodiscard]] std::chrono::system_clock::time_point toSysTime() const;
    /// Seconds since the GPS epoch (converts to GPS time).
    [[nodiscard]] double gpsSeconds() const;
    /// Week number and time of week in the given GNSS system's scale.
    [[nodiscard]] GnssWeekTow gnssWeekTow(GnssSystem system) const;

    /// Add SI seconds. UTC is advanced through TAI so leap seconds are counted;
    /// other scales advance their Julian Date directly.
    [[nodiscard]] Time operator+(double seconds) const;
    [[nodiscard]] Time operator-(double seconds) const { return *this + (-seconds); }
    Time& operator+=(double seconds) { return *this = *this + seconds; }
    Time& operator-=(double seconds) { return *this = *this + (-seconds); }

    /// Elapsed seconds a - b, measured in a's scale (b is converted first; UTC
    /// differences are taken in TAI so leap seconds are counted).
    friend double operator-(const Time& a, const Time& b);
    /// Ordering after converting b to a's scale.
    friend std::partial_ordering operator<=>(const Time& a, const Time& b);
    friend bool operator==(const Time& a, const Time& b) { return (a <=> b) == 0; }

private:
    JulianDate jd_{constants::jdJ2000, 0.0};
    TimeScale scale_{TimeScale::TT};
};

/// TDB - TT in seconds (Fairhead & Bretagnon 1990, full 787-term series plus
/// topocentric terms; ~3 ns over 1950-2050). `tt` is TT (or TDB; the difference
/// is negligible). For the geocentric value leave the observer terms at zero.
/// @param ut1Fraction  UT1 fraction of day (only needed for topocentric terms)
/// @param elong        observer east longitude, rad
/// @param uKm          observer distance from Earth's spin axis, km
/// @param vKm          observer distance north of the equatorial plane, km
[[nodiscard]] double tdbMinusTt(const JulianDate& tt, double ut1Fraction = 0.0, double elong = 0.0,
                                double uKm = 0.0, double vKm = 0.0) noexcept;

/// Resolve a truncated (rolled-over) GNSS week number, e.g. 10-bit GPS LNAV
/// weeks, to the full week closest to `reference`.
[[nodiscard]] int resolveWeekRollover(int truncatedWeek, int bits, GnssSystem system, const Time& reference);

} // namespace astro
