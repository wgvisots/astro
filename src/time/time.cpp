/**********************************************************************
 * time.cpp
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

#include <algorithm>
#include <array>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstdio>

#include "../core/enum_io.h"
#include "../data/fairhead_tdb.h"
#include "astro/core/error.h"
#include "astro/frames/eop.h"
#include "astro/time/leap_seconds.h"

namespace astro {

using constants::secondsPerDay;

namespace {

constexpr double kJdBdtEpoch        = 2453736.5; // 2006-01-01 00:00:00 BDT
constexpr int    kGalileoWeekOffset = 1024;      // GST week 0 == GPS week 1024

// --------------------------------------------------------------- UTC helpers
//
// UTC to and from TAI follows ERFA's eraUtctai and eraTaiutc, including the "quasi-JD" day fraction
// on a day with a leap second and the rate changes before 1972. See THIRD_PARTY_NOTICES.md.

struct UtcDay {
    int    year, month, day;
    double dat0;  // TAI-UTC at 0h
    double dlod;  // pre-1972 rate: extra SI seconds per UTC day
    double dleap; // leap second at the end of this day (s), normally 0
};

UtcDay utcDay(double jd1Midnight) {
    UtcDay     u{};
    const auto jdn = static_cast<long long>(std::floor(jd1Midnight + 0.5));
    dateFromJulianDayNumber(jdn, u.year, u.month, u.day);
    int y2, m2, d2;
    dateFromJulianDayNumber(jdn + 1, y2, m2, d2);
    const auto table   = leapSecondTable();
    u.dat0             = table->taiMinusUtc(u.year, u.month, u.day, 0.0);
    const double dat12 = table->taiMinusUtc(u.year, u.month, u.day, 0.5);
    const double dat24 = table->taiMinusUtc(y2, m2, d2, 0.0);
    u.dlod             = 2.0 * (dat12 - u.dat0);
    u.dleap            = dat24 - (2.0 * dat12 - u.dat0);
    return u;
}

JulianDate utcToTai(const JulianDate& utcIn) {
    const JulianDate utc  = utcIn.normalized();
    const UtcDay     u    = utcDay(utc.jd1);
    double           fd   = utc.jd2;
    fd                   *= (secondsPerDay + u.dleap) / secondsPerDay; // undo quasi-JD leap-day scaling
    fd                   *= (secondsPerDay + u.dlod) / secondsPerDay;  // pre-1972 UTC seconds -> SI seconds
    return JulianDate{utc.jd1, fd + u.dat0 / secondsPerDay}.normalized();
}

JulianDate taiToUtc(const JulianDate& tai) {
    JulianDate u = tai.normalized();
    for (int i = 0; i < 4; ++i) {
        const JulianDate g = utcToTai(u);
        u                  = JulianDate{u.jd1, u.jd2 + ((tai.jd1 - g.jd1) + (tai.jd2 - g.jd2))}.normalized();
    }
    return u;
}

double utcDayLength(double jd1Midnight) {
    return secondsPerDay + utcDay(jd1Midnight).dleap;
}

// --------------------------------------------------------------- relativistic scales

JulianDate add(const JulianDate& jd, double days) {
    return JulianDate{jd.jd1, jd.jd2 + days}.normalized();
}

// Days since T0 = 1977-01-01T00:00:32.184 (as TT/TCG/TDB/TCB JD 2443144.5003725).
double daysSinceT0(const JulianDate& jd) {
    return (jd.jd1 - 2443144.5) + (jd.jd2 - 0.0003725);
}

JulianDate ttToTcg(const JulianDate& tt) {
    constexpr double elg = constants::LG / (1.0 - constants::LG);
    return add(tt, daysSinceT0(tt) * elg);
}
JulianDate tcgToTt(const JulianDate& tcg) {
    return add(tcg, -daysSinceT0(tcg) * constants::LG);
}
JulianDate tdbToTcb(const JulianDate& tdb) {
    constexpr double elb  = constants::LB / (1.0 - constants::LB);
    constexpr double tdb0 = constants::TDB0 / secondsPerDay;
    return add(tdb, -tdb0 / (1.0 - constants::LB) + daysSinceT0(tdb) * elb);
}
JulianDate tcbToTdb(const JulianDate& tcb) {
    constexpr double tdb0 = constants::TDB0 / secondsPerDay;
    return add(tcb, -daysSinceT0(tcb) * constants::LB + tdb0);
}
JulianDate ttToTdb(const JulianDate& tt) {
    return add(tt, tdbMinusTt(tt) / secondsPerDay);
}
JulianDate tdbToTt(const JulianDate& tdb) {
    return add(tdb, -tdbMinusTt(tdb) / secondsPerDay);
}

const EopProvider& requireEop(const EopProvider* eop) {
    if (eop == nullptr) {
        throw Error("UT1 conversion requires Earth orientation data (EopProvider)");
    }
    return *eop;
}

JulianDate taiToUt1(const JulianDate& tai, const EopProvider& eop) {
    const JulianDate utc   = taiToUtc(tai);
    const UtcDay     u     = utcDay(utc.jd1);
    const auto       table = leapSecondTable();
    const double     dat   = table->taiMinusUtc(u.year, u.month, u.day, utc.jd2);
    const double     dut1  = eop.at(utc.mjd()).dut1;
    return add(tai, (dut1 - dat) / secondsPerDay);
}

JulianDate toTai(const JulianDate& jd, TimeScale scale, const EopProvider* eop) {
    switch (scale) {
    case TimeScale::TAI: return jd;
    case TimeScale::UTC: return utcToTai(jd);
    case TimeScale::TT:  return add(jd, -constants::ttMinusTai / secondsPerDay);
    case TimeScale::GPS:
    case TimeScale::GST: return add(jd, constants::taiMinusGps / secondsPerDay);
    case TimeScale::BDT: return add(jd, constants::taiMinusBdt / secondsPerDay);
    case TimeScale::TCG: return add(tcgToTt(jd), -constants::ttMinusTai / secondsPerDay);
    case TimeScale::TDB: return add(tdbToTt(jd), -constants::ttMinusTai / secondsPerDay);
    case TimeScale::TCB: return add(tdbToTt(tcbToTdb(jd)), -constants::ttMinusTai / secondsPerDay);
    case TimeScale::UT1: {
        const EopProvider& e   = requireEop(eop);
        JulianDate         tai = add(jd, 37.0 / secondsPerDay);
        for (int i = 0; i < 4; ++i) {
            const JulianDate u = taiToUt1(tai, e);
            tai                = add(tai, jd.daysSince(u));
        }
        return tai;
    }
    }
    throw Error("unknown time scale");
}

JulianDate fromTai(const JulianDate& tai, TimeScale scale, const EopProvider* eop) {
    switch (scale) {
    case TimeScale::TAI: return tai;
    case TimeScale::UTC: return taiToUtc(tai);
    case TimeScale::TT:  return add(tai, constants::ttMinusTai / secondsPerDay);
    case TimeScale::GPS:
    case TimeScale::GST: return add(tai, -constants::taiMinusGps / secondsPerDay);
    case TimeScale::BDT: return add(tai, -constants::taiMinusBdt / secondsPerDay);
    case TimeScale::TCG: return ttToTcg(add(tai, constants::ttMinusTai / secondsPerDay));
    case TimeScale::TDB: return ttToTdb(add(tai, constants::ttMinusTai / secondsPerDay));
    case TimeScale::TCB: return tdbToTcb(ttToTdb(add(tai, constants::ttMinusTai / secondsPerDay)));
    case TimeScale::UT1: return taiToUt1(tai, requireEop(eop));
    }
    throw Error("unknown time scale");
}

double gnssEpochJd(GnssSystem system) {
    return system == GnssSystem::BeiDou ? kJdBdtEpoch : constants::jdGpsEpoch;
}

TimeScale gnssScale(GnssSystem system) {
    switch (system) {
    case GnssSystem::GPS:     return TimeScale::GPS;
    case GnssSystem::Galileo: return TimeScale::GST;
    case GnssSystem::BeiDou:  return TimeScale::BDT;
    }
    return TimeScale::GPS;
}

// --------------------------------------------------------------- ISO-8601 parsing

struct Cursor {
    std::string_view          s;
    std::size_t               p{0};
    [[nodiscard]] bool        done() const { return p >= s.size(); }
    [[nodiscard]] char        peek() const { return done() ? '\0' : s[p]; }
    [[nodiscard]] std::size_t digitsAhead() const {
        std::size_t n = 0;
        while (p + n < s.size() && std::isdigit(static_cast<unsigned char>(s[p + n]))) {
            ++n;
        }
        return n;
    }
    int integer(std::size_t ndigits) {
        if (digitsAhead() < ndigits) {
            fail();
        }
        int v = 0;
        for (std::size_t i = 0; i < ndigits; ++i) {
            v = v * 10 + (s[p++] - '0');
        }
        return v;
    }
    void expect(char c) {
        if (peek() != c) {
            fail();
        }
        ++p;
    }
    [[noreturn]] void fail() const { throw Error("malformed ISO-8601 time: '" + std::string(s) + "'"); }
};

} // namespace

// ------------------------------------------------------------------ names

std::string_view toString(TimeScale scale) noexcept {
    return detail::enumName(scale, "?");
}

TimeScale timeScaleFromString(std::string_view name) {
    std::string up;
    for (char c : name) {
        up.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
    }
    static constexpr std::array all{TimeScale::UTC,
                                    TimeScale::TAI,
                                    TimeScale::TT,
                                    TimeScale::TDB,
                                    TimeScale::TCG,
                                    TimeScale::TCB,
                                    TimeScale::UT1,
                                    TimeScale::GPS,
                                    TimeScale::GST,
                                    TimeScale::BDT};
    for (TimeScale s : all) {
        if (up == toString(s)) {
            return s;
        }
    }
    if (up == "GPST") {
        return TimeScale::GPS;
    }
    throw Error("unknown time scale '" + std::string(name) + "'");
}

// ------------------------------------------------------------------ Time

Time::Time(const JulianDate& jd, TimeScale scale)
    : jd_(jd.normalized())
    , scale_(scale) {}

Time Time::fromJd(double jd1, double jd2, TimeScale scale) {
    return Time(JulianDate{jd1, jd2}, scale);
}

Time Time::fromMjd(double mjd, TimeScale scale) {
    return Time(JulianDate::fromMjd(mjd), scale);
}

Time Time::fromCalendar(TimeScale scale, int year, int month, int day, int hour, int minute, double second) {
    if (month < 1 || month > 12 || day < 1 || day > daysInMonth(year, month)) {
        throw Error("Time::fromCalendar: invalid date");
    }
    if (hour < 0 || hour > 23 || minute < 0 || minute > 59 || second < 0.0) {
        throw Error("Time::fromCalendar: invalid time of day");
    }
    const double jd1 = static_cast<double>(julianDayNumber(year, month, day)) - 0.5;
    const double sod = hour * 3600.0 + minute * 60.0 + second;
    if (scale != TimeScale::UTC) {
        if (second >= 60.0) {
            throw Error("Time::fromCalendar: seconds must be < 60");
        }
        return Time(JulianDate{jd1, sod / secondsPerDay}, scale);
    }
    const double dleap    = utcDay(jd1).dleap;
    const double secLimit = (hour == 23 && minute == 59) ? 60.0 + std::max(dleap, 0.0) : 60.0;
    if (second >= secLimit) {
        throw Error("Time::fromCalendar: seconds out of range for this UTC day");
    }
    return Time(JulianDate{jd1, sod / (secondsPerDay + dleap)}, scale);
}

Time Time::fromCalendar(TimeScale scale, const CalendarDate& c) {
    return fromCalendar(scale, c.year, c.month, c.day, c.hour, c.minute, c.second);
}

Time Time::fromUnix(long long wholeSeconds, double fraction) {
    const double fWhole  = std::floor(fraction);
    wholeSeconds        += static_cast<long long>(fWhole);
    fraction            -= fWhole;
    long long days       = wholeSeconds / 86400;
    long long sodWhole   = wholeSeconds % 86400;
    if (sodWhole < 0) {
        sodWhole += 86400;
        --days;
    }
    const double jd1 = constants::jdUnixEpoch + static_cast<double>(days);
    const double sod = static_cast<double>(sodWhole) + fraction;
    return Time(JulianDate{jd1, sod / utcDayLength(jd1)}, TimeScale::UTC);
}

Time Time::fromUnix(double unixSeconds) {
    const double whole = std::floor(unixSeconds);
    return fromUnix(static_cast<long long>(whole), unixSeconds - whole);
}

Time Time::fromSysTime(std::chrono::system_clock::time_point tp) {
    using namespace std::chrono;
    const auto ns    = duration_cast<nanoseconds>(tp.time_since_epoch()).count();
    long long  whole = ns / 1000000000LL;
    long long  rem   = ns % 1000000000LL;
    if (rem < 0) {
        rem += 1000000000LL;
        --whole;
    }
    return fromUnix(whole, static_cast<double>(rem) * 1e-9);
}

Time Time::now() {
    return fromSysTime(std::chrono::system_clock::now());
}

Time Time::fromGpsSeconds(double seconds) {
    return Time(JulianDate{constants::jdGpsEpoch, 0.0}.addSeconds(seconds), TimeScale::GPS);
}

Time Time::fromGnssWeekTow(GnssSystem system, int week, double tow) {
    const int        w  = system == GnssSystem::Galileo ? week + kGalileoWeekOffset : week;
    const JulianDate jd = JulianDate{gnssEpochJd(system) + 7.0 * w, 0.0}.addSeconds(tow);
    return Time(jd, gnssScale(system));
}

Time Time::fromIso8601(std::string_view text, TimeScale defaultScale) {
    while (!text.empty() && std::isspace(static_cast<unsigned char>(text.front()))) {
        text.remove_prefix(1);
    }
    while (!text.empty() && std::isspace(static_cast<unsigned char>(text.back()))) {
        text.remove_suffix(1);
    }
    Cursor    c{text};
    const int year = c.integer(4);
    c.expect('-');
    int month = 1, day = 1;
    if (c.digitsAhead() == 3) {
        monthDayFromDayOfYear(year, c.integer(3), month, day);
    } else {
        month = c.integer(2);
        c.expect('-');
        day = c.integer(2);
    }
    int    hour = 0, minute = 0;
    double second = 0.0;
    if (c.peek() == 'T'
        || (c.peek() == ' ' && c.p + 1 < text.size() && std::isdigit(static_cast<unsigned char>(text[c.p + 1])))) {
        ++c.p;
        hour = c.integer(2);
        c.expect(':');
        minute = c.integer(2);
        if (c.peek() == ':') {
            ++c.p;
            const char* b = text.data() + c.p;
            const char* e = text.data() + text.size();
            if (c.digitsAhead() < 2) {
                c.fail();
            }
            auto [ptr, ec] = std::from_chars(b, e, second);
            if (ec != std::errc{}) {
                c.fail();
            }
            c.p = static_cast<std::size_t>(ptr - text.data());
        }
    }
    TimeScale scale         = defaultScale;
    double    offsetSeconds = 0.0;
    while (c.peek() == ' ') {
        ++c.p;
    }
    if (!c.done()) {
        const char z = c.peek();
        if ((z == 'Z' || z == 'z') && c.p + 1 == text.size()) {
            scale = TimeScale::UTC;
            ++c.p;
        } else if (z == '+' || z == '-') {
            ++c.p;
            const int oh = c.integer(2);
            int       om = 0;
            if (c.peek() == ':') {
                ++c.p;
            }
            if (!c.done()) {
                om = c.integer(2);
            }
            offsetSeconds = (z == '-' ? -1.0 : 1.0) * (oh * 3600.0 + om * 60.0);
            scale         = TimeScale::UTC;
        } else {
            scale = timeScaleFromString(text.substr(c.p));
            c.p   = text.size();
        }
    }
    if (!c.done()) {
        c.fail();
    }
    Time t = fromCalendar(scale, year, month, day, hour, minute, second);
    return offsetSeconds != 0.0 ? t + (-offsetSeconds) : t;
}

Time Time::to(TimeScale target, const EopProvider* eop) const {
    if (target == scale_) {
        return *this;
    }
    return Time(fromTai(toTai(jd_, scale_, eop), target, eop), target);
}

CalendarDate Time::calendar(int decimals) const {
    if (scale_ != TimeScale::UTC) {
        return jd_.toCalendar(decimals);
    }
    const double dayLen = utcDayLength(jd_.jd1);
    const double scale  = std::pow(10.0, decimals);
    double       ticks  = std::round(jd_.jd2 * dayLen * scale);
    auto         jdn    = static_cast<long long>(std::floor(jd_.jd1 + 0.5));
    if (ticks >= std::round(dayLen * scale)) {
        ticks = 0.0;
        ++jdn;
    }
    CalendarDate cal;
    dateFromJulianDayNumber(jdn, cal.year, cal.month, cal.day);
    const double tPerHour = 3600.0 * scale, tPerMin = 60.0 * scale;
    cal.hour    = std::min(23, static_cast<int>(std::floor(ticks / tPerHour)));
    ticks      -= cal.hour * tPerHour;
    cal.minute  = std::min(59, static_cast<int>(std::floor(ticks / tPerMin)));
    ticks      -= cal.minute * tPerMin;
    cal.second  = ticks / scale;
    return cal;
}

std::string Time::iso8601(int decimals, bool withScale) const {
    const CalendarDate c = calendar(decimals);
    char               sec[32];
    if (decimals > 0) {
        std::snprintf(sec, sizeof(sec), "%0*.*f", decimals + 3, decimals, c.second);
    } else {
        std::snprintf(sec, sizeof(sec), "%02d", static_cast<int>(c.second));
    }
    char buf[96];
    std::snprintf(buf, sizeof(buf), "%04d-%02d-%02dT%02d:%02d:%s", c.year, c.month, c.day, c.hour, c.minute, sec);
    std::string out(buf);
    if (withScale) {
        if (scale_ == TimeScale::UTC) {
            out += 'Z';
        } else {
            out += ' ';
            out += toString(scale_);
        }
    }
    return out;
}

double Time::unixSeconds() const {
    const Time   u   = to(TimeScale::UTC);
    const double sod = u.jd_.jd2 * utcDayLength(u.jd_.jd1);
    return (u.jd_.jd1 - constants::jdUnixEpoch) * secondsPerDay + sod;
}

std::chrono::system_clock::time_point Time::toSysTime() const {
    using namespace std::chrono;
    const Time   u        = to(TimeScale::UTC);
    const double sod      = u.jd_.jd2 * utcDayLength(u.jd_.jd1);
    const auto   days     = static_cast<long long>(std::llround(u.jd_.jd1 - constants::jdUnixEpoch));
    const double sodWhole = std::floor(sod);
    const auto   ns       = (days * 86400LL + static_cast<long long>(sodWhole)) * 1000000000LL
                  + static_cast<long long>(std::llround((sod - sodWhole) * 1e9));
    return system_clock::time_point(duration_cast<system_clock::duration>(nanoseconds(ns)));
}

double Time::gpsSeconds() const {
    return to(TimeScale::GPS).jd_.secondsSince(JulianDate{constants::jdGpsEpoch, 0.0});
}

GnssWeekTow Time::gnssWeekTow(GnssSystem system) const {
    const Time   t    = to(gnssScale(system));
    const double days = t.jd_.jd1 - gnssEpochJd(system); // exact integer
    const double week = std::floor(days / 7.0);
    GnssWeekTow  out;
    out.week = static_cast<int>(week);
    out.tow  = (days - 7.0 * week) * secondsPerDay + t.jd_.jd2 * secondsPerDay;
    if (system == GnssSystem::Galileo) {
        out.week -= kGalileoWeekOffset;
    }
    return out;
}

Time Time::operator+(double seconds) const {
    if (scale_ == TimeScale::UTC) {
        const JulianDate tai = utcToTai(jd_).addSeconds(seconds);
        return Time(taiToUtc(tai), TimeScale::UTC);
    }
    return Time(jd_.addSeconds(seconds), scale_);
}

double operator-(const Time& a, const Time& b) {
    if (a.scale_ == TimeScale::UTC) {
        return utcToTai(a.jd_).secondsSince(toTai(b.jd_, b.scale_, nullptr));
    }
    const Time bb = b.to(a.scale_);
    return a.jd_.secondsSince(bb.jd_);
}

std::partial_ordering operator<=>(const Time& a, const Time& b) {
    return (a - b) <=> 0.0;
}

// ------------------------------------------------------------------ TDB - TT
//
// Follows ERFA's eraDtdb (same series, topocentric terms and adjustments). See THIRD_PARTY_NOTICES.md.

double tdbMinusTt(const JulianDate& tt, double ut1Fraction, double elong, double uKm, double vKm) noexcept {
    using constants::degToRad;
    const double t = tt.daysSinceJ2000() / constants::daysPerJulianMillennium;

    // Topocentric terms (Moyer 1981, Murray 1983).
    double wt = 0.0;
    if (uKm != 0.0 || vKm != 0.0) {
        const double tsol  = std::fmod(ut1Fraction, 1.0) * constants::twoPi + elong;
        const double w     = t / 3600.0;
        const double elsun = std::fmod(280.46645683 + 1296027711.03429 * w, 360.0) * degToRad;
        const double emsun = std::fmod(357.52910918 + 1295965810.481 * w, 360.0) * degToRad;
        const double d     = std::fmod(297.85019547 + 16029616012.090 * w, 360.0) * degToRad;
        const double elj   = std::fmod(34.35151874 + 109306899.89453 * w, 360.0) * degToRad;
        const double els   = std::fmod(50.07744430 + 44046398.47038 * w, 360.0) * degToRad;
        wt = 0.00029e-10 * uKm * std::sin(tsol + elsun - els) + 0.00100e-10 * uKm * std::sin(tsol - 2.0 * emsun)
           + 0.00133e-10 * uKm * std::sin(tsol - d) + 0.00133e-10 * uKm * std::sin(tsol + elsun - elj)
           - 0.00229e-10 * uKm * std::sin(tsol + 2.0 * elsun + emsun) - 0.02200e-10 * vKm * std::cos(elsun + emsun)
           + 0.05312e-10 * uKm * std::sin(tsol - emsun) - 0.13677e-10 * uKm * std::sin(tsol + 2.0 * elsun)
           - 1.31840e-10 * vKm * std::cos(elsun) + 3.17679e-10 * uKm * std::sin(tsol);
    }

    // Fairhead & Bretagnon series, summed smallest terms first.
    const auto& f   = detail::kFairheadTdb;
    auto        sum = [&](int first, int last) {
        double s = 0.0;
        for (int j = last; j >= first; --j) {
            s += f[j][0] * std::sin(f[j][1] * t + f[j][2]);
        }
        return s;
    };
    const double w0 = sum(0, 473), w1 = sum(474, 678), w2 = sum(679, 763), w3 = sum(764, 783), w4 = sum(784, 786);
    const double wf = t * (t * (t * (t * w4 + w3) + w2) + w1) + w0;

    // Adjustments to use JPL planetary masses instead of IAU.
    const double wj = 0.00065e-6 * std::sin(6069.776754 * t + 4.021194)
                    + 0.00033e-6 * std::sin(213.299095 * t + 5.543132)
                    - 0.00196e-6 * std::sin(6208.294251 * t + 5.696701)
                    - 0.00173e-6 * std::sin(74.781599 * t + 2.435900) + 0.03638e-6 * t * t;
    return wt + wf + wj;
}

int resolveWeekRollover(int truncatedWeek, int bits, GnssSystem system, const Time& reference) {
    if (bits <= 0 || bits > 30) {
        throw Error("resolveWeekRollover: bits out of range");
    }
    const int    mod = 1 << bits;
    const int    ref = reference.gnssWeekTow(system).week;
    const double k   = std::round(static_cast<double>(ref - truncatedWeek) / mod);
    return truncatedWeek + static_cast<int>(k) * mod;
}

std::ostream& operator<<(std::ostream& os, TimeScale value) {
    return detail::writeEnum(os, "TimeScale", value);
}

std::ostream& operator<<(std::ostream& os, GnssSystem value) {
    return detail::writeEnum(os, "GnssSystem", value);
}

} // namespace astro
