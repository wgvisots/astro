/**********************************************************************
 * angle.cpp
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

#include "astro/core/angle.h"

#include <array>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <system_error>

#include "astro/core/error.h"

namespace astro {

using constants::pi;
using constants::twoPi;

double wrapTwoPi(double rad) noexcept {
    double r = std::fmod(rad, twoPi);
    if (r < 0.0) {
        r += twoPi;
    }
    return r >= twoPi ? 0.0 : r;
}

// wrapPi and wrap180 take the remainder first and then move it by at most one turn. Both steps are
// exact: std::fmod has no rounding error, and the turn is added to (or subtracted from) a value within a
// factor of two of it, which is exact by Sterbenz's lemma. A value already in range is therefore
// returned unchanged. (Shifting by half a turn, wrapping to [0, turn) and shifting back rounds twice:
// it moved 0.1 by 8e-17, sent 1e-20 to 0 and sent the double just below pi to -pi.)
double wrapPi(double rad) noexcept {
    double r = std::fmod(rad, twoPi); // in (-2 pi, 2 pi), with the sign of rad
    if (r >= pi) {
        r -= twoPi;
    } else if (r < -pi) {
        r += twoPi;
    }
    return r;
}

double wrap360(double deg) noexcept {
    double r = std::fmod(deg, 360.0);
    if (r < 0.0) {
        r += 360.0;
    }
    return r >= 360.0 ? 0.0 : r;
}

double wrap180(double deg) noexcept {
    double r = std::fmod(deg, 360.0); // in (-360, 360), with the sign of deg
    if (r >= 180.0) {
        r -= 360.0;
    } else if (r < -180.0) {
        r += 360.0;
    }
    return r;
}

// The sign belongs to the whole angle, as in "-5 deg 30 min": it comes from the first argument or from
// `negative`, and the three parts are then magnitudes. (Dropping the sign of the first argument, as this
// did, turned dmsToRad(-5, 30, 0) into +5.5 degrees.)
double hmsToRad(int h, int m, double s, bool negative) noexcept {
    const double v = (std::abs(h) + std::abs(m) / 60.0 + std::abs(s) / 3600.0) * 15.0 * constants::degToRad;
    return (negative || h < 0) ? -v : v;
}

double dmsToRad(int d, int m, double s, bool negative) noexcept {
    const double v = (std::abs(d) + std::abs(m) / 60.0 + std::abs(s) / 3600.0) * constants::degToRad;
    return (negative || d < 0) ? -v : v;
}

namespace {

Sexagesimal decompose(double units, int decimals) {
    Sexagesimal out;
    out.negative       = units < 0.0;
    const double scale = std::pow(10.0, decimals);
    // Total in rounded units of the last printed second digit.
    const double totalTicks     = std::round(std::abs(units) * 3600.0 * scale);
    const double ticksPerMinute = 60.0 * scale;
    const double ticksPerWhole  = 3600.0 * scale;
    const double whole          = std::floor(totalTicks / ticksPerWhole);
    // The whole part must fit Sexagesimal::whole, an int. Converting a not-a-number (NaN), an infinity
    // or a value beyond the largest int is undefined behaviour, so those are reported in `seconds`
    // instead, as documented at radToHms in angle.h: NaN for a NaN angle, infinity for an angle that is
    // infinite or too large. The comparison is written so that a NaN fails it.
    if (!(whole <= static_cast<double>(std::numeric_limits<int>::max()))) {
        out.seconds =
            std::isnan(whole) ? std::numeric_limits<double>::quiet_NaN() : std::numeric_limits<double>::infinity();
        return out;
    }
    const double remWhole = totalTicks - whole * ticksPerWhole;
    const double minutes  = std::floor(remWhole / ticksPerMinute);
    const double secTicks = remWhole - minutes * ticksPerMinute;
    out.whole             = static_cast<int>(whole);
    out.minutes           = static_cast<int>(minutes);
    out.seconds           = secTicks / scale;
    if (out.whole == 0 && out.minutes == 0 && secTicks == 0.0) {
        out.negative = false;
    }
    return out;
}

// What the first field of a sexagesimal string counts: it decides which unit marks are allowed.
enum class Whole {
    Hours,
    Degrees
};

// A unit mark found in the text: the field it labels and the number of bytes it takes.
struct UnitMark {
    std::size_t field{0};  // 0 = hours or degrees, 1 = minutes, 2 = seconds
    std::size_t length{0}; // 0 if there is no unit mark here
    bool        allowed{true};
};

// The unit mark at the start of `rest`, if there is one. Hours take 'h', 'm', 's'. Degrees take 'd' or
// the degree sign, then 'm' or an apostrophe, then 's' or a double quote. A mark of the other quantity
// (degrees in an hours string, say) is found but flagged as not allowed, so that it is reported as an
// error instead of being read as something else.
UnitMark unitMarkAt(std::string_view rest, Whole whole) {
    if (rest.empty()) {
        return {};
    }
    const bool hours = whole == Whole::Hours;
    const auto c     = static_cast<unsigned char>(rest[0]);
    switch (c) {
    case 'h':
    case 'H':  return {0, 1, hours};
    case 'd':
    case 'D':  return {0, 1, !hours};
    case 0xB0: return {0, 1, !hours}; // degree sign in ISO 8859-1 (Latin-1)
    case 0xC2:                        // degree sign in UTF-8 is the two bytes C2 B0
        return rest.size() > 1 && static_cast<unsigned char>(rest[1]) == 0xB0 ? UnitMark{0, 2, !hours} : UnitMark{};
    case 'm':
    case 'M':  return {1, 1, true};
    case '\'': return {1, 1, !hours};
    case 's':
    case 'S':  return {2, 1, true};
    case '"':  return {2, 1, !hours};
    default:   return {};
    }
}

// Throws astro::Error naming the public function (parseHms or parseDms), what is wrong and the text.
[[noreturn]] void failSexagesimal(Whole whole, const char* what, std::string_view text) {
    const char* function = whole == Whole::Hours ? "parseHms" : "parseDms";
    throw Error(std::string(function) + ": " + what + " (got '" + std::string(text) + "')");
}

// Parses the grammar documented at parseHms in angle.h: an optional sign, then one to three numbers
// with exactly one separator between each pair. Returns the value in whole units (hours or degrees).
//
// The numbers are checked character by character before they are converted, so only plain decimal
// digits get through. (Handing each field to std::from_chars and skipping any run of separator
// characters, as this did, let through "nan", "inf", exponents, empty fields such as "5::30", a fraction
// in any field, so that "12.5:30" was 13, and any unit letter anywhere.)
double parseSexagesimal(std::string_view text, Whole whole) {
    std::size_t pos        = 0;
    const auto  atEnd      = [&] { return pos == text.size(); };
    const auto  skipBlanks = [&] {
        const std::size_t start = pos;
        while (!atEnd() && (text[pos] == ' ' || text[pos] == '\t')) {
            ++pos;
        }
        return pos != start;
    };
    const auto skipDigits = [&] {
        const std::size_t start = pos;
        while (!atEnd() && text[pos] >= '0' && text[pos] <= '9') {
            ++pos;
        }
        return pos - start;
    };

    skipBlanks();
    bool negative = false;
    if (!atEnd() && (text[pos] == '-' || text[pos] == '+')) {
        negative = text[pos] == '-';
        ++pos;
        skipBlanks();
    }

    std::array<double, 3> fields{};
    std::size_t           count       = 0;
    bool                  hadFraction = false;
    while (!atEnd()) {
        // A number: digits, with an optional fraction.
        const std::size_t begin  = pos;
        std::size_t       digits = skipDigits();
        const bool        point  = !atEnd() && text[pos] == '.';
        if (point) {
            ++pos;
            digits += skipDigits();
        }
        if (digits == 0) {
            failSexagesimal(whole, "malformed sexagesimal string", text);
        }
        if (count == fields.size()) {
            failSexagesimal(whole, "sexagesimal string must have 1-3 fields", text);
        }
        if (hadFraction) {
            failSexagesimal(whole, "only the last field of a sexagesimal string may have a fraction", text);
        }
        hadFraction = point;
        const std::from_chars_result converted =
            std::from_chars(text.data() + begin, text.data() + pos, fields[count], std::chars_format::fixed);
        if (converted.ec != std::errc{}) {
            failSexagesimal(whole, "malformed sexagesimal string", text); // more digits than a double can hold
        }
        if (count > 0 && fields[count] >= 60.0) {
            failSexagesimal(whole, "sexagesimal minutes/seconds must be < 60", text);
        }

        // What follows the number: its unit mark, or one ':' or ',', or blanks alone. Blanks may
        // surround a mark or a separator, but two separators in a row would be an empty field.
        const bool     blanks = skipBlanks();
        const UnitMark mark   = unitMarkAt(text.substr(pos), whole);
        if (mark.length > 0) {
            if (!mark.allowed || mark.field != count) {
                failSexagesimal(whole, "sexagesimal unit mark is of the wrong kind or out of place", text);
            }
            pos += mark.length;
            skipBlanks();
        } else if (!atEnd() && (text[pos] == ':' || text[pos] == ',')) {
            ++pos;
            skipBlanks();
            if (atEnd()) {
                failSexagesimal(whole, "malformed sexagesimal string", text); // nothing after the separator
            }
        } else if (!atEnd() && !blanks) {
            failSexagesimal(whole, "malformed sexagesimal string", text); // a character that cannot follow a number
        }
        ++count;
    }
    if (count == 0) {
        failSexagesimal(whole, "sexagesimal string must have 1-3 fields", text);
    }
    const double v = fields[0] + fields[1] / 60.0 + fields[2] / 3600.0;
    return negative ? -v : v;
}

std::string format(const Sexagesimal& s, int decimals, char sep, bool forceSign, int wholeWidth) {
    const char* sign = s.negative ? "-" : (forceSign ? "+" : "");
    // An angle that is not finite, or too large to decompose, has no digits: see decompose() above.
    if (std::isnan(s.seconds)) {
        return "nan";
    }
    if (std::isinf(s.seconds)) {
        return std::string(sign) + "inf";
    }
    char      secBuf[32];
    const int secWidth = decimals > 0 ? decimals + 3 : 2;
    std::snprintf(secBuf, sizeof(secBuf), "%0*.*f", secWidth, decimals, s.seconds);
    char buf[96];
    std::snprintf(buf, sizeof(buf), "%s%0*d%c%02d%c%s", sign, wholeWidth, s.whole, sep, s.minutes, sep, secBuf);
    return buf;
}

} // namespace

Sexagesimal radToHms(double rad, int decimals) noexcept {
    return decompose(rad * constants::radToDeg / 15.0, decimals);
}

Sexagesimal radToDms(double rad, int decimals) noexcept {
    return decompose(rad * constants::radToDeg, decimals);
}

// A number with hundreds of digits is a number, but not an angle a double can hold: refuse it rather
// than return infinity.
static double finiteAngle(double radians, std::string_view text, Whole whole) {
    if (!std::isfinite(radians)) {
        failSexagesimal(whole, "sexagesimal value is too large", text);
    }
    return radians;
}

double parseHms(std::string_view text) {
    return finiteAngle(parseSexagesimal(text, Whole::Hours) * (15.0 * constants::degToRad), text, Whole::Hours);
}

double parseDms(std::string_view text) {
    return finiteAngle(parseSexagesimal(text, Whole::Degrees) * constants::degToRad, text, Whole::Degrees);
}

std::string formatHms(double rad, int decimals, char sep) {
    return format(radToHms(rad, decimals), decimals, sep, false, 2);
}

std::string formatDms(double rad, int decimals, char sep) {
    return format(radToDms(rad, decimals), decimals, sep, true, 2);
}

} // namespace astro
