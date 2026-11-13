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

#include <cctype>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <vector>

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
    const double remWhole       = totalTicks - whole * ticksPerWhole;
    const double minutes        = std::floor(remWhole / ticksPerMinute);
    const double secTicks       = remWhole - minutes * ticksPerMinute;
    out.whole                   = static_cast<int>(whole);
    out.minutes                 = static_cast<int>(minutes);
    out.seconds                 = secTicks / scale;
    if (out.whole == 0 && out.minutes == 0 && secTicks == 0.0) {
        out.negative = false;
    }
    return out;
}

bool isSeparator(char c) {
    return c == ' ' || c == ':' || c == 'h' || c == 'H' || c == 'd' || c == 'D' || c == 'm' || c == 'M' || c == 's'
        || c == 'S' || c == '\'' || c == '"' || c == '\t' || c == ',';
}

// Parses one-to-three numeric fields with an optional sign. Returns value in "whole" units.
double parseSexagesimal(std::string_view text) {
    std::size_t pos = 0;
    // Strip UTF-8 degree sign (0xC2 0xB0) by treating its bytes as separators.
    auto skipSeps = [&]() {
        while (pos < text.size()) {
            const auto uc = static_cast<unsigned char>(text[pos]);
            if (isSeparator(text[pos]) || uc == 0xC2 || uc == 0xB0) {
                ++pos;
            } else {
                break;
            }
        }
    };
    while (pos < text.size() && std::isspace(static_cast<unsigned char>(text[pos]))) {
        ++pos;
    }
    bool negative = false;
    if (pos < text.size() && (text[pos] == '-' || text[pos] == '+')) {
        negative = text[pos] == '-';
        ++pos;
    }
    std::vector<double> fields;
    while (pos < text.size()) {
        skipSeps();
        if (pos >= text.size()) {
            break;
        }
        double      v     = 0.0;
        const char* begin = text.data() + pos;
        const char* end   = text.data() + text.size();
        auto [ptr, ec]    = std::from_chars(begin, end, v);
        if (ec != std::errc{} || ptr == begin || v < 0.0) {
            throw Error("malformed sexagesimal string: '" + std::string(text) + "'");
        }
        fields.push_back(v);
        pos = static_cast<std::size_t>(ptr - text.data());
    }
    if (fields.empty() || fields.size() > 3) {
        throw Error("sexagesimal string must have 1-3 fields: '" + std::string(text) + "'");
    }
    for (std::size_t i = 1; i < fields.size(); ++i) {
        if (fields[i] >= 60.0) {
            throw Error("sexagesimal minutes/seconds must be < 60: '" + std::string(text) + "'");
        }
    }
    double v = fields[0];
    if (fields.size() > 1) {
        v += fields[1] / 60.0;
    }
    if (fields.size() > 2) {
        v += fields[2] / 3600.0;
    }
    return negative ? -v : v;
}

std::string format(const Sexagesimal& s, int decimals, char sep, bool forceSign, int wholeWidth) {
    char      secBuf[32];
    const int secWidth = decimals > 0 ? decimals + 3 : 2;
    std::snprintf(secBuf, sizeof(secBuf), "%0*.*f", secWidth, decimals, s.seconds);
    char        buf[96];
    const char* sign = s.negative ? "-" : (forceSign ? "+" : "");
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

double parseHms(std::string_view text) {
    return parseSexagesimal(text) * 15.0 * constants::degToRad;
}

double parseDms(std::string_view text) {
    return parseSexagesimal(text) * constants::degToRad;
}

std::string formatHms(double rad, int decimals, char sep) {
    return format(radToHms(rad, decimals), decimals, sep, false, 2);
}

std::string formatDms(double rad, int decimals, char sep) {
    return format(radToDms(rad, decimals), decimals, sep, true, 2);
}

} // namespace astro
