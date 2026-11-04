/**********************************************************************
 * tle.cpp
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

#include "astro/orbit/tle.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <sstream>

#include "astro/core/angle.h"
#include "astro/core/error.h"

namespace astro {

namespace {

std::string trim(std::string_view s) {
    std::size_t b = 0, e = s.size();
    while (b < e && std::isspace(static_cast<unsigned char>(s[b]))) {
        ++b;
    }
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) {
        --e;
    }
    return std::string(s.substr(b, e - b));
}

// Columns are 1-based inclusive as in the TLE specification.
std::string_view cols(std::string_view line, std::size_t first, std::size_t last) {
    if (line.size() < first) {
        return {};
    }
    return line.substr(first - 1, std::min(last, line.size()) - first + 1);
}

double toDouble(std::string_view s, const char* what) {
    const std::string t = trim(s);
    if (t.empty()) {
        return 0.0;
    }
    try {
        std::size_t  pos = 0;
        const double v   = std::stod(t, &pos);
        if (pos != t.size()) {
            throw std::invalid_argument("trailing");
        }
        return v;
    } catch (const std::exception&) {
        throw Error(std::string("TLE: bad ") + what + " field '" + t + "'");
    }
}

int toInt(std::string_view s, const char* what) {
    const std::string t = trim(s);
    if (t.empty()) {
        return 0;
    }
    try {
        std::size_t pos = 0;
        const int   v   = std::stoi(t, &pos);
        if (pos != t.size()) {
            throw std::invalid_argument("trailing");
        }
        return v;
    } catch (const std::exception&) {
        throw Error(std::string("TLE: bad ") + what + " field '" + t + "'");
    }
}

// "-12345-6" style assumed-decimal exponent field: -0.12345e-6.
double parseExpField(std::string_view s, const char* what) {
    std::string t = trim(s);
    if (t.empty()) {
        return 0.0;
    }
    double sign = 1.0;
    if (t[0] == '-' || t[0] == '+') {
        sign = t[0] == '-' ? -1.0 : 1.0;
        t.erase(0, 1);
    }
    const auto epos = t.find_last_of("+-");
    if (epos == std::string::npos || epos == 0) {
        return sign * toDouble("0." + t, what);
    }
    const double mant = toDouble("0." + trim(t.substr(0, epos)), what);
    const int    ex   = toInt(t.substr(epos), what);
    return sign * mant * std::pow(10.0, ex);
}

std::string formatExpField(double v) {
    // Format as " 12345-6" (8 chars): sign, 5-digit mantissa, signed exponent digit.
    if (v == 0.0) {
        return " 00000-0";
    }
    const char sign = v < 0.0 ? '-' : ' ';
    double     a    = std::abs(v);
    int        ex   = static_cast<int>(std::floor(std::log10(a))) + 1;
    long       mant = std::lround(a / std::pow(10.0, ex) * 1e5);
    if (mant >= 100000) {
        mant /= 10;
        ++ex;
    }
    char buf[16];
    std::snprintf(buf, sizeof(buf), "%c%05ld%c%d", sign, mant % 100000, ex < 0 ? '-' : '+', std::abs(ex) % 10);
    return buf;
}

std::string withChecksum(std::string line) {
    line.resize(68, ' ');
    line += static_cast<char>('0' + tleChecksum(line));
    return line;
}

} // namespace

int tleChecksum(std::string_view line) noexcept {
    int sum = 0;
    for (std::size_t i = 0; i < line.size() && i < 68; ++i) {
        const char c = line[i];
        if (c >= '0' && c <= '9') {
            sum += c - '0';
        }
        if (c == '-') {
            sum += 1;
        }
    }
    return sum % 10;
}

int decodeAlpha5(std::string_view s) {
    const std::string t = trim(s);
    if (t.empty()) {
        throw Error("TLE: empty catalog number");
    }
    if (std::isalpha(static_cast<unsigned char>(t[0]))) {
        const char c = static_cast<char>(std::toupper(static_cast<unsigned char>(t[0])));
        if (c == 'I' || c == 'O') {
            throw Error("TLE: invalid Alpha-5 letter");
        }
        int v = c - 'A' + 10;
        if (c > 'I') {
            --v;
        }
        if (c > 'O') {
            --v;
        }
        return v * 10000 + toInt(t.substr(1), "catalog number");
    }
    return toInt(t, "catalog number");
}

std::string encodeAlpha5(int satnum) {
    char buf[8];
    if (satnum < 100000) {
        std::snprintf(buf, sizeof(buf), "%05d", std::max(0, satnum) % 100000);
        return buf;
    }
    int lead = satnum / 10000;
    if (lead > 33) {
        throw Error("encodeAlpha5: catalog number too large");
    }
    char c = static_cast<char>('A' + lead - 10);
    if (c >= 'I') {
        ++c;
    }
    if (c >= 'O') {
        ++c;
    }
    std::snprintf(buf, sizeof(buf), "%c%04d", c, satnum % 10000);
    return buf;
}

Tle Tle::parse(std::string_view l1, std::string_view l2, std::string_view name, bool verifyChecksum) {
    std::string a = trim(l1), b = trim(l2);
    if (a.size() < 64 || a[0] != '1') {
        throw Error("TLE: line 1 malformed");
    }
    if (b.size() < 63 || b[0] != '2') {
        throw Error("TLE: line 2 malformed");
    }
    if (verifyChecksum) {
        if (a.size() >= 69 && std::isdigit(static_cast<unsigned char>(a[68])) && tleChecksum(a) != a[68] - '0') {
            throw Error("TLE: line 1 checksum mismatch");
        }
        if (b.size() >= 69 && std::isdigit(static_cast<unsigned char>(b[68])) && tleChecksum(b) != b[68] - '0') {
            throw Error("TLE: line 2 checksum mismatch");
        }
    }
    Tle t;
    t.name = trim(name);
    if (!t.name.empty() && t.name[0] == '0' && t.name.size() > 1 && t.name[1] == ' ') {
        t.name = trim(t.name.substr(2));
    }
    t.satnum = decodeAlpha5(cols(a, 3, 7));
    if (decodeAlpha5(cols(b, 3, 7)) != t.satnum) {
        throw Error("TLE: catalog numbers on lines 1 and 2 differ");
    }
    t.classification  = a.size() >= 8 ? a[7] : 'U';
    t.intlDesignator  = trim(cols(a, 10, 17));
    const int    yy   = toInt(cols(a, 19, 20), "epoch year");
    const double day  = toDouble(cols(a, 21, 32), "epoch day");
    const int    year = yy < 57 ? 2000 + yy : 1900 + yy;
    t.epoch = Time(JulianDate{static_cast<double>(julianDayNumber(year, 1, 1)) - 0.5, day - 1.0}, TimeScale::UTC);
    t.meanMotionDot    = toDouble(cols(a, 34, 43), "n-dot");
    t.meanMotionDdot   = parseExpField(cols(a, 45, 52), "n-ddot");
    t.bstar            = parseExpField(cols(a, 54, 61), "B*");
    t.ephemerisType    = toInt(cols(a, 63, 63), "ephemeris type");
    t.elementSetNumber = toInt(cols(a, 65, 68), "element set number");
    t.inclination      = deg2rad(toDouble(cols(b, 9, 16), "inclination"));
    t.raan             = deg2rad(toDouble(cols(b, 18, 25), "RAAN"));
    t.eccentricity     = toDouble("0." + trim(cols(b, 27, 33)), "eccentricity");
    t.argPerigee       = deg2rad(toDouble(cols(b, 35, 42), "argument of perigee"));
    t.meanAnomaly      = deg2rad(toDouble(cols(b, 44, 51), "mean anomaly"));
    t.meanMotion       = toDouble(cols(b, 53, 63), "mean motion");
    t.revolutionNumber = toInt(cols(b, 64, 68), "revolution number");
    return t;
}

std::vector<Tle> Tle::parseMany(std::string_view text, bool verifyChecksum) {
    std::vector<std::string> lines;
    std::istringstream       in{std::string(text)};
    std::string              line;
    while (std::getline(in, line)) {
        const std::string t = trim(line);
        if (!t.empty() && t[0] != '#') {
            lines.push_back(line);
        }
    }
    std::vector<Tle> out;
    std::string      pendingName;
    for (std::size_t i = 0; i < lines.size(); ++i) {
        const std::string t = trim(lines[i]);
        if (t.size() >= 2 && t[0] == '1' && t[1] == ' ' && i + 1 < lines.size() && trim(lines[i + 1])[0] == '2') {
            out.push_back(parse(lines[i], lines[i + 1], pendingName, verifyChecksum));
            pendingName.clear();
            ++i;
        } else {
            pendingName = t;
        }
    }
    return out;
}

std::string Tle::line1() const {
    // Epoch: year and fractional day of year.
    const CalendarDate c = epoch.to(TimeScale::UTC).calendar(9);
    const JulianDate   jan1{static_cast<double>(julianDayNumber(c.year, 1, 1)) - 0.5, 0.0};
    const double       day = epoch.to(TimeScale::UTC).jd().daysSince(jan1) + 1.0;
    char               nd[16];
    std::snprintf(nd, sizeof(nd), "%.8f", meanMotionDot);
    std::string ndot = nd;
    // " .00002182" / "-.00002182" style.
    if (ndot.rfind("-0.", 0) == 0) {
        ndot = "-" + ndot.substr(2);
    } else if (ndot.rfind("0.", 0) == 0) {
        ndot = " " + ndot.substr(1);
    }
    char buf[96];
    std::snprintf(buf,
                  sizeof(buf),
                  "1 %s%c %-8s %02d%012.8f %10s %8s %8s %1d %4d",
                  encodeAlpha5(satnum).c_str(),
                  classification,
                  intlDesignator.c_str(),
                  c.year % 100,
                  day,
                  ndot.c_str(),
                  formatExpField(meanMotionDdot).c_str(),
                  formatExpField(bstar).c_str(),
                  ephemerisType,
                  elementSetNumber % 10000);
    return withChecksum(buf);
}

std::string Tle::line2() const {
    char ecc[16];
    std::snprintf(ecc, sizeof(ecc), "%07ld", std::lround(eccentricity * 1e7));
    char buf[96];
    std::snprintf(buf,
                  sizeof(buf),
                  "2 %s %8.4f %8.4f %s %8.4f %8.4f %11.8f%5d",
                  encodeAlpha5(satnum).c_str(),
                  rad2deg(inclination),
                  rad2deg(wrapTwoPi(raan)),
                  ecc,
                  rad2deg(wrapTwoPi(argPerigee)),
                  rad2deg(wrapTwoPi(meanAnomaly)),
                  meanMotion,
                  revolutionNumber % 100000);
    return withChecksum(buf);
}

} // namespace astro
