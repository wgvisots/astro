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
#include <charconv>
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

// Columns `first` to `last` of `line`, counted from 1 and inclusive as in the TLE specification. The result always
// has the full width of the field: columns past the end of a short line are blank.
std::string cols(std::string_view line, std::size_t first, std::size_t last) {
    std::string field(last - first + 1, ' ');
    if (line.size() >= first) {
        const std::string_view present = line.substr(first - 1, field.size());
        field.replace(0, present.size(), present);
    }
    return field;
}

[[noreturn]] void badField(const char* what, std::string_view field) {
    throw Error(std::string("TLE: bad ") + what + " field '" + std::string(field) + "'");
}

bool isDigit(char c) {
    return c >= '0' && c <= '9';
}

// `s` without the blanks at its two ends.
std::string_view stripBlanks(std::string_view s) {
    while (!s.empty() && s.front() == ' ') {
        s.remove_prefix(1);
    }
    while (!s.empty() && s.back() == ' ') {
        s.remove_suffix(1);
    }
    return s;
}

// The number fields are read as a fixed-format (Fortran) read would take them: digits, at most one sign and one
// decimal point, and blanks, with an all-blank field read as zero. The C library's conversions (std::stod,
// std::stoi) are not used to decide what a number is: they also take "nan", "inf", hexadecimal and exponents,
// none of which a TLE field can hold, and they follow the program's locale.

// A field with an explicit decimal point: blanks, an optional sign, digits with at most one decimal point, blanks.
double toDouble(std::string_view s, const char* what) {
    const std::string_view t = stripBlanks(s);
    if (t.empty()) {
        return 0.0;
    }
    const bool        negative = t[0] == '-';
    const std::size_t start    = (negative || t[0] == '+') ? 1 : 0;
    std::size_t       digits = 0, points = 0;
    for (std::size_t i = start; i < t.size(); ++i) {
        if (isDigit(t[i])) {
            ++digits;
        } else if (t[i] == '.') {
            ++points;
        } else {
            badField(what, s);
        }
    }
    if (digits == 0 || points > 1) {
        badField(what, s);
    }
    double      v        = 0.0;
    const char* end      = t.data() + t.size();
    const auto [ptr, ec] = std::from_chars(t.data() + start, end, v, std::chars_format::fixed);
    if (ec != std::errc{} || ptr != end) {
        badField(what, s);
    }
    return negative ? -v : v;
}

// A whole-number field: blanks, digits, blanks. None of these fields is signed, and none is wider than five
// columns, so the value cannot overflow.
int toInt(std::string_view s, const char* what) {
    const std::string_view t = stripBlanks(s);
    if (t.size() > 5) {
        badField(what, s);
    }
    int v = 0;
    for (const char c : t) {
        if (!isDigit(c)) {
            badField(what, s);
        }
        v = 10 * v + (c - '0');
    }
    return v;
}

// A field of digits with the decimal point implied before the first column: "0006703" is 0.0006703. A blank column
// is a zero digit wherever it is, as in Vallado's reference reader (twoline2rv fills these columns with zeros
// before it reads the number). Taking the blanks for padding instead would move the digits up and make the value
// ten times too large for every leading blank. Returns false if a column holds anything but a digit or a blank.
bool toImpliedDecimal(std::string_view s, double& value) {
    double digits = 0.0, scale = 1.0;
    for (const char c : s) {
        if (c != ' ' && !isDigit(c)) {
            return false;
        }
        digits  = 10.0 * digits + (c == ' ' ? 0 : c - '0');
        scale  *= 10.0;
    }
    // Both numbers are exact, so the quotient is the double nearest to the decimal number, as a decimal read of
    // "0.NNNNNNN" would give.
    value = digits / scale;
    return true;
}

// The "+NNNNN-N" fields (second derivative of the mean motion, and B*): a sign, five mantissa digits with the
// decimal point implied before them, the sign of the exponent and one exponent digit, so "-11606-4" is
// -0.11606e-4. A blank sign is a plus and a blank digit is a zero. The value is formed as the reference reader
// forms it, mantissa times pow(10, exponent), so that SGP4 starts from the same B* bit for bit.
double parseExpField(std::string_view s, const char* what) {
    const auto isSign   = [](char c) { return c == ' ' || c == '+' || c == '-'; };
    double     mantissa = 0.0;
    if (s.size() != 8 || !isSign(s[0]) || !isSign(s[6]) || (s[7] != ' ' && !isDigit(s[7]))
        || !toImpliedDecimal(s.substr(1, 5), mantissa)) {
        badField(what, s);
    }
    const int    exponent = s[7] == ' ' ? 0 : s[7] - '0';
    const double value    = mantissa * std::pow(10.0, s[6] == '-' ? -exponent : exponent);
    return s[0] == '-' ? -value : value;
}

// The epoch of an element set is a day of the year and a fraction of an 86400 s day, counted from 0h UTC, whatever
// the length of that day. The sets are produced from, and read back into, a Julian date that knows nothing of leap
// seconds (in the reference code, twoline2rv turns the fraction into hours, minutes and seconds of a 24 h day).
// astro::Time holds UTC differently: as a quasi-Julian date whose day fraction is the time of day over the length
// of *that* day, which is 86401 s when the day ends with a leap second (see time.h). The two fractions differ by
// up to 1 s on such a day, 7.7 km along the track in low Earth orbit.
//
// This is the factor between them for the UTC day that starts at the Julian date `midnight`: (86400 s) / (length of
// the day). 12:00:00 is 43200 s into any day, so it is twice the quasi-Julian fraction of noon. It is exactly 1 on
// an ordinary day. (Sgp4 needs the same factor and has a copy of this function.)
double tleDayPerUtcDay(double midnight) {
    int year = 0, month = 0, day = 0;
    dateFromJulianDayNumber(std::llround(midnight + 0.5), year, month, day);
    return 2.0 * Time::fromCalendar(TimeScale::UTC, year, month, day, 12).jd().jd2;
}

// The epoch from the two-digit year (columns 19-20) and the day of the year with its fraction (columns 21-32).
Time parseEpoch(std::string_view yearField, std::string_view dayField) {
    const int    yy   = toInt(yearField, "epoch year");
    const double day  = toDouble(dayField, "epoch day");
    const int    year = yy < 57 ? 2000 + yy : 1900 + yy;
    // Days are numbered from 1, so the field runs from 1.0 to just under 366.0 (367.0 in a leap year).
    const int    daysInYear = isLeapYear(year) ? 366 : 365;
    const double wholeDay   = std::floor(day);
    if (!(wholeDay >= 1.0 && wholeDay <= daysInYear)) {
        throw Error("TLE: bad epoch day field '" + std::string(dayField) + "' (" + std::to_string(year)
                    + " has days 1 to " + std::to_string(daysInYear) + ")");
    }
    const double midnight = static_cast<double>(julianDayNumber(year, 1, 1)) - 0.5 + (wholeDay - 1.0);
    return Time(JulianDate{midnight, (day - wholeDay) * tleDayPerUtcDay(midnight)}, TimeScale::UTC);
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
    t.classification   = a.size() >= 8 ? a[7] : 'U';
    t.intlDesignator   = trim(cols(a, 10, 17));
    t.epoch            = parseEpoch(cols(a, 19, 20), cols(a, 21, 32));
    t.meanMotionDot    = toDouble(cols(a, 34, 43), "n-dot");
    t.meanMotionDdot   = parseExpField(cols(a, 45, 52), "n-ddot");
    t.bstar            = parseExpField(cols(a, 54, 61), "B*");
    t.ephemerisType    = toInt(cols(a, 63, 63), "ephemeris type");
    t.elementSetNumber = toInt(cols(a, 65, 68), "element set number");
    t.inclination      = deg2rad(toDouble(cols(b, 9, 16), "inclination"));
    t.raan             = deg2rad(toDouble(cols(b, 18, 25), "RAAN"));
    if (!toImpliedDecimal(cols(b, 27, 33), t.eccentricity)) {
        badField("eccentricity", cols(b, 27, 33));
    }
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
    // Epoch: the year, and the day of the year with the time of day as a fraction of 86400 s (see tleDayPerUtcDay).
    const JulianDate utc  = epoch.to(TimeScale::UTC).jd();
    int              year = 0, month = 0, dayOfMonth = 0;
    dateFromJulianDayNumber(std::llround(utc.jd1 + 0.5), year, month, dayOfMonth);
    const double day = dayOfYear(year, month, dayOfMonth) + utc.jd2 / tleDayPerUtcDay(utc.jd1);
    char         nd[16];
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
                  year % 100,
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
