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

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstdlib>
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

// Writes a "+NNNNN-N" field (see parseExpField) as the value nearest to `v` that the field can hold: five mantissa
// digits and an exponent from -9 to +9. Below 0.1e-9 the exponent stays at -9 and the mantissa gets leading zeros,
// down to the smallest step, 0.00001e-9 = 1e-14; a value nearer to zero than half a step is written as zero.
// Throws astro::Error naming `what` if the value is not finite or rounds beyond the largest one, 0.99999e9.
std::string formatExpField(double v, const char* what) {
    if (!std::isfinite(v)) {
        throw Error(std::string("Tle::line1: ") + what + " is not finite");
    }
    const double a    = std::abs(v);
    const char*  zero = " 00000-0";
    if (a == 0.0) {
        return zero;
    }
    // The C library rounds correctly to the digits asked for, here five significant ones: "d.dddde-XX". Taking the
    // exponent from a logarithm instead can be off by one next to a power of ten.
    char text[32];
    std::snprintf(text, sizeof(text), "%.4e", a);
    long mantissa = 10000 * (text[0] - '0') + std::strtol(text + 2, nullptr, 10);
    int  exponent = static_cast<int>(std::strtol(text + 7, nullptr, 10)) + 1; // of 0.ddddd, not of d.dddd
    if (exponent > 9) {
        std::snprintf(text, sizeof(text), "%g", v);
        throw Error(std::string("Tle::line1: ") + what + " " + text
                    + " is too large for its field (the largest value is 0.99999e9)");
    }
    if (exponent < -9) {
        // Fourteen decimals, of which the first nine are zeros and the last five are the mantissa.
        std::snprintf(text, sizeof(text), "%.14f", a);
        mantissa = std::strtol(text + 11, nullptr, 10);
        exponent = -9;
        if (mantissa == 0) {
            return zero;
        }
    }
    std::snprintf(text,
                  sizeof(text),
                  "%c%05ld%c%d",
                  v < 0.0 ? '-' : ' ',
                  mantissa,
                  exponent < 0 ? '-' : '+',
                  std::abs(exponent));
    return text;
}

// The writers below refuse a value its field cannot hold. Writing it anyway would push the later fields out of
// their columns, or wrap the value into another one, and the line would still carry a valid checksum.

[[noreturn]] void cannotWrite(const char* function, const char* what, double value, const char* why) {
    char text[32];
    std::snprintf(text, sizeof(text), "%.10g", value);
    throw Error(std::string(function) + ": " + what + " " + text + " " + why);
}

// A number with `decimals` decimals, right-justified in `width` columns. Throws astro::Error naming `what` if it is
// not finite or needs more columns than the field has.
std::string fixedField(double v, int width, int decimals, const char* function, const char* what) {
    if (!std::isfinite(v)) {
        cannotWrite(function, what, v, "is not finite");
    }
    char buf[32];
    if (std::snprintf(buf, sizeof(buf), "%*.*f", width, decimals, v) != width) {
        cannotWrite(function, what, v, "does not fit its field");
    }
    return buf;
}

// An angle of line 2 other than the inclination: degrees in [0, 360) with four decimals. An angle that rounds up to
// 360.0000 is written as 0.0000, because the top of the range is not in it.
std::string angleField(double rad, const char* what) {
    if (!std::isfinite(rad)) {
        cannotWrite("Tle::line2", what, rad, "is not finite");
    }
    const std::string text = fixedField(rad2deg(wrapTwoPi(rad)), 8, 4, "Tle::line2", what);
    return text == "360.0000" ? "  0.0000" : text;
}

// The first derivative of the mean motion: a minus sign or a blank, then the decimal point and eight decimals. The
// field has no column for a digit before the point.
std::string ndotField(double v) {
    if (!std::isfinite(v)) {
        cannotWrite("Tle::line1", "n-dot", v, "is not finite");
    }
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.8f", v);
    const std::string text = buf; // "0.00002182" or "-0.00002182"
    if (text.rfind("-0.", 0) == 0) {
        return "-" + text.substr(2);
    }
    if (text.rfind("0.", 0) == 0) {
        return " " + text.substr(1);
    }
    cannotWrite("Tle::line1", "n-dot", v, "does not fit its field (the magnitude must round to less than 1)");
}

// The epoch field: the two-digit year, the day of the year, and the time of day as a fraction of 86400 s to eight
// decimals (see tleDayPerUtcDay).
std::string epochField(const Time& epoch) {
    if (!std::isfinite(epoch.jd().jd1 + epoch.jd().jd2)) {
        throw Error("Tle::line1: epoch is not finite");
    }
    const JulianDate utc  = epoch.to(TimeScale::UTC).jd();
    int              year = 0, month = 0, dayOfMonth = 0;
    dateFromJulianDayNumber(std::llround(utc.jd1 + 0.5), year, month, dayOfMonth);
    const double perUtcDay = tleDayPerUtcDay(utc.jd1);
    const double fraction  = utc.jd2 / perUtcDay; // of 86400 s; it reaches 1 and goes on during a leap second
    int          day       = dayOfYear(year, month, dayOfMonth);
    long long    units     = std::llround(fraction * 1e8); // of 1e-8 day (0.864 ms), the resolution of the field
    if (units >= 100000000) {
        // The fraction rounds up to the end of an 86400 s day, which the field cannot say ("366.00000000" is not a
        // time on day 365). Normally the end of the day is 0h of the next day, and that is the epoch to write. On
        // a day with a leap second 0h is a second further on.
        const double toNextDay = 1.0 / perUtcDay - fraction; // in days of 86400 s
        if (toNextDay * 1e8 <= 0.5) {
            units = 0;
            if (++day > (isLeapYear(year) ? 366 : 365)) {
                day = 1;
                ++year;
            }
        } else if (fraction < 1.0) {
            units = 99999999; // the last half step before a leap second: this is the nearest epoch the field holds
        } else {
            throw Error("Tle::line1: epoch " + epoch.to(TimeScale::UTC).iso8601(3)
                        + " is inside a leap second, which a fraction of an 86400 s day cannot name");
        }
    }
    if (year < 1957 || year > 2056) {
        throw Error("Tle::line1: epoch year " + std::to_string(year)
                    + " does not fit its two-digit field (1957 to 2056)");
    }
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%02d%03d.%08lld", year % 100, day, units);
    return buf;
}

// Whether a character can stand in a line: printable ASCII, the blank included.
bool isPrintable(char c) {
    return c >= ' ' && c <= '~';
}

// Adds the checksum to the 68 columns of a line. Every field was written to its exact width, so the length is a
// check on this file and not on the caller: a line cut or padded to 68 columns would be the malformed line with a
// valid checksum that the field checks exist to prevent.
std::string withChecksum(std::string line) {
    if (line.size() != 68) {
        throw Error("Tle: a line of " + std::to_string(line.size()) + " columns was formatted (internal error)");
    }
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
    if (satnum < 0 || satnum > 339999) {
        throw Error("encodeAlpha5: catalog number " + std::to_string(satnum)
                    + " cannot be written in 5 characters (0 to 339999)");
    }
    char buf[8];
    if (satnum < 100000) {
        std::snprintf(buf, sizeof(buf), "%05d", satnum);
        return buf;
    }
    const int lead = satnum / 10000;
    char      c    = static_cast<char>('A' + lead - 10);
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
    if (!isPrintable(classification)) {
        throw Error("Tle::line1: classification is not a printable character (code "
                    + std::to_string(static_cast<int>(classification)) + ")");
    }
    if (intlDesignator.size() > 8 || !std::all_of(intlDesignator.begin(), intlDesignator.end(), isPrintable)) {
        throw Error("Tle::line1: international designator '" + intlDesignator
                    + "' does not fit its field (at most 8 printable characters)");
    }
    if (ephemerisType < 0 || ephemerisType > 9) {
        throw Error("Tle::line1: ephemeris type " + std::to_string(ephemerisType) + " is not one digit");
    }
    if (elementSetNumber < 0) {
        throw Error("Tle::line1: element set number " + std::to_string(elementSetNumber) + " is negative");
    }
    char buf[96];
    std::snprintf(buf,
                  sizeof(buf),
                  "1 %s%c %-8s %s %s %s %s %1d %4d",
                  encodeAlpha5(satnum).c_str(),
                  classification,
                  intlDesignator.c_str(),
                  epochField(epoch).c_str(),
                  ndotField(meanMotionDot).c_str(),
                  formatExpField(meanMotionDdot, "n-ddot").c_str(),
                  formatExpField(bstar, "B*").c_str(),
                  ephemerisType,
                  elementSetNumber % 10000);
    return withChecksum(buf);
}

std::string Tle::line2() const {
    // Seven digits after an implied decimal point: 0 to 0.9999999.
    const double eccentricityDigits = std::isfinite(eccentricity) ? std::round(eccentricity * 1e7) : -1.0;
    if (!(eccentricityDigits >= 0.0 && eccentricityDigits <= 9999999.0)) {
        cannotWrite("Tle::line2", "eccentricity", eccentricity, "does not fit its field (0 to 0.9999999)");
    }
    if (revolutionNumber < 0) {
        throw Error("Tle::line2: revolution number " + std::to_string(revolutionNumber) + " is negative");
    }
    char buf[96];
    std::snprintf(buf,
                  sizeof(buf),
                  "2 %s %s %s %07ld %s %s %s%5d",
                  encodeAlpha5(satnum).c_str(),
                  fixedField(rad2deg(inclination), 8, 4, "Tle::line2", "inclination").c_str(),
                  angleField(raan, "RAAN").c_str(),
                  static_cast<long>(eccentricityDigits),
                  angleField(argPerigee, "argument of perigee").c_str(),
                  angleField(meanAnomaly, "mean anomaly").c_str(),
                  fixedField(meanMotion, 11, 8, "Tle::line2", "mean motion").c_str(),
                  revolutionNumber % 100000);
    return withChecksum(buf);
}

} // namespace astro
