/**********************************************************************
 * leap_seconds.cpp
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

#include "astro/time/leap_seconds.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <fstream>
#include <sstream>

#include "astro/core/error.h"
#include "astro/time/julian_date.h"

namespace astro {

namespace {

// TAI-UTC history (USNO tai-utc.dat / IERS Bulletin C).
const std::vector<LeapSecondEntry>& builtinEntries() {
    static const std::vector<LeapSecondEntry> e = {
        {1960, 1, 1.4178180, 37300.0, 0.0012960},
        {1961, 1, 1.4228180, 37300.0, 0.0012960},
        {1961, 8, 1.3728180, 37300.0, 0.0012960},
        {1962, 1, 1.8458580, 37665.0, 0.0011232},
        {1963, 11, 1.9458580, 37665.0, 0.0011232},
        {1964, 1, 3.2401300, 38761.0, 0.0012960},
        {1964, 4, 3.3401300, 38761.0, 0.0012960},
        {1964, 9, 3.4401300, 38761.0, 0.0012960},
        {1965, 1, 3.5401300, 38761.0, 0.0012960},
        {1965, 3, 3.6401300, 38761.0, 0.0012960},
        {1965, 7, 3.7401300, 38761.0, 0.0012960},
        {1965, 9, 3.8401300, 38761.0, 0.0012960},
        {1966, 1, 4.3131700, 39126.0, 0.0025920},
        {1968, 2, 4.2131700, 39126.0, 0.0025920},
        {1972, 1, 10.0},
        {1972, 7, 11.0},
        {1973, 1, 12.0},
        {1974, 1, 13.0},
        {1975, 1, 14.0},
        {1976, 1, 15.0},
        {1977, 1, 16.0},
        {1978, 1, 17.0},
        {1979, 1, 18.0},
        {1980, 1, 19.0},
        {1981, 7, 20.0},
        {1982, 7, 21.0},
        {1983, 7, 22.0},
        {1985, 7, 23.0},
        {1988, 1, 24.0},
        {1990, 1, 25.0},
        {1991, 1, 26.0},
        {1992, 7, 27.0},
        {1993, 7, 28.0},
        {1994, 7, 29.0},
        {1996, 1, 30.0},
        {1997, 7, 31.0},
        {1999, 1, 32.0},
        {2006, 1, 33.0},
        {2009, 1, 34.0},
        {2012, 7, 35.0},
        {2015, 7, 36.0},
        {2017, 1, 37.0},
    };
    return e;
}

constexpr long long monthKey(int y, int m) noexcept {
    return static_cast<long long>(y) * 12 + (m - 1);
}

double mjdOfDate(int y, int m, int d) {
    return static_cast<double>(julianDayNumber(y, m, d)) - 2400001.0;
}

const char* kMonths[12] = {"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};

// The three file formats, as named in error messages.
constexpr const char* kTaiUtcDat       = "tai-utc.dat";
constexpr const char* kLeapSecondDat   = "Leap_Second.dat";
constexpr const char* kLeapSecondsList = "leap-seconds.list";

constexpr std::string_view kWhiteSpace = " \t\r\n\f\v";

/// `text` without the white space around it (which includes the '\r' of a Windows line ending).
std::string_view trimmed(std::string_view text) {
    const auto first = text.find_first_not_of(kWhiteSpace);
    if (first == std::string_view::npos) {
        return {};
    }
    return text.substr(first, text.find_last_not_of(kWhiteSpace) - first + 1);
}

/// The data on one line of a leap-second file: the text before any '#' comment. Empty for a blank line and
/// for a comment line.
std::string_view dataOf(std::string_view line) {
    return trimmed(line.substr(0, line.find('#')));
}

/// True once nothing but white space is left to read.
bool atEnd(std::istream& in) {
    return in.eof() || (in >> std::ws).eof();
}

/// Reports a line that cannot be used, by its number and text. (A line that is skipped instead goes unnoticed,
/// and TAI-UTC is then wrong by a second from that line's date on.)
[[noreturn]] void failLine(const char* format, int lineNumber, std::string_view line, const char* why) {
    throw Error(std::string(format) + ": line " + std::to_string(lineNumber) + ": " + why + ": '"
                + std::string(trimmed(line)) + "'");
}

// The built-in history before the first entry read from a file, followed by all entries of the file.
// Leap_Second.dat and leap-seconds.list start in 1972 and tai-utc.dat in 1961, so each needs the earlier
// years from here; whatever the file itself states is kept.
std::vector<LeapSecondEntry> withEarlierHistory(const std::vector<LeapSecondEntry>& parsed) {
    long long first = monthKey(parsed.front().year, parsed.front().month);
    for (const auto& e : parsed) {
        first = std::min(first, monthKey(e.year, e.month));
    }
    std::vector<LeapSecondEntry> out;
    for (const auto& e : builtinEntries()) {
        if (monthKey(e.year, e.month) < first) {
            out.push_back(e);
        }
    }
    out.insert(out.end(), parsed.begin(), parsed.end());
    return out;
}

/// The number that starts at `pos` of a tai-utc.dat line, after blanks. It must end the line or be followed
/// by white space, by ')' or by the unit letter 'S'; any other character ("1O.0") means a damaged number.
double numberAt(const std::string& line, std::size_t pos, int lineNumber) {
    double      value = 0.0;
    std::size_t used  = 0;
    try {
        value = std::stod(line.substr(pos), &used);
    } catch (const std::exception&) { // std::invalid_argument, std::out_of_range
        failLine(kTaiUtcDat, lineNumber, line, "malformed number");
    }
    const std::size_t next = pos + used;
    if (next < line.size() && kWhiteSpace.find(line[next]) == std::string_view::npos && line[next] != ')'
        && line[next] != 'S') {
        failLine(kTaiUtcDat, lineNumber, line, "malformed number");
    }
    return value;
}

LeapSecondTable parseTaiUtcDat(std::string_view content) {
    // " 1961 JAN  1 =JD 2437300.5  TAI-UTC=   1.4228180 S + (MJD - 37300.) X 0.001296 S"
    std::vector<LeapSecondEntry> entries;
    std::istringstream           in{std::string(content)};
    std::string                  line;
    int                          lineNumber = 0;
    while (std::getline(in, line)) {
        ++lineNumber;
        if (dataOf(line).empty()) {
            continue;
        }
        std::istringstream ls(line);
        int                year = 0, day = 0;
        std::string        mon;
        ls >> year >> mon >> day;
        int month = 0;
        for (int i = 0; i < 12; ++i) {
            if (mon == kMonths[i]) {
                month = i + 1;
            }
        }
        const auto pJd   = line.find("=JD");
        const auto pOff  = line.find("TAI-UTC=");
        const auto pRef  = line.find("MJD -");
        const auto pRate = line.find(" X ");
        if (!ls || month == 0 || pJd == std::string::npos || pOff == std::string::npos || pRef == std::string::npos
            || pRate == std::string::npos) {
            failLine(kTaiUtcDat,
                     lineNumber,
                     line,
                     "expected 'year MON day =JD <JD>  TAI-UTC= <s> S + (MJD - <MJD>) X <s/day> S'");
        }
        // A step takes effect at 0h on the first of a month, and the JD repeats the date. A line that
        // says otherwise is damaged; using it would put the step on the wrong day.
        if (day != 1 || numberAt(line, pJd + 3, lineNumber) != mjdOfDate(year, month, day) + constants::mjdOffset) {
            failLine(kTaiUtcDat, lineNumber, line, "the date is not the first of a month, or is not the JD given");
        }
        LeapSecondEntry e{year,
                          month,
                          numberAt(line, pOff + 8, lineNumber),
                          numberAt(line, pRef + 5, lineNumber),
                          numberAt(line, pRate + 3, lineNumber)};
        entries.push_back(e);
    }
    if (entries.empty()) {
        throw Error("tai-utc.dat: no entries found");
    }
    return LeapSecondTable(withEarlierHistory(entries));
}

LeapSecondTable parseLeapSecondDat(std::string_view content) {
    // IERS Leap_Second.dat:
    //   #  File expires on 28 June 2026
    //   41317.0    1  1 1972       10          (MJD, day, month, year, TAI-UTC)
    std::vector<LeapSecondEntry> entries;
    std::optional<double>        expires;
    std::istringstream           in{std::string(content)};
    std::string                  line;
    int                          lineNumber = 0;
    while (std::getline(in, line)) {
        ++lineNumber;
        const auto pExp = line.find("File expires on");
        if (pExp != std::string::npos) {
            std::istringstream es(line.substr(pExp + 15));
            int                d = 0, y = 0;
            std::string        monName;
            es >> d >> monName >> y;
            static const char* full[12] = {"January",
                                           "February",
                                           "March",
                                           "April",
                                           "May",
                                           "June",
                                           "July",
                                           "August",
                                           "September",
                                           "October",
                                           "November",
                                           "December"};
            for (int i = 0; i < 12; ++i) {
                if (monName == full[i]) {
                    expires = mjdOfDate(y, i + 1, d);
                }
            }
            continue;
        }
        const std::string_view data = dataOf(line);
        if (data.empty()) {
            continue;
        }
        std::istringstream ls{std::string(data)};
        double             mjd = 0.0, offset = 0.0;
        int                d = 0, m = 0, y = 0;
        if (!(ls >> mjd >> d >> m >> y >> offset) || !atEnd(ls)) {
            failLine(kLeapSecondDat, lineNumber, line, "expected 'MJD day month year TAI-UTC'");
        }
        // A step takes effect at 0h on the first of a month, and the MJD repeats the date. A line that
        // says otherwise is damaged; using it would put the step on the wrong day.
        if (m < 1 || m > 12 || d != 1 || mjd != mjdOfDate(y, m, d)) {
            failLine(kLeapSecondDat, lineNumber, line, "the date is not the first of a month, or is not the MJD given");
        }
        entries.push_back({y, m, offset});
    }
    if (entries.empty()) {
        throw Error("Leap_Second.dat: no entries found");
    }
    return LeapSecondTable(withEarlierHistory(entries), expires);
}

constexpr double kMaxNtpSeconds = 1.0e11; // ~3000 years after 1900

LeapSecondTable parseLeapSecondsList(std::string_view content) {
    // IETF leap-seconds.list: "<NTP seconds> <TAI-UTC> # comment", "#@ <NTP expiry>".
    constexpr double             ntpEpochMjd = 15020.0; // 1900-01-01
    std::vector<LeapSecondEntry> entries;
    std::optional<double>        expires;
    std::istringstream           in{std::string(content)};
    std::string                  line;
    int                          lineNumber = 0;
    while (std::getline(in, line)) {
        ++lineNumber;
        if (line.rfind("#@", 0) == 0) {
            std::istringstream es(line.substr(2));
            double             ntp = 0.0;
            if (!(es >> ntp) || !atEnd(es)) {
                failLine(kLeapSecondsList, lineNumber, line, "expected '#@ <NTP seconds>'");
            }
            expires = ntpEpochMjd + ntp / 86400.0;
            continue;
        }
        const std::string_view data = dataOf(line);
        if (data.empty()) {
            continue;
        }
        std::istringstream ls{std::string(data)};
        double             ntp = 0.0, offset = 0.0;
        if (!(ls >> ntp >> offset) || !atEnd(ls)) {
            failLine(kLeapSecondsList, lineNumber, line, "expected '<NTP seconds> <TAI-UTC>'");
        }
        if (!(ntp >= 0.0 && ntp < kMaxNtpSeconds)) {
            failLine(kLeapSecondsList, lineNumber, line, "NTP time out of range");
        }
        const double    days = ntp / 86400.0;
        const long long mjd  = static_cast<long long>(std::llround(ntpEpochMjd + days));
        int             y = 0, m = 0, d = 0;
        dateFromJulianDayNumber(mjd + 2400001, y, m, d);
        // A step takes effect at 0h on the first of a month; the table cannot hold any other instant.
        if (days != std::floor(days) || d != 1) {
            failLine(kLeapSecondsList, lineNumber, line, "the NTP time is not 0h on the first of a month");
        }
        entries.push_back({y, m, offset});
    }
    if (entries.empty()) {
        throw Error("leap-seconds.list: no entries found");
    }
    return LeapSecondTable(withEarlierHistory(entries), expires);
}

std::atomic<std::shared_ptr<const LeapSecondTable>>& globalTable() {
    static std::atomic<std::shared_ptr<const LeapSecondTable>> table{std::make_shared<const LeapSecondTable>()};
    return table;
}

} // namespace

LeapSecondTable::LeapSecondTable()
    : entries_(builtinEntries()) {}

LeapSecondTable::LeapSecondTable(std::vector<LeapSecondEntry> entries, std::optional<double> expiresMjd)
    : entries_(std::move(entries))
    , expiresMjd_(expiresMjd) {
    if (entries_.empty()) {
        throw Error("LeapSecondTable: empty table");
    }
    for (const auto& e : entries_) {
        if (e.year < 1960 || e.year > 9999 || e.month < 1 || e.month > 12 || !std::isfinite(e.offset)
            || !std::isfinite(e.mjdRef) || !std::isfinite(e.rate)) {
            throw Error("LeapSecondTable: invalid entry for " + std::to_string(e.year) + "-" + std::to_string(e.month));
        }
    }
    if (expiresMjd_ && !std::isfinite(*expiresMjd_)) {
        throw Error("LeapSecondTable: invalid expiration date");
    }
    std::stable_sort(entries_.begin(), entries_.end(), [](const LeapSecondEntry& a, const LeapSecondEntry& b) {
        return monthKey(a.year, a.month) < monthKey(b.year, b.month);
    });
}

LeapSecondTable LeapSecondTable::builtin() {
    return LeapSecondTable();
}

LeapSecondTable LeapSecondTable::parse(std::string_view content) {
    // The format is recognised from the shape of the first data line, never from the comments: those are
    // free text, and the published files use each other's words in them.
    //   tai-utc.dat        " 1961 JAN  1 =JD 2437300.5  TAI-UTC=   1.4228180 S + (MJD - 37300.) X 0.001296 S"
    //   leap-seconds.list  "2272060800      10      # 1 Jan 1972"     two fields
    //   Leap_Second.dat    "    41317.0    1  1 1972       10"        five fields
    std::istringstream in{std::string(content)};
    std::string        line;
    int                lineNumber = 0;
    while (std::getline(in, line)) {
        ++lineNumber;
        const std::string_view data = dataOf(line);
        if (data.empty()) {
            continue;
        }
        if (data.find("TAI-UTC=") != std::string_view::npos) {
            return parseTaiUtcDat(content);
        }
        std::istringstream fields{std::string(data)};
        std::string        field;
        int                count = 0;
        while (fields >> field) {
            ++count;
        }
        if (count == 2) {
            return parseLeapSecondsList(content);
        }
        if (count == 5) {
            return parseLeapSecondDat(content);
        }
        failLine("LeapSecondTable",
                 lineNumber,
                 line,
                 "not a line of Leap_Second.dat, leap-seconds.list or tai-utc.dat");
    }
    throw Error("LeapSecondTable: no entries found");
}

LeapSecondTable LeapSecondTable::load(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw Error("LeapSecondTable: cannot open " + path);
    }
    std::stringstream ss;
    ss << in.rdbuf();
    return parse(ss.str());
}

double LeapSecondTable::taiMinusUtc(int year, int month, int day, double fractionOfDay) const {
    const long long key = monthKey(year, month);
    if (key < monthKey(entries_.front().year, entries_.front().month)) {
        throw Error("TAI-UTC undefined before " + std::to_string(entries_.front().year) + " (UTC did not exist)");
    }
    // Last entry whose start month <= requested month.
    auto it = std::upper_bound(entries_.begin(), entries_.end(), key, [](long long k, const LeapSecondEntry& e) {
        return k < monthKey(e.year, e.month);
    });
    const LeapSecondEntry& e = *(it - 1);
    double                 v = e.offset;
    if (e.rate != 0.0) {
        v += (mjdOfDate(year, month, day) + fractionOfDay - e.mjdRef) * e.rate;
    }
    return v;
}

double LeapSecondTable::taiMinusUtc(double mjdUtc) const {
    if (!(std::abs(mjdUtc) < 1.0e9)) {
        throw Error("TAI-UTC: MJD out of range");
    }
    const double whole = std::floor(mjdUtc);
    int          y = 0, m = 0, d = 0;
    dateFromJulianDayNumber(static_cast<long long>(whole) + 2400001, y, m, d);
    return taiMinusUtc(y, m, d, mjdUtc - whole);
}

std::shared_ptr<const LeapSecondTable> leapSecondTable() {
    return globalTable().load();
}

void setLeapSecondTable(LeapSecondTable table) {
    globalTable().store(std::make_shared<const LeapSecondTable>(std::move(table)));
}

} // namespace astro
