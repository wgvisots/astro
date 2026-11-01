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
        {1972, 1, 10.0}, {1972, 7, 11.0}, {1973, 1, 12.0}, {1974, 1, 13.0},
        {1975, 1, 14.0}, {1976, 1, 15.0}, {1977, 1, 16.0}, {1978, 1, 17.0},
        {1979, 1, 18.0}, {1980, 1, 19.0}, {1981, 7, 20.0}, {1982, 7, 21.0},
        {1983, 7, 22.0}, {1985, 7, 23.0}, {1988, 1, 24.0}, {1990, 1, 25.0},
        {1991, 1, 26.0}, {1992, 7, 27.0}, {1993, 7, 28.0}, {1994, 7, 29.0},
        {1996, 1, 30.0}, {1997, 7, 31.0}, {1999, 1, 32.0}, {2006, 1, 33.0},
        {2009, 1, 34.0}, {2012, 7, 35.0}, {2015, 7, 36.0}, {2017, 1, 37.0},
    };
    return e;
}

constexpr int monthKey(int y, int m) noexcept { return y * 12 + (m - 1); }

double mjdOfDate(int y, int m, int d) {
    return static_cast<double>(julianDayNumber(y, m, d)) - 2400001.0;
}

const char* kMonths[12] = {"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};

// Merge parsed post-1972 entries with the built-in pre-1972 history.
std::vector<LeapSecondEntry> withPre1972(std::vector<LeapSecondEntry> parsed) {
    std::vector<LeapSecondEntry> out;
    for (const auto& e : builtinEntries())
        if (e.year < 1972) out.push_back(e);
    for (const auto& e : parsed)
        if (e.year >= 1972) out.push_back(e);
    return out;
}

LeapSecondTable parseTaiUtcDat(std::string_view content) {
    // " 1961 JAN  1 =JD 2437300.5  TAI-UTC=   1.4228180 S + (MJD - 37300.) X 0.001296 S"
    std::vector<LeapSecondEntry> entries;
    std::istringstream in{std::string(content)};
    std::string line;
    while (std::getline(in, line)) {
        if (line.find("TAI-UTC=") == std::string::npos) continue;
        std::istringstream ls(line);
        int year = 0, day = 0;
        std::string mon;
        ls >> year >> mon >> day;
        int month = 0;
        for (int i = 0; i < 12; ++i)
            if (mon == kMonths[i]) month = i + 1;
        const auto pOff = line.find("TAI-UTC=");
        const auto pRef = line.find("MJD -");
        const auto pRate = line.find(" X ");
        if (month == 0 || pRef == std::string::npos || pRate == std::string::npos)
            throw Error("tai-utc.dat: malformed line: " + line);
        LeapSecondEntry e{year, month, std::stod(line.substr(pOff + 8)), std::stod(line.substr(pRef + 5)),
                          std::stod(line.substr(pRate + 3))};
        entries.push_back(e);
    }
    if (entries.empty()) throw Error("tai-utc.dat: no entries found");
    return LeapSecondTable(std::move(entries));
}

LeapSecondTable parseLeapSecondDat(std::string_view content) {
    // IERS Leap_Second.dat:
    //   #  File expires on 28 June 2026
    //   41317.0    1  1 1972       10
    std::vector<LeapSecondEntry> entries;
    std::optional<double> expires;
    std::istringstream in{std::string(content)};
    std::string line;
    while (std::getline(in, line)) {
        const auto pExp = line.find("File expires on");
        if (pExp != std::string::npos) {
            std::istringstream es(line.substr(pExp + 15));
            int d = 0, y = 0;
            std::string monName;
            es >> d >> monName >> y;
            static const char* full[12] = {"January", "February", "March",     "April",   "May",      "June",
                                           "July",    "August",   "September", "October", "November", "December"};
            for (int i = 0; i < 12; ++i)
                if (monName == full[i]) expires = mjdOfDate(y, i + 1, d);
            continue;
        }
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ls(line);
        double mjd = 0.0, offset = 0.0;
        int d = 0, m = 0, y = 0;
        if (ls >> mjd >> d >> m >> y >> offset) entries.push_back({y, m, offset});
    }
    if (entries.empty()) throw Error("Leap_Second.dat: no entries found");
    return LeapSecondTable(withPre1972(std::move(entries)), expires);
}

LeapSecondTable parseLeapSecondsList(std::string_view content) {
    // IETF leap-seconds.list: "<NTP seconds> <TAI-UTC> # comment", "#@ <NTP expiry>".
    constexpr double ntpEpochMjd = 15020.0;  // 1900-01-01
    std::vector<LeapSecondEntry> entries;
    std::optional<double> expires;
    std::istringstream in{std::string(content)};
    std::string line;
    while (std::getline(in, line)) {
        if (line.rfind("#@", 0) == 0) {
            expires = ntpEpochMjd + std::stod(line.substr(2)) / 86400.0;
            continue;
        }
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ls(line);
        double ntp = 0.0, offset = 0.0;
        if (!(ls >> ntp >> offset)) continue;
        const long long mjd = static_cast<long long>(std::llround(ntpEpochMjd + ntp / 86400.0));
        int y = 0, m = 0, d = 0;
        dateFromJulianDayNumber(mjd + 2400001, y, m, d);
        entries.push_back({y, m, offset});
    }
    if (entries.empty()) throw Error("leap-seconds.list: no entries found");
    return LeapSecondTable(withPre1972(std::move(entries)), expires);
}

std::atomic<std::shared_ptr<const LeapSecondTable>>& globalTable() {
    static std::atomic<std::shared_ptr<const LeapSecondTable>> table{std::make_shared<const LeapSecondTable>()};
    return table;
}

} // namespace

LeapSecondTable::LeapSecondTable() : entries_(builtinEntries()) {}

LeapSecondTable::LeapSecondTable(std::vector<LeapSecondEntry> entries, std::optional<double> expiresMjd)
    : entries_(std::move(entries)), expiresMjd_(expiresMjd) {
    if (entries_.empty()) throw Error("LeapSecondTable: empty table");
    std::stable_sort(entries_.begin(), entries_.end(), [](const LeapSecondEntry& a, const LeapSecondEntry& b) {
        return monthKey(a.year, a.month) < monthKey(b.year, b.month);
    });
}

LeapSecondTable LeapSecondTable::builtin() { return LeapSecondTable(); }

LeapSecondTable LeapSecondTable::parse(std::string_view content) {
    if (content.find("TAI-UTC=") != std::string_view::npos) return parseTaiUtcDat(content);
    if (content.find("#@") != std::string_view::npos || content.find("NTP") != std::string_view::npos)
        return parseLeapSecondsList(content);
    return parseLeapSecondDat(content);
}

LeapSecondTable LeapSecondTable::load(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw Error("LeapSecondTable: cannot open " + path);
    std::stringstream ss;
    ss << in.rdbuf();
    return parse(ss.str());
}

double LeapSecondTable::taiMinusUtc(int year, int month, int day, double fractionOfDay) const {
    const int key = monthKey(year, month);
    if (key < monthKey(entries_.front().year, entries_.front().month))
        throw Error("TAI-UTC undefined before " + std::to_string(entries_.front().year) + " (UTC did not exist)");
    // Last entry whose start month <= requested month.
    auto it = std::upper_bound(entries_.begin(), entries_.end(), key, [](int k, const LeapSecondEntry& e) {
        return k < monthKey(e.year, e.month);
    });
    const LeapSecondEntry& e = *(it - 1);
    double v = e.offset;
    if (e.rate != 0.0) v += (mjdOfDate(year, month, day) + fractionOfDay - e.mjdRef) * e.rate;
    return v;
}

double LeapSecondTable::taiMinusUtc(double mjdUtc) const {
    const double whole = std::floor(mjdUtc);
    int y = 0, m = 0, d = 0;
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
