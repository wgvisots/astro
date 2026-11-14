/**********************************************************************
 * eop.cpp
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

#include "astro/frames/eop.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string_view>

#include "../core/enum_io.h"
#include "astro/core/constants.h"
#include "astro/core/error.h"
#include "astro/time/leap_seconds.h"

namespace astro {

using constants::arcsecToRad;
using constants::masToRad;

namespace {

std::string readFile(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw Error("EopTable: cannot open " + path);
    }
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

// A number as the IERS files write it: an optional sign, then digits with at most one decimal point,
// and nothing else. Returns false for anything else. (std::stod is not used to decide what a number
// is: it also takes "nan", "inf", hexadecimal, exponents and a number followed by other text, and it
// follows the program's locale.)
bool plainNumber(std::string_view s, double& out) {
    std::size_t i      = (!s.empty() && (s[0] == '+' || s[0] == '-')) ? 1 : 0;
    std::size_t digits = 0, points = 0;
    for (std::size_t k = i; k < s.size(); ++k) {
        if (s[k] >= '0' && s[k] <= '9') {
            ++digits;
        } else if (s[k] == '.') {
            ++points;
        } else {
            return false;
        }
    }
    if (digits == 0 || points > 1) {
        return false;
    }
    double      value    = 0.0;
    const char* first    = s.data() + i;
    const char* last     = s.data() + s.size();
    const auto [ptr, ec] = std::from_chars(first, last, value, std::chars_format::fixed);
    if (ec != std::errc{} || ptr != last) {
        return false;
    }
    out = s[0] == '-' ? -value : value;
    return true;
}

// Reports a data line that cannot be used, by its number and the field at fault. (A line that is
// skipped instead goes unnoticed, and the table then interpolates across the gap.)
[[noreturn]] void failField(const char* format, int lineNumber, const char* what, std::string_view text) {
    throw Error(std::string("EopTable: ") + format + " line " + std::to_string(lineNumber) + ": malformed " + what
                + " field '" + std::string(text) + "'");
}

// The fixed-width field [begin, end) (0-based) of a finals2000A line. Returns false if it is blank
// (the prediction rows leave fields out); throws if it holds anything but a plain number.
bool field(const std::string& line, std::size_t begin, std::size_t end, double& out, int lineNumber, const char* what) {
    if (line.size() <= begin) {
        return false;
    }
    std::string_view s = std::string_view(line).substr(begin, end - begin);
    while (!s.empty() && (s.front() == ' ' || s.front() == '\t')) {
        s.remove_prefix(1);
    }
    while (!s.empty() && (s.back() == ' ' || s.back() == '\t' || s.back() == '\r')) {
        s.remove_suffix(1);
    }
    if (s.empty()) {
        return false;
    }
    if (!plainNumber(s, out)) {
        failField("finals2000A", lineNumber, what, s);
    }
    return true;
}

} // namespace

EopTable::EopTable(std::vector<EopRecord> records, EopInterpolation interp)
    : records_(std::move(records))
    , interp_(interp) {
    if (records_.size() < 2) {
        throw Error("EopTable: need at least two records");
    }
    // An epoch that is not a number has no place in the order: sorting it is undefined behaviour.
    // A value that is not a number would come back out of every interpolation that touches it.
    for (const EopRecord& r : records_) {
        if (!std::isfinite(r.mjd)) {
            throw Error("EopTable: a record's epoch is not finite");
        }
        const EopValues& v = r.values;
        if (!std::isfinite(v.xp + v.yp + v.dut1 + v.lod + v.dX + v.dY)) {
            throw Error("EopTable: the record at MJD " + std::to_string(r.mjd) + " has a value that is not finite");
        }
    }
    std::sort(records_.begin(), records_.end(), [](const EopRecord& a, const EopRecord& b) { return a.mjd < b.mjd; });
    // Two records at one epoch would put a zero in the denominators of the interpolation weights in
    // at(). IERS files have one row per day, so a repeated epoch is a fault in the input: refuse it.
    const auto repeat = std::adjacent_find(records_.begin(),
                                           records_.end(),
                                           [](const EopRecord& a, const EopRecord& b) { return a.mjd == b.mjd; });
    if (repeat != records_.end()) {
        throw Error("EopTable: duplicate epoch MJD " + std::to_string(repeat->mjd));
    }
    const auto leap = leapSecondTable();
    ut1MinusTai_.reserve(records_.size());
    for (const auto& r : records_) {
        ut1MinusTai_.push_back(r.values.dut1 - leap->taiMinusUtc(r.mjd));
    }
}

EopTable EopTable::parseFinals2000A(std::string_view content, EopInterpolation interp) {
    std::vector<EopRecord> recs;
    std::istringstream     in{std::string(content)};
    std::string            line;
    int                    lineNumber = 0;
    while (std::getline(in, line)) {
        ++lineNumber;
        EopRecord r;
        double    xp = 0, yp = 0, dut1 = 0;
        // A row without a date, polar motion or UT1 - UTC is not a data row yet (the file runs on
        // past the last prediction with dates only). A field that is there must be a number.
        const bool hasMjd  = field(line, 7, 15, r.mjd, lineNumber, "MJD");
        const bool hasXp   = field(line, 18, 27, xp, lineNumber, "polar motion x");
        const bool hasYp   = field(line, 37, 46, yp, lineNumber, "polar motion y");
        const bool hasDut1 = field(line, 58, 68, dut1, lineNumber, "UT1-UTC");
        double     lod = 0, dX = 0, dY = 0;
        field(line, 79, 86, lod, lineNumber, "LOD");
        field(line, 97, 106, dX, lineNumber, "dX");
        field(line, 116, 125, dY, lineNumber, "dY");
        if (!hasMjd || !hasXp || !hasYp || !hasDut1) {
            continue;
        }
        r.values = {xp * arcsecToRad, yp * arcsecToRad, dut1, lod * 1e-3, dX * masToRad, dY * masToRad};
        recs.push_back(r);
    }
    if (recs.empty()) {
        throw Error("EopTable: no finals2000A records found");
    }
    return EopTable(std::move(recs), interp);
}

EopTable EopTable::loadFinals2000A(const std::string& path, EopInterpolation interp) {
    return parseFinals2000A(readFile(path), interp);
}

EopTable EopTable::parseC04(std::string_view content, EopInterpolation interp) {
    std::vector<EopRecord> recs;
    std::istringstream     in{std::string(content)};
    std::string            line;
    int                    lineNumber = 0;
    while (std::getline(in, line)) {
        ++lineNumber;
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::istringstream       ls(line);
        std::vector<std::string> f;
        std::string              tok;
        while (ls >> tok) {
            f.push_back(tok);
        }
        if (f.size() < 10) {
            continue;
        }
        // Data lines start with a 4-digit year; anything else is part of the header.
        if (f[0].size() != 4 || f[0].find_first_not_of("0123456789") != std::string::npos) {
            continue;
        }
        // From here on the line is a data line, and each number it must hold is read or the line is
        // reported: a damaged line that was skipped would leave a silent gap in the table.
        const auto number = [&](std::size_t column, const char* what) {
            double value = 0.0;
            if (!plainNumber(f[column], value)) {
                failField("C04", lineNumber, what, f[column]);
            }
            return value;
        };
        EopRecord r;
        // The "20" format has an hour column (1-2 digits) before the MJD.
        const bool format20 = f[3].size() <= 2 && f[4].find('.') != std::string::npos;
        if (format20) { // YYYY MM DD HH MJD x y UT1-UTC dX dY xrt yrt LOD
            if (f.size() < 13) {
                failField("C04", lineNumber, "LOD", "");
            }
            r.mjd    = number(4, "MJD");
            r.values = {number(5, "polar motion x") * arcsecToRad,
                        number(6, "polar motion y") * arcsecToRad,
                        number(7, "UT1-UTC"),
                        number(12, "LOD"),
                        number(8, "dX") * arcsecToRad,
                        number(9, "dY") * arcsecToRad};
        } else { // YYYY MM DD MJD x y UT1-UTC LOD dX dY
            r.mjd    = number(3, "MJD");
            r.values = {number(4, "polar motion x") * arcsecToRad,
                        number(5, "polar motion y") * arcsecToRad,
                        number(6, "UT1-UTC"),
                        number(7, "LOD"),
                        number(8, "dX") * arcsecToRad,
                        number(9, "dY") * arcsecToRad};
        }
        recs.push_back(r);
    }
    if (recs.empty()) {
        throw Error("EopTable: no C04 records found");
    }
    return EopTable(std::move(recs), interp);
}

EopTable EopTable::loadC04(const std::string& path, EopInterpolation interp) {
    return parseC04(readFile(path), interp);
}

EopValues EopTable::at(double mjd) const {
    if (mjd < firstMjd() || mjd > lastMjd()) {
        throw Error("EopTable: MJD " + std::to_string(mjd) + " outside table range [" + std::to_string(firstMjd())
                    + ", " + std::to_string(lastMjd()) + "]");
    }
    const auto it =
        std::upper_bound(records_.begin(), records_.end(), mjd, [](double m, const EopRecord& r) { return m < r.mjd; });
    std::size_t i1 = static_cast<std::size_t>(it - records_.begin());
    if (i1 >= records_.size()) {
        i1 = records_.size() - 1;
    }
    const std::size_t i0 = i1 - 1;

    // Interpolation weights over a set of record indices.
    std::size_t idx[4];
    double      w[4];
    std::size_t n = 0;
    if (interp_ == EopInterpolation::Lagrange4 && records_.size() >= 4) {
        std::size_t start = (i0 >= 1) ? i0 - 1 : 0;
        if (start + 4 > records_.size()) {
            start = records_.size() - 4;
        }
        for (std::size_t k = 0; k < 4; ++k) {
            idx[k] = start + k;
        }
        n = 4;
        for (std::size_t k = 0; k < 4; ++k) {
            double wk = 1.0;
            for (std::size_t m = 0; m < 4; ++m) {
                if (m != k) {
                    wk *= (mjd - records_[idx[m]].mjd) / (records_[idx[k]].mjd - records_[idx[m]].mjd);
                }
            }
            w[k] = wk;
        }
    } else {
        const double t = (mjd - records_[i0].mjd) / (records_[i1].mjd - records_[i0].mjd);
        idx[0]         = i0;
        idx[1]         = i1;
        w[0]           = 1.0 - t;
        w[1]           = t;
        n              = 2;
    }
    EopValues v;
    double    ut1Tai = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        const EopValues& r  = records_[idx[k]].values;
        v.xp               += w[k] * r.xp;
        v.yp               += w[k] * r.yp;
        v.lod              += w[k] * r.lod;
        v.dX               += w[k] * r.dX;
        v.dY               += w[k] * r.dY;
        ut1Tai             += w[k] * ut1MinusTai_[idx[k]];
    }
    v.dut1 = ut1Tai + leapSecondTable()->taiMinusUtc(mjd);
    // Finite records can still interpolate to something that is not a number: epochs a rounding
    // error apart give weights of 1e300. Refuse that rather than return it.
    if (!std::isfinite(v.xp + v.yp + v.dut1 + v.lod + v.dX + v.dY)) {
        throw Error("EopTable: the interpolation at MJD " + std::to_string(mjd)
                    + " is not finite (records too close together, or values too large)");
    }
    return v;
}

std::ostream& operator<<(std::ostream& os, EopInterpolation value) {
    return detail::writeEnum(os, "EopInterpolation", value);
}

} // namespace astro
