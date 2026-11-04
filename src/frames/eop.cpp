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
#include <cmath>
#include <fstream>
#include <sstream>

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

// Parse a fixed-width field [begin, end) (0-based); returns false if blank.
bool field(const std::string& line, std::size_t begin, std::size_t end, double& out) {
    if (line.size() <= begin) {
        return false;
    }
    const std::string s = line.substr(begin, end - begin);
    if (s.find_first_not_of(" \t\r") == std::string::npos) {
        return false;
    }
    try {
        out = std::stod(s);
    } catch (const std::exception&) {
        return false;
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
    std::sort(records_.begin(), records_.end(), [](const EopRecord& a, const EopRecord& b) { return a.mjd < b.mjd; });
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
    while (std::getline(in, line)) {
        EopRecord r;
        double    xp = 0, yp = 0, dut1 = 0;
        if (!field(line, 7, 15, r.mjd) || !field(line, 18, 27, xp) || !field(line, 37, 46, yp)
            || !field(line, 58, 68, dut1)) {
            continue;
        }
        double lod = 0, dX = 0, dY = 0;
        field(line, 79, 86, lod);
        field(line, 97, 106, dX);
        field(line, 116, 125, dY);
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
    while (std::getline(in, line)) {
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
        // Data lines start with a 4-digit year.
        if (f[0].size() != 4 || f[0].find_first_not_of("0123456789") != std::string::npos) {
            continue;
        }
        try {
            EopRecord r;
            // "20" format has an hour column (1-2 digits) before the MJD.
            const bool format20 = f[3].size() <= 2 && f[4].find('.') != std::string::npos;
            if (format20) { // YYYY MM DD HH MJD x y UT1-UTC dX dY xrt yrt LOD
                if (f.size() < 13) {
                    continue;
                }
                r.mjd    = std::stod(f[4]);
                r.values = {std::stod(f[5]) * arcsecToRad,
                            std::stod(f[6]) * arcsecToRad,
                            std::stod(f[7]),
                            std::stod(f[12]),
                            std::stod(f[8]) * arcsecToRad,
                            std::stod(f[9]) * arcsecToRad};
            } else { // YYYY MM DD MJD x y UT1-UTC LOD dX dY
                r.mjd    = std::stod(f[3]);
                r.values = {std::stod(f[4]) * arcsecToRad,
                            std::stod(f[5]) * arcsecToRad,
                            std::stod(f[6]),
                            std::stod(f[7]),
                            std::stod(f[8]) * arcsecToRad,
                            std::stod(f[9]) * arcsecToRad};
            }
            recs.push_back(r);
        } catch (const std::exception&) {
            continue;
        }
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
    return v;
}

std::ostream& operator<<(std::ostream& os, EopInterpolation value) {
    return detail::writeEnum(os, "EopInterpolation", value);
}

} // namespace astro
