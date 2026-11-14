/**********************************************************************
 * rotation_series.cpp
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

#include "astro/frames/rotation_series.h"

#include <algorithm>
#include <cmath>

#include "astro/core/constants.h"
#include "astro/core/error.h"
#include "astro/core/rotation.h"

namespace astro {

namespace {

/// ERA rate per UT1 second: 2 pi * 1.00273781191135448 / 86400 (IERS 2010 eq. 5.15).
constexpr double kEraPerUt1Second = constants::twoPi * 1.00273781191135448 / constants::secondsPerDay;

} // namespace

ItrsRotationSeries::ItrsRotationSeries(const Time&         start,
                                       double              duration,
                                       const EopProvider&  eop,
                                       const FrameOptions& options,
                                       double              nodeSpacing)
    : start_(start)
    , duration_(duration) {
    if (!(duration >= 0.0) || !std::isfinite(duration)) {
        throw Error("ItrsRotationSeries: duration must be finite and non-negative");
    }
    if (!(nodeSpacing > 0.0) || !std::isfinite(nodeSpacing)) {
        throw Error("ItrsRotationSeries: node spacing must be finite and positive");
    }
    // At least one segment, so the rotation rate always comes from two nodes.
    const double span = std::max(duration, 1.0);
    const auto   n    = static_cast<std::size_t>(std::ceil(span / nodeSpacing));
    h_                = span / static_cast<double>(n);

    std::vector<Time> ut1;
    for (std::size_t k = 0; k <= n; ++k) {
        const FrameTransformer ft(start + h_ * static_cast<double>(k), eop, options);
        const Mat3&            cirs = ft.gcrsTo(Frame::CIRS);
        const Mat3&            tirs = ft.gcrsTo(Frame::TIRS);
        q_.push_back(cirs);
        w_.push_back(ft.gcrsTo(Frame::ITRS) * tirs.transpose());
        era_.push_back(ft.earthRotationAngle());
        ut1.push_back(ft.ut1());
    }
    for (std::size_t k = 0; k < n; ++k) {
        // ERA is linear in UT1, so the UT1 elapsed over the segment gives the rotation rate,
        // including the length-of-day variation.
        rate_.push_back(kEraPerUt1Second * (ut1[k + 1] - ut1[k]) / h_);
    }
}

ItrsRotationSeries::Parts ItrsRotationSeries::at(double dt) const {
    if (!(dt >= -1.0 && dt <= duration_ + 1.0)) {
        throw Error("ItrsRotationSeries: time outside the window");
    }
    const double u    = dt / h_;
    const double last = static_cast<double>(rate_.size() - 1);
    const auto   k    = static_cast<std::size_t>(std::clamp(std::floor(u), 0.0, last));
    const double w    = u - static_cast<double>(k);
    const double tau  = dt - h_ * static_cast<double>(k);
    Parts        p;
    p.Q    = q_[k] * (1.0 - w) + q_[k + 1] * w;
    p.W    = w_[k] * (1.0 - w) + w_[k + 1] * w;
    p.era  = era_[k] + rate_[k] * tau;
    p.rate = rate_[k];
    return p;
}

Mat3 ItrsRotationSeries::gcrsToItrsMatrix(double dt) const {
    const Parts p = at(dt);
    return p.W * R3(p.era) * p.Q;
}

Vec3 ItrsRotationSeries::gcrsToItrs(const Vec3& r, double dt) const {
    const Parts p = at(dt);
    return p.W * (R3(p.era) * (p.Q * r));
}

StateVector ItrsRotationSeries::gcrsToItrs(const StateVector& s, double dt) const {
    // TIRS position and velocity relative to the rotating Earth: v_t = R3 Q v - omega x r_t, with
    // omega = rate * z (the CIP); then polar motion.
    const Parts p  = at(dt);
    const Mat3  R  = R3(p.era);
    const Vec3  rt = R * (p.Q * s.r);
    const Vec3  vt = R * (p.Q * s.v) - cross(Vec3{0.0, 0.0, p.rate}, rt);
    return {p.W * rt, p.W * vt};
}

StateVector ItrsRotationSeries::itrsToGcrs(const StateVector& s, double dt) const {
    const Parts p  = at(dt);
    const Mat3  Rt = R3(p.era).transpose();
    const Mat3  Wt = p.W.transpose();
    const Mat3  Qt = p.Q.transpose();
    const Vec3  rt = Wt * s.r;
    const Vec3  vt = Wt * s.v + cross(Vec3{0.0, 0.0, p.rate}, rt);
    return {Qt * (Rt * rt), Qt * (Rt * vt)};
}

double ItrsRotationSeries::rotationRate(double dt) const {
    return at(dt).rate;
}

} // namespace astro
