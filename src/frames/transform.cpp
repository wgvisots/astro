/**********************************************************************
 * transform.cpp
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

#include "astro/frames/transform.h"

#include "../core/enum_io.h"
#include "astro/core/angle.h"
#include "astro/core/rotation.h"
#include "astro/frames/earth_rotation.h"

namespace astro {

namespace {
constexpr std::size_t idx(Frame f) noexcept {
    return static_cast<std::size_t>(f);
}
} // namespace

std::string_view toString(Frame f) noexcept {
    return detail::enumName(f, "?");
}

bool isEarthFixed(Frame f) noexcept {
    return f == Frame::TIRS || f == Frame::PEF || f == Frame::ITRS;
}

FrameTransformer::FrameTransformer(const Time& t, const EopProvider& eop, const FrameOptions& options) {
    const Time utc = t.to(TimeScale::UTC, eop);
    eop_           = eop.at(utc.mjd());
    if (!options.applyPoleOffsets) {
        eop_.dX = 0.0;
        eop_.dY = 0.0;
    }
    tt_            = t.to(TimeScale::TT, eop);
    ut1_           = t.to(TimeScale::UT1, eop);
    const double T = tt_.centuriesSinceJ2000();

    // CIO-based chain.
    const CipXY xyModel = cipXY(T, options.model);
    s_                  = cioLocatorS06(T, xyModel);
    xy_                 = {xyModel.x + eop_.dX, xyModel.y + eop_.dY};
    const Mat3 c2i      = celestialToIntermediateMatrix(xy_, s_);
    era_                = astro::earthRotationAngle(ut1_.jd());
    const Mat3 rEra     = R3(era_);
    polar_              = polarMotionMatrix(eop_.xp, eop_.yp, tioLocatorSp(T));
    omega_              = kEarthRotationRateEra * (1.0 - eop_.lod / constants::secondsPerDay);

    // Equinox-based chain, with the pole offsets mapped to nutation corrections.
    const Nutation          dd = poleOffsetsToNutation(T, eop_.dX, eop_.dY);
    const FukushimaWilliams fw = fukushimaWilliams06(T);
    const Nutation          n  = nutation(T, options.model);
    const Mat3 npb = fwToMatrix(fw.gammaBar, fw.phiBar, fw.psiBar + n.dpsi + dd.dpsi, fw.epsA + n.deps + dd.deps);
    const Mat3 bp  = fwToMatrix(fw.gammaBar, fw.phiBar, fw.psiBar, fw.epsA);
    gast_          = wrapTwoPi(era_ - equationOfOrigins(npb, cioLocatorS06(T, cipXYFromMatrix(npb))));
    gmst_          = gmst06(ut1_.jd(), tt_.jd());
    gmst82_        = astro::gmst82(ut1_.jd());

    m_[idx(Frame::GCRS)]  = Mat3::identity();
    m_[idx(Frame::J2000)] = frameBiasMatrix06();
    m_[idx(Frame::MOD)]   = bp;
    m_[idx(Frame::TOD)]   = npb;
    m_[idx(Frame::CIRS)]  = c2i;
    m_[idx(Frame::TIRS)]  = rEra * c2i;
    m_[idx(Frame::PEF)]   = R3(gast_) * npb;
    m_[idx(Frame::ITRS)]  = polar_ * m_[idx(Frame::TIRS)];
    m_[idx(Frame::TEME)]  = R3(-gmst82_) * m_[idx(Frame::PEF)];
}

const Mat3& FrameTransformer::gcrsTo(Frame f) const noexcept {
    return m_[idx(f)];
}

Mat3 FrameTransformer::matrix(Frame from, Frame to) const noexcept {
    return m_[idx(to)] * m_[idx(from)].transpose();
}

Vec3 FrameTransformer::angularVelocity(Frame f) const noexcept {
    switch (f) {
    case Frame::TIRS:
    case Frame::PEF:  return {0.0, 0.0, omega_};
    case Frame::ITRS: return polar_ * Vec3{0.0, 0.0, omega_};
    default:          return {};
    }
}

Vec3 FrameTransformer::transform(const Vec3& r, Frame from, Frame to) const noexcept {
    return matrix(from, to) * r;
}

StateVector FrameTransformer::transform(const StateVector& s, Frame from, Frame to) const noexcept {
    // Inertial velocity expressed in `from`, then rotated, then made relative to `to`.
    const Vec3  vInertialFrom = s.v + cross(angularVelocity(from), s.r);
    const Mat3  M             = matrix(from, to);
    StateVector out;
    out.r = M * s.r;
    out.v = M * vInertialFrom - cross(angularVelocity(to), out.r);
    return out;
}

Mat3 gcrsToItrsMatrix(const Time& t, const EopProvider& eop, const FrameOptions& options) {
    return FrameTransformer(t, eop, options).gcrsTo(Frame::ITRS);
}

Vec3 gcrsToItrs(const Vec3& r, const Time& t, const EopProvider& eop) {
    return FrameTransformer(t, eop).transform(r, Frame::GCRS, Frame::ITRS);
}
Vec3 itrsToGcrs(const Vec3& r, const Time& t, const EopProvider& eop) {
    return FrameTransformer(t, eop).transform(r, Frame::ITRS, Frame::GCRS);
}
StateVector gcrsToItrs(const StateVector& s, const Time& t, const EopProvider& eop) {
    return FrameTransformer(t, eop).transform(s, Frame::GCRS, Frame::ITRS);
}
StateVector itrsToGcrs(const StateVector& s, const Time& t, const EopProvider& eop) {
    return FrameTransformer(t, eop).transform(s, Frame::ITRS, Frame::GCRS);
}
StateVector temeToItrs(const StateVector& s, const Time& t, const EopProvider& eop) {
    return FrameTransformer(t, eop).transform(s, Frame::TEME, Frame::ITRS);
}
StateVector itrsToTeme(const StateVector& s, const Time& t, const EopProvider& eop) {
    return FrameTransformer(t, eop).transform(s, Frame::ITRS, Frame::TEME);
}
StateVector temeToGcrs(const StateVector& s, const Time& t, const EopProvider& eop) {
    return FrameTransformer(t, eop).transform(s, Frame::TEME, Frame::GCRS);
}
StateVector gcrsToTeme(const StateVector& s, const Time& t, const EopProvider& eop) {
    return FrameTransformer(t, eop).transform(s, Frame::GCRS, Frame::TEME);
}

std::ostream& operator<<(std::ostream& os, Frame value) {
    return detail::writeEnum(os, "Frame", value);
}

} // namespace astro
