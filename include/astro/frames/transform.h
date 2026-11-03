/**********************************************************************
 * transform.h
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

#pragma once

#include <iosfwd>
#include <string_view>

#include <boost/describe/enum.hpp>

#include "astro/core/vec.h"
#include "astro/frames/eop.h"
#include "astro/frames/precession_nutation.h"
#include "astro/time/time.h"

/// Celestial <-> terrestrial reference frame transformations.
///
/// The primary chain is the IERS 2010 CIO-based transformation
///     [ITRS] = W^T(xp, yp, s') R3(ERA) C(X+dX, Y+dY, s) [GCRS]
/// with the classical equinox-based frames (J2000/EME2000, MOD, TOD, PEF) and
/// TEME (the SGP4 output frame) available for interoperability.
namespace astro {

enum class Frame {
    GCRS,  ///< Geocentric Celestial Reference System ("ECI")
    J2000, ///< mean equator and equinox of J2000.0 (EME2000), frame bias applied
    MOD,   ///< mean equator and equinox of date
    TOD,   ///< true equator and equinox of date
    TEME,  ///< true equator, mean equinox (SGP4 / two-line elements)
    CIRS,  ///< Celestial Intermediate Reference System
    TIRS,  ///< Terrestrial Intermediate Reference System (Earth-fixed, before polar motion)
    PEF,   ///< pseudo Earth-fixed (equinox-based counterpart of TIRS)
    ITRS,  ///< International Terrestrial Reference System ("ECEF")
};

BOOST_DESCRIBE_ENUM(Frame, GCRS, J2000, MOD, TOD, TEME, CIRS, TIRS, PEF, ITRS)

/// Writes "Frame::<enumerator>", or "Frame::INVALID_VALUE" for any other value.
std::ostream& operator<<(std::ostream& os, Frame value);

[[nodiscard]] std::string_view toString(Frame f) noexcept;
/// True for frames that rotate with the Earth (TIRS, PEF, ITRS).
[[nodiscard]] bool isEarthFixed(Frame f) noexcept;

struct FrameOptions {
    PrecessionNutationModel model{PrecessionNutationModel::IAU2006_2000A};
    /// Apply the EOP celestial pole offsets dX, dY.
    bool applyPoleOffsets{true};
};

/// All frame rotations at one epoch. Construct once per epoch and reuse; the
/// constructor evaluates the full precession-nutation series.
class FrameTransformer {
public:
    /// @param t    epoch in any time scale (converted internally to TT, UT1, UTC)
    /// @param eop  Earth orientation data (e.g. EopTable, or ConstantEop{} for none)
    FrameTransformer(const Time& t, const EopProvider& eop, const FrameOptions& options = {});

    /// Rotation from GCRS to frame f: r_f = gcrsTo(f) * r_gcrs.
    [[nodiscard]] const Mat3& gcrsTo(Frame f) const noexcept;
    /// Rotation between two frames: r_to = matrix(from, to) * r_from.
    [[nodiscard]] Mat3 matrix(Frame from, Frame to) const noexcept;
    /// Earth's angular velocity vector expressed in frame f (zero for non-rotating frames).
    [[nodiscard]] Vec3 angularVelocity(Frame f) const noexcept;

    [[nodiscard]] Vec3 transform(const Vec3& r, Frame from, Frame to) const noexcept;
    /// Position and velocity, including the Earth-rotation (omega x r) term for Earth-fixed frames.
    [[nodiscard]] StateVector transform(const StateVector& s, Frame from, Frame to) const noexcept;

    [[nodiscard]] const Time&      tt() const noexcept { return tt_; }
    [[nodiscard]] const Time&      ut1() const noexcept { return ut1_; }
    [[nodiscard]] const EopValues& eop() const noexcept { return eop_; }
    [[nodiscard]] const CipXY&     cip() const noexcept { return xy_; }
    [[nodiscard]] double           cioLocator() const noexcept { return s_; }
    [[nodiscard]] double           earthRotationAngle() const noexcept { return era_; }
    [[nodiscard]] double           gast() const noexcept { return gast_; }
    [[nodiscard]] double           gmst() const noexcept { return gmst_; }
    [[nodiscard]] double           gmst82() const noexcept { return gmst82_; }

private:
    Time      tt_;
    Time      ut1_;
    EopValues eop_;
    CipXY     xy_;
    double    s_{0.0};
    double    era_{0.0};
    double    gast_{0.0};
    double    gmst_{0.0};
    double    gmst82_{0.0};
    double    omega_{0.0};
    Mat3      polar_;
    Mat3      m_[9];
};

/// Rotation from GCRS to the celestial (non-rotating) frame f at t: r_f = M r_gcrs. Identity for
/// GCRS. Uses IAU 2006/2000B precession-nutation (about 1 mas) and no Earth-orientation data, which
/// suits conversions applied in both directions with the same model. Throws astro::Error for the
/// Earth-fixed frames (TIRS, PEF, ITRS).
[[nodiscard]] Mat3 gcrsToCelestialMatrix(Frame f, const Time& t);

/// Rotation from the celestial frame `from` to a frame whose z axis is the Celestial Intermediate
/// Pole (Earth's true spin axis) at t, for aligning a zonal gravity field (GravityField's
/// toBodyFrame). Exactly the identity when the z axis of `from` already is the true pole (TOD, TEME,
/// CIRS). Same model and restrictions as gcrsToCelestialMatrix.
[[nodiscard]] Mat3 poleOfDateMatrix(Frame from, const Time& t);

/// GCRS -> ITRS rotation matrix ("ECI to ECEF").
[[nodiscard]] Mat3 gcrsToItrsMatrix(const Time& t, const EopProvider& eop, const FrameOptions& options = {});

[[nodiscard]] Vec3        gcrsToItrs(const Vec3& r, const Time& t, const EopProvider& eop);
[[nodiscard]] Vec3        itrsToGcrs(const Vec3& r, const Time& t, const EopProvider& eop);
[[nodiscard]] StateVector gcrsToItrs(const StateVector& s, const Time& t, const EopProvider& eop);
[[nodiscard]] StateVector itrsToGcrs(const StateVector& s, const Time& t, const EopProvider& eop);
[[nodiscard]] StateVector temeToItrs(const StateVector& s, const Time& t, const EopProvider& eop);
[[nodiscard]] StateVector itrsToTeme(const StateVector& s, const Time& t, const EopProvider& eop);
[[nodiscard]] StateVector temeToGcrs(const StateVector& s, const Time& t, const EopProvider& eop);
[[nodiscard]] StateVector gcrsToTeme(const StateVector& s, const Time& t, const EopProvider& eop);

/// Aliases using the common names.
[[nodiscard]] inline Vec3 eciToEcef(const Vec3& r, const Time& t, const EopProvider& eop) {
    return gcrsToItrs(r, t, eop);
}
[[nodiscard]] inline Vec3 ecefToEci(const Vec3& r, const Time& t, const EopProvider& eop) {
    return itrsToGcrs(r, t, eop);
}
[[nodiscard]] inline StateVector eciToEcef(const StateVector& s, const Time& t, const EopProvider& eop) {
    return gcrsToItrs(s, t, eop);
}
[[nodiscard]] inline StateVector ecefToEci(const StateVector& s, const Time& t, const EopProvider& eop) {
    return itrsToGcrs(s, t, eop);
}

} // namespace astro
