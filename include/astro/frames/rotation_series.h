/**********************************************************************
 * rotation_series.h
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

#include <vector>

#include "astro/core/vec.h"
#include "astro/frames/eop.h"
#include "astro/frames/transform.h"
#include "astro/time/time.h"

/// \file
/// GCRS <-> ITRS rotations over a time window, precomputed for fast evaluation at many times.
namespace astro {

/// GCRS <-> ITRS ("ECI <-> ECEF") rotations for any time in a window, at a small fraction of the
/// cost of a FrameTransformer per time.
///
/// The rotation is split as in the IERS 2010 conventions, ITRS = W R3(ERA) Q GCRS:
///
/// | Part | What it is                         | Changes                   | Here                            |
/// |------|------------------------------------|---------------------------|---------------------------------|
/// | Q    | precession-nutation (GCRS -> CIRS) | ~1e-7 rad per hour        | exact at nodes, linear between  |
/// | ERA  | Earth Rotation Angle (from UT1)    | one turn per sidereal day | exact at every time             |
/// | W    | polar motion (TIRS -> ITRS)        | ~1e-8 rad per hour        | exact at nodes, linear between  |
///
/// Nodes are at most `nodeSpacing` seconds apart (1 h by default). Inside the window the rotation
/// agrees with FrameTransformer to a few 1e-11 rad (0.2 mm at the Earth's surface; see the tests).
/// The nodes cost one FrameTransformer each: about 125 us with IAU 2006/2000A, 6 us with 2000B.
/// Build one series per window and share it: it is immutable and thread-safe.
class ItrsRotationSeries {
public:
    /// Rotations for start + dt, 0 <= dt <= duration (seconds, as for Time::operator+).
    /// Throws astro::Error if duration is negative or not finite, nodeSpacing is not positive, or
    /// the EOP do not cover the window.
    ItrsRotationSeries(const Time&         start,
                       double              duration,
                       const EopProvider&  eop,
                       const FrameOptions& options     = {},
                       double              nodeSpacing = 3600.0);

    /// Start of the window.
    [[nodiscard]] const Time& start() const noexcept { return start_; }
    /// Length of the window (s).
    [[nodiscard]] double duration() const noexcept { return duration_; }

    /// Rotation from GCRS to ITRS at start + dt: r_itrs = M r_gcrs. Throws astro::Error if dt is
    /// outside the window by more than 1 s.
    [[nodiscard]] Mat3 gcrsToItrsMatrix(double dt) const;
    /// GCRS position to ITRS at start + dt.
    [[nodiscard]] Vec3 gcrsToItrs(double dt, const Vec3& r) const;
    /// GCRS position and velocity to ITRS at start + dt. The velocity is relative to the rotating
    /// Earth (includes -omega x r), as for FrameTransformer::transform.
    [[nodiscard]] StateVector gcrsToItrs(double dt, const StateVector& s) const;
    /// ITRS position and velocity to GCRS at start + dt (inverse of gcrsToItrs).
    [[nodiscard]] StateVector itrsToGcrs(double dt, const StateVector& s) const;
    /// Earth's rotation rate at start + dt (rad/s). It comes from the UT1 elapsed between nodes, so
    /// it is exactly the rate of the rotation above and includes the length-of-day variation.
    /// (FrameTransformer takes it from the EOP length of day instead; the two agree when the EOP
    /// are consistent, as published series are.)
    [[nodiscard]] double rotationRate(double dt) const;

private:
    struct Parts {
        Mat3   Q;         // GCRS -> CIRS
        Mat3   W;         // TIRS -> ITRS
        double era{0.0};  // Earth Rotation Angle
        double rate{0.0}; // d(ERA)/dt (rad/s)
    };
    [[nodiscard]] Parts at(double dt) const;

    Time                start_;
    double              duration_{0.0};
    double              h_{1.0}; // node spacing actually used
    std::vector<Mat3>   q_;      // Q at the nodes
    std::vector<Mat3>   w_;      // W at the nodes
    std::vector<double> era_;    // ERA at the nodes
    std::vector<double> rate_;   // ERA rate over each segment (rad/s)
};

} // namespace astro
