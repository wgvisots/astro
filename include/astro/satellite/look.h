/**********************************************************************
 * look.h
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
#include "astro/geodesy/geodetic.h"
#include "astro/geodesy/local.h"
#include "astro/orbit/trajectory.h"
#include "astro/time/time.h"

/// \file
/// Ground-station geometry: Doppler shift and pass prediction. The look angles and range rate
/// from a site to a satellite are `ecefToLookAngles` (astro/geodesy/local.h).
namespace astro {

/// Received-frequency shift (Hz) for transmit frequency `f` and range rate:
/// first order, -f * rangeRate / c.
[[nodiscard]] double dopplerShift(double frequencyHz, double rangeRate) noexcept;
/// Exact relativistic shift for purely radial motion: f (sqrt((1-b)/(1+b)) - 1), b = rangeRate / c.
/// At a range rate of +c the result is -f, at -c it is infinite, and beyond the speed of light
/// it is NaN.
[[nodiscard]] double dopplerShiftRelativistic(double frequencyHz, double rangeRate) noexcept;

/// A visibility pass of a satellite over a site.
struct Pass {
    Time   aos;                       ///< acquisition of signal (rise above the mask)
    Time   tca;                       ///< time of maximum elevation
    Time   los;                       ///< loss of signal
    double maxElevation{0.0};         ///< elevation at tca (rad)
    double aosAzimuth{0.0};           ///< azimuth at aos (rad)
    double losAzimuth{0.0};           ///< azimuth at los (rad)
    bool   startsBeforeWindow{false}; ///< already above the mask at the search start
    bool   endsAfterWindow{false};    ///< still above the mask at the search end
};

/// Find passes of satellite `sat` (a Trajectory in Frame::ITRS, such as
/// `sgp4Trajectory(sgp4, eop, Frame::ITRS)`) above `minElevation` in [start, end].
/// The elevation is sampled every `step` seconds (passes shorter than ~step can be
/// missed) and events are refined with Brent's method.
///
/// A pass in progress at `start` or at `end` is cut to the window: its `aos` is
/// `start` (startsBeforeWindow) or its `los` is `end` (endsAfterWindow), and `tca`
/// and `maxElevation` are the time and the value of the highest elevation within
/// the window. That is the cut edge itself when the pass is already descending at
/// `start`, or still rising at `end`.
///
/// An empty window (end == start) returns no pass. Throws astro::Error if `sat` is
/// not in Frame::ITRS, if `end` is before `start`, if either time or `minElevation`
/// is NaN, if `step` is not positive and finite, or if the window holds more than
/// 2e9 steps. Also throws if `end - start` does (one of the two in UT1, or a UTC
/// time before 1960; see Time). An exception thrown by `sat` passes through.
[[nodiscard]] std::vector<Pass> findPasses(const Trajectory& sat,
                                           const Geodetic&   site,
                                           const Time&       start,
                                           const Time&       end,
                                           double            minElevation = 0.0,
                                           double            step         = 60.0,
                                           const Ellipsoid&  ell          = Ellipsoid::wgs84());

} // namespace astro
