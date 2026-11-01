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

#include <functional>
#include <vector>

#include "astro/core/vec.h"
#include "astro/frames/eop.h"
#include "astro/geodesy/geodetic.h"
#include "astro/orbit/sgp4.h"
#include "astro/time/time.h"

/// Ground-station geometry: look angles, range rate, Doppler, and pass prediction.
namespace astro {

struct LookAngles {
    double az{0.0};         ///< azimuth, rad [0, 2pi), clockwise from north
    double el{0.0};         ///< elevation, rad
    double range{0.0};      ///< m
    double rangeRate{0.0};  ///< m/s, positive receding
};

/// Look angles from a site to a satellite whose state is given in the ITRS (ECEF).
[[nodiscard]] LookAngles lookAngles(const StateVector& satItrs, const Geodetic& site,
                                    const Ellipsoid& ell = Ellipsoid::wgs84()) noexcept;

/// Received-frequency shift (Hz) for transmit frequency `f` and range rate:
/// first order, -f * rangeRate / c.
[[nodiscard]] double dopplerShift(double frequencyHz, double rangeRate) noexcept;
/// Exact relativistic shift for purely radial motion: f (sqrt((1-b)/(1+b)) - 1), b = rangeRate / c.
[[nodiscard]] double dopplerShiftRelativistic(double frequencyHz, double rangeRate) noexcept;

/// A visibility pass of a satellite over a site.
struct Pass {
    Time aos;               ///< acquisition of signal (rise above the mask)
    Time tca;               ///< time of maximum elevation
    Time los;               ///< loss of signal
    double maxElevation{0.0};
    double aosAzimuth{0.0};
    double losAzimuth{0.0};
    bool startsBeforeWindow{false};  ///< already above the mask at the search start
    bool endsAfterWindow{false};     ///< still above the mask at the search end
};

/// Satellite ITRS state as a function of time.
using ItrsStateFunction = std::function<StateVector(const Time&)>;

/// Find passes above `minElevation` in [start, end]. The elevation is sampled every
/// `step` seconds (passes shorter than ~step can be missed) and events are refined
/// with Brent's method.
[[nodiscard]] std::vector<Pass> findPasses(const ItrsStateFunction& sat, const Geodetic& site, const Time& start,
                                           const Time& end, double minElevation = 0.0, double step = 60.0,
                                           const Ellipsoid& ell = Ellipsoid::wgs84());

/// ITRS state function for an SGP4 satellite (TEME -> ITRS with the given EOP;
/// the IAU 2000B nutation model is used for speed).
[[nodiscard]] ItrsStateFunction sgp4ItrsFunction(const Sgp4& sgp4, const EopProvider& eop);

} // namespace astro
