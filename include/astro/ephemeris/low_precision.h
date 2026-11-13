/**********************************************************************
 * low_precision.h
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

#include "astro/ephemeris/ephemeris.h"

/// \file
/// Analytical low-precision Sun and Moon positions; no data files needed.
namespace astro {

/// Geocentric geometric position of the Sun, GCRS axes, meters.
/// Meeus (1998) ch. 25 series, precessed to the GCRS and corrected for the
/// Earth's motion about the Earth-Moon barycenter. Measured against DE440 over
/// 1950-2050: better than 30 arcsec in direction and 6e-5 relative in distance.
/// Throws astro::Error if `t` cannot be converted to TT: a UTC time before
/// 1960, or a UT1 time (there is no Earth orientation data here).
[[nodiscard]] Vec3 sunPositionLowPrecision(const Time& t);

/// Geocentric geometric position of the Moon, GCRS axes, meters.
/// Montenbruck & Gill (2000) truncated lunar series. Measured against DE440
/// over 1950-2050: better than 6 arcmin in direction and 0.14% in distance.
/// Throws astro::Error as sunPositionLowPrecision does.
[[nodiscard]] Vec3 moonPositionLowPrecision(const Time& t);

/// Ephemeris backed by the analytical series. Supports Sun, Moon and Earth
/// relative to each other (velocities by numerical differentiation).
class LowPrecisionEphemeris final : public Ephemeris {
public:
    /// State of `target` relative to `center` at time t, from the series above.
    /// Throws astro::Error for a body other than Sun, Moon and Earth, and as
    /// sunPositionLowPrecision does for `t`.
    [[nodiscard]] StateVector state(Body target, Body center, const Time& t) const override;
};

} // namespace astro
