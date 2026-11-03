/**********************************************************************
 * ephemeris.h
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
#include "astro/time/time.h"

/// \file
/// The Ephemeris interface and the Body identifiers (NAIF ids).
namespace astro {

/// Solar-system bodies, numbered with their NAIF integer IDs.
enum class Body : int {
    SolarSystemBarycenter = 0,
    MercuryBarycenter     = 1,
    VenusBarycenter       = 2,
    EarthMoonBarycenter   = 3,
    MarsBarycenter        = 4,
    JupiterBarycenter     = 5,
    SaturnBarycenter      = 6,
    UranusBarycenter      = 7,
    NeptuneBarycenter     = 8,
    PlutoBarycenter       = 9,
    Sun                   = 10,
    Mercury               = 199,
    Venus                 = 299,
    Moon                  = 301,
    Earth                 = 399,
    Mars                  = 499,
};

BOOST_DESCRIBE_ENUM(Body,
                    SolarSystemBarycenter,
                    MercuryBarycenter,
                    VenusBarycenter,
                    EarthMoonBarycenter,
                    MarsBarycenter,
                    JupiterBarycenter,
                    SaturnBarycenter,
                    UranusBarycenter,
                    NeptuneBarycenter,
                    PlutoBarycenter,
                    Sun,
                    Mercury,
                    Venus,
                    Moon,
                    Earth,
                    Mars)

/// Writes "Body::<enumerator>", or "Body::INVALID_VALUE" for any other value.
std::ostream& operator<<(std::ostream& os, Body value);

/// Enumerator name of `b` (e.g. "Moon"), or "?" for any other value.
[[nodiscard]] std::string_view toString(Body b) noexcept;
/// NAIF integer ID of `b`.
[[nodiscard]] constexpr int naifId(Body b) noexcept {
    return static_cast<int>(b);
}

/// Source of solar-system body states.
///
/// States are returned in the ICRF/GCRS axes (for a geocentric `center`, this
/// is the GCRS), in meters and meters/second, geometric (no light time or
/// aberration).
class Ephemeris {
public:
    virtual ~Ephemeris() = default;
    /// State of `target` relative to `center` at time t (any scale; converted to TDB).
    [[nodiscard]] virtual StateVector state(Body target, Body center, const Time& t) const = 0;
    /// Position (m) of `target` relative to `center` at time t: the position part of state().
    [[nodiscard]] Vec3 position(Body target, Body center, const Time& t) const { return state(target, center, t).r; }
};

} // namespace astro
