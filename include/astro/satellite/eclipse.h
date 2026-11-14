/**********************************************************************
 * eclipse.h
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
#include <iosfwd>
#include <vector>

#include <boost/describe/enum.hpp>

#include "astro/core/constants.h"
#include "astro/core/vec.h"
#include "astro/time/time.h"

/// \file
/// Earth shadow (eclipse) geometry and the solar beta angle. Positions are
/// geocentric, in any common frame, in meters.
namespace astro {

/// Earth shadow model for shadow() and findEclipses().
enum class ShadowModel {
    Cylindrical, ///< umbra only, Earth's shadow as a cylinder
    Conical,     ///< umbra and penumbra with the finite solar disk
};

BOOST_DESCRIBE_ENUM(ShadowModel, Cylindrical, Conical)

/// Writes "ShadowModel::<enumerator>", or "ShadowModel::INVALID_VALUE" for any other value.
std::ostream& operator<<(std::ostream& os, ShadowModel value);

/// Where a satellite is relative to the Earth's shadow.
enum class ShadowState {
    Sunlit,   ///< the whole solar disk is visible
    Penumbra, ///< the Earth hides part of the solar disk
    Umbra,    ///< the Earth hides the whole solar disk
    Annular   ///< the Earth's disk lies inside the solar disk (a ring of Sun stays visible)
};

BOOST_DESCRIBE_ENUM(ShadowState, Sunlit, Penumbra, Umbra, Annular)

/// Writes "ShadowState::<enumerator>", or "ShadowState::INVALID_VALUE" for any other value.
std::ostream& operator<<(std::ostream& os, ShadowState value);

/// Result of shadow().
struct ShadowResult {
    ShadowState state{ShadowState::Sunlit}; ///< shadow region (only Sunlit or Umbra for the cylindrical model)
    double      illumination{1.0};          ///< fraction of the solar disk visible, 0..1
};

/// Shadow condition of a satellite at `sat` with the Sun at `sun` (geocentric).
/// The conical model follows Montenbruck & Gill (2000), sec. 3.4.2, with a
/// spherical Earth of radius `earthRadius`. A satellite at the Earth's centre or
/// at the Sun's position has no shadow geometry: the illumination is then NaN.
[[nodiscard]] ShadowResult shadow(const Vec3& sat,
                                  const Vec3& sun,
                                  ShadowModel model       = ShadowModel::Conical,
                                  double      earthRadius = constants::wgs84A,
                                  double      sunRadius   = constants::sunRadius) noexcept;

/// Solar beta angle (rad): elevation of the Sun above the orbit plane.
[[nodiscard]] double betaAngle(const StateVector& sat, const Vec3& sun) noexcept;

/// One eclipse found by findEclipses().
struct EclipseInterval {
    Time penumbraEntry; ///< start of the eclipse (shadow entry for the cylindrical model)
    Time umbraEntry;    ///< equals penumbraEntry for the cylindrical model or a penumbral-only eclipse
    Time umbraExit;     ///< equals penumbraExit for the cylindrical model or a penumbral-only eclipse
    Time penumbraExit;  ///< end of the eclipse
    bool umbral{true};  ///< false for a penumbral-only (grazing) eclipse
};

/// Find eclipse intervals in [start, end] given geocentric satellite and Sun
/// position functions (same frame). Sampled every `step` s and refined with Brent.
/// Eclipses already in progress at `start` or `end` are clipped to the window. An
/// eclipse (or a sunlit gap) shorter than `step` can fall between two samples and be
/// missed. The Sun's radius is the nominal one (constants::sunRadius).
/// An empty window (end == start) returns no eclipse. Throws astro::Error if `end`
/// is before `start`, if either time is NaN, if `step` is not positive and finite,
/// or if the window holds more than 2e9 steps. Also throws if `end - start` does
/// (one of the two in UT1, or a UTC time before 1960; see Time). An exception
/// thrown by `satPosition` or `sunPosition` passes through.
[[nodiscard]] std::vector<EclipseInterval> findEclipses(const std::function<Vec3(const Time&)>& satPosition,
                                                        const std::function<Vec3(const Time&)>& sunPosition,
                                                        const Time&                             start,
                                                        const Time&                             end,
                                                        double                                  step = 30.0,
                                                        ShadowModel model       = ShadowModel::Conical,
                                                        double      earthRadius = constants::wgs84A);

} // namespace astro
