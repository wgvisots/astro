/**********************************************************************
 * rise_set.h
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

#include "astro/astrometry/apparent.h"

/// \file
/// Rise, set, transit and twilight times, and lunar phase.
namespace astro {

enum class SkyEventType {
    Rise,
    Set,
    Transit
};

BOOST_DESCRIBE_ENUM(SkyEventType, Rise, Set, Transit)

/// Writes "SkyEventType::<enumerator>", or "SkyEventType::INVALID_VALUE" for any other value.
std::ostream& operator<<(std::ostream& os, SkyEventType value);

struct SkyEvent {
    SkyEventType type;
    Time         time;
    double       elevation; ///< elevation function value at the event (rad)
};

/// Find upward/downward crossings of `threshold` by `elevation(t)` and its local
/// maxima (transits) in [start, end]. The function is sampled every `step`
/// seconds and events are refined with Brent's method to ~1 ms; features
/// shorter than `step` can be missed.
[[nodiscard]] std::vector<SkyEvent> findElevationEvents(const std::function<double(const Time&)>& elevation,
                                                        const Time&                               start,
                                                        const Time&                               end,
                                                        double                                    threshold,
                                                        double                                    step = 600.0);

/// Standard horizon for sunrise/sunset: upper limb with 34' refraction (-0.8333 deg).
inline constexpr double kSunriseElevation = -0.8333 * constants::degToRad;

enum class Twilight {
    Civil,
    Nautical,
    Astronomical
};

BOOST_DESCRIBE_ENUM(Twilight, Civil, Nautical, Astronomical)

/// Writes "Twilight::<enumerator>", or "Twilight::INVALID_VALUE" for any other value.
std::ostream& operator<<(std::ostream& os, Twilight value);

/// Sun elevation defining the twilight (-6, -12 or -18 deg).
[[nodiscard]] double twilightElevation(Twilight tw) noexcept;

/// Geometric topocentric elevation of a body's center (no refraction).
[[nodiscard]] double
bodyElevation(Body body, const Time& t, const Geodetic& site, const EopProvider& eop, const Ephemeris& eph);

/// Sunrise/sunset (and solar transits) in [start, end] for a horizon elevation.
[[nodiscard]] std::vector<SkyEvent> sunRiseSet(const Geodetic&    site,
                                               const Time&        start,
                                               const Time&        end,
                                               const EopProvider& eop,
                                               const Ephemeris&   eph,
                                               double             horizon = kSunriseElevation);
/// Twilight begin (Rise) / end (Set) events.
[[nodiscard]] std::vector<SkyEvent> twilight(Twilight           tw,
                                             const Geodetic&    site,
                                             const Time&        start,
                                             const Time&        end,
                                             const EopProvider& eop,
                                             const Ephemeris&   eph);
/// Moonrise/moonset of the upper limb with 34' refraction (topocentric, so parallax is included).
[[nodiscard]] std::vector<SkyEvent>
moonRiseSet(const Geodetic& site, const Time& start, const Time& end, const EopProvider& eop, const Ephemeris& eph);

struct MoonPhase {
    double phaseAngle{0.0};          ///< Sun-Moon-Earth angle (rad), 0 = full
    double illuminatedFraction{0.0}; ///< 0..1
    double elongation{0.0};          ///< Moon-Earth-Sun angle (rad)
    double age{0.0};                 ///< ecliptic longitude Moon - Sun, [0, 2pi): 0 new, pi/2 first quarter
    bool   waxing{false};
};
[[nodiscard]] MoonPhase moonPhase(const Time& t, const Ephemeris& eph);

} // namespace astro
