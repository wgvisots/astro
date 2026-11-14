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

/// Kind of event found by findElevationEvents().
enum class SkyEventType {
    Rise,   ///< upward crossing of the threshold elevation
    Set,    ///< downward crossing of the threshold elevation
    Transit ///< local maximum of elevation
};

BOOST_DESCRIBE_ENUM(SkyEventType, Rise, Set, Transit)

/// Writes "SkyEventType::<enumerator>", or "SkyEventType::INVALID_VALUE" for any other value.
std::ostream& operator<<(std::ostream& os, SkyEventType value);

/// A rise, set or transit event.
struct SkyEvent {
    SkyEventType type;      ///< kind of event
    Time         time;      ///< time of the event, in the time scale of the search start
    double       elevation; ///< elevation function value at the event (rad)
};

/// Find upward/downward crossings of `threshold` by `elevation(t)` and its local
/// maxima (transits) in [start, end]. The function is sampled every `step`
/// seconds and events are refined with Brent's method to ~1 ms; features
/// shorter than `step` can be missed. `elevation` is only called for times in
/// [start, end].
///
/// A transit is a maximum strictly inside the window, wherever it falls,
/// including the first and the last sample interval. A window edge is never a
/// transit: where the elevation is highest at `start` or at `end` itself, it
/// cannot be told, from inside the window, whether it turns there or is still
/// falling (or rising) through the edge. So with windows laid end to end, each
/// transit is found in exactly one of them, unless it falls on a shared edge.
///
/// An empty window (end == start) returns no event. A step longer than the
/// window is allowed: the window is always sampled at its two ends and at least
/// once between. Throws astro::Error if `end` is before `start`, if either time
/// or `threshold` is NaN, if `step` is not positive and finite, or if the window
/// holds more than 2e9 steps. Also throws if `end - start` does (one of the two
/// in UT1, or a UTC time before 1960; see Time). An exception thrown by
/// `elevation` passes through.
[[nodiscard]] std::vector<SkyEvent> findElevationEvents(const std::function<double(const Time&)>& elevation,
                                                        const Time&                               start,
                                                        const Time&                               end,
                                                        double                                    threshold,
                                                        double                                    step = 600.0);

/// Twilight definitions, by the Sun's elevation.
enum class Twilight {
    Civil,       ///< Sun 6 deg below the horizon
    Nautical,    ///< Sun 12 deg below the horizon
    Astronomical ///< Sun 18 deg below the horizon
};

BOOST_DESCRIBE_ENUM(Twilight, Civil, Nautical, Astronomical)

/// Writes "Twilight::<enumerator>", or "Twilight::INVALID_VALUE" for any other value.
std::ostream& operator<<(std::ostream& os, Twilight value);

/// Sun elevation defining the twilight (-6, -12 or -18 deg).
[[nodiscard]] double twilightElevation(Twilight tw) noexcept;

/// Geometric topocentric elevation of a body's center (no refraction). Throws
/// astro::Error if `eop` has no data at `t` or `t` is before 1960 (see
/// FrameTransformer), and if `eph` has no data for `body` at `t`.
[[nodiscard]] double
bodyElevation(Body body, const Time& t, const Geodetic& site, const EopProvider& eop, const Ephemeris& eph);

/// Sunrise/sunset (and solar transits) in [start, end] for a horizon elevation.
/// The window is treated as by findElevationEvents(): empty gives no event, and
/// `end` before `start` throws astro::Error. The same holds for twilightEvents() and
/// moonRiseSet(). All three also throw as bodyElevation() does, for any time
/// in the window, and sunRiseSet() for a `horizon` that is NaN.
[[nodiscard]] std::vector<SkyEvent> sunRiseSet(const Geodetic&    site,
                                               const Time&        start,
                                               const Time&        end,
                                               const EopProvider& eop,
                                               const Ephemeris&   eph,
                                               double             horizon = constants::sunriseElevation);
/// Twilight begin (Rise) / end (Set) events. Throws astro::Error as
/// sunRiseSet() does.
[[nodiscard]] std::vector<SkyEvent> twilightEvents(Twilight           tw,
                                                   const Geodetic&    site,
                                                   const Time&        start,
                                                   const Time&        end,
                                                   const EopProvider& eop,
                                                   const Ephemeris&   eph);
/// Moonrise/moonset of the upper limb with 34' refraction (topocentric, so parallax is included).
/// Like sunRiseSet(), it also returns the transits.
/// Throws astro::Error as sunRiseSet() does.
[[nodiscard]] std::vector<SkyEvent>
moonRiseSet(const Geodetic& site, const Time& start, const Time& end, const EopProvider& eop, const Ephemeris& eph);

/// Geocentric lunar phase quantities.
struct MoonPhase {
    double phaseAngle{0.0};          ///< Sun-Moon-Earth angle (rad), 0 = full
    double illuminatedFraction{0.0}; ///< 0..1
    double elongation{0.0};          ///< Moon-Earth-Sun angle (rad)
    double age{0.0};                 ///< ecliptic longitude Moon - Sun, [0, 2pi): 0 new, pi/2 first quarter
    bool   waxing{false};            ///< true while `age` < pi (new to full)
};
/// Geocentric phase of the Moon at `t`, from the Moon and Sun positions in `eph`.
/// Throws astro::Error if `eph` has no data for the Moon or the Sun at `t`.
[[nodiscard]] MoonPhase moonPhase(const Time& t, const Ephemeris& eph);

} // namespace astro
