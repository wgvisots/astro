/**********************************************************************
 * rise_set.cpp
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

#include "astro/astrometry/rise_set.h"

#include <algorithm>
#include <cmath>

#include "../core/enum_io.h"
#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/core/rootfind.h"
#include "astro/frames/precession_nutation.h"

namespace astro {

std::vector<SkyEvent> findElevationEvents(const std::function<double(const Time&)>& elevation,
                                          const Time&                               start,
                                          const Time&                               end,
                                          double                                    threshold,
                                          double                                    step) {
    const double span = end - start;
    if (!std::isfinite(span) || span < 0.0) {
        throw Error("findElevationEvents: start and end must be finite, with end not before start");
    }
    if (!(step > 0.0) || !std::isfinite(step)) {
        throw Error("findElevationEvents: step must be positive and finite");
    }
    if (std::isnan(threshold)) {
        throw Error("findElevationEvents: threshold must not be NaN");
    }
    if (!(span / step <= 2e9)) { // the count of steps is kept in an int
        throw Error("findElevationEvents: the window holds more than 2e9 steps");
    }
    std::vector<SkyEvent> events;
    if (span == 0.0) {
        return events; // an empty window holds no event
    }
    const int    n = std::max(2, static_cast<int>(std::ceil(span / step)));
    const double h = span / n;
    auto         f = [&](double s) { return elevation(start + s); };

    std::vector<double> x(static_cast<std::size_t>(n) + 1), y(x.size());
    for (int i = 0; i <= n; ++i) {
        x[static_cast<std::size_t>(i)] = i * h;
        y[static_cast<std::size_t>(i)] = f(i * h);
    }
    // The highest elevation between `lo` and `hi` seconds after start. The minimizer works on the
    // bracket mapped to [1, 2], where a tolerance means the same whether it is taken as relative or
    // as absolute. 2e-8 is just above the square root of machine epsilon, the finest a minimizer can
    // resolve, and is 0.05 ms of a 20 minute bracket.
    auto highest = [&](double lo, double hi) {
        auto       s = [&](double u) { return lo + (u - 1.0) * (hi - lo); };
        const auto m = brentMinimize([&](double u) { return -f(s(u)); }, 1.0, 2.0, 2e-8);
        return SkyEvent{SkyEventType::Transit, start + s(m.x), -m.fx};
    };
    for (std::size_t i = 1; i < x.size(); ++i) {
        const double a = y[i - 1] - threshold, b = y[i] - threshold;
        if ((a < 0.0) != (b < 0.0)) {
            const double root = brentRoot([&](double s) { return f(s) - threshold; }, x[i - 1], x[i], 1e-3);
            events.push_back({a < 0.0 ? SkyEventType::Rise : SkyEventType::Set, start + root, threshold});
        }
        // Local maximum bracketed by samples i-1, i, i+1.
        if (i + 1 < x.size() && y[i] >= y[i - 1] && y[i] > y[i + 1]) {
            events.push_back(highest(x[i - 1], x[i + 1]));
        }
    }
    // The first and the last sample interval have a sample on one side only. There, a window edge
    // that is higher than the sample next to it means one of two things: the elevation is still
    // falling from (or rising to) the edge, or it has a maximum between the two, nearer the edge.
    // The highest point of the interval tells them apart: it is a transit only if it is higher than
    // the edge. So a maximum exactly at the edge is not one.
    for (const std::size_t edge : {std::size_t{0}, x.size() - 1}) {
        const std::size_t inner = edge == 0 ? 1 : edge - 1;
        if (y[edge] > y[inner]) {
            const SkyEvent top = highest(std::min(x[edge], x[inner]), std::max(x[edge], x[inner]));
            if (top.elevation > y[edge]) {
                events.push_back(top);
            }
        }
    }
    std::sort(events.begin(), events.end(), [](const SkyEvent& p, const SkyEvent& q) { return p.time < q.time; });
    return events;
}

double twilightElevation(Twilight tw) noexcept {
    switch (tw) {
    case Twilight::Civil:        return -6.0 * constants::degToRad;
    case Twilight::Nautical:     return -12.0 * constants::degToRad;
    case Twilight::Astronomical: return -18.0 * constants::degToRad;
    }
    return 0.0;
}

double bodyElevation(Body body, const Time& t, const Geodetic& site, const EopProvider& eop, const Ephemeris& eph) {
    // The 2000B model (1 mas) is ample here and ~20x faster.
    const FrameTransformer ft(t, eop, {.model = PrecessionNutationModel::IAU2006_2000B});
    return topocentricFromGcrs(eph.position(body, Body::Earth, t), ft, site).el;
}

std::vector<SkyEvent> sunRiseSet(const Geodetic&    site,
                                 const Time&        start,
                                 const Time&        end,
                                 const EopProvider& eop,
                                 const Ephemeris&   eph,
                                 double             horizon) {
    return findElevationEvents([&](const Time& t) { return bodyElevation(Body::Sun, t, site, eop, eph); },
                               start,
                               end,
                               horizon);
}

std::vector<SkyEvent> twilight(Twilight           tw,
                               const Geodetic&    site,
                               const Time&        start,
                               const Time&        end,
                               const EopProvider& eop,
                               const Ephemeris&   eph) {
    auto ev = sunRiseSet(site, start, end, eop, eph, twilightElevation(tw));
    std::erase_if(ev, [](const SkyEvent& e) { return e.type == SkyEventType::Transit; });
    return ev;
}

std::vector<SkyEvent>
moonRiseSet(const Geodetic& site, const Time& start, const Time& end, const EopProvider& eop, const Ephemeris& eph) {
    constexpr double refr = 34.0 / 60.0 * constants::degToRad;
    auto             f    = [&](const Time& t) {
        const FrameTransformer ft(t, eop, {.model = PrecessionNutationModel::IAU2006_2000B});
        const Aer              a = topocentricFromGcrs(eph.position(Body::Moon, Body::Earth, t), ft, site);
        const double           semidiameter = std::asin(constants::moonRadius / a.range);
        return a.el + semidiameter + refr; // upper limb, refracted
    };
    return findElevationEvents(f, start, end, 0.0, 600.0);
}

MoonPhase moonPhase(const Time& t, const Ephemeris& eph) {
    const Vec3 moon = eph.position(Body::Moon, Body::Earth, t);
    const Vec3 sun  = eph.position(Body::Sun, Body::Earth, t);
    MoonPhase  p;
    p.elongation          = angleBetween(moon, sun);
    p.phaseAngle          = angleBetween(-moon, sun - moon);
    p.illuminatedFraction = 0.5 * (1.0 + std::cos(p.phaseAngle));
    const double eps      = meanObliquity06(0.0);
    const double lonMoon  = equatorialToEcliptic(moon, eps).ra;
    const double lonSun   = equatorialToEcliptic(sun, eps).ra;
    p.age                 = wrapTwoPi(lonMoon - lonSun);
    p.waxing              = p.age < constants::pi;
    return p;
}

std::ostream& operator<<(std::ostream& os, SkyEventType value) {
    return detail::writeEnum(os, "SkyEventType", value);
}

std::ostream& operator<<(std::ostream& os, Twilight value) {
    return detail::writeEnum(os, "Twilight", value);
}

} // namespace astro
