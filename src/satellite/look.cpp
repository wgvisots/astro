/**********************************************************************
 * look.cpp
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

#include "astro/satellite/look.h"

#include <cmath>

#include "../orbit/trajectory_detail.h"
#include "astro/astrometry/rise_set.h"
#include "astro/core/error.h"
#include "astro/frames/transform.h"
#include "astro/geodesy/local.h"

namespace astro {

double dopplerShift(double f, double rangeRate) noexcept {
    return -f * rangeRate / constants::c;
}

double dopplerShiftRelativistic(double f, double rangeRate) noexcept {
    const double b = rangeRate / constants::c;
    return f * (std::sqrt((1.0 - b) / (1.0 + b)) - 1.0);
}

std::vector<Pass> findPasses(const Trajectory& sat,
                             const Geodetic&   site,
                             const Time&       start,
                             const Time&       end,
                             double            minElevation,
                             double            step,
                             const Ellipsoid&  ell) {
    detail::requireFrame(sat, Frame::ITRS, "findPasses", "the satellite trajectory");
    const double span = end - start;
    if (!std::isfinite(span) || span < 0.0) {
        throw Error("findPasses: start and end must be finite, with end not before start");
    }
    if (!(step > 0.0) || !std::isfinite(step)) {
        throw Error("findPasses: step must be positive and finite");
    }
    if (std::isnan(minElevation)) {
        throw Error("findPasses: minElevation must not be NaN");
    }
    if (!(span / step <= 2e9)) { // the limit of findElevationEvents
        throw Error("findPasses: the window holds more than 2e9 steps");
    }
    std::vector<Pass> passes;
    if (span == 0.0) {
        return passes; // an empty window holds no pass
    }
    auto         elevation = [&](const Time& t) { return ecefToLookAngles(sat.state(t), site, ell).el; };
    const auto   events    = findElevationEvents(elevation, start, end, minElevation, step);
    const double elStart   = elevation(start);
    bool         open      = elStart > minElevation;
    Pass         cur;
    if (open) {
        // A pass cut by the start of the window. Its maximum within the window is the elevation at
        // the start, unless a transit inside the window is higher (the pass was still rising).
        cur.aos                = start;
        cur.startsBeforeWindow = true;
        cur.aosAzimuth         = ecefToLookAngles(sat.state(start), site, ell).az;
        cur.maxElevation       = elStart;
        cur.tca                = start;
    }
    for (const SkyEvent& e : events) {
        switch (e.type) {
        case SkyEventType::Rise:
            cur              = Pass{};
            cur.aos          = e.time;
            cur.aosAzimuth   = ecefToLookAngles(sat.state(e.time), site, ell).az;
            cur.maxElevation = -1e9;
            open             = true;
            break;
        case SkyEventType::Transit:
            if (open && e.elevation > cur.maxElevation) {
                cur.maxElevation = e.elevation;
                cur.tca          = e.time;
            }
            break;
        case SkyEventType::Set:
            if (open) {
                cur.los        = e.time;
                cur.losAzimuth = ecefToLookAngles(sat.state(e.time), site, ell).az;
                if (cur.maxElevation < minElevation) { // no interior maximum sampled: use the endpoints
                    cur.maxElevation = minElevation;
                    cur.tca          = cur.aos + 0.5 * (cur.los - cur.aos);
                }
                passes.push_back(cur);
                open = false;
            }
            break;
        }
    }
    if (open) {
        // A pass cut by the end of the window: likewise, the elevation at the end is its maximum
        // within the window unless a transit (or the start) was higher.
        cur.los             = end;
        cur.endsAfterWindow = true;
        cur.losAzimuth      = ecefToLookAngles(sat.state(end), site, ell).az;
        const double elEnd  = elevation(end);
        if (elEnd > cur.maxElevation) {
            cur.maxElevation = elEnd;
            cur.tca          = end;
        }
        passes.push_back(cur);
    }
    return passes;
}

} // namespace astro
