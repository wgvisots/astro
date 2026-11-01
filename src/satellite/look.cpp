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

#include "astro/astrometry/rise_set.h"
#include "astro/frames/transform.h"
#include "astro/geodesy/local.h"

namespace astro {

LookAngles lookAngles(const StateVector& sat, const Geodetic& site, const Ellipsoid& ell) noexcept {
    const Vec3 rho = sat.r - geodeticToEcef(site, ell);
    const Aer a = enuToAer(ecefToEnuMatrix(site.lat, site.lon) * rho);
    return {a.az, a.el, a.range, dot(rho, sat.v) / a.range};
}

double dopplerShift(double f, double rangeRate) noexcept {
    return -f * rangeRate / constants::c;
}

double dopplerShiftRelativistic(double f, double rangeRate) noexcept {
    const double b = rangeRate / constants::c;
    return f * (std::sqrt((1.0 - b) / (1.0 + b)) - 1.0);
}

std::vector<Pass> findPasses(const ItrsStateFunction& sat, const Geodetic& site, const Time& start, const Time& end,
                             double minElevation, double step, const Ellipsoid& ell) {
    auto elevation = [&](const Time& t) { return lookAngles(sat(t), site, ell).el; };
    const auto events = findElevationEvents(elevation, start, end, minElevation, step);
    std::vector<Pass> passes;
    bool open = elevation(start) > minElevation;
    Pass cur;
    if (open) {
        cur.aos = start;
        cur.startsBeforeWindow = true;
        cur.aosAzimuth = lookAngles(sat(start), site, ell).az;
        cur.maxElevation = -1e9;
    }
    for (const SkyEvent& e : events) {
        switch (e.type) {
            case SkyEventType::Rise:
                cur = Pass{};
                cur.aos = e.time;
                cur.aosAzimuth = lookAngles(sat(e.time), site, ell).az;
                cur.maxElevation = -1e9;
                open = true;
                break;
            case SkyEventType::Transit:
                if (open && e.elevation > cur.maxElevation) {
                    cur.maxElevation = e.elevation;
                    cur.tca = e.time;
                }
                break;
            case SkyEventType::Set:
                if (open) {
                    cur.los = e.time;
                    cur.losAzimuth = lookAngles(sat(e.time), site, ell).az;
                    if (cur.maxElevation < minElevation) {  // no interior maximum sampled: use the endpoints
                        cur.maxElevation = minElevation;
                        cur.tca = cur.aos + 0.5 * (cur.los - cur.aos);
                    }
                    passes.push_back(cur);
                    open = false;
                }
                break;
        }
    }
    if (open) {
        cur.los = end;
        cur.endsAfterWindow = true;
        cur.losAzimuth = lookAngles(sat(end), site, ell).az;
        const double elEnd = elevation(end);
        if (elEnd > cur.maxElevation) {
            cur.maxElevation = elEnd;
            cur.tca = end;
        }
        if (cur.startsBeforeWindow) {
            const double elStart = elevation(start);
            if (elStart > cur.maxElevation) {
                cur.maxElevation = elStart;
                cur.tca = start;
            }
        }
        passes.push_back(cur);
    }
    return passes;
}

ItrsStateFunction sgp4ItrsFunction(const Sgp4& sgp4, const EopProvider& eop) {
    return [&sgp4, &eop](const Time& t) {
        const FrameTransformer ft(t, eop, {.model = PrecessionNutationModel::IAU2006_2000B});
        return ft.transform(sgp4.propagate(t), Frame::TEME, Frame::ITRS);
    };
}

} // namespace astro
