/**********************************************************************
 * sgp4_passes.cpp
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

// Two-line elements -> SGP4 -> passes over a ground station, with look angles and Doppler.

#include <iomanip>
#include <iostream>

#include "astro/core/angle.h"
#include "astro/frames/eop.h"
#include "astro/frames/transform.h"
#include "astro/orbit/sgp4.h"
#include "astro/orbit/tle.h"
#include "astro/orbit/trajectory.h"
#include "astro/satellite/look.h"

using namespace astro;

int main() {
    //! [sgp4]
    // ISS element set (the example from the Wikipedia article on two-line elements).
    const Tle         tle = Tle::parse("1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
                               "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
                               "ISS (ZARYA)");
    const Sgp4        sgp4(tle);
    const ConstantEop eop;
    const Time        t    = sgp4.epoch() + 3600.0;
    const StateVector teme = sgp4.propagate(t); // TEME, m and m/s
    const StateVector gcrs = temeToGcrs(teme, t, eop);
    std::cout << tle.name << " at " << t.iso8601(0) << "\n  TEME " << teme.r << " m\n  GCRS " << gcrs.r << " m\n";
    //! [sgp4]

    //! [passes]
    // Passes above 10 deg over a station during the day after the epoch. findPasses needs the
    // satellite in the Earth-fixed ITRS: sgp4Trajectory converts the SGP4 output (TEME) to it.
    const Geodetic   station{degToRad(38.8895), degToRad(-77.0353), 100.0};
    const Trajectory iss    = sgp4Trajectory(sgp4, eop, Frame::ITRS);
    const auto       passes = findPasses(iss, station, sgp4.epoch(), sgp4.epoch() + 86400.0, degToRad(10.0));
    for (const Pass& p : passes) {
        std::cout << "  AOS " << p.aos.iso8601(0) << "  max el " << std::fixed << std::setprecision(1)
                  << radToDeg(p.maxElevation) << " deg  LOS " << p.los.iso8601(0) << '\n';
    }
    //! [passes]

    //! [look]
    // Look angles and the Doppler shift of a 437 MHz downlink at the first pass's closest approach.
    if (!passes.empty()) {
        const LookAngles look = ecefToLookAngles(iss.state(passes[0].tca), station);
        std::cout << "  at maximum elevation: range " << look.range / 1e3 << " km, range rate " << look.rangeRate
                  << " m/s, Doppler " << dopplerShift(437e6, look.rangeRate) << " Hz\n";
    }
    //! [look]
    return 0;
}
