/**********************************************************************
 * impact.cpp
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

// Ballistic impact on the WGS84 ellipsoid: where, when, and the impact error ellipse.

#include <cmath>
#include <iomanip>
#include <iostream>
#include <optional>

#include "astro/core/angle.h"
#include "astro/geodesy/local.h"
#include "astro/orbit/impact.h"

using namespace astro;

int main() {
    // A launch from the ellipsoid at 30 N, 100 W: 5 km/s relative to the ground, 35 deg above the
    // horizon, heading north-east. Built here in ITRS and converted to GCRS; normally the state
    // comes from your tracker in an ECI frame.
    const Time        epoch = Time::fromIso8601("2026-03-01T12:00:00Z");
    const ConstantEop eop; // in real use: EopTable::loadFinals2000A("finals2000A.all")
    const Geodetic    site{degToRad(30.0), degToRad(-100.0), 0.0};
    const Mat3        enu = ecefToEnuMatrix(site.lat, site.lon);
    const double      V = 5000.0, gamma = degToRad(35.0), azimuth = degToRad(45.0);
    const Vec3        vEnu{V * std::cos(gamma) * std::sin(azimuth),
                    V * std::cos(gamma) * std::cos(azimuth),
                    V * std::sin(gamma)};
    const StateVector s0 = itrsToGcrs(StateVector{geodeticToEcef(site), enu.transpose() * vEnu}, epoch, eop);

    //! [impact]
    // 1-sigma: 100 m per axis in position, 0.1 m/s per axis in velocity (GCRS, m^2 and m^2/s^2).
    Mat6 P0{};
    for (std::size_t i = 0; i < 3; ++i) {
        P0(i, i)         = 100.0 * 100.0;
        P0(i + 3, i + 3) = 0.1 * 0.1;
    }
    const std::optional<Impact> impact = predictImpact(epoch, s0, P0, Frame::GCRS, eop);
    if (!impact) {
        std::cout << "no impact within ImpactOptions::maxFlightTime\n";
        return 1;
    }
    const ErrorEllipse e1  = impact->uncertainty->ellipse; // 1-sigma: 39 % of impacts
    const ErrorEllipse e95 = e1.scaled(0.95);              // 95 % of impacts
    //! [impact]

    //! [print]
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "impact time      " << impact->time.iso8601() << '\n';
    std::cout << "flight time      " << impact->flightTime << " s\n";
    std::cout << "latitude         " << radToDeg(impact->point.lat) << " deg\n";
    std::cout << "longitude        " << radToDeg(impact->point.lon) << " deg\n";
    std::cout << "impact angle     " << radToDeg(impact->angle) << " deg below horizontal\n";
    std::cout << "speed            " << impact->speed << " m/s (relative to the ground)\n";
    std::cout << std::setprecision(1);
    std::cout << "1-sigma ellipse  " << e1.semiMajor << " x " << e1.semiMinor << " m, rotation "
              << radToDeg(e1.rotation) << " deg (clockwise from north)\n";
    std::cout << "95 % ellipse     " << e95.semiMajor << " x " << e95.semiMinor << " m\n";
    std::cout << std::setprecision(3) << "time sigma       " << impact->uncertainty->timeSigma << " s\n";
    //! [print]

    //! [shared]
    // Many predictions in the same time window: build the Earth rotations once and share them.
    const ItrsRotationSeries rotations(epoch, 3 * 3600.0, eop, ImpactOptions{}.frameOptions);
    ImpactOptions            options;
    options.rotations = &rotations;
    for (const double speed : {4000.0, 4500.0, 5000.0}) {
        const Vec3        v = vEnu * (speed / V);
        const StateVector s = itrsToGcrs(StateVector{geodeticToEcef(site), enu.transpose() * v}, epoch, eop);
        if (const auto hit = predictImpact(epoch, s, Frame::GCRS, eop, options)) {
            std::cout << std::setprecision(0) << speed << std::setprecision(3) << " m/s -> latitude "
                      << radToDeg(hit->point.lat) << ", longitude " << radToDeg(hit->point.lon) << " deg\n";
        }
    }
    //! [shared]
    return 0;
}
