/**********************************************************************
 * geodesy.cpp
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

// Geodesy: WGS84 geodetic <-> ECEF, local ENU and azimuth/elevation/range, geodesic distance.

#include <iomanip>
#include <iostream>

#include "astro/core/angle.h"
#include "astro/geodesy/geodesic.h"
#include "astro/geodesy/geodetic.h"
#include "astro/geodesy/local.h"

using namespace astro;

int main() {
    std::cout << std::fixed << std::setprecision(3);

    //! [ecef]
    // WGS84 latitude, longitude (radians) and height above the ellipsoid (m) -> ECEF (m), and back.
    const Geodetic site{degToRad(38.8895), degToRad(-77.0353), 100.0};
    const Vec3     ecef = geodeticToEcef(site);
    const Geodetic back = ecefToGeodetic(ecef);
    std::cout << "ECEF " << ecef << " m\n";
    std::cout << "lat " << radToDeg(back.lat) << " deg, lon " << radToDeg(back.lon) << " deg, h " << back.h << " m\n";
    //! [ecef]

    //! [aer]
    // Where is a point 10 km east, 5 km north and 3 km up, as seen from the site?
    const Vec3 target = enuToEcef({10e3, 5e3, 3e3}, site);
    const Aer  aer    = enuToAer(ecefToEnu(target, site));
    std::cout << "az " << radToDeg(aer.az) << " deg, el " << radToDeg(aer.el) << " deg, range " << aer.range << " m\n";
    //! [aer]

    //! [geodesic]
    // Ellipsoidal distance and initial azimuth from Washington to Paris (Vincenty).
    const GeodesicInverseResult g = geodesicInverse(site.lat, site.lon, degToRad(48.8566), degToRad(2.3522));
    std::cout << "distance " << g.distance / 1e3 << " km, initial azimuth " << radToDeg(g.azimuth1) << " deg\n";
    //! [geodesic]
    return 0;
}
