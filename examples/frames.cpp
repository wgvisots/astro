/**********************************************************************
 * frames.cpp
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

// Reference frames: ECEF (ITRS) <-> ECI (GCRS), TEME, true-of-date, and Earth orientation data.
//
//   astro_example_frames [finals2000A.all]

#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>

#include "astro/core/angle.h"
#include "astro/frames/eop.h"
#include "astro/frames/transform.h"
#include "astro/geodesy/geodetic.h"

using namespace astro;

int main(int argc, char** argv) {
    //! [eop]
    // Earth orientation: an IERS finals2000A file for full accuracy, or ConstantEop{} (UT1 = UTC, no
    // polar motion: good to ~1 s of Earth rotation, ~400 m at the equator).
    std::unique_ptr<EopProvider> eop = std::make_unique<ConstantEop>();
    if (argc > 1) {
        eop = std::make_unique<EopTable>(EopTable::loadFinals2000A(argv[1]));
    }
    //! [eop]

    //! [ecef_eci]
    std::cout << std::fixed << std::setprecision(3);
    const Time        t = Time::fromIso8601("2024-03-20T12:00:00Z");
    const StateVector ecef{geodeticToEcef({degToRad(38.8895), degToRad(-77.0353), 100.0}), {}};
    const StateVector eci  = itrsToGcrs(ecef, t, *eop); // includes the omega x r velocity
    const StateVector back = gcrsToItrs(eci, t, *eop);
    std::cout << "GCRS r " << eci.r << " m, v " << eci.v << " m/s\n";
    std::cout << "round trip error " << std::scientific << norm(back.r - ecef.r) << std::fixed << " m\n";
    //! [ecef_eci]

    //! [transformer]
    // Many vectors at one epoch: build the rotations once.
    const FrameTransformer ft(t, *eop);
    const Vec3             teme = ft.transform(eci.r, Frame::GCRS, Frame::TEME);
    const Vec3             tod  = ft.transform(eci.r, Frame::GCRS, Frame::TOD);
    std::cout << "TEME " << teme << "\nTOD  " << tod << '\n';
    //! [transformer]

    //! [pole]
    // The Earth's true pole is not GCRS z (0.134 deg apart in 2024, growing ~0.0056 deg a year with
    // precession). Zonal gravity must be applied about the true pole.
    const Mat3 toPole = poleOfDateMatrix(Frame::GCRS, t);
    std::cout << std::setprecision(7) << "true pole in GCRS: " << Vec3{toPole(2, 0), toPole(2, 1), toPole(2, 2)}
              << ", tilt " << std::setprecision(4) << radToDeg(std::acos(toPole(2, 2))) << " deg\n";
    //! [pole]
    return 0;
}
