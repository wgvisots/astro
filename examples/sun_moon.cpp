/**********************************************************************
 * sun_moon.cpp
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

// Sun and Moon: low-precision ephemeris (no files) or a JPL SPK kernel; sunrise and sunset; Moon
// phase; whether a satellite is in Earth's shadow.
//
//   astro_example_sun_moon [de440s.bsp]

#include <iomanip>
#include <iostream>
#include <memory>

#include "astro/astrometry/rise_set.h"
#include "astro/core/angle.h"
#include "astro/ephemeris/low_precision.h"
#include "astro/ephemeris/spk.h"
#include "astro/frames/eop.h"
#include "astro/satellite/eclipse.h"

using namespace astro;

int main(int argc, char** argv) {
    //! [ephemeris]
    // Low precision (Sun < 30", Moon < 6') needs no files; a JPL kernel gives the full DE440 accuracy.
    std::shared_ptr<const Ephemeris> eph = std::make_shared<LowPrecisionEphemeris>();
    if (argc > 1) {
        eph = std::make_shared<SpkEphemeris>(argv[1]);
    }
    const Time t   = Time::fromIso8601("2024-02-01T00:00:00Z");
    const Vec3 sun = eph->position(Body::Sun, Body::Earth, t); // geocentric, m
    std::cout << "Sun distance " << std::fixed << std::setprecision(6) << norm(sun) / constants::au << " au\n";
    //! [ephemeris]

    //! [rise_set]
    const Geodetic site{deg2rad(51.4779), 0.0, 45.0}; // Greenwich
    for (const SkyEvent& e : sunRiseSet(site, t, t + 86400.0, ConstantEop{}, *eph)) {
        const char* what = e.type == SkyEventType::Rise ? "sunrise "
                         : e.type == SkyEventType::Set  ? "sunset  "
                                                        : "noon    ";
        std::cout << what << e.time.iso8601(0) << '\n';
    }
    const MoonPhase phase = moonPhase(t, *eph);
    std::cout << "Moon " << std::setprecision(0) << phase.illuminatedFraction * 100 << "% illuminated, "
              << (phase.waxing ? "waxing" : "waning") << '\n';
    //! [rise_set]

    //! [shadow]
    // A satellite 7000 km from the Earth's center, directly behind it as seen from the Sun.
    const ShadowResult s = shadow(-normalized(sun) * 7000e3, sun);
    std::cout << "satellite: " << s.state << ", illumination " << s.illumination << '\n';
    //! [shadow]
    return 0;
}
