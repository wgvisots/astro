/**********************************************************************
 * low_precision_ephemeris.cpp
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

// Low-precision ephemeris study (docs/accuracy.md): the geocentric Sun and Moon of
// LowPrecisionEphemeris against a JPL kernel, every 6 hours over 1950-2050 (TT). Reports the worst
// error in direction and the worst relative error in distance, and when each occurred.
//
// The committed test kernel spans only January-March 2024, so this needs a full kernel:
//
//   scripts/fetch_iers_data.sh --de440s -o data
//   astro_study_low_precision data/de440s.bsp [first year] [last year]   (default 1950 2050)

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "astro/core/angle.h"
#include "astro/ephemeris/low_precision.h"
#include "astro/ephemeris/spk.h"
#include "astro/time/time.h"

using namespace astro;

namespace {

struct Worst {
    double direction{0.0}; // rad
    double distance{0.0};  // relative
    double directionYear{0.0};
    double distanceYear{0.0};
};

void report(const char* body, const Worst& w) {
    std::printf("| %-4s | %8.1f\" (%.2f') in %.1f | %.2e in %.1f |\n",
                body,
                radToDeg(w.direction) * 3600.0,
                radToDeg(w.direction) * 60.0,
                w.directionYear,
                w.distance,
                w.distanceYear);
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: astro_study_low_precision <de440s.bsp> [first year] [last year]\n");
        return 2;
    }
    const int first = argc > 2 ? std::atoi(argv[2]) : 1950;
    const int last  = argc > 3 ? std::atoi(argv[3]) : 2050;

    const SpkEphemeris          de(argv[1]);
    const LowPrecisionEphemeris lp;
    const Time                  start = Time::fromCalendar(TimeScale::TT, {first, 1, 1, 0, 0, 0.0});
    const Time                  end   = Time::fromCalendar(TimeScale::TT, {last, 1, 1, 0, 0, 0.0});
    const double                step  = 6.0 * 3600.0;

    Worst       sun;
    Worst       moon;
    std::size_t samples = 0;
    for (Time t = start; t <= end; t += step, ++samples) {
        const double year = 2000.0 + t.jd().daysSinceJ2000() / 365.25;
        for (auto [body, worst] : {
                 std::pair{Body::Sun,  &sun },
                 std::pair{Body::Moon, &moon}
        }) {
            const Vec3   truth = de.state(body, Body::Earth, t).r;
            const Vec3   model = lp.state(body, Body::Earth, t).r;
            const double dir   = angleBetween(truth, model);
            const double dist  = std::abs(norm(model) - norm(truth)) / norm(truth);
            if (dir > worst->direction) {
                worst->direction     = dir;
                worst->directionYear = year;
            }
            if (dist > worst->distance) {
                worst->distance     = dist;
                worst->distanceYear = year;
            }
        }
    }

    std::printf("Low-precision Sun and Moon against %s, %d-%d, every 6 h (%zu epochs)\n\n",
                argv[1],
                first,
                last,
                samples);
    std::printf("| Body | Worst direction error | Worst relative distance error |\n|---|---|---|\n");
    report("Sun", sun);
    report("Moon", moon);
    return 0;
}
