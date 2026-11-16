/**********************************************************************
 * data_files.cpp
 *
 * Copyright 2024-2026 W. Greg Visotski
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

// External data files: load whichever of the leap-second table, the Earth orientation data and a JPL
// ephemeris are present in a directory, fall back to what is built in for the others, and check
// how long the loaded data remains valid.
//
//   astro_example_data_files [data directory]     (as written by scripts/fetch_iers_data.sh)

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>

#include "astro/ephemeris/low_precision.h"
#include "astro/ephemeris/spk.h"
#include "astro/frames/eop.h"
#include "astro/time/leap_seconds.h"
#include "astro/time/time.h"

using namespace astro;
namespace fs = std::filesystem;

int main(int argc, char** argv) {
    const fs::path dir = argc > 1 ? fs::path(argv[1]) : fs::path("data");

    //! [startup]
    // Once at start-up, before any UTC conversion: a newer leap-second table, if there is one.
    if (fs::exists(dir / "Leap_Second.dat")) {
        setLeapSecondTable(LeapSecondTable::load((dir / "Leap_Second.dat").string()));
    }

    // Earth orientation: the IERS file, or UT1 = UTC with no polar motion (about 400 m at the equator).
    std::shared_ptr<const EopProvider> eop = std::make_shared<ConstantEop>();
    if (fs::exists(dir / "finals2000A.all")) {
        eop = std::make_shared<EopTable>(EopTable::loadFinals2000A((dir / "finals2000A.all").string()));
    }

    // Sun, Moon and planets: a JPL kernel, or the built-in low-precision Sun and Moon.
    std::shared_ptr<const Ephemeris> eph = std::make_shared<LowPrecisionEphemeris>();
    if (fs::exists(dir / "de440s.bsp")) {
        eph = std::make_shared<SpkEphemeris>((dir / "de440s.bsp").string());
    }
    //! [startup]

    //! [validity]
    // How long is the loaded data good for? Check at start-up, not when a conversion fails.
    const double today = Time::now().mjd(); // MJD (UTC)
    if (const auto expires = leapSecondTable()->expiresMjd()) {
        std::cout << "leap-second table expires in " << std::fixed << std::setprecision(0) << (*expires - today)
                  << " days\n";
    } else {
        std::cout << "built-in leap-second table: no expiry date; it ends with the 2017 leap second\n";
    }
    if (const auto* table = dynamic_cast<const EopTable*>(eop.get())) {
        std::cout << "EOP data (with predictions) run out in " << std::fixed << std::setprecision(0)
                  << (table->lastMjd() - today) << " days\n";
    } else {
        std::cout << "no EOP file: UT1 = UTC, no polar motion\n";
    }
    //! [validity]
    return 0;
}
