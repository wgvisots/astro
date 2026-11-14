/**********************************************************************
 * tle.cpp
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

// Fuzz target: two-line element sets, and SGP4 on whatever parses.

#include <cmath>
#include <string>

#include "astro/orbit/sgp4.h"
#include "astro/orbit/tle.h"
#include "common.h"

using namespace astro;

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size) {
    const std::string_view text = fuzz::asText(data, size);
    for (bool verifyChecksum : {true, false}) {
        try {
            for (const Tle& tle : Tle::parseMany(text, verifyChecksum)) {
                // Whatever parses must either be written back as two 69-column lines that parse again
                // to the same lines, or be refused with astro::Error.
                try {
                    const std::string l1 = tle.line1(), l2 = tle.line2();
                    fuzz::require(l1.size() == 69 && l2.size() == 69);
                    const Tle again = Tle::parse(l1, l2, tle.name, true);
                    fuzz::require(again.line1() == l1 && again.line2() == l2);
                } catch (const Error&) {
                }
                // SGP4 must give a finite state with code 0 or 6, an error code, or astro::Error.
                try {
                    const Sgp4 sgp4(tle);
                    for (double dt : {0.0, 60.0, -86400.0, 864000.0}) { // seconds
                        StateVector s;
                        const int   code = sgp4.propagate(dt, s);
                        if (code == 0) {
                            fuzz::require(std::isfinite(s.r.x + s.r.y + s.r.z + s.v.x + s.v.y + s.v.z));
                        }
                    }
                } catch (const Error&) {
                }
            }
        } catch (const Error&) {
        }
    }
    return 0;
}
