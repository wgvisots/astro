/**********************************************************************
 * spk.cpp
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

// Fuzz target: SPK ephemeris kernels (.bsp). The reader takes a file, so each input is written to
// a temporary file first.

#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

#include "astro/ephemeris/spk.h"
#include "common.h"

using namespace astro;

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size) {
    static const std::string path = "/tmp/astro_fuzz_spk_" + std::to_string(getpid()) + ".bsp";
    std::FILE*               file = std::fopen(path.c_str(), "wb");
    if (file == nullptr) {
        return 0;
    }
    if (size > 0) {
        std::fwrite(data, 1, size, file);
    }
    std::fclose(file);
    try {
        const SpkEphemeris spk(path);
        // Ask every segment (a hostile file may list many: the first few are enough) for a state at
        // its ends and its middle. Each must be finite, or refused with astro::Error.
        std::size_t asked = 0;
        for (const SpkSegmentInfo& seg : spk.segments()) {
            if (++asked > 8) {
                break;
            }
            for (double f : {0.0, 0.5, 1.0}) {
                try {
                    const StateVector s =
                        spk.stateById(seg.target, seg.center, seg.startEt + f * (seg.endEt - seg.startEt));
                    fuzz::require(std::isfinite(s.r.x + s.r.y + s.r.z + s.v.x + s.v.y + s.v.z));
                } catch (const Error&) {
                }
            }
        }
    } catch (const Error&) {
    }
    return 0;
}
