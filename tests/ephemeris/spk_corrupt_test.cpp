/**********************************************************************
 * spk_corrupt_test.cpp
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

// The SPK reader on a kernel whose data are damaged: a state must be finite or refused, never a
// not-a-number returned as a position. (Found by fuzzing the reader; see fuzz/spk.cpp.)

#include <unistd.h>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "astro/core/error.h"
#include "astro/ephemeris/spk.h"

using namespace astro;

namespace {

std::vector<char> readAll(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

/// A file of its own for each call, removed when it goes out of scope.
struct TempFile {
    std::string path;
    explicit TempFile(const std::vector<char>& bytes) {
        static int counter = 0;
        path               = (std::filesystem::temp_directory_path()
                / ("astro_spk_corrupt_" + std::to_string(getpid()) + "_" + std::to_string(counter++) + ".bsp"))
                   .string();
        std::ofstream(path, std::ios::binary).write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
    }
    ~TempFile() { std::remove(path.c_str()); }
};

bool isFinite(const StateVector& s) {
    return std::isfinite(s.r.x + s.r.y + s.r.z + s.v.x + s.v.y + s.v.z);
}

/// A 32-bit integer of the file at `offset`, in the file's byte order.
std::int32_t int32At(const std::vector<char>& bytes, std::size_t offset, bool swap) {
    unsigned char b[4];
    std::memcpy(b, bytes.data() + offset, 4);
    if (swap) {
        std::swap(b[0], b[3]);
        std::swap(b[1], b[2]);
    }
    std::int32_t value = 0;
    std::memcpy(&value, b, 4);
    return value;
}

} // namespace

TEST(Spk, DamagedCoefficientsAreRefusedNotReturnedAsNaN) {
    // Where the first segment's data begin is read from the file itself (NAIF's "DAF Required
    // Reading"): bytes 76-79 of the file record give the number of the first summary record; a
    // summary record holds three control numbers and then, for each segment, two doubles and six
    // integers, the fifth being the address (in 8-byte words, from 1) of the segment's first datum.
    // A data record starts with two time tags, which the reader already checked; the coefficients
    // after them were used as they came. Overwrite the first ten with NaN and ask for states.
    for (const char* name : {"de440s_2024q1_le.bsp", "de440s_2024q1_be.bsp"}) {
        std::vector<char> bytes = readAll(std::string(ASTRO_TEST_DATA_DIR) + "/" + name);
        ASSERT_GT(bytes.size(), 4u * 1024u) << name;
        const bool fileIsBig = std::string(bytes.data() + 88, 8) == "BIG-IEEE";
        const bool swap      = fileIsBig != (std::endian::native == std::endian::big);
        const auto summary   = static_cast<std::size_t>(int32At(bytes, 76, swap) - 1) * 1024;
        const auto firstWord = static_cast<std::size_t>(int32At(bytes, summary + 24 + 16 + 4 * 4, swap));
        const auto data      = (firstWord - 1) * 8;
        ASSERT_LT(data + 16 + 80, bytes.size()) << name;
        const double nan = std::numeric_limits<double>::quiet_NaN();
        for (std::size_t k = 0; k < 10; ++k) {
            // A NaN whose eight bytes are stored in either order is still a NaN to the reader: its
            // exponent bits are all ones at whichever end they are read from.
            std::memcpy(bytes.data() + data + 16 + 8 * k, &nan, 8);
            if (swap) {
                std::reverse(bytes.begin() + static_cast<std::ptrdiff_t>(data + 16 + 8 * k),
                             bytes.begin() + static_cast<std::ptrdiff_t>(data + 24 + 8 * k));
            }
        }
        const TempFile     file(bytes);
        const SpkEphemeris spk(file.path);
        int                refused = 0, finite = 0;
        for (const SpkSegmentInfo& seg : spk.segments()) {
            for (double f : {0.0, 0.5, 1.0}) {
                try {
                    const StateVector s =
                        spk.stateById(seg.target, seg.center, seg.startEt + f * (seg.endEt - seg.startEt));
                    EXPECT_TRUE(isFinite(s))
                        << name << ": body " << seg.target << " returned a state that is not finite";
                    ++finite;
                } catch (const Error& e) {
                    EXPECT_NE(std::string(e.what()).find("corrupt"), std::string::npos) << e.what();
                    ++refused;
                }
            }
        }
        EXPECT_GT(refused, 0) << name << ": the damage was not noticed";
        EXPECT_GT(finite, 0) << name << ": the undamaged records must still work";
    }
}
