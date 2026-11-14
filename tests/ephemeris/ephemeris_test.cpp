/**********************************************************************
 * ephemeris_test.cpp
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

#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <random>
#include <string>

#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/ephemeris/low_precision.h"
#include "astro/ephemeris/spk.h"
#include "test_util.h"

using namespace astro;
using astro::test::dataPath;
using astro::test::readCsv;

namespace {
Time tdb(double jd) {
    return Time::fromJd(jd, 0.0, TimeScale::TDB);
}

/// Contents of a file as bytes.
std::string fileBytes(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

/// A kernel written to the temporary directory for one test, removed again at the end of it. The name
/// has a random part because several builds of these tests may run at once on one machine.
class TempKernel {
public:
    explicit TempKernel(const std::string& bytes)
        : path_((std::filesystem::temp_directory_path()
                 / ("astro_spk_" + std::to_string(std::random_device{}()) + ".bsp"))
                    .string()) {
        std::ofstream out(path_, std::ios::binary);
        out.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
    }
    ~TempKernel() { std::filesystem::remove(path_); }
    TempKernel(const TempKernel&)            = delete;
    TempKernel& operator=(const TempKernel&) = delete;

    [[nodiscard]] const std::string& path() const { return path_; }

private:
    std::string path_;
};

/// Offset in a little-endian test kernel of the reference-frame id of segment `i`. A summary record
/// is three control doubles and then the summaries; each summary is two doubles (the coverage) and
/// six 32-bit integers: target, center, frame, type, first address, last address.
std::size_t frameIdOffsetLe(const std::string& bytes, std::size_t i) {
    const auto        b = [&](std::size_t k) { return static_cast<std::size_t>(static_cast<unsigned char>(bytes[k])); };
    const std::size_t fward = b(76) | (b(77) << 8) | (b(78) << 16) | (b(79) << 24); // first summary record
    return (fward - 1) * 1024 + 24 + i * 40 + 16 + 8;
}

/// Every state of spk_reference.csv (jplephem, from the full DE440s) for which `use` is true, read
/// from `spk`, with the tolerances of SpkTest.MatchesJplephem. Returns the number compared.
template <class Use>
int expectReferenceStates(const SpkEphemeris& spk, Use use) {
    const auto csv = readCsv(dataPath("spk_reference.csv"));
    int        n   = 0;
    for (const auto& r : csv.rows) {
        const double jd     = std::stod(r[0]);
        const int    center = std::stoi(r[1]), target = std::stoi(r[2]);
        if (!use(center, target)) {
            continue;
        }
        const Vec3        p{std::stod(r[3]) * 1e3, std::stod(r[4]) * 1e3, std::stod(r[5]) * 1e3};
        const Vec3        v{std::stod(r[6]) * 1e3, std::stod(r[7]) * 1e3, std::stod(r[8]) * 1e3};
        const StateVector s = spk.stateById(target, center, (jd - constants::jdJ2000) * 86400.0);
        EXPECT_LT(norm(s.r - p), 1e-12 * norm(p) + 1e-3) << center << "->" << target << " at " << jd;
        EXPECT_LT(norm(s.v - v), 1e-6) << center << "->" << target << " at " << jd;
        ++n;
    }
    return n;
}

/// The message of the astro::Error that `f` throws ("" if it throws nothing).
template <class F>
std::string errorMessage(F f) {
    try {
        f();
    } catch (const Error& e) {
        return e.what();
    }
    return {};
}

bool contains(const std::string& text, const std::string& part) {
    return text.find(part) != std::string::npos;
}
} // namespace

class SpkTest : public ::testing::TestWithParam<const char*> {};

TEST_P(SpkTest, MatchesJplephem) {
    const SpkEphemeris spk(dataPath(GetParam()));
    EXPECT_EQ(spk.segments().size(), 5u);
    EXPECT_EQ(spk.bigEndian(), std::string(GetParam()).find("_be") != std::string::npos);
    const auto csv = readCsv(dataPath("spk_reference.csv"));
    ASSERT_EQ(csv.rows.size(), 125u);
    for (const auto& r : csv.rows) {
        const double      jd     = std::stod(r[0]);
        const int         center = std::stoi(r[1]), target = std::stoi(r[2]);
        const Vec3        p{std::stod(r[3]) * 1e3, std::stod(r[4]) * 1e3, std::stod(r[5]) * 1e3};
        const Vec3        v{std::stod(r[6]) * 1e3, std::stod(r[7]) * 1e3, std::stod(r[8]) * 1e3};
        const StateVector s = spk.stateById(target, center, (jd - constants::jdJ2000) * 86400.0);
        EXPECT_LT(norm(s.r - p), 1e-12 * norm(p) + 1e-3) << center << "->" << target << " at " << jd;
        EXPECT_LT(norm(s.v - v), 1e-6) << center << "->" << target << " at " << jd;
    }
}

INSTANTIATE_TEST_SUITE_P(BothEndians, SpkTest, ::testing::Values("de440s_2024q1_le.bsp", "de440s_2024q1_be.bsp"));

TEST(Spk, ChainingAndErrors) {
    const SpkEphemeris spk(dataPath("de440s_2024q1_le.bsp"));
    const Time         t  = tdb(2460330.5);
    const double       et = ephemerisTime(t);
    // Moon relative to Earth = (EMB->Moon) - (EMB->Earth).
    const StateVector me = spk.state(Body::Moon, Body::Earth, t);
    const StateVector a = spk.stateById(301, 3, et), b = spk.stateById(399, 3, et);
    test::expectVecNear(me.r, a.r - b.r, 1e-6);
    EXPECT_GT(norm(me.r), 3.5e8);
    EXPECT_LT(norm(me.r), 4.1e8);
    // Sun relative to Earth goes through the SSB; reversing the pair negates the state.
    const StateVector se = spk.state(Body::Sun, Body::Earth, t);
    const StateVector es = spk.state(Body::Earth, Body::Sun, t);
    test::expectVecNear(se.r, -es.r, 1e-3);
    test::expectVecNear(se.v, -es.v, 1e-9);
    EXPECT_NEAR(norm(se.r) / constants::au, 0.984, 0.01); // January: near perihelion
    // Type-3 segment (Mars barycenter).
    EXPECT_GT(norm(spk.state(Body::MarsBarycenter, Body::Sun, t).r), 1.3 * constants::au);
    EXPECT_THROW((void)spk.state(Body::Venus, Body::Earth, t), Error);             // not in file
    EXPECT_THROW((void)spk.state(Body::Moon, Body::Earth, tdb(2460500.5)), Error); // out of range
    EXPECT_THROW(SpkEphemeris(dataPath("spk_reference.csv")), Error);
}

TEST(Spk, CorruptFilesThrow) {
    // Malformed kernels must fail with astro::Error, never read out of bounds or hang.
    std::string bytes;
    {
        std::ifstream in(dataPath("de440s_2024q1_le.bsp"), std::ios::binary);
        bytes.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    }
    ASSERT_GT(bytes.size(), 3 * 1024u);
    std::int32_t fward = 0;
    std::memcpy(&fward, bytes.data() + 76, 4);
    const std::size_t summary = static_cast<std::size_t>(fward - 1) * 1024;
    const std::string path =
        (std::filesystem::temp_directory_path()
         / ("astro_corrupt_" + std::to_string(::testing::UnitTest::GetInstance()->random_seed()) + ".bsp"))
            .string();
    auto expectThrows = [&](std::string b, const char* what) {
        {
            std::ofstream out(path, std::ios::binary);
            out.write(b.data(), static_cast<std::streamsize>(b.size()));
        }
        EXPECT_THROW(
            {
                const SpkEphemeris spk(path);
                for (const auto& s : spk.segments()) {
                    (void)spk.stateById(s.target, s.center, 0.5 * (s.startEt + s.endEt));
                }
            },
            Error)
            << what;
    };
    auto withDouble = [&](std::size_t offset, double v) {
        std::string b = bytes;
        std::memcpy(b.data() + offset, &v, 8);
        return b;
    };
    expectThrows(withDouble(summary + 16, 1e6), "summary count too large");
    expectThrows(withDouble(summary + 16, std::nan("")), "summary count NaN");
    expectThrows(withDouble(summary, static_cast<double>(fward)), "summary chain loops");
    expectThrows(withDouble(summary + 24, std::nan("")), "segment start NaN");
    expectThrows(bytes.substr(0, bytes.size() / 2), "truncated");
    // Segment 0's directory: RSIZE = 0.
    std::int32_t ints[6];
    std::memcpy(ints, bytes.data() + summary + 24 + 16, sizeof(ints));
    expectThrows(withDouble(static_cast<std::size_t>(ints[5] - 2) * 8, 0.0), "RSIZE zero");
    std::filesystem::remove(path);
}

TEST(Spk, SegmentInAnotherFrameIsRefused) {
    // A segment's states are in the reference frame that its summary names (NAIF SPK Required Reading,
    // "Segments--The Fundamental SPK Building Blocks"). The reader returns them as ICRF axes, which is
    // right for frame 1 ("J2000") only. Here the Sun segment of the test kernel is relabelled as frame
    // 17 (ECLIPJ2000, tilted 23.4 deg): it must be refused, never returned as if it were frame 1.
    std::string bytes = fileBytes(dataPath("de440s_2024q1_le.bsp"));
    ASSERT_GT(bytes.size(), 3 * 1024u);
    bytes[frameIdOffsetLe(bytes, 0)] = 17;
    const TempKernel   kernel(bytes);
    const SpkEphemeris spk(kernel.path());
    ASSERT_EQ(spk.segments()[0].target, 10);
    ASSERT_EQ(spk.segments()[0].center, 0);
    EXPECT_EQ(spk.segments()[0].frame, 17);

    // Asked for directly: the error names the segment (as the file names it) and its frame.
    const Time        t      = tdb(2460330.5);
    const std::string direct = errorMessage([&] { (void)spk.stateById(10, 0, ephemerisTime(t)); });
    EXPECT_TRUE(contains(direct, "SEG 0->10 type 2")) << direct;
    EXPECT_TRUE(contains(direct, "frame 17")) << direct;
    // As a link in a chain: the Sun relative to the Earth goes through the solar-system barycenter.
    EXPECT_THROW((void)spk.state(Body::Sun, Body::Earth, t), Error);
    EXPECT_THROW((void)spk.state(Body::Earth, Body::Sun, t), Error);
    // States that do not need that segment are still the kernel's (4 other segments, 25 epochs each).
    EXPECT_EQ(expectReferenceStates(spk, [](int, int target) { return target != 10; }), 100);
}

TEST(Spk, OlderNaifDafFileRecordIsAccepted) {
    // Kernels written by early SPICE Toolkits (the de405.bsp that NAIF distributes is one) begin with
    // the id word "NAIF/DAF" instead of "DAF/SPK ", and their file record ends after the FREE pointer:
    // bytes 88 onward are zero, so there is no binary-format word ("BIG-IEEE" / "LTL-IEEE") and no FTP
    // test string. The byte order then has to be found from ND = 2 and NI = 6, as NAIF's own reader and
    // jplephem find it. The rest of the file is laid out the same, so the states must be the kernel's.
    for (const char* name : {"de440s_2024q1_le.bsp", "de440s_2024q1_be.bsp"}) {
        const bool  big   = std::string(name).find("_be") != std::string::npos;
        std::string bytes = fileBytes(dataPath(name));
        ASSERT_GT(bytes.size(), 3 * 1024u);
        const auto all = [](int, int) { return true; };
        bytes.replace(0, 8, "NAIF/DAF");
        { // Only the id word is the older one.
            const TempKernel   kernel(bytes);
            const SpkEphemeris spk(kernel.path());
            EXPECT_EQ(spk.bigEndian(), big) << name;
            EXPECT_EQ(expectReferenceStates(spk, all), 125) << name;
        }
        std::fill(bytes.begin() + 88, bytes.begin() + 1024, '\0');
        { // The whole file record as in de405.bsp.
            const TempKernel   kernel(bytes);
            const SpkEphemeris spk(kernel.path());
            EXPECT_EQ(spk.bigEndian(), big) << name;
            EXPECT_EQ(expectReferenceStates(spk, all), 125) << name;
        }
        std::fill(bytes.begin() + 8, bytes.begin() + 16, '\0');
        { // ND and NI are 2 and 6 in neither byte order: not an SPK file.
            const TempKernel kernel(bytes);
            EXPECT_THROW(SpkEphemeris(kernel.path()), Error) << name;
        }
    }
}

TEST(Spk, MarsIsNotReplacedByItsBarycenter) {
    // DE430 and DE440 (and this fixture, cut from DE440s) have the Mars system barycenter (NAIF id 4)
    // but no segment for Mars itself (499). The reader must not hand back one for the other: the error
    // has to say that body 499 is missing and that the kernel has Body::MarsBarycenter.
    const SpkEphemeris spk(dataPath("de440s_2024q1_le.bsp"));
    const Time         t    = tdb(2460330.5);
    const std::string  mars = errorMessage([&] { (void)spk.state(Body::Mars, Body::Earth, t); });
    EXPECT_TRUE(contains(mars, "no segment for body 499")) << mars;
    EXPECT_TRUE(contains(mars, "body 4 (Body::MarsBarycenter)")) << mars;
    EXPECT_NO_THROW((void)spk.state(Body::MarsBarycenter, Body::Earth, t));
    // Venus (299): the fixture has neither Venus nor its barycenter (2), so there is nothing to offer.
    const std::string venus = errorMessage([&] { (void)spk.state(Body::Venus, Body::Earth, t); });
    EXPECT_TRUE(contains(venus, "no segment for body 299")) << venus;
    EXPECT_FALSE(contains(venus, "Barycenter")) << venus;
    // A body the kernel does have, asked for outside its coverage: the message is about the epoch.
    const std::string moon = errorMessage([&] { (void)spk.state(Body::Moon, Body::Earth, tdb(2460500.5)); });
    EXPECT_TRUE(contains(moon, "no segment for body 301 at ET")) << moon;
    EXPECT_FALSE(contains(moon, "Barycenter")) << moon;
}

TEST(LowPrecision, SunAndMoonAgainstDe440) {
    const SpkEphemeris          spk(dataPath("de440s_2024q1_le.bsp"));
    const LowPrecisionEphemeris lp;
    for (double jd = 2460311.0; jd < 2460370.0; jd += 2.3) {
        const Time t      = tdb(jd);
        const Vec3 sunRef = spk.position(Body::Sun, Body::Earth, t);
        const Vec3 sun    = sunPositionLowPrecision(t);
        EXPECT_LT(rad2deg(angleBetween(sun, sunRef)), 0.01) << jd;
        EXPECT_NEAR(norm(sun) / norm(sunRef), 1.0, 6e-5) << jd;
        const Vec3 moonRef = spk.position(Body::Moon, Body::Earth, t);
        const Vec3 moon    = moonPositionLowPrecision(t);
        EXPECT_LT(rad2deg(angleBetween(moon, moonRef)), 0.1) << jd;
        EXPECT_NEAR(norm(moon) / norm(moonRef), 1.0, 2e-3) << jd;
        // Ephemeris interface + numerical velocity.
        const StateVector m = lp.state(Body::Moon, Body::Earth, t);
        EXPECT_LT(norm(m.v - spk.state(Body::Moon, Body::Earth, t).v), 30.0);
        test::expectVecNear(lp.position(Body::Earth, Body::Sun, t), -sun, 1e-2);
    }
    EXPECT_THROW((void)lp.state(Body::Mars, Body::Earth, tdb(2460311.0)), Error);
}
