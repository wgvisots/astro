/**********************************************************************
 * pole_test.cpp
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

#include <gtest/gtest.h>

#include <cmath>

#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/frames/eop.h"
#include "astro/frames/transform.h"

using namespace astro;

namespace {

const Time kT = Time::fromIso8601("2026-11-02T00:00:00Z");

Vec3 row(const Mat3& m, std::size_t i) {
    return {m(i, 0), m(i, 1), m(i, 2)};
}

} // namespace

TEST(PoleOfDate, GcrsPoleIsTheCip) {
    const Mat3 M = poleOfDateMatrix(Frame::GCRS, kT);
    // The new z axis, expressed in GCRS, is the CIP from the full IAU 2006/2000A model.
    const FrameTransformer full(kT, ConstantEop{});
    EXPECT_LT(angleBetween(row(M, 2), row(full.gcrsTo(Frame::CIRS), 2)), 1e-8); // 2000B: ~1 mas
    // About 0.15 deg from GCRS z in late 2026 (precession since J2000).
    EXPECT_NEAR(rad2deg(std::acos(M(2, 2))), 0.1503, 5e-4);
}

TEST(PoleOfDate, OfDateFramesAreAlreadyAligned) {
    for (Frame f : {Frame::TOD, Frame::TEME, Frame::CIRS}) {
        EXPECT_EQ(poleOfDateMatrix(f, kT), Mat3::identity()) << f;
    }
    // J2000 differs from GCRS by the frame bias only (~20 mas).
    const Vec3 pj = row(poleOfDateMatrix(Frame::J2000, kT), 2);
    const Vec3 pg = row(poleOfDateMatrix(Frame::GCRS, kT), 2);
    EXPECT_LT(std::abs(angleBetween(pj, Vec3{0, 0, 1}) - angleBetween(pg, Vec3{0, 0, 1})), 2e-7);
    // Earth-fixed frames are rejected.
    for (Frame f : {Frame::TIRS, Frame::PEF, Frame::ITRS}) {
        EXPECT_THROW((void)poleOfDateMatrix(f, kT), Error) << f;
        EXPECT_THROW((void)gcrsToCelestialMatrix(f, kT), Error) << f;
    }
}

TEST(PoleOfDate, CelestialMatricesMatchTheFullModel) {
    const FrameTransformer full(kT, ConstantEop{});
    EXPECT_EQ(gcrsToCelestialMatrix(Frame::GCRS, kT), Mat3::identity());
    for (Frame f : {Frame::J2000, Frame::MOD, Frame::TOD, Frame::TEME, Frame::CIRS}) {
        const Mat3 d = gcrsToCelestialMatrix(f, kT) - full.gcrsTo(f);
        double     m = 0.0;
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                m = std::max(m, std::abs(d(i, j)));
            }
        }
        EXPECT_LT(m, 1e-8) << f; // IAU 2000B vs 2000A, ~1 mas
    }
}
