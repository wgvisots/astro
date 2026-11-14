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
#include <string>

#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/frames/eop.h"
#include "astro/frames/precession_nutation.h"
#include "astro/frames/transform.h"
#include "test_util.h"

using namespace astro;
using astro::test::expectMatNear;

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

TEST(CelestialMatrix, NeedsOnlyTtBefore1960) {
    // Precession and nutation are functions of TT alone, so the celestial frames they define exist at
    // any epoch that has a TT: also before 1960, where UTC (and with it the leap-second table) does
    // not. The truth is each matrix built from TT centuries directly with the same IAU 2006/2000B
    // model (precession_nutation.h, itself checked against ERFA in PrecessionNutation.MatchesErfa).
    constexpr auto   model = PrecessionNutationModel::IAU2006_2000B;
    constexpr double tol   = 1e-15; // same formulas and arguments: rounding only (elements are <= 1)
    for (int year : {1900, 1950}) {
        const std::string what = std::to_string(year);
        const Time        t    = Time::fromCalendar(TimeScale::TT, year, 1, 1, 12);
        const double      T    = t.centuriesSinceJ2000();
        const Mat3        bias = frameBiasMatrix06();
        const Mat3        mod  = biasPrecessionMatrix06(T);
        const CipXY       xy   = cipXY(T, model);
        const Mat3        cirs = celestialToIntermediateMatrix(xy, cioLocatorS06(T, xy));

        EXPECT_EQ(gcrsToCelestialMatrix(Frame::GCRS, t), Mat3::identity()) << what;
        expectMatNear(gcrsToCelestialMatrix(Frame::J2000, t), bias, tol, (what + " J2000").c_str());
        expectMatNear(gcrsToCelestialMatrix(Frame::MOD, t), mod, tol, (what + " MOD").c_str());
        expectMatNear(gcrsToCelestialMatrix(Frame::TOD, t), npbMatrix06(T, model), tol, (what + " TOD").c_str());
        expectMatNear(gcrsToCelestialMatrix(Frame::CIRS, t), cirs, tol, (what + " CIRS").c_str());

        // The pole of date is the z axis of CIRS: from frame f the rotation is CIRS * f^T.
        expectMatNear(poleOfDateMatrix(Frame::GCRS, t), cirs, tol, (what + " pole from GCRS").c_str());
        expectMatNear(poleOfDateMatrix(Frame::J2000, t),
                      cirs * bias.transpose(),
                      tol,
                      (what + " pole from J2000").c_str());
        expectMatNear(poleOfDateMatrix(Frame::MOD, t), cirs * mod.transpose(), tol, (what + " pole from MOD").c_str());
        for (Frame f : {Frame::TOD, Frame::TEME, Frame::CIRS}) {
            EXPECT_EQ(poleOfDateMatrix(f, t), Mat3::identity()) << what << " " << f;
        }

        // The same instant in the other scales that are tied to TT without UTC.
        for (TimeScale scale : {TimeScale::TAI, TimeScale::TDB, TimeScale::GPS}) {
            expectMatNear(gcrsToCelestialMatrix(Frame::MOD, t.to(scale)),
                          mod,
                          tol,
                          (what + " MOD, other scale").c_str());
        }

        // TEME is defined with GMST 1982, a function of UT1, which without Earth-orientation data is
        // taken as UTC: that frame does need an epoch with a UTC. So does an epoch given in UTC,
        // except for the frames that are fixed (GCRS, J2000), which do not look at the epoch.
        EXPECT_THROW((void)gcrsToCelestialMatrix(Frame::TEME, t), Error) << what;
        const Time utc(t.jd(), TimeScale::UTC); // the same Julian date, labelled UTC
        EXPECT_THROW((void)gcrsToCelestialMatrix(Frame::MOD, utc), Error) << what;
        EXPECT_THROW((void)poleOfDateMatrix(Frame::GCRS, utc), Error) << what;
        EXPECT_EQ(gcrsToCelestialMatrix(Frame::J2000, utc), gcrsToCelestialMatrix(Frame::J2000, t)) << what;
    }
}

TEST(CelestialMatrix, SameAsAFrameTransformerWithTheSameModel) {
    // Where UTC exists, the matrices are those of a FrameTransformer with IAU 2006/2000B, no pole
    // offsets and no Earth-orientation data: the route through TT alone changes no result.
    const FrameTransformer ft(kT,
                              ConstantEop{},
                              {.model = PrecessionNutationModel::IAU2006_2000B, .applyPoleOffsets = false});
    constexpr double       tol = 1e-15; // same formulas and arguments: rounding only
    for (Frame f : {Frame::GCRS, Frame::J2000, Frame::MOD, Frame::TOD, Frame::TEME, Frame::CIRS}) {
        expectMatNear(gcrsToCelestialMatrix(f, kT), ft.gcrsTo(f), tol, std::string(toString(f)).c_str());
    }
    for (Frame f : {Frame::GCRS, Frame::J2000, Frame::MOD}) { // the others already have the pole as z
        expectMatNear(poleOfDateMatrix(f, kT), ft.matrix(f, Frame::CIRS), tol, std::string(toString(f)).c_str());
    }
}

TEST(CelestialMatrix, MatchesErfaFrom1900To2100) {
    // ERFA's matrices at 60 epochs from 1902 to 2097 (tests/data/frames.csv), 22 of them before 1960.
    const auto csv = test::readCsv(test::dataPath("frames.csv"));
    ASSERT_GE(csv.rows.size(), 50u);
    int before1960 = 0;
    for (const auto& r : csv.rows) {
        const auto v      = [&](const std::string& n) { return std::stod(r[csv.col(n)]); };
        const auto matrix = [&](const std::string& prefix) {
            Mat3 m;
            for (std::size_t i = 0; i < 3; ++i) {
                for (std::size_t j = 0; j < 3; ++j) {
                    m(i, j) = v(prefix + std::to_string(i) + std::to_string(j));
                }
            }
            return m;
        };
        const Time t(JulianDate{v("tt1"), v("tt2")}, TimeScale::TT);
        before1960 += t.mjd() < 36934.0 ? 1 : 0; // 1960-01-01
        // Frame bias and precession are the same IAU 2006 model on both sides: rounding only (the
        // tolerances of PrecessionNutation.MatchesErfa).
        expectMatNear(gcrsToCelestialMatrix(Frame::J2000, t), matrix("rb"), 1e-15, "J2000");
        expectMatNear(gcrsToCelestialMatrix(Frame::MOD, t), matrix("rbp"), 1e-14, "MOD");
        // True of date and CIRS carry the nutation: ERFA's are IAU 2000A, this function uses 2000B,
        // which is within about 1 mas of 2000A over 1995-2050 and a little looser over this wider
        // span (the largest element difference at these epochs is 1.46 mas).
        constexpr double tol2000B = 2.0 * constants::masToRad;
        const Mat3       cirs     = celestialToIntermediateMatrix({v("x"), v("y")}, v("s"));
        expectMatNear(gcrsToCelestialMatrix(Frame::TOD, t), matrix("npb"), tol2000B, "TOD");
        expectMatNear(gcrsToCelestialMatrix(Frame::CIRS, t), cirs, tol2000B, "CIRS");
    }
    EXPECT_GE(before1960, 20); // the rows that used to throw
}
