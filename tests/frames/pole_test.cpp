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

/// The frame bias from GCRS to J2000 (EME2000) built from its published angles, IERS Conventions
/// (2010) eq. 5.22: B = R1(-eta0) R2(xi0) R3(d alpha0) with xi0 = -16.6170 mas, eta0 = -6.8192 mas
/// (eq. 5.21) and d alpha0 = -14.6 mas (eq. 5.23). The rotations are written out here.
Mat3 publishedFrameBias() {
    const double xi0 = -16.6170 * constants::masToRad, eta0 = -6.8192 * constants::masToRad;
    const double da0 = -14.6 * constants::masToRad;
    const auto   r1 = [](double a) { return Mat3{1, 0, 0, 0, std::cos(a), std::sin(a), 0, -std::sin(a), std::cos(a)}; };
    const auto   r2 = [](double a) { return Mat3{std::cos(a), 0, -std::sin(a), 0, 1, 0, std::sin(a), 0, std::cos(a)}; };
    const auto   r3 = [](double a) { return Mat3{std::cos(a), std::sin(a), 0, -std::sin(a), std::cos(a), 0, 0, 0, 1}; };
    return r1(-eta0) * r2(xi0) * r3(da0);
}

/// Compares a GCRS -> J2000 matrix with publishedFrameBias(). The elements carrying xi0 and eta0
/// ((0,2), (1,2) and their transposes) agree to 0.2 uas (measured: the IAU 2006 bias, taken from the
/// Fukushima-Williams angles at J2000, against eq. 5.22); tolerance 1 uas. Those carrying d alpha0
/// ((0,1), (1,0)) are compared to the printed precision of d alpha0, 0.05 mas. The diagonal is
/// 1 - 3e-15 and is compared to rounding.
void expectPublishedFrameBias(const Mat3& m, const char* what) {
    const Mat3   b   = publishedFrameBias();
    const double uas = constants::uasToRad, mas = constants::masToRad;
    // clang-format off
    const double tol[3][3] = {{1e-15,      0.05 * mas, 1.0 * uas},
                              {0.05 * mas, 1e-15,      1.0 * uas},
                              {1.0 * uas,  1.0 * uas,  1e-15    }};
    // clang-format on
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            EXPECT_NEAR(m(i, j), b(i, j), tol[i][j]) << what << " (" << i << "," << j << ")";
        }
    }
}

} // namespace

TEST(FrameBias, MatchesThePublishedAngles) {
    // The bias is 23 mas in all; leaving it out (an identity matrix) misses by 14.6 to 16.6 mas in the
    // off-diagonal elements, far beyond these tolerances.
    expectPublishedFrameBias(frameBiasMatrix06(), "frameBiasMatrix06");
    expectPublishedFrameBias(gcrsToCelestialMatrix(Frame::J2000, kT), "gcrsToCelestialMatrix(J2000)");
    expectPublishedFrameBias(FrameTransformer(kT, ConstantEop{}).gcrsTo(Frame::J2000), "FrameTransformer");
}

TEST(PoleOfDate, GcrsPoleIsTheCip) {
    const Mat3 M = poleOfDateMatrix(Frame::GCRS, kT);
    // The new z axis, expressed in GCRS, is the CIP from the full IAU 2006/2000A model.
    const FrameTransformer full(kT, ConstantEop{});
    EXPECT_LT(angleBetween(row(M, 2), row(full.gcrsTo(Frame::CIRS), 2)), 1e-8); // 2000B: ~1 mas
    // Its angle from GCRS z, from the polynomial part of the CIP coordinates (IERS Conventions 2010,
    // eq. 5.16) and the largest nutation term, of the Moon's node Omega (18.6 years; table 5.3a:
    // 17.2064161" in longitude, 9.2052331" in obliquity), X += sin(eps0) dpsi, Y += deps:
    const double t    = kT.centuriesSinceJ2000();
    const double as   = constants::arcsecToRad;
    const double node = deg2rad(125.04455501) - 6962890.5431 * as * t; // eq. 5.43
    const double eps0 = 84381.406 * as;
    const double X =
        (-0.016617 + 2004.191898 * t - 0.4297829 * t * t) * as - std::sin(eps0) * 17.2064161 * as * std::sin(node);
    const double Y = (-0.006951 - 0.025896 * t - 22.4072747 * t * t) * as + 9.2052331 * as * std::cos(node);
    // 541" (0.1503 deg) in late 2026. The terms left out move it by 0.47" here (measured); the next
    // largest nutation terms reach 0.5" in X and 0.6" in Y. Tolerance 1.5".
    EXPECT_NEAR(std::acos(M(2, 2)), std::asin(std::hypot(X, Y)), 1.5 * as);
}

TEST(PoleOfDate, OfDateFramesAreAlreadyAligned) {
    for (Frame f : {Frame::TOD, Frame::TEME, Frame::CIRS}) {
        EXPECT_EQ(poleOfDateMatrix(f, kT), Mat3::identity()) << f;
    }
    // J2000 differs from GCRS by the frame bias only: the pole expressed in J2000 is the published bias
    // applied to the pole expressed in GCRS. The bias moves the pole by 18 mas; the tolerance is that
    // of the bias's xi0 and eta0 elements (expectPublishedFrameBias()), 1 uas.
    const Vec3 pj = row(poleOfDateMatrix(Frame::J2000, kT), 2);
    const Vec3 pg = row(poleOfDateMatrix(Frame::GCRS, kT), 2);
    test::expectVecNear(pj, publishedFrameBias() * pg, 1.0 * constants::uasToRad, "pole in J2000");
    EXPECT_GT(norm(pj - pg), 15.0 * constants::masToRad); // the bias is there
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
