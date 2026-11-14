/**********************************************************************
 * frames_test.cpp
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

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/core/rotation.h"
#include "astro/frames/earth_rotation.h"
#include "astro/frames/eop.h"
#include "astro/frames/precession_nutation.h"
#include "astro/frames/transform.h"
#include "test_util.h"

using namespace astro;
using astro::test::dataPath;
using astro::test::expectMatNear;
using astro::test::expectOrthonormal;
using astro::test::expectVecNear;
using astro::test::readCsv;

namespace {
constexpr double kUas = constants::uasToRad;

Mat3 matFrom(const test::Csv& csv, const std::vector<std::string>& r, const std::string& prefix) {
    Mat3 m;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            m(i, j) = std::stod(r[csv.col(prefix + std::to_string(i) + std::to_string(j))]);
        }
    }
    return m;
}
} // namespace

TEST(PrecessionNutation, MatchesErfa) {
    const auto csv = readCsv(dataPath("frames.csv"));
    ASSERT_GE(csv.rows.size(), 50u);
    for (const auto& r : csv.rows) {
        auto             v = [&](const char* n) { return std::stod(r[csv.col(n)]); };
        const JulianDate tt{v("tt1"), v("tt2")};
        const JulianDate ut1{v("tt1"), v("ut12")};
        const double     t  = tt.centuriesSinceJ2000();
        const CipXY      xy = cipXY06(t);
        EXPECT_NEAR(xy.x, v("x"), 0.5 * kUas) << "X at t=" << t;
        EXPECT_NEAR(xy.y, v("y"), 0.5 * kUas) << "Y at t=" << t;
        EXPECT_NEAR(cioLocatorS06(t, {v("x"), v("y")}), v("s"), 0.05 * kUas) << "s";
        EXPECT_NEAR(tioLocatorSp(t), v("sp"), 1e-18);
        EXPECT_NEAR(earthRotationAngle(ut1), v("era"), 1e-12);
        EXPECT_NEAR(gmst06(ut1, tt), v("gmst06"), 1e-12);
        EXPECT_NEAR(gmst82(ut1), v("gmst82"), 1e-12);
        // The IERS R06 nutation table and SOFA's nut06a apply the IAU 2006 adjustment differently: they agree
        // to ~0.1 uas near J2000 and drift apart by a few uas at +/-100 years (GAST inherits this via EE).
        EXPECT_NEAR(wrapPi(gast06(ut1, tt) - v("gst06a")), 0.0, 5.0 * kUas);
        EXPECT_NEAR(meanObliquity06(t), v("obl06"), 1e-15);
        const Nutation na = nutationIau2000A(t);
        EXPECT_NEAR(na.dpsi, v("dpsi06a"), 5.0 * kUas) << "dpsi";
        EXPECT_NEAR(na.deps, v("deps06a"), 1.0 * kUas) << "deps";
        const Nutation nb = nutationIau2000B(t);
        EXPECT_NEAR(nb.dpsi, v("dpsi00b"), 1e-15);
        EXPECT_NEAR(nb.deps, v("deps00b"), 1e-15);
        // 2000B is a ~1 mas approximation of 2000A (over 1995-2050; looser over 1900-2100).
        EXPECT_NEAR(nb.dpsi, na.dpsi, 5.0 * constants::masToRad);

        expectMatNear(frameBiasMatrix06(), matFrom(csv, r, "rb"), 1e-15, "bias");
        expectMatNear(precessionMatrix06(t), matFrom(csv, r, "rp"), 1e-14, "precession");
        expectMatNear(biasPrecessionMatrix06(t), matFrom(csv, r, "rbp"), 1e-14, "bias-precession");
        expectMatNear(npbMatrix06(t), matFrom(csv, r, "npb"), 5.0 * kUas, "npb");

        // Full CIO-based GCRS -> ITRS with polar motion and pole offsets.
        const CipXY xyc{xy.x + v("dx"), xy.y + v("dy")};
        const Mat3  c2t = polarMotionMatrix(v("xp"), v("yp"), tioLocatorSp(t)) * R3(earthRotationAngle(ut1))
                       * celestialToIntermediateMatrix(xyc, cioLocatorS06(t, xy));
        expectMatNear(c2t, matFrom(csv, r, "c2t"), 1.0 * kUas, "c2t");
        expectOrthonormal(c2t);
    }
}

TEST(PrecessionNutation, BuildingBlocksMatchErfa) {
    // The parts not covered by MatchesErfa, each against ERFA's values in frames.csv:
    // - the Fukushima-Williams angles through fwToMatrix give ERFA's bias-precession matrix, and their
    //   obliquity is ERFA's obl06;
    // - nutationMatrix(obliquity, nutation) times that matrix is ERFA's NPB matrix (as ERFA's pn06
    //   composes it);
    // - the equation of the origins from ERFA's NPB matrix and s is ERA - GST (ERFA: gst06a = era - eo).
    const auto csv = readCsv(dataPath("frames.csv"));
    for (const auto& r : csv.rows) {
        auto                    v  = [&](const char* n) { return std::stod(r[csv.col(n)]); };
        const double            t  = JulianDate{v("tt1"), v("tt2")}.centuriesSinceJ2000();
        const FukushimaWilliams fw = fukushimaWilliams06(t);
        SCOPED_TRACE(testing::Message() << "t = " << t);
        // The tolerances of the same quantities in MatchesErfa.
        expectMatNear(fwToMatrix(fw.gammaBar, fw.phiBar, fw.psiBar, fw.epsA), matFrom(csv, r, "rbp"), 1e-14, "fw");
        EXPECT_NEAR(fw.epsA, v("obl06"), 1e-15);
        const Mat3 n = nutationMatrix(v("obl06"), {v("dpsi06a"), v("deps06a")});
        expectMatNear(n * matFrom(csv, r, "rbp"), matFrom(csv, r, "npb"), 1e-14, "nutation matrix");
        // Rounding in elements of the matrix (1e-16) and in angles near 2 pi: measured 3.7e-14.
        EXPECT_NEAR(wrapPi(equationOfOrigins(matFrom(csv, r, "npb"), v("s")) - (v("era") - v("gst06a"))), 0.0, 1e-13);
    }
}

TEST(PrecessionNutation, FwToMatrixIsTheProductOfItsFourRotations) {
    // IERS Conventions (2010) eq. 5.40: R1(-eps) R3(-psi) R1(phi) R3(gamma), with the rotations
    // written out (Rk(a) turns the axes by +a about axis k).
    const double g = 0.3, p = 0.4, s = -0.2, e = 0.41;
    const auto   r1 = [](double a) { return Mat3{1, 0, 0, 0, std::cos(a), std::sin(a), 0, -std::sin(a), std::cos(a)}; };
    const auto   r3 = [](double a) { return Mat3{std::cos(a), std::sin(a), 0, -std::sin(a), std::cos(a), 0, 0, 0, 1}; };
    expectMatNear(fwToMatrix(g, p, s, e), r1(-e) * r3(-s) * r1(p) * r3(g), 1e-15);
}

TEST(PrecessionNutation, FundamentalArgumentsArePublishedPolynomials) {
    // IERS Conventions (2010) eq. 5.43 (Delaunay arguments; constant terms in degrees, the others in
    // arcseconds) and 5.44 (planetary mean longitudes and general precession, radians), at J2000 and
    // two other epochs.
    const double as = constants::arcsecToRad;
    for (double t : {0.0, 0.27, -0.83}) {
        SCOPED_TRACE(testing::Message() << "t = " << t);
        const double t2 = t * t, t3 = t2 * t, t4 = t3 * t;
        // clang-format off
        const double expected[14] = {
            (134.96340251 * 3600.0 + 1717915923.2178 * t + 31.8792 * t2 + 0.051635 * t3 - 0.00024470 * t4) * as, // l
            (357.52910918 * 3600.0 +  129596581.0481 * t -  0.5532 * t2 + 0.000136 * t3 - 0.00001149 * t4) * as, // l'
            ( 93.27209062 * 3600.0 + 1739527262.8478 * t - 12.7512 * t2 - 0.001037 * t3 + 0.00000417 * t4) * as, // F
            (297.85019547 * 3600.0 + 1602961601.2090 * t -  6.3706 * t2 + 0.006593 * t3 - 0.00003169 * t4) * as, // D
            (125.04455501 * 3600.0 -    6962890.5431 * t +  7.4722 * t2 + 0.007702 * t3 - 0.00005939 * t4) * as, // Om
            4.402608842 + 2608.7903141574 * t, // Mercury
            3.176146697 + 1021.3285546211 * t, // Venus
            1.753470314 +  628.3075849991 * t, // Earth
            6.203480913 +  334.0612426700 * t, // Mars
            0.599546497 +   52.9690962641 * t, // Jupiter
            0.874016757 +   21.3299104960 * t, // Saturn
            5.481293872 +    7.4781598567 * t, // Uranus
            5.311886287 +    3.8133035638 * t, // Neptune
            0.02438175 * t + 0.00000538691 * t2, // general precession in longitude
        };
        // clang-format on
        const std::array<double, 14> f = fundamentalArguments(t);
        for (std::size_t k = 0; k < 14; ++k) {
            // The Delaunay arguments are evaluated in arcseconds, up to 1.4e9" here: 2.4e-7" (1.2e-12
            // rad) of rounding (measured 1.5e-12).
            EXPECT_NEAR(wrapPi(f[k] - expected[k]), 0.0, 5e-12) << "argument " << k;
        }
    }
}

TEST(PrecessionNutation, CipOfEitherModelMatchesErfa) {
    // cipXY(t, model) against ERFA: for 2000A, ERFA's X, Y; for 2000B, the z row of the NPB matrix that
    // ERFA's 2000B nutation (nut00b) gives with ERFA's bias-precession matrix and obliquity, composed as
    // in BuildingBlocksMatchErfa.
    const auto csv = readCsv(dataPath("frames.csv"));
    for (const auto& r : csv.rows) {
        auto         v = [&](const char* n) { return std::stod(r[csv.col(n)]); };
        const double t = JulianDate{v("tt1"), v("tt2")}.centuriesSinceJ2000();
        SCOPED_TRACE(testing::Message() << "t = " << t);
        const CipXY a = cipXY(t, PrecessionNutationModel::IAU2006_2000A);
        EXPECT_NEAR(a.x, v("x"), 0.5 * kUas); // as for cipXY06 in MatchesErfa
        EXPECT_NEAR(a.y, v("y"), 0.5 * kUas);
        const Mat3  npbB = nutationMatrix(v("obl06"), {v("dpsi00b"), v("deps00b")}) * matFrom(csv, r, "rbp");
        const CipXY b    = cipXY(t, PrecessionNutationModel::IAU2006_2000B);
        // Rounding of matrix elements.
        EXPECT_NEAR(b.x, npbB(2, 0), 1e-14);
        EXPECT_NEAR(b.y, npbB(2, 1), 1e-14);
    }
}

TEST(FrameTransformer, WithEarthOrientationDataMatchesErfa) {
    // ERFA's GCRS -> ITRS matrix with polar motion, UT1 and the pole offsets dX, dY (tests/data/
    // frames.csv), against a FrameTransformer given the same values through an EopProvider: this checks
    // how the transformer applies each parameter, and its accessors. Epochs from 1960 on (UTC needed).
    const auto csv  = readCsv(dataPath("frames.csv"));
    int        rows = 0;
    for (const auto& r : csv.rows) {
        auto       v = [&](const char* n) { return std::stod(r[csv.col(n)]); };
        const Time t(JulianDate{v("tt1"), v("tt2")}, TimeScale::TT);
        if (t.mjd() < 36934.0 + 1.0) {
            continue; // before 1960, where UTC (and with it this route) does not exist
        }
        ++rows;
        SCOPED_TRACE(testing::Message() << "MJD " << t.mjd());
        // UT1 - UTC from ERFA's UT1 (Julian date tt1 + ut12) and the library's UTC.
        const double           dut1 = JulianDate{v("tt1"), v("ut12")}.secondsSince(t.to(TimeScale::UTC).jd());
        const EopValues        e{.xp = v("xp"), .yp = v("yp"), .dut1 = dut1, .lod = 0.0, .dX = v("dx"), .dY = v("dy")};
        const FrameTransformer ft(t, ConstantEop(e));
        // 1 uas: the tolerance of the same matrix built from the parts in PrecessionNutation.MatchesErfa.
        expectMatNear(ft.gcrsTo(Frame::ITRS), matFrom(csv, r, "c2t"), 1.0 * kUas, "c2t");
        // The CIP includes the offsets; s is that of the model's CIP (0.05 uas as in MatchesErfa).
        EXPECT_NEAR(ft.cip().x, v("x") + v("dx"), 0.5 * kUas);
        EXPECT_NEAR(ft.cip().y, v("y") + v("dy"), 0.5 * kUas);
        EXPECT_NEAR(ft.cioLocator(), v("s"), 0.05 * kUas);
        // Earth rotation from UT1: 1e-12 rad as in MatchesErfa, plus the rounding of dut1 through two
        // Julian dates (1e-11 s, 7e-16 rad).
        EXPECT_NEAR(ft.earthRotationAngle(), v("era"), 2e-12);
        EXPECT_NEAR(wrapPi(ft.gmst82() - v("gmst82")), 0.0, 2e-12);
        EXPECT_NEAR(wrapPi(ft.gmst() - v("gmst06")), 0.0, 2e-12);
        EXPECT_EQ(ft.eop().xp, e.xp);
        EXPECT_EQ(ft.eop().dX, e.dX);
        EXPECT_EQ(ft.tt().scale(), TimeScale::TT);
        EXPECT_NEAR(ft.tt().jd().secondsSince(t.jd()), 0.0, 1e-9);
        EXPECT_NEAR(ft.ut1().jd().secondsSince(JulianDate{v("tt1"), v("ut12")}), 0.0, 1e-9);

        // Without the pole offsets: the CIP is the model's, and the offsets read back as zero.
        const FrameTransformer noOffsets(t, ConstantEop(e), {.applyPoleOffsets = false});
        EXPECT_NEAR(noOffsets.cip().x, v("x"), 0.5 * kUas);
        EXPECT_NEAR(noOffsets.cip().y, v("y"), 0.5 * kUas);
        EXPECT_EQ(noOffsets.eop().dX, 0.0);
        EXPECT_EQ(noOffsets.eop().dY, 0.0);
        // ERFA's GAST (gst06a) has no pole offsets either (with them it moves by up to 0.5 mas here).
        EXPECT_NEAR(wrapPi(noOffsets.gast() - v("gst06a")), 0.0, 5.0 * kUas); // as in MatchesErfa
    }
    EXPECT_GE(rows, 30);
}

TEST(FrameTransformer, AngularVelocityIsAlongTheCipAtTheRateOfTheEarthRotationAngle) {
    // The Earth turns about the CIP at d(ERA)/dt = 2 pi 1.00273781191135448 rad per UT1 day (IERS
    // Conventions 2010 eq. 5.15), shortened by the excess length of day: omega (1 - LOD / 86400). In
    // TIRS and PEF the CIP is z; in ITRS it is (xp, -yp, 1) to first order in the polar motion.
    const Time             t  = Time::fromIso8601("2021-07-04T12:34:56Z");
    const double           xp = arcsec2rad(0.2), yp = arcsec2rad(0.4), lod = 0.0015;
    const ConstantEop      eop(EopValues{.xp = xp, .yp = yp, .dut1 = -0.15, .lod = lod});
    const double           omega = constants::twoPi * 1.00273781191135448 / 86400.0 * (1.0 - lod / 86400.0);
    const FrameTransformer ft(t, eop);
    // Rounding: a few units in the last place of 7.3e-5 rad/s (1.4e-20; measured 1.4e-20).
    expectVecNear(ft.angularVelocity(Frame::TIRS), {0.0, 0.0, omega}, 5e-20, "TIRS");
    expectVecNear(ft.angularVelocity(Frame::PEF), {0.0, 0.0, omega}, 5e-20, "PEF");
    // Second order in the polar motion: (xp^2 + yp^2) / 2 = 2.3e-12 of omega, 1.7e-16 rad/s.
    expectVecNear(ft.angularVelocity(Frame::ITRS), {omega * xp, -omega * yp, omega}, 5e-16, "ITRS");
    // The celestial frames do not rotate.
    for (Frame f : {Frame::GCRS, Frame::J2000, Frame::MOD, Frame::TOD, Frame::TEME, Frame::CIRS}) {
        EXPECT_EQ(ft.angularVelocity(f), Vec3{}) << f;
    }
    // So a point fixed in the ITRS moves at omega times its distance from the CIP, eastward.
    const StateVector fixed{
        {6378137.0, 0, 0},
        {}
    };
    const StateVector tirs = ft.transform(fixed, Frame::ITRS, Frame::TIRS);
    expectVecNear(tirs.v, {0.0, 0.0, 0.0}, 1e-12, "fixed in TIRS"); // it rotates with the frame
    const StateVector gcrs = ft.transform(fixed, Frame::ITRS, Frame::GCRS);
    const Vec3        axis = ft.gcrsTo(Frame::TIRS).row(2); // the CIP, in GCRS
    expectVecNear(gcrs.v, omega * cross(axis, gcrs.r), 1e-9, "inertial velocity");
}

TEST(FrameTransformer, EarthFixedFramesAreTirsPefAndItrs) {
    for (Frame f : {Frame::TIRS, Frame::PEF, Frame::ITRS}) {
        EXPECT_TRUE(isEarthFixed(f)) << f;
    }
    for (Frame f : {Frame::GCRS, Frame::J2000, Frame::MOD, Frame::TOD, Frame::TEME, Frame::CIRS}) {
        EXPECT_FALSE(isEarthFixed(f)) << f;
    }
}

TEST(FrameTransformer, EpochWithoutEarthOrientationDataOrBefore1960Throws) {
    // A table covering ten days in 2023, asked for an epoch a year later, and an epoch before 1960
    // (no UTC): the constructor and the one-off conversions throw.
    std::vector<EopRecord> records;
    for (int d = 0; d < 10; ++d) {
        records.push_back({60000.0 + d, {}});
    }
    const EopTable    table(records);
    const Time        inside = Time::fromMjd(60004.5, TimeScale::UTC);
    const Time        later  = Time::fromMjd(60400.0, TimeScale::UTC);
    const Time        early  = Time::fromCalendar(TimeScale::TT, 1955, 1, 1, 0);
    const StateVector s{
        {7000e3, 0,     0},
        {0,      7.5e3, 0}
    };
    EXPECT_NO_THROW(FrameTransformer(inside, table));
    EXPECT_THROW(FrameTransformer(later, table), Error);
    EXPECT_THROW((void)gcrsToItrs(s, later, table), Error);
    EXPECT_THROW((void)temeToItrs(s, later, table), Error);
    EXPECT_THROW(FrameTransformer(early, ConstantEop{}), Error);
    EXPECT_THROW((void)gcrsToItrsMatrix(early, ConstantEop{}), Error);
}

TEST(PrecessionNutation, PoleOffsetsToNutationRoundTrip) {
    const double   t  = 0.2;
    const double   dX = 0.3 * constants::masToRad, dY = -0.2 * constants::masToRad;
    const Nutation n = poleOffsetsToNutation(t, dX, dY);
    // Applying the nutation correction should move the CIP by (dX, dY).
    const CipXY a = cipXYFromMatrix(npbMatrix06(t));
    const CipXY b = cipXYFromMatrix(npbMatrix06(t, PrecessionNutationModel::IAU2006_2000A, n.dpsi, n.deps));
    EXPECT_NEAR(b.x - a.x, dX, 1e-3 * constants::masToRad);
    EXPECT_NEAR(b.y - a.y, dY, 1e-3 * constants::masToRad);
}

// Vallado, Fundamentals of Astrodynamics and Applications (4th ed.), Example 3-15,
// and Vallado et al. (2006) "Revisiting Spacetrack Report #3", Appendix C.
class ValladoExample : public ::testing::Test {
protected:
    const Time        utc = Time::fromCalendar(TimeScale::UTC, 2004, 4, 6, 7, 51, 28.386009);
    const ConstantEop eop{
        EopValues{.xp   = arcsec2rad(-0.140682),
                  .yp   = arcsec2rad(0.333309),
                  .dut1 = -0.4399619,
                  .lod  = 0.0015563,
                  .dX   = arcsec2rad(-0.000205),
                  .dY   = arcsec2rad(-0.000136)}
    };
    const StateVector itrf{
        {-1033.4793830e3, 7901.2952754e3, 6380.3565958e3},
        {-3.225636520e3,  -2.872451450e3, 5.531924446e3 }
    };
};

TEST_F(ValladoExample, ItrfToGcrf) {
    const StateVector gcrf = itrsToGcrs(itrf, utc, eop);
    // Vallado's values use IAU 2000; IAU 2006 differs at the centimetre level: measured 1.0 cm and
    // 4.5e-6 m/s per component. Tolerances 2 cm and 1e-5 m/s, below the 1.05e-5 m/s that the excess
    // length of day (LOD) contributes to the velocity, so that leaving LOD out is seen.
    expectVecNear(gcrf.r, {5102.508959e3, 6123.011403e3, 6378.136925e3}, 0.02, "r GCRF");
    expectVecNear(gcrf.v, {-4.743220156e3, 0.790536497e3, 5.533755724e3}, 1e-5, "v GCRF");
    const StateVector back = gcrsToItrs(gcrf, utc, eop);
    expectVecNear(back.r, itrf.r, 1e-6, "round trip r");
    expectVecNear(back.v, itrf.v, 1e-9, "round trip v");
}

TEST_F(ValladoExample, TemeAndPef) {
    const FrameTransformer ft(utc, eop);
    const StateVector      teme = ft.transform(itrf, Frame::ITRS, Frame::TEME);
    // Measured 6.7 mm and 4.7e-6 m/s per component (TEME's GMST 1982 against Vallado's); tolerances
    // 2 cm and 1e-5 m/s. LOD contributes 1.1e-5 m/s here.
    expectVecNear(teme.r, {5094.18016210e3, 6127.64465950e3, 6380.34453270e3}, 0.02, "r TEME");
    expectVecNear(teme.v, {-4.746131487e3, 0.785818041e3, 5.531931288e3}, 1e-5, "v TEME");
    const StateVector pef = ft.transform(itrf, Frame::ITRS, Frame::PEF);
    // PEF is ITRS without polar motion: no precession-nutation model enters, and the agreement is
    // 0.1 mm and 3e-7 m/s (measured). Tolerances 0.5 mm and 1e-6 m/s.
    expectVecNear(pef.r, {-1033.47503130e3, 7901.30558560e3, 6380.34453270e3}, 5e-4, "r PEF");
    expectVecNear(pef.v, {-3.2256327470e3, -2.8724425110e3, 5.5319312880e3}, 1e-6, "v PEF");
    // TIRS and PEF differ by the TIO locator s' (about 1e-10 rad here): 3.5e-5 m (measured).
    expectVecNear(ft.transform(itrf.r, Frame::ITRS, Frame::TIRS), pef.r, 1e-4, "TIRS vs PEF");
    // J2000 (EME2000) differs from GCRS by the frame bias (~23 mas -> ~0.7 m here).
    const Vec3 j2000 = ft.transform(itrf.r, Frame::ITRS, Frame::J2000);
    const Vec3 gcrs  = ft.transform(itrf.r, Frame::ITRS, Frame::GCRS);
    EXPECT_GT(norm(j2000 - gcrs), 0.3);
    EXPECT_LT(norm(j2000 - gcrs), 1.5);
}

TEST_F(ValladoExample, OneOffConversionsBetweenTemeItrsAndGcrs) {
    // The one-off functions against Vallado's published states in the three frames (the same values
    // as above). Each tolerance covers the two published states it connects: GCRS 1.0 cm and
    // 4.5e-6 m/s, TEME 6.7 mm and 4.7e-6 m/s from ITRS (measured above), so 2 cm and 1e-5 m/s
    // throughout.
    const StateVector gcrf{
        {5102.508959e3,  6123.011403e3, 6378.136925e3},
        {-4.743220156e3, 0.790536497e3, 5.533755724e3}
    };
    const StateVector teme{
        {5094.18016210e3, 6127.64465950e3, 6380.34453270e3},
        {-4.746131487e3,  0.785818041e3,   5.531931288e3  }
    };
    const auto expectState = [](const StateVector& got, const StateVector& want, const char* what) {
        expectVecNear(got.r, want.r, 0.02, what);
        expectVecNear(got.v, want.v, 1e-5, what);
    };
    expectState(temeToItrs(teme, utc, eop), itrf, "TEME to ITRS");
    expectState(itrsToTeme(itrf, utc, eop), teme, "ITRS to TEME");
    expectState(temeToGcrs(teme, utc, eop), gcrf, "TEME to GCRS");
    expectState(gcrsToTeme(gcrf, utc, eop), teme, "GCRS to TEME");
    expectState(eciToEcef(gcrf, utc, eop), itrf, "ECI to ECEF");
    expectState(ecefToEci(itrf, utc, eop), gcrf, "ECEF to ECI");
    // The position-only forms have no Earth-rotation term and give the same positions.
    expectVecNear(eciToEcef(gcrf.r, utc, eop), itrf.r, 0.02, "ECI to ECEF, position");
    expectVecNear(ecefToEci(itrf.r, utc, eop), gcrf.r, 0.02, "ECEF to ECI, position");
}

TEST(FrameTransformer, ConsistencyAndRoundTrips) {
    const Time             t = Time::fromIso8601("2021-07-04T12:34:56Z");
    const ConstantEop      eop(EopValues{.xp = arcsec2rad(0.2), .yp = arcsec2rad(0.4), .dut1 = -0.15, .lod = 0.001});
    const FrameTransformer ft(t, eop);
    const StateVector      s{
             {7000e3, -1200e3, 300e3},
             {1.0e3,  7.2e3,   0.5e3}
    };
    const Frame all[] = {Frame::GCRS,
                         Frame::J2000,
                         Frame::MOD,
                         Frame::TOD,
                         Frame::TEME,
                         Frame::CIRS,
                         Frame::TIRS,
                         Frame::PEF,
                         Frame::ITRS};
    for (Frame a : all) {
        expectOrthonormal(ft.gcrsTo(a), 1e-14);
        for (Frame b : all) {
            const StateVector x = ft.transform(s, Frame::GCRS, a);
            const StateVector y = ft.transform(ft.transform(x, a, b), b, Frame::GCRS);
            expectVecNear(y.r, s.r, 1e-6);
            expectVecNear(y.v, s.v, 1e-9);
        }
    }
    // An ITRS-fixed point has inertial velocity omega x r.
    const StateVector fixed{
        {6378137.0, 0, 0},
        {}
    };
    const StateVector inertial = ft.transform(fixed, Frame::ITRS, Frame::GCRS);
    EXPECT_NEAR(norm(inertial.v), kEarthRotationRateEra * 6378137.0, 0.01);
    // The CIO and equinox routes agree (TIRS vs PEF) to well under a milliarcsecond.
    const Mat3 d = ft.gcrsTo(Frame::TIRS) * ft.gcrsTo(Frame::PEF).transpose();
    EXPECT_LT(std::abs(d(0, 1)), 0.1 * constants::masToRad);
    // GAST - GMST is the equation of the equinoxes (< 1.2 s of time).
    EXPECT_LT(std::abs(wrapPi(ft.gast() - ft.gmst())), deg2rad(1.2 * 15.0 / 3600.0));
    // Without Earth-orientation data (ConstantEop{}: UT1 taken as UTC, no polar motion) the transformer
    // is still built.
    EXPECT_NO_THROW(FrameTransformer(t, ConstantEop{}));
}

TEST(FrameTransformer, Model2000BWithinMilliarcsecond) {
    // IAU 2000B is documented as a 1 mas approximation of 2000A for the CIO-based (ITRS) chain over
    // 1995-2050. Measured here: 0.27 mas.
    const Time        t = Time::fromIso8601("2030-01-01T00:00:00Z");
    const ConstantEop eop;
    const Mat3        a     = gcrsToItrsMatrix(t, eop);
    const Mat3        b     = gcrsToItrsMatrix(t, eop, {.model = PrecessionNutationModel::IAU2006_2000B});
    const Mat3        d     = a * b.transpose();
    const double      angle = std::sqrt(d(0, 1) * d(0, 1) + d(0, 2) * d(0, 2) + d(1, 2) * d(1, 2));
    EXPECT_LT(angle, 1.0 * constants::masToRad);
    EXPECT_GT(angle, 0.01 * constants::masToRad); // the option does change the model
}

namespace {
std::string finalsLine(int    y,
                       int    m,
                       int    d,
                       double mjd,
                       double xp,
                       double yp,
                       double dut1,
                       double lodMs,
                       double dXmas,
                       double dYmas) {
    char buf[200];
    std::snprintf(buf,
                  sizeof(buf),
                  "%2d%02d%02d %8.2f %c %9.6f%9.6f %9.6f%9.6f  %c%10.7f%10.7f %7.4f%7.4f  %c %9.3f%9.3f %9.3f%9.3f",
                  y % 100,
                  m,
                  d,
                  mjd,
                  'I',
                  xp,
                  0.00003,
                  yp,
                  0.00003,
                  'I',
                  dut1,
                  0.00001,
                  lodMs,
                  0.007,
                  'I',
                  dXmas,
                  0.1,
                  dYmas,
                  0.1);
    return buf;
}
} // namespace

TEST(EopTable, Finals2000AParsingAndLeapSecondInterpolation) {
    // Around the 2016-12-31 leap second: UT1-UTC jumps by +1 s.
    std::string content;
    content += finalsLine(2016, 12, 29, 57751.0, 0.030, 0.280, -0.4050000, 1.0, 0.10, -0.10) + "\n";
    content += finalsLine(2016, 12, 30, 57752.0, 0.029, 0.279, -0.4060000, 1.0, 0.11, -0.11) + "\n";
    content += finalsLine(2016, 12, 31, 57753.0, 0.028, 0.278, -0.4070000, 1.0, 0.12, -0.12) + "\n";
    content += finalsLine(2017, 1, 1, 57754.0, 0.027, 0.277, 0.5920000, 1.0, 0.13, -0.13) + "\n";
    content += finalsLine(2017, 1, 2, 57755.0, 0.026, 0.276, 0.5910000, 1.0, 0.14, -0.14) + "\n";
    content += "170103 57756.00 P  \n"; // prediction rows without data are skipped
    for (auto interp : {EopInterpolation::Linear, EopInterpolation::Lagrange4}) {
        const EopTable table = EopTable::parseFinals2000A(content, interp);
        ASSERT_EQ(table.records().size(), 5u);
        const EopValues a = table.at(57752.0);
        EXPECT_NEAR(rad2arcsec(a.xp), 0.029, 1e-12);
        EXPECT_NEAR(a.dut1, -0.406, 1e-12);
        EXPECT_NEAR(a.lod, 0.001, 1e-15);
        EXPECT_NEAR(rad2mas(a.dX), 0.11, 1e-12);
        // Mid-day before the leap: smooth, no 0.5 s artifact.
        EXPECT_NEAR(table.at(57753.5).dut1, -0.4075, 1e-4);
        // After the leap.
        EXPECT_NEAR(table.at(57754.5).dut1, 0.5915, 1e-4);
        EXPECT_THROW((void)table.at(57760.0), Error);
    }
}

TEST(EopTable, C04Formats) {
    const char*    c04_14 = "# header\n"
                            "2020   1   1  58849   0.076614   0.282383  -0.1772561   0.0003791   0.000098  -0.000233\n"
                            "2020   1   2  58850   0.074996   0.282597  -0.1775922   0.0002929   0.000099  -0.000227\n";
    const EopTable a      = EopTable::parseC04(c04_14);
    EXPECT_NEAR(rad2arcsec(a.at(58849.0).xp), 0.076614, 1e-12);
    EXPECT_NEAR(a.at(58849.0).dut1, -0.1772561, 1e-12);
    EXPECT_NEAR(a.at(58849.0).lod, 0.0003791, 1e-15);
    const char* c04_20 =
        "# YR MM DD HH MJD x y UT1-UTC dX dY xrt yrt LOD\n"
        "2020  1  1  0 58849.00   0.076614   0.282383  -0.1772561   0.000098  -0.000233   0.0  0.0   0.0003791\n"
        "2020  1  2  0 58850.00   0.074996   0.282597  -0.1775922   0.000099  -0.000227   0.0  0.0   0.0002929\n";
    const EopTable b = EopTable::parseC04(c04_20);
    EXPECT_NEAR(rad2arcsec(b.at(58850.0).yp), 0.282597, 1e-12);
    EXPECT_NEAR(b.at(58850.0).lod, 0.0002929, 1e-15);
    EXPECT_NEAR(rad2arcsec(b.at(58850.0).dX), 0.000099, 1e-12);
}

TEST(EopTable, EpochThatIsNotFiniteIsRejected) {
    // Such a record used to reach std::sort, whose comparison must be a strict order (it is not for
    // NaN: undefined behaviour), and was refused only later, with a message about leap seconds.
    for (double bad : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()}) {
        EopRecord a, b, c;
        a.mjd = 59000.0;
        b.mjd = bad;
        c.mjd = 59002.0;
        try {
            const EopTable table({a, b, c});
            FAIL() << "accepted an epoch of " << bad;
        } catch (const Error& e) {
            EXPECT_NE(std::string(e.what()).find("not finite"), std::string::npos) << e.what();
        }
    }
}

TEST(EopTable, DuplicateEpochIsRejected) {
    // Two records at one epoch put a zero in the denominators of the interpolation weights: the table
    // used to be accepted and at() then returned infinity or NaN. Real IERS files have one row per
    // day, so a repeated epoch is an error in the input and is refused when the table is built.
    const auto rec = [](double mjd, double xpArcsec) {
        EopRecord r;
        r.mjd       = mjd;
        r.values.xp = arcsec2rad(xpArcsec);
        return r;
    };
    const auto message = [](const std::vector<EopRecord>& records, EopInterpolation interp) -> std::string {
        try {
            const EopTable table(records, interp);
            return "accepted; xp(59000.5) = " + std::to_string(table.at(59000.5).xp);
        } catch (const Error& e) {
            return e.what();
        }
    };
    // Records may come in any order, so the two at MJD 59001 need not be next to each other.
    const std::vector<EopRecord> repeated{rec(59001.0, 2.0),
                                          rec(59003.0, 4.0),
                                          rec(59000.0, 1.0),
                                          rec(59001.0, 2.0),
                                          rec(59002.0, 3.0)};
    for (auto interp : {EopInterpolation::Linear, EopInterpolation::Lagrange4}) {
        const std::string what = message(repeated, interp);
        EXPECT_NE(what.find("duplicate epoch"), std::string::npos) << what;
        EXPECT_NE(what.find("59001"), std::string::npos) << what; // names the epoch
    }
    // The same through the file readers.
    std::string finals;
    finals += finalsLine(2020, 5, 31, 59000.0, 0.10, 0.40, -0.2400000, 0.5, 0.10, -0.10) + "\n";
    finals += finalsLine(2020, 6, 1, 59001.0, 0.11, 0.41, -0.2405000, 0.5, 0.11, -0.11) + "\n";
    finals += finalsLine(2020, 6, 1, 59001.0, 0.11, 0.41, -0.2405000, 0.5, 0.11, -0.11) + "\n";
    finals += finalsLine(2020, 6, 2, 59002.0, 0.12, 0.42, -0.2410000, 0.5, 0.12, -0.12) + "\n";
    EXPECT_THROW((void)EopTable::parseFinals2000A(finals), Error);
    const char* c04 = "2020   5  31  59000   0.100000   0.400000  -0.2400000   0.0005000   0.000100  -0.000100\n"
                      "2020   6   1  59001   0.110000   0.410000  -0.2405000   0.0005000   0.000110  -0.000110\n"
                      "2020   6   1  59001   0.110000   0.410000  -0.2405000   0.0005000   0.000110  -0.000110\n"
                      "2020   6   2  59002   0.120000   0.420000  -0.2410000   0.0005000   0.000120  -0.000120\n";
    EXPECT_THROW((void)EopTable::parseC04(c04), Error);

    // Without the repeat, the same records in the same (unsorted) order are accepted and sorted. The
    // values lie on a straight line (xp = MJD - 58999 arcsec), which both interpolations reproduce.
    std::vector<EopRecord> distinct = repeated;
    distinct.erase(distinct.begin());
    for (auto interp : {EopInterpolation::Linear, EopInterpolation::Lagrange4}) {
        const EopTable table(distinct, interp);
        ASSERT_EQ(table.records().size(), 4u);
        EXPECT_EQ(table.firstMjd(), 59000.0);
        EXPECT_EQ(table.lastMjd(), 59003.0);
        for (double mjd : {59000.0, 59000.5, 59001.0, 59001.75, 59002.25, 59003.0}) {
            EXPECT_NEAR(rad2arcsec(table.at(mjd).xp), mjd - 58999.0, 1e-12) << mjd; // rounding only
        }
    }
}

namespace {

/// Daily EOP records from MJD 60000 to 60007 (2023, no leap second) whose values are polynomials of
/// u = MJD - 60003.4: polar motion x, UT1-UTC, LOD, dX and dY cubic, polar motion y quartic.
struct PolynomialEop {
    static constexpr double c4 = 1e-3; ///< coefficient of u^4 in y, arcsec

    static double u(double mjd) { return mjd - 60003.4; }
    static double xp(double mjd) {
        return arcsec2rad(0.1 + 0.01 * u(mjd) - 0.002 * std::pow(u(mjd), 2) + 3e-4 * std::pow(u(mjd), 3));
    }
    static double yp(double mjd) { return arcsec2rad(0.3 + 0.02 * u(mjd) + c4 * std::pow(u(mjd), 4)); }
    static double dut1(double mjd) { return -0.1 + 1e-3 * u(mjd) - 2e-4 * std::pow(u(mjd), 3); }
    static double lod(double mjd) { return 1e-3 + 1e-4 * std::pow(u(mjd), 2) - 1e-5 * std::pow(u(mjd), 3); }
    static double dX(double mjd) { return mas2rad(0.2 - 0.03 * u(mjd) + 0.004 * std::pow(u(mjd), 3)); }
    static double dY(double mjd) { return mas2rad(-0.1 + 0.01 * std::pow(u(mjd), 2) - 0.002 * std::pow(u(mjd), 3)); }

    static std::vector<EopRecord> records() {
        std::vector<EopRecord> out;
        for (int d = 0; d <= 7; ++d) {
            const double m = 60000.0 + d;
            out.push_back({
                m,
                {xp(m), yp(m), dut1(m), lod(m), dX(m), dY(m)}
            });
        }
        return out;
    }
};

} // namespace

TEST(EopTable, Lagrange4ReproducesACubicExactly) {
    // Four-point Lagrange interpolation is exact for a polynomial of degree three: at nodes, in the
    // first and last intervals (where the four points cannot straddle the epoch) and inside.
    const EopTable table(PolynomialEop::records(), EopInterpolation::Lagrange4);
    for (double mjd : {60000.0, 60000.3, 60001.5, 60003.7, 60005.0, 60006.8, 60007.0}) {
        SCOPED_TRACE(testing::Message() << "MJD " << mjd);
        const EopValues v = table.at(mjd);
        // Rounding only: a few units in the last place of each value, summed with weights below 2 (and
        // UT1-UTC passes through UT1-TAI, 37 s). Measured 1.1e-22 rad, 5.3e-15 s, 4.3e-19 s, 2.1e-25 rad.
        EXPECT_NEAR(v.xp, PolynomialEop::xp(mjd), 1e-21);
        EXPECT_NEAR(v.dut1, PolynomialEop::dut1(mjd), 5e-14);
        EXPECT_NEAR(v.lod, PolynomialEop::lod(mjd), 5e-18);
        EXPECT_NEAR(v.dX, PolynomialEop::dX(mjd), 2e-24);
        EXPECT_NEAR(v.dY, PolynomialEop::dY(mjd), 2e-24);
    }
}

TEST(EopTable, Lagrange4ErrorOnAQuarticIsTheNodePolynomial) {
    // For f(u) = ... + c4 u^4 the interpolation error of a cubic through four nodes x_k is exactly
    // f''''/4! prod(x - x_k) = c4 prod(x - x_k). The nodes are the two records on each side of the
    // epoch, or the first (last) four records in the first (last) interval.
    const EopTable table(PolynomialEop::records(), EopInterpolation::Lagrange4);
    for (double mjd : {60000.3, 60001.5, 60003.7, 60004.2, 60006.8}) {
        SCOPED_TRACE(testing::Message() << "MJD " << mjd);
        const double lower          = std::clamp(std::floor(mjd) - 1.0, 60000.0, 60004.0); // first of the four nodes
        double       nodePolynomial = 1.0;
        for (int k = 0; k < 4; ++k) {
            nodePolynomial *= mjd - (lower + k);
        }
        const double expected = PolynomialEop::yp(mjd) - arcsec2rad(PolynomialEop::c4 * nodePolynomial);
        // The error itself is 1.7e-9 to 4.7e-9 rad (0.35 to 0.96 mas) at these epochs; rounding is
        // 4.2e-22 rad (measured).
        EXPECT_NEAR(table.at(mjd).yp, expected, 4e-21);
        EXPECT_GT(std::abs(PolynomialEop::yp(mjd) - expected), 1e-9); // the error term is not negligible
    }
}

TEST(EopTable, FilesAreReadLikeTheirContents) {
    // loadFinals2000A and loadC04 read a file and parse it as parseFinals2000A and parseC04 parse the
    // same text; a file that cannot be opened throws.
    std::string finals;
    finals                += finalsLine(2020, 5, 31, 59000.0, 0.10, 0.40, -0.2400000, 0.5, 0.10, -0.10) + "\n";
    finals                += finalsLine(2020, 6, 1, 59001.0, 0.11, 0.41, -0.2405000, 0.6, 0.11, -0.11) + "\n";
    finals                += finalsLine(2020, 6, 2, 59002.0, 0.12, 0.42, -0.2410000, 0.7, 0.12, -0.12) + "\n";
    const std::string c04  = "2020   5  31  59000   0.100000   0.400000  -0.2400000   0.0005000   0.000100  -0.000100\n"
                             "2020   6   1  59001   0.110000   0.410000  -0.2405000   0.0006000   0.000110  -0.000110\n";
    const auto        write = [](const std::string& name, const std::string& text) {
        const std::string path = testing::TempDir() + name;
        std::FILE*        f    = std::fopen(path.c_str(), "w");
        EXPECT_NE(f, nullptr) << path;
        if (f != nullptr) {
            std::fputs(text.c_str(), f);
            std::fclose(f);
        }
        return path;
    };
    const EopTable fromFile = EopTable::loadFinals2000A(write("astro_finals_test.txt", finals));
    const EopTable fromText = EopTable::parseFinals2000A(finals);
    ASSERT_EQ(fromFile.records().size(), 3u);
    EXPECT_EQ(fromFile.firstMjd(), 59000.0);
    EXPECT_EQ(fromFile.lastMjd(), 59002.0);
    for (double mjd : {59000.0, 59001.3, 59002.0}) {
        EXPECT_EQ(fromFile.at(mjd).xp, fromText.at(mjd).xp) << mjd;
        EXPECT_EQ(fromFile.at(mjd).dut1, fromText.at(mjd).dut1) << mjd;
    }
    EXPECT_NEAR(rad2arcsec(fromFile.at(59001.0).xp), 0.11, 1e-12); // the value in the file
    const EopTable c04FromFile = EopTable::loadC04(write("astro_c04_test.txt", c04));
    ASSERT_EQ(c04FromFile.records().size(), 2u);
    EXPECT_NEAR(c04FromFile.at(59001.0).lod, 0.0006, 1e-15);
    EXPECT_EQ(c04FromFile.at(59000.5).yp, EopTable::parseC04(c04).at(59000.5).yp);
    EXPECT_THROW((void)EopTable::loadFinals2000A(testing::TempDir() + "no_such_finals_file.txt"), Error);
    EXPECT_THROW((void)EopTable::loadC04(testing::TempDir() + "no_such_c04_file.txt"), Error);
}

TEST(EopTable, TooFewRecordsOrAnEpochOutsideTheTableThrows) {
    // Documented: fewer than two records (or no usable rows in a file), an epoch before 1960, and at()
    // outside [firstMjd, lastMjd] or at NaN.
    EopRecord a, b;
    a.mjd = 59000.0;
    b.mjd = 59001.0;
    EXPECT_THROW(EopTable({}), Error);
    EXPECT_THROW(EopTable({a}), Error);
    EopRecord old;
    old.mjd = 36000.0; // 1957
    EXPECT_THROW(EopTable({old, a}), Error);
    EXPECT_THROW((void)EopTable::parseFinals2000A("# nothing but a comment\n"), Error);
    EXPECT_THROW((void)EopTable::parseC04("# nothing but a comment\n"), Error);
    const EopTable table({a, b});
    EXPECT_NO_THROW((void)table.at(59000.0));
    EXPECT_NO_THROW((void)table.at(59001.0));
    EXPECT_THROW((void)table.at(58999.999), Error);
    EXPECT_THROW((void)table.at(59001.001), Error);
    EXPECT_THROW((void)table.at(std::numeric_limits<double>::quiet_NaN()), Error);
}

TEST(EopTable, LinearInterpolationIsTheChord) {
    const EopTable table(PolynomialEop::records(), EopInterpolation::Linear);
    for (double mjd : {60000.3, 60003.7, 60006.8}) {
        const double x0 = std::floor(mjd), t = mjd - x0;
        EXPECT_NEAR(table.at(mjd).yp, (1.0 - t) * PolynomialEop::yp(x0) + t * PolynomialEop::yp(x0 + 1.0), 4e-21)
            << mjd;
        EXPECT_NEAR(table.at(mjd).dut1, (1.0 - t) * PolynomialEop::dut1(x0) + t * PolynomialEop::dut1(x0 + 1.0), 5e-14)
            << mjd;
    }
}

TEST(FrameTransformer, InvalidFrameValueGivesNaN) {
    // A Frame that is not one of the enumerators (a bad cast, an uninitialised value) used to index
    // past the end of the transformer's table of matrices. gcrsTo() and matrix() are noexcept, so they
    // answer with NaN, which cannot be taken for a rotation; the functions that may throw, throw.
    const Time             t = Time::fromIso8601("2021-07-04T12:34:56Z");
    const FrameTransformer ft(t, ConstantEop{});
    const auto             allNaN = [](const Mat3& m) {
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                if (!std::isnan(m(i, j))) {
                    return false;
                }
            }
        }
        return true;
    };
    const auto        allNaNVec = [](const Vec3& v) { return std::isnan(v.x) && std::isnan(v.y) && std::isnan(v.z); };
    const StateVector s{
        {7000e3, -1200e3, 300e3},
        {1.0e3,  7.2e3,   0.5e3}
    };
    for (int raw : {9, 99, -1}) { // 9 is one past the last enumerator (ITRS = 8)
        const Frame bad = static_cast<Frame>(raw);
        EXPECT_TRUE(allNaN(ft.gcrsTo(bad))) << raw;
        EXPECT_TRUE(allNaN(ft.matrix(bad, Frame::GCRS))) << raw;
        EXPECT_TRUE(allNaN(ft.matrix(Frame::ITRS, bad))) << raw;
        EXPECT_TRUE(allNaNVec(ft.transform(s.r, bad, Frame::ITRS))) << raw;
        EXPECT_TRUE(allNaNVec(ft.transform(s.r, Frame::GCRS, bad))) << raw;
        for (const StateVector& out : {ft.transform(s, bad, Frame::ITRS), ft.transform(s, Frame::GCRS, bad)}) {
            EXPECT_TRUE(allNaNVec(out.r)) << raw;
            EXPECT_TRUE(allNaNVec(out.v)) << raw;
        }
        EXPECT_THROW((void)gcrsToCelestialMatrix(bad, t), Error) << raw;
        EXPECT_THROW((void)poleOfDateMatrix(bad, t), Error) << raw;
    }
}
