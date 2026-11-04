/**********************************************************************
 * impact_test.cpp
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
#include <memory>
#include <vector>

#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/core/rootfind.h"
#include "astro/geodesy/local.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/impact.h"
#include "test_util.h"

using namespace astro;

namespace {

const Time kEpoch = Time::fromIso8601("2026-03-01T12:00:00Z");

ConstantEop someEop() {
    EopValues v;
    v.xp   = 0.10 * constants::arcsecToRad;
    v.yp   = 0.35 * constants::arcsecToRad;
    v.dut1 = 0.3;
    v.lod  = 1.0e-3;
    return ConstantEop(v);
}

/// GCRS state of a launch from `site` with ground-relative speed V at flight-path angle gamma and
/// azimuth az.
StateVector launch(const Geodetic& site, double V, double gamma, double az, const EopProvider& eop) {
    const Mat3 enu = ecefToEnuMatrix(site.lat, site.lon);
    const Vec3 v{V * std::cos(gamma) * std::sin(az), V * std::cos(gamma) * std::cos(az), V * std::sin(gamma)};
    return itrsToGcrs(StateVector{geodeticToEcef(site), enu.transpose() * v}, kEpoch, eop);
}

const Geodetic kSite{deg2rad(30.0), deg2rad(-100.0), 0.0};

/// Brute-force reference: tight tolerance, the full IAU 2006/2000A transformation at every
/// evaluation, Brent's method on the propagated height.
struct Reference {
    double dt;
    Vec3   itrs;
};
Reference referenceImpact(const StateVector& s0, double guess, const EopProvider& eop) {
    const ImpactOptions o;
    AdaptiveOptions     tight;
    tight.relTol = 1e-13;
    tight.absTol = 1e-9;
    const NumericalPropagator prop(
        kEpoch,
        s0,
        {std::make_shared<GravityField>(o.J, o.mu, o.re, poleOfDateMatrix(Frame::GCRS, kEpoch))},
        tight);
    const auto height = [&](double dt) { return ecefToGeodetic(gcrsToItrs(prop.propagate(dt).r, kEpoch + dt, eop)).h; };
    const double dt   = brentRoot(height, guess - 2.0, guess + 2.0, 1e-9);
    return {dt, gcrsToItrs(prop.propagate(dt).r, kEpoch + dt, eop)};
}

} // namespace

TEST(ErrorEllipse, AxesAndRotationFromCovariance) {
    // C = a^2 u u^T + b^2 w w^T with u the major-axis direction at azimuth alpha (east = sin, north = cos).
    for (double alphaDeg : {0.0, 30.0, 90.0, 120.0, 179.0}) {
        const double       a = 300.0, b = 100.0, alpha = deg2rad(alphaDeg);
        const double       ue = std::sin(alpha), un = std::cos(alpha), we = std::cos(alpha), wn = -std::sin(alpha);
        const ErrorEllipse e = errorEllipse(a * a * ue * ue + b * b * we * we,
                                            a * a * un * un + b * b * wn * wn,
                                            a * a * ue * un + b * b * we * wn);
        EXPECT_NEAR(e.semiMajor, a, 1e-9) << alphaDeg;
        EXPECT_NEAR(e.semiMinor, b, 1e-9) << alphaDeg;
        EXPECT_NEAR(e.rotation, alpha, 1e-12) << alphaDeg;
    }
    EXPECT_EQ(errorEllipse(4.0, 4.0, 0.0).rotation, 0.0); // a circle
    EXPECT_NEAR(errorEllipse(4.0, 4.0, 0.0).semiMinor, 2.0, 1e-15);
    EXPECT_NEAR(ellipseScale(0.5), 1.1774100225154747, 1e-12);
    EXPECT_NEAR(ellipseScale(0.95), 2.4477468306808161, 1e-12);
    const ErrorEllipse e95 = ErrorEllipse{10.0, 5.0, 0.4}.scaled(0.95);
    EXPECT_NEAR(e95.semiMajor, 24.477468306808161, 1e-10);
    EXPECT_NEAR(e95.semiMinor, 12.238734153404080, 1e-10);
    EXPECT_EQ(e95.rotation, 0.4);
    EXPECT_THROW((void)ellipseScale(0.0), Error);
    EXPECT_THROW((void)ellipseScale(1.0), Error);
}

TEST(Impact, MatchesBruteForceReference) {
    // Launches from the ellipsoid surface: 400 s to 40 min of flight. The default options (IAU
    // 2000B, interpolated Earth rotation, relTol 1e-9) agree with the brute-force reference to
    // about a centimetre.
    const ConstantEop eop = someEop();
    for (const double V : {3000.0, 5000.0, 7000.0}) {
        const StateVector s0  = launch(kSite, V, deg2rad(35.0), deg2rad(45.0), eop);
        const auto        imp = predictImpact(kEpoch, s0, Frame::GCRS, eop);
        ASSERT_TRUE(imp.has_value()) << V;
        const Reference ref = referenceImpact(s0, imp->flightTime, eop);
        EXPECT_NEAR(imp->flightTime, ref.dt, 1e-5) << V;
        EXPECT_LT(norm(imp->itrs.r - ref.itrs), 0.02) << V;
        EXPECT_NEAR(imp->point.h, 0.0, 1e-3) << V;
        EXPECT_NEAR(imp->time - kEpoch, imp->flightTime, 1e-6);
        // Ballistic symmetry over a nearly spherical Earth: impact angle and speed close to launch.
        EXPECT_NEAR(rad2deg(imp->angle), 35.0, 0.5) << V;
        EXPECT_NEAR(imp->speed / V, 1.0, 0.01) << V;
        // The ECI and ECEF results describe the same point.
        test::expectVecNear(gcrsToItrs(imp->state.r, imp->time, eop), imp->itrs.r, 0.05);
    }
}

TEST(Impact, FramesSurfaceHeightAndSharedRotations) {
    const ConstantEop eop  = someEop();
    const StateVector s0   = launch(kSite, 5000.0, deg2rad(40.0), deg2rad(300.0), eop);
    const auto        gcrs = predictImpact(kEpoch, s0, Frame::GCRS, eop);
    ASSERT_TRUE(gcrs.has_value());
    EXPECT_GT(gcrs->flightTime, 600.0); // the launch point (h = 0 at the epoch) is not the impact

    // The same physical state given in other celestial frames gives the same impact, and the
    // impact state comes back in the input frame of the impact time.
    for (Frame f : {Frame::J2000, Frame::TOD, Frame::TEME}) {
        const Mat3 M = gcrsToCelestialMatrix(f, kEpoch);
        const auto i = predictImpact(kEpoch, StateVector{M * s0.r, M * s0.v}, f, eop);
        ASSERT_TRUE(i.has_value());
        EXPECT_NEAR(i->flightTime, gcrs->flightTime, 1e-7) << f;
        test::expectVecNear(i->itrs.r, gcrs->itrs.r, 1e-4);
        test::expectVecNear(i->state.r, gcrsToCelestialMatrix(f, i->time) * gcrs->state.r, 1e-4);
    }

    // A raised surface is reached earlier, at that height.
    ImpactOptions high;
    high.surfaceHeight = 1500.0;
    const auto hi      = predictImpact(kEpoch, s0, Frame::GCRS, eop, high);
    ASSERT_TRUE(hi.has_value());
    EXPECT_LT(hi->flightTime, gcrs->flightTime);
    EXPECT_NEAR(hi->point.h, 1500.0, 1e-3);

    // A shared rotation series gives the same answer; one that does not cover the flight throws.
    const ItrsRotationSeries series(kEpoch - 600.0, 3 * 3600.0, eop, ImpactOptions{}.frameOptions);
    ImpactOptions            shared;
    shared.rotations = &series;
    const auto sh    = predictImpact(kEpoch, s0, Frame::GCRS, eop, shared);
    ASSERT_TRUE(sh.has_value());
    EXPECT_NEAR(sh->flightTime, gcrs->flightTime, 1e-6);
    test::expectVecNear(sh->itrs.r, gcrs->itrs.r, 1e-3);
    const ItrsRotationSeries shortSeries(kEpoch, 600.0, eop);
    shared.rotations = &shortSeries;
    EXPECT_THROW((void)predictImpact(kEpoch, s0, Frame::GCRS, eop, shared), Error);
}

TEST(Impact, NoImpactAndInvalidInput) {
    const ConstantEop eop;
    // A 400 km circular orbit never comes down.
    const StateVector orbit = elementsToState({6778e3, 0.0, deg2rad(51.6), 0.3, 0.5, 0.7});
    EXPECT_FALSE(predictImpact(kEpoch, orbit, Frame::GCRS, eop).has_value());
    EXPECT_THROW((void)predictImpact(kEpoch, orbit, Frame::ITRS, eop), Error);
    ImpactOptions bad;
    bad.maxFlightTime = -1.0;
    EXPECT_THROW((void)predictImpact(kEpoch, orbit, Frame::GCRS, eop, bad), Error);
    bad               = {};
    bad.surfaceHeight = std::nan("");
    EXPECT_THROW((void)predictImpact(kEpoch, orbit, Frame::GCRS, eop, bad), Error);
}

TEST(Impact, UncertaintyMatchesFiniteDifferences) {
    // The impact ellipse is the linear map of the initial covariance: check the whole chain
    // (state transition matrix, time-of-impact correction, frames) against central differences of
    // the impact point, with a correlated covariance given in TOD.
    const ConstantEop eop = someEop();
    const Mat3        M   = gcrsToCelestialMatrix(Frame::TOD, kEpoch);
    const StateVector g0  = launch(kSite, 6000.0, deg2rad(30.0), deg2rad(70.0), eop);
    const StateVector s0{M * g0.r, M * g0.v};
    Mat6              L{}; // P0 = L L^T, correlated
    const double      sr[3] = {80.0, 120.0, 60.0}, sv[3] = {0.05, 0.08, 0.03};
    for (std::size_t i = 0; i < 3; ++i) {
        L[i][i]         = sr[i];
        L[i + 3][i + 3] = sv[i];
        L[i + 3][i]     = 0.3 * sv[i]; // velocity correlated with position
    }
    L[1][0]       = 40.0;
    const Mat6 P0 = multiply(L, transpose(L));

    const auto nominal = predictImpact(kEpoch, s0, P0, Frame::TOD, eop);
    ASSERT_TRUE(nominal.has_value());
    ASSERT_TRUE(nominal->uncertainty.has_value());
    const ImpactUncertainty& u = *nominal->uncertainty;

    // Jacobian of (east, north, time) with respect to the initial state.
    const Mat3 enu = ecefToEnuMatrix(nominal->point.lat, nominal->point.lon);
    double     J[3][6];
    for (std::size_t k = 0; k < 6; ++k) {
        const double h  = k < 3 ? 1.0 : 1e-3;
        StateVector  sp = s0, sm = s0;
        (k < 3 ? sp.r[k] : sp.v[k - 3]) += h;
        (k < 3 ? sm.r[k] : sm.v[k - 3]) -= h;
        const auto ip                    = predictImpact(kEpoch, sp, Frame::TOD, eop);
        const auto im                    = predictImpact(kEpoch, sm, Frame::TOD, eop);
        ASSERT_TRUE(ip && im);
        const Vec3 d = enu * (ip->itrs.r - im->itrs.r) / (2.0 * h);
        J[0][k]      = d.x;
        J[1][k]      = d.y;
        J[2][k]      = (ip->flightTime - im->flightTime) / (2.0 * h);
    }
    double C[3][3] = {};
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            for (std::size_t a = 0; a < 6; ++a) {
                for (std::size_t b = 0; b < 6; ++b) {
                    C[i][j] += J[i][a] * P0[a][b] * J[j][b];
                }
            }
        }
    }
    const double scale = u.ellipse.semiMajor * u.ellipse.semiMajor;
    EXPECT_NEAR(u.eastVariance, C[0][0], 1e-4 * scale);
    EXPECT_NEAR(u.northVariance, C[1][1], 1e-4 * scale);
    EXPECT_NEAR(u.eastNorthCovariance, C[0][1], 1e-4 * scale);
    EXPECT_NEAR(u.timeSigma, std::sqrt(C[2][2]), 1e-4 * u.timeSigma);
    EXPECT_GT(u.ellipse.semiMajor, u.ellipse.semiMinor);

    // The covariance at impact is reported in the input frame, and process noise only adds to it.
    EXPECT_GT(u.covariance[0][0], 0.0);
    ImpactOptions noisy;
    noisy.accelerationNoisePsd = 1e-6;
    const auto n               = predictImpact(kEpoch, s0, P0, Frame::TOD, eop, noisy);
    ASSERT_TRUE(n && n->uncertainty);
    EXPECT_GT(n->uncertainty->ellipse.semiMajor, u.ellipse.semiMajor);
    EXPECT_GE(n->uncertainty->ellipse.semiMinor, u.ellipse.semiMinor);
}
