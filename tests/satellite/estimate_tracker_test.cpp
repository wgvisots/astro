/**********************************************************************
 * estimate_tracker_test.cpp
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

#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <random>
#include <vector>

#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/core/mat6.h"
#include "astro/frames/eop.h"
#include "astro/frames/transform.h"
#include "astro/orbit/elements.h"
#include "astro/satellite/estimate_tracker.h"
#include "astro/satellite/orbital_frames.h"
#include "test_util.h"

using namespace astro;

namespace {

constexpr double km = 1000.0;

const Time kEpoch = Time::fromIso8601("2026-03-01T00:00:00Z");

StateVector leo() {
    return elementsToState({6778 * km, 0.0005, degToRad(51.6), 0.3, 0.5, 0.7});
}

Mat6 diagonal(double positionVar, double velocityVar) {
    Mat6 P{};
    for (std::size_t i = 0; i < 3; ++i) {
        P(i, i)         = positionVar;
        P(i + 3, i + 3) = velocityVar;
    }
    return P;
}

/// "Truth": the tracker's own force model with a tight tolerance.
StateVector truthAt(const StateVector& s0, double dt) {
    TrackerOptions tight;
    tight.integrator.relTol = 1e-13;
    tight.integrator.absTol = 1e-10;
    return trackerPropagator({kEpoch, s0, Mat6{}}, tight).propagate(dt);
}

} // namespace

TEST(EstimateTracker, ChiSquareSurvival) {
    EXPECT_NEAR(chiSquareSurvival6(16.8119), 0.01, 1e-6);
    EXPECT_NEAR(chiSquareSurvival6(22.4577), 0.001, 1e-7);
    EXPECT_NEAR(chiSquareSurvival3(11.3449), 0.01, 1e-6);
    EXPECT_NEAR(chiSquareSurvival3(7.8147), 0.05, 1e-6);
    EXPECT_EQ(chiSquareSurvival6(0.0), 1.0);
    EXPECT_EQ(chiSquareSurvival3(-1.0), 1.0);
}

TEST(EstimateTracker, ChiSquareSurvivalOfAHugeStatisticIsZero) {
    // P(X >= x) is erfc(sqrt(x/2)) + sqrt(2x/pi) exp(-x/2) for 3 degrees of freedom and
    // exp(-x/2) (1 + x/2 + x^2/8) for 6. Evaluated in extended precision, both are still representable
    // doubles at x = 1400 ...
    for (long double x : {100.0L, 700.0L, 1400.0L}) {
        const long double q3 = std::erfc(std::sqrt(0.5L * x))
                             + std::sqrt(2.0L * x / std::numbers::pi_v<long double>) * std::exp(-0.5L * x);
        const long double q6 = std::exp(-0.5L * x) * (1.0L + 0.5L * x + 0.125L * x * x);
        EXPECT_NEAR(chiSquareSurvival3(static_cast<double>(x)) / static_cast<double>(q3), 1.0, 1e-12);
        EXPECT_NEAR(chiSquareSurvival6(static_cast<double>(x)) / static_cast<double>(q6), 1.0, 1e-12);
    }
    // ... and from x = 1e4 on they are below 1e-2000, far under the smallest double: the value is 0,
    // all the way to infinity.
    const double inf = std::numeric_limits<double>::infinity();
    for (double x : {1e4, 1e10, 1e154, 1e155, 1e300, std::numeric_limits<double>::max(), inf}) {
        EXPECT_EQ(chiSquareSurvival3(x), 0.0) << x;
        EXPECT_EQ(chiSquareSurvival6(x), 0.0) << x;
    }
}

TEST(EstimateTracker, OrbitTurnedOverIsInfinitelyFarWithZeroPValue) {
    // The same position with the velocity reversed: a retrograde orbit, 128.4 deg from the predicted
    // plane. The tracker calls that infinitely far in elements; the chance of it is zero, not NaN.
    const StateVector   s0 = leo();
    const StateEstimate prior{kEpoch, s0, diagonal(100.0, 1e-4)};
    const StateVector   truth = truthAt(s0, 60.0);
    const StateEstimate next{
        kEpoch + 60.0,
        {truth.r, -truth.v},
        diagonal(100.0, 1e-4)
    };
    const auto c = checkConsistency(prior, next);
    EXPECT_TRUE(c.valid);
    EXPECT_EQ(c.elementsMahalanobis2, std::numeric_limits<double>::infinity());
    EXPECT_EQ(c.elementsPValue, 0.0);
    EXPECT_FALSE(c.consistent);
}

TEST(EstimateTracker, NewEstimateOnAnOpenOrbitIsInfinitelyFarInElements) {
    // The same position at one and a half times the speed: a hyperbola, which has no equinoctial
    // elements. Against an elliptic prediction that is infinitely far, with p-value 0; the other two
    // tests still give their (finite, enormous) distances.
    const StateVector   s0 = leo();
    const StateEstimate prior{kEpoch, s0, diagonal(100.0, 1e-4)};
    const StateVector   truth = truthAt(s0, 60.0);
    const StateEstimate next{
        kEpoch + 60.0,
        {truth.r, 1.5 * truth.v},
        diagonal(100.0, 1e-4)
    };
    const auto c = checkConsistency(prior, next);
    EXPECT_TRUE(c.valid);
    EXPECT_EQ(c.elementsMahalanobis2, std::numeric_limits<double>::infinity());
    EXPECT_EQ(c.elementsPValue, 0.0);
    EXPECT_FALSE(c.consistent);
    EXPECT_TRUE(std::isfinite(c.mahalanobis2));
    EXPECT_GT(c.mahalanobis2, 1e6);
    // A prediction that is itself on an open orbit has no elements statistic: it is NaN, the other
    // two statistics are unaffected, and nothing throws. With the (default) elements test the
    // delivery counts as not consistent; with the position test a delivery that follows the
    // hyperbola is accepted.
    const StateVector   open{s0.r, 1.5 * s0.v};
    const StateEstimate openPrior{kEpoch, open, diagonal(100.0, 1e-4)};
    const StateEstimate onTheHyperbola{kEpoch + 60.0, truthAt(open, 60.0), diagonal(100.0, 1e-4)};
    const auto          e = checkConsistency(openPrior, onTheHyperbola);
    EXPECT_TRUE(e.valid);
    EXPECT_TRUE(std::isnan(e.elementsMahalanobis2));
    EXPECT_TRUE(std::isnan(e.elementsPValue));
    EXPECT_TRUE(std::isnan(e.semiMajorAxisResidual));
    EXPECT_LT(e.mahalanobis2, 1.0); // the new estimate is the prediction itself
    EXPECT_GT(e.positionPValue, 0.9);
    EXPECT_FALSE(e.consistent);
    TrackerOptions byPosition;
    byPosition.test = ConsistencyTest::Position;
    EXPECT_TRUE(checkConsistency(openPrior, onTheHyperbola, byPosition).consistent);
    // A tracker on such an object keeps taking deliveries.
    EstimateTracker tracker(byPosition);
    (void)tracker.update(openPrior);
    const auto update = tracker.update(onTheHyperbola);
    ASSERT_TRUE(update.has_value());
    EXPECT_TRUE(update->consistent);
    EXPECT_TRUE(update->replaced);
}

TEST(EstimateTracker, Mahalanobis) {
    const Mat6                  D = diagonal(4.0, 1e-4);
    const std::array<double, 6> d{2.0, 0.0, -4.0, 0.01, 0.0, 0.0};
    EXPECT_NEAR(*mahalanobis2(d, D), 1.0 + 4.0 + 1.0, 1e-12);
    EXPECT_NEAR(*mahalanobis2(d, D, 3), 5.0, 1e-12);
    // Correlated S = A A^T and d = S x give d^T S^-1 d = x^T S x.
    Mat6 A{};
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            A(i, j) = (i == j ? 2.0 : 0.3) * (i < 3 ? 10.0 : 0.01);
        }
    }
    const Mat6                  S = A * A.transpose();
    const std::array<double, 6> x{0.1, -0.2, 0.05, 30.0, -10.0, 4.0};
    std::array<double, 6>       Sx{};
    double                      xSx = 0.0;
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            Sx[i] += S(i, j) * x[j];
        }
    }
    for (std::size_t i = 0; i < 6; ++i) {
        xSx += x[i] * Sx[i];
    }
    EXPECT_NEAR(*mahalanobis2(Sx, S) / xSx, 1.0, 1e-10);
    EXPECT_FALSE(mahalanobis2(d, Mat6{}).has_value());
    Mat6 bad  = D;
    bad(0, 1) = bad(1, 0) = 5.0; // |rho| > 1
    EXPECT_FALSE(mahalanobis2(d, bad).has_value());
}

TEST(EstimateTracker, ConsistentDeliveryPasses) {
    const StateVector   s0 = leo();
    const StateEstimate prior{kEpoch, s0, diagonal(100.0, 1e-4)};
    const StateEstimate next{kEpoch + 1800.0, truthAt(s0, 1800.0), diagonal(25.0, 1e-5)};
    const auto          c = checkConsistency(prior, next);
    EXPECT_TRUE(c.valid);
    EXPECT_TRUE(c.consistent);
    EXPECT_TRUE(c.replaced);
    EXPECT_DOUBLE_EQ(c.dt, 1800.0);
    EXPECT_LT(norm(c.positionResidualRic), 1.0); // only the integrator error at relTol 1e-9
    EXPECT_LT(c.mahalanobis2, 1e-2);
    EXPECT_LT(c.elementsMahalanobis2, 1e-2);
    EXPECT_LT(std::abs(c.semiMajorAxisResidual), 0.1);
    // The in-track uncertainty grows fastest.
    EXPECT_GT(c.positionSigmaRic.y, c.positionSigmaRic.x);
    EXPECT_GT(c.positionSigmaRic.y, c.positionSigmaRic.z);
}

TEST(EstimateTracker, DistanceIsChiSquareForSampledPriors) {
    // Priors drawn from N(truth, P0): the predicted-minus-truth distance must follow chi-square with
    // 6 degrees of freedom (mean 6, variance 12) if the covariance propagation is right.
    const StateVector                s0 = leo();
    const Mat6                       P0 = diagonal(400.0, 4e-4); // 20 m, 2 cm/s
    const double                     dt = 2700.0;
    const StateEstimate              truth{kEpoch + dt, truthAt(s0, dt), diagonal(1e-6, 1e-12)};
    std::mt19937_64                  rng(12345);
    std::normal_distribution<double> n01(0.0, 1.0);
    const int                        samples = 400;
    double                           sum     = 0.0;
    int                              tail    = 0;
    for (int k = 0; k < samples; ++k) {
        StateVector s  = s0;
        s.r           += Vec3{n01(rng), n01(rng), n01(rng)} * 20.0;
        s.v           += Vec3{n01(rng), n01(rng), n01(rng)} * 0.02;
        const auto c   = checkConsistency({kEpoch, s, P0}, truth);
        ASSERT_TRUE(c.valid);
        sum  += c.mahalanobis2;
        tail += c.pValue < 0.05 ? 1 : 0;
    }
    EXPECT_NEAR(sum / samples, 6.0, 0.6); // standard error of the mean ~0.17
    EXPECT_NEAR(static_cast<double>(tail) / samples, 0.05, 0.035);
}

TEST(EstimateTracker, PositionTestHoldsOverLongGaps) {
    // 4 h after a 100 m / 0.1 m/s prior the Cartesian linear covariance no longer describes the full
    // 6-D error, but its position marginal still does: position d^2 stays chi-square with 3 dof.
    const StateVector                s0 = leo();
    const Mat6                       P0 = diagonal(1e4, 1e-2);
    const double                     dt = 14400.0;
    const StateEstimate              truth{kEpoch + dt, truthAt(s0, dt), diagonal(1e-6, 1e-12)};
    std::mt19937_64                  rng(99);
    std::normal_distribution<double> n01(0.0, 1.0);
    const int                        samples = 300;
    double                           sum     = 0.0;
    for (int k = 0; k < samples; ++k) {
        StateVector s  = s0;
        s.r           += Vec3{n01(rng), n01(rng), n01(rng)} * 100.0;
        s.v           += Vec3{n01(rng), n01(rng), n01(rng)} * 0.1;
        sum           += checkConsistency({kEpoch, s, P0}, truth).positionMahalanobis2;
    }
    EXPECT_NEAR(sum / samples, 3.0, 0.45); // standard error ~0.14
}

TEST(EstimateTracker, ElementTestHoldsOverLongGaps) {
    // Same setting: the full 6-D distance in elements stays chi-square with 6 dof, while the
    // Cartesian 6-D distance does not (measured mean ~78 here).
    const StateVector                s0 = leo();
    const Mat6                       P0 = diagonal(1e4, 1e-2);
    const double                     dt = 14400.0;
    const StateEstimate              truth{kEpoch + dt, truthAt(s0, dt), diagonal(1e-6, 1e-12)};
    std::mt19937_64                  rng(98);
    std::normal_distribution<double> n01(0.0, 1.0);
    const int                        samples  = 300;
    double                           elements = 0.0, cartesian = 0.0;
    for (int k = 0; k < samples; ++k) {
        StateVector s  = s0;
        s.r           += Vec3{n01(rng), n01(rng), n01(rng)} * 100.0;
        s.v           += Vec3{n01(rng), n01(rng), n01(rng)} * 0.1;
        const auto c   = checkConsistency({kEpoch, s, P0}, truth);
        elements      += c.elementsMahalanobis2;
        cartesian     += c.mahalanobis2;
    }
    EXPECT_NEAR(elements / samples, 6.0, 0.9); // standard error ~0.2
    EXPECT_GT(cartesian / samples, 20.0);
}

TEST(EstimateTracker, ManoeuvreIsFlaggedInTrack) {
    // A 0.5 m/s along-track burn 15 minutes after the prior.
    const StateVector s0    = leo();
    StateVector       burn  = truthAt(s0, 900.0);
    burn.v                 += normalized(burn.v) * 0.5;
    TrackerOptions tight;
    tight.integrator.relTol   = 1e-13;
    tight.integrator.absTol   = 1e-10;
    const StateVector   after = trackerPropagator({kEpoch + 900.0, burn, Mat6{}}, tight).propagate(900.0);
    const StateEstimate prior{kEpoch, s0, diagonal(100.0, 1e-4)};
    const StateEstimate next{kEpoch + 1800.0, after, diagonal(25.0, 1e-5)};
    const auto          c = checkConsistency(prior, next);
    ASSERT_TRUE(c.valid);
    EXPECT_FALSE(c.consistent);
    EXPECT_LT(c.pValue, 1e-6);
    // Mostly along-track and radial (the in-plane response to a tangential burn), not cross-track.
    EXPECT_GT(std::abs(c.positionResidualRic.y), 100.0);
    EXPECT_LT(std::abs(c.positionResidualRic.z), 1.0);
}

TEST(EstimateTracker, FreshBurnNeedsAVelocityTest) {
    // A 0.5 m/s burn 20 s before the delivery has barely moved the position: the element test
    // (default) and the Cartesian 6-D test see it, the position test does not.
    const StateVector s0    = leo();
    StateVector       burn  = truthAt(s0, 1780.0);
    burn.v                 += normalized(burn.v) * 0.5;
    TrackerOptions tight;
    tight.integrator.relTol   = 1e-13;
    tight.integrator.absTol   = 1e-10;
    const StateVector   after = trackerPropagator({kEpoch + 1780.0, burn, Mat6{}}, tight).propagate(20.0);
    const StateEstimate prior{kEpoch, s0, diagonal(100.0, 1e-4)};
    const StateEstimate next{kEpoch + 1800.0, after, diagonal(25.0, 1e-5)};
    TrackerOptions      position, full;
    position.test = ConsistencyTest::Position;
    full.test     = ConsistencyTest::PositionVelocity;
    EXPECT_EQ(TrackerOptions{}.test, ConsistencyTest::Elements);
    EXPECT_FALSE(checkConsistency(prior, next).consistent);
    EXPECT_FALSE(checkConsistency(prior, next, full).consistent);
    EXPECT_TRUE(checkConsistency(prior, next, position).consistent);
}

TEST(EstimateTracker, SemiMajorAxisResidualOfABurn) {
    // An in-track burn dv changes a by about 2 a^2 v dv / mu (vis-viva).
    const StateVector s0    = leo();
    StateVector       burn  = truthAt(s0, 1500.0);
    const double      dv    = 0.3;
    const double      a     = stateToElements(burn).a;
    const double      da    = 2.0 * a * a * norm(burn.v) * dv / constants::gmEarth;
    burn.v                 += normalized(burn.v) * dv;
    TrackerOptions tight;
    tight.integrator.relTol = 1e-13;
    tight.integrator.absTol = 1e-10;
    const StateVector after = trackerPropagator({kEpoch + 1500.0, burn, Mat6{}}, tight).propagate(300.0);
    const auto        c =
        checkConsistency({kEpoch, s0, diagonal(100.0, 1e-4)}, {kEpoch + 1800.0, after, diagonal(25.0, 1e-5)});
    EXPECT_NEAR(c.semiMajorAxisResidual / da, 1.0, 0.05);
    EXPECT_GT(c.semiMajorAxisSigma, 0.0);
    EXPECT_GT(c.semiMajorAxisResidual / c.semiMajorAxisSigma, 10.0);
    EXPECT_FALSE(c.consistent);
}

TEST(EstimateTracker, ElementTestIsMirrorInvariant) {
    // Retrograde orbits are mirrored (180 deg about x) before the element conversion. Mirroring the
    // whole problem by hand gives a prograde orbit; with point-mass gravity in GCRS (no pole, no
    // turning frame) the mirror is an exact symmetry, so the distance must agree.
    const StateVector s0  = elementsToState({7000 * km, 0.002, degToRad(140.0), 0.5, 1.0, 2.0});
    StateVector       sp  = s0;
    sp.r                 += Vec3{30.0, -10.0, 5.0};
    sp.v                 += Vec3{0.01, 0.02, -0.01};
    Mat6 P                = diagonal(400.0, 4e-4);
    P(0, 4) = P(4, 0) = 0.05;
    const StateEstimate prior{kEpoch, sp, P};
    const StateEstimate next{kEpoch + 3600.0, truthAt(s0, 3600.0), diagonal(1e-6, 1e-12)};
    const Mat3          M{1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0};
    auto                mirror = [&](const StateEstimate& e) {
        return StateEstimate{
            e.epoch,
                           {M * e.state.r, M * e.state.v},
            rotateCovariance(e.covariance, M)
        };
    };
    TrackerOptions pointMass;
    pointMass.body.J = {};
    const auto c     = checkConsistency(prior, next, pointMass);
    const auto cm    = checkConsistency(mirror(prior), mirror(next), pointMass);
    ASSERT_TRUE(c.valid);
    ASSERT_TRUE(cm.valid);
    EXPECT_GT(c.elementsMahalanobis2, 1.0);
    EXPECT_NEAR(c.elementsMahalanobis2 / cm.elementsMahalanobis2, 1.0, 1e-6);
    EXPECT_NEAR(c.mahalanobis2 / cm.mahalanobis2, 1.0, 1e-6);
}

TEST(EstimateTracker, LatestWins) {
    EstimateTracker tracker;
    EXPECT_FALSE(tracker.hasEstimate());
    EXPECT_THROW((void)tracker.current(), Error);
    EXPECT_THROW((void)tracker.predict(kEpoch), Error);

    const StateVector s0 = leo();
    EXPECT_FALSE(tracker.update({kEpoch, s0, diagonal(100.0, 1e-4)}).has_value());
    ASSERT_TRUE(tracker.hasEstimate());

    const StateEstimate next{kEpoch + 600.0, truthAt(s0, 600.0), diagonal(50.0, 5e-5)};
    const auto          c = tracker.update(next);
    ASSERT_TRUE(c.has_value());
    EXPECT_TRUE(c->consistent);
    EXPECT_TRUE(c->replaced);
    EXPECT_EQ(tracker.current().epoch, next.epoch);

    // An older delivery is checked but does not replace the current estimate.
    const auto stale = tracker.update({kEpoch + 300.0, truthAt(s0, 300.0), diagonal(50.0, 5e-5)});
    ASSERT_TRUE(stale.has_value());
    EXPECT_FALSE(stale->replaced);
    EXPECT_TRUE(stale->consistent);
    EXPECT_EQ(tracker.current().epoch, next.epoch);

    // Predictions start from the current estimate.
    const auto grid = tracker.predictGrid(0.0, 60.0, 31);
    EXPECT_EQ(grid.front().state, next.state);
    EXPECT_EQ(grid.front().covariance, next.covariance);
    const StateCovariance p = tracker.predict(next.epoch + 1800.0);
    astro::test::expectVecNear(p.state.r, grid.back().state.r, 1e-6);
    EXPECT_NEAR(p.covariance(1, 1) / grid.back().covariance(1, 1), 1.0, 1e-6);
}

TEST(EstimateTracker, InvalidCovarianceIsReported) {
    const StateVector s0 = leo();
    const auto        c  = checkConsistency({kEpoch, s0, Mat6{}}, {kEpoch + 60.0, truthAt(s0, 60.0), Mat6{}});
    EXPECT_FALSE(c.valid);
    EXPECT_FALSE(c.consistent);
    EXPECT_TRUE(std::isnan(c.mahalanobis2));
}

TEST(EstimateTracker, NonFiniteDeliveryIsRejected) {
    // Bad data upstream must not pass the check or replace the estimate.
    EstimateTracker   tracker;
    const StateVector s0 = leo();
    (void)tracker.update({kEpoch, s0, diagonal(100.0, 1e-4)});
    const double  nan = std::nan("");
    StateEstimate bad{kEpoch + 60.0, truthAt(s0, 60.0), diagonal(100.0, 1e-4)};
    bad.state.r.y = nan;
    EXPECT_THROW((void)tracker.update(bad), Error);
    bad                  = {kEpoch + 60.0, truthAt(s0, 60.0), diagonal(100.0, 1e-4)};
    bad.covariance(4, 4) = std::numeric_limits<double>::infinity();
    EXPECT_THROW((void)tracker.update(bad), Error);
    EXPECT_THROW((void)checkConsistency({kEpoch, s0, diagonal(100.0, 1e-4)}, bad), Error);
    EXPECT_EQ(tracker.current().epoch, kEpoch);
    EXPECT_EQ(tracker.current().state, s0);

    EXPECT_TRUE(std::isnan(chiSquareSurvival3(nan)));
    EXPECT_TRUE(std::isnan(chiSquareSurvival6(nan)));
    const std::array<double, 6> d{1.0, nan, 0.0, 0.0, 0.0, 0.0};
    EXPECT_FALSE(mahalanobis2(d, diagonal(1.0, 1.0)).has_value());
}

TEST(EstimateTracker, TrueOfDateDeliveries) {
    // Truth in GCRS; deliveries in TOD of their own epoch, converted with the full IAU 2006/2000A
    // model. The tracker must reproduce TOD-of-date at each output time (freezing TOD at the
    // estimate's epoch would be off by ~0.85 m after 4 h).
    const Time        t0 = Time::fromIso8601("2026-11-02T00:00:00Z");
    const StateVector sg = leo();
    TrackerOptions    tight;
    tight.integrator.relTol = 1e-13;
    tight.integrator.absTol = 1e-10;
    const auto truthG       = [&](double dt) { return trackerPropagator({t0, sg, Mat6{}}, tight).propagate(dt); };
    const auto toTod        = [&](const StateVector& s, double dt) {
        const Mat3 B = FrameTransformer(t0 + dt, ConstantEop{}).gcrsTo(Frame::TOD);
        return StateVector{B * s.r, B * s.v};
    };
    const Mat6 P = diagonal(100.0, 1e-4);

    TrackerOptions todOptions;
    todOptions.frame = Frame::TOD;
    EstimateTracker tod(todOptions), gcrs;
    (void)tod.update({t0, toTod(sg, 0.0), P});
    (void)gcrs.update({t0, sg, P});
    const auto ct = tod.update({t0 + 3600.0, toTod(truthG(3600.0), 3600.0), P});
    const auto cg = gcrs.update({t0 + 3600.0, truthG(3600.0), P});
    ASSERT_TRUE(ct && cg);
    EXPECT_TRUE(ct->consistent);
    EXPECT_LT(norm(ct->positionResidualRic), 0.1);
    // The check does not depend on the delivery frame.
    astro::test::expectVecNear(ct->positionResidualRic, cg->positionResidualRic, 0.05);
    EXPECT_NEAR(ct->positionSigmaRic.y / cg->positionSigmaRic.y, 1.0, 1e-9);

    // Predictions over 4 h, in TOD of each output time: the GCRS tracker's predictions converted with
    // the full model (which isolates the frame handling from the integrator's ~0.2 m at relTol 1e-9).
    const auto pt     = tod.predictGrid(0.0, 600.0, 25);
    const auto pg     = gcrs.predictGrid(0.0, 600.0, 25);
    double     ofDate = 0.0, frozen = 0.0, truth = 0.0;
    for (std::size_t k = 0; k < pt.size(); ++k) {
        const double dt = 3600.0 + 600.0 * static_cast<double>(k);
        ofDate          = std::max(ofDate, norm(pt[k].state.r - toTod(pg[k].state, dt).r));
        frozen          = std::max(frozen, norm(pt[k].state.r - toTod(pg[k].state, 3600.0).r));
        truth           = std::max(truth, norm(pt[k].state.r - toTod(truthG(dt), dt).r));
        // Covariances are the GCRS ones rotated into TOD of date.
        const Mat6 rot = rotateCovariance(pg[k].covariance, gcrsToCelestialMatrix(Frame::TOD, t0 + dt));
        EXPECT_NEAR(pt[k].covariance(1, 1) / rot(1, 1), 1.0, 1e-6);
    }
    EXPECT_LT(ofDate, 0.05);
    EXPECT_GT(frozen, 0.5);
    EXPECT_LT(truth, 0.5); // integrator error at relTol 1e-9
    // predict() agrees with the grid.
    astro::test::expectVecNear(tod.predict(t0 + 3600.0 + 7200.0).state.r, pt[12].state.r, 1e-3);
}

TEST(EstimateTracker, RejectsOptionsThatCannotWork) {
    // These were accepted, and then every update() threw from the propagation: a tracker that could
    // never take a delivery. The constructor refuses them.
    const auto rejected = [](auto&& set) {
        TrackerOptions options;
        set(options);
        try {
            const EstimateTracker tracker(options);
            return false;
        } catch (const Error&) {
            return true;
        }
    };
    const double nan = std::numeric_limits<double>::quiet_NaN();
    EXPECT_TRUE(rejected([](TrackerOptions& o) { o.accelerationNoisePsd = -1e-9; }));
    EXPECT_TRUE(rejected([&](TrackerOptions& o) { o.accelerationNoisePsd = nan; }));
    EXPECT_TRUE(rejected([&](TrackerOptions& o) { o.alarmProbability = nan; }));
    EXPECT_TRUE(rejected([](TrackerOptions& o) { o.alarmProbability = -0.1; }));
    EXPECT_TRUE(rejected([](TrackerOptions& o) { o.alarmProbability = 1.5; }));
    EXPECT_TRUE(rejected([](TrackerOptions& o) { o.body.mu = 0.0; }));
    EXPECT_TRUE(rejected([&](TrackerOptions& o) { o.body.re = nan; }));
    // The ends of the range are meaningful: 0 never flags a delivery, 1 flags every one.
    EXPECT_FALSE(rejected([](TrackerOptions& o) { o.alarmProbability = 0.0; }));
    EXPECT_FALSE(rejected([](TrackerOptions& o) { o.alarmProbability = 1.0; }));
    EXPECT_FALSE(rejected([](TrackerOptions& o) { o.accelerationNoisePsd = 0.0; }));
}

TEST(EstimateTracker, RejectsEarthFixedFrames) {
    TrackerOptions itrs;
    itrs.frame = Frame::ITRS;
    EXPECT_THROW(EstimateTracker{itrs}, Error);
    const StateEstimate e{kEpoch, leo(), diagonal(1.0, 1e-6)};
    EXPECT_THROW((void)checkConsistency(e, e, itrs), Error);
}

namespace {

/// R P R^T for a 3x3 block of a covariance, written out.
Mat3 rotatedBlock(const Mat6& P, std::size_t row, std::size_t col, const Mat3& R) {
    Mat3 out;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            double s = 0.0;
            for (std::size_t k = 0; k < 3; ++k) {
                for (std::size_t l = 0; l < 3; ++l) {
                    s += R(i, k) * P(row + k, col + l) * R(j, l);
                }
            }
            out(i, j) = s;
        }
    }
    return out;
}

} // namespace

TEST(EstimateTracker, ResidualsAndSigmasAreThoseOfThePredictionInRic) {
    // A delivery placed at known offsets from the prediction, along the prediction's radial,
    // in-track and cross-track axes: the residuals are those offsets, and the sigmas the square roots
    // of the diagonal of the predicted covariance turned into the same axes.
    const StateVector   s0 = leo();
    const StateEstimate prior{kEpoch, s0, diagonal(100.0, 1e-4)};
    const Time          t1 = kEpoch + 1200.0;
    const auto          c0 = checkConsistency(prior, {t1, truthAt(s0, 1200.0), diagonal(25.0, 1e-5)});
    const Mat3          R  = ricMatrix(c0.predicted);
    const Vec3          dr{1.0, -2.0, 3.0}, dv{0.01, -0.02, 0.03};
    const StateEstimate next{
        t1,
        {c0.predicted.r + R.transpose() * dr, c0.predicted.v + R.transpose() * dv},
        diagonal(25.0, 1e-5)
    };
    const auto c = checkConsistency(prior, next);
    // The prediction is the same (same prior, same epoch): residuals to rounding in coordinates of
    // 7e6 m and velocities of 7.7e3 m/s (measured 3.3e-10 m; below 1e-11 m/s).
    test::expectVecNear(c.positionResidualRic, dr, 2e-9, "position residual");
    test::expectVecNear(c.velocityResidualRic, dv, 1e-11, "velocity residual");
    const Mat3 Prr = rotatedBlock(c.predictedCovariance, 0, 0, R);
    const Mat3 Pvv = rotatedBlock(c.predictedCovariance, 3, 3, R);
    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_NEAR(c.positionSigmaRic[i], std::sqrt(Prr(i, i)), 1e-14 * std::sqrt(Prr(i, i))) << i;
        EXPECT_NEAR(c.velocitySigmaRic[i], std::sqrt(Pvv(i, i)), 1e-14 * std::sqrt(Pvv(i, i))) << i;
    }
    // The in-track uncertainty has grown most over 20 minutes.
    EXPECT_GT(c.positionSigmaRic.y, c.positionSigmaRic.z);
}

TEST(EstimateTracker, AlarmProbabilityIsTheThresholdOnThePValue) {
    // The selected test's p-value against TrackerOptions::alarmProbability: consistent at or above it.
    const StateVector   s0 = leo();
    const StateEstimate prior{kEpoch, s0, diagonal(100.0, 1e-4)};
    StateVector         moved  = truthAt(s0, 1200.0);
    moved.r                   += Vec3{20.0, -15.0, 10.0}; // a few sigma: a p-value well inside (0, 1)
    const StateEstimate next{kEpoch + 1200.0, moved, diagonal(25.0, 1e-5)};
    const double        p = checkConsistency(prior, next).elementsPValue;
    ASSERT_GT(p, 1e-12);
    ASSERT_LT(p, 0.5);
    TrackerOptions below, above;
    below.alarmProbability = 0.5 * p;
    above.alarmProbability = 2.0 * p;
    EXPECT_TRUE(checkConsistency(prior, next, below).consistent);
    EXPECT_FALSE(checkConsistency(prior, next, above).consistent);
}

TEST(EstimateTracker, AccelerationNoiseAddsTheWhiteNoiseCovariance) {
    // White acceleration noise of spectral density q adds, over dt, q dt^3/3 to each position
    // variance, q dt^2/2 to each position-velocity covariance and q dt to each velocity variance (a
    // free particle; gravity gradients change this by about (n dt)^2 = 4e-3 of it over 60 s; measured
    // 1.5e-3). Tolerance 1e-2.
    const StateVector s0 = leo();
    const double      q = 1e-6, dt = 60.0;
    TrackerOptions    noisy;
    noisy.accelerationNoisePsd = q;
    const StateEstimate prior{kEpoch, s0, diagonal(100.0, 1e-4)};
    const StateEstimate next{kEpoch + dt, truthAt(s0, dt), diagonal(25.0, 1e-5)};
    const Mat6          quiet = checkConsistency(prior, next).predictedCovariance;
    const Mat6          noise = checkConsistency(prior, next, noisy).predictedCovariance;
    for (std::size_t i = 0; i < 3; ++i) {
        SCOPED_TRACE(testing::Message() << "axis " << i);
        EXPECT_NEAR((noise(i, i) - quiet(i, i)) / (q * dt * dt * dt / 3.0), 1.0, 1e-2);
        EXPECT_NEAR((noise(i, i + 3) - quiet(i, i + 3)) / (q * dt * dt / 2.0), 1.0, 1e-2);
        EXPECT_NEAR((noise(i + 3, i + 3) - quiet(i + 3, i + 3)) / (q * dt), 1.0, 1e-2);
    }
}

TEST(EstimateTracker, PredictGridWithoutAnEstimateOrWithAZeroStepThrows) {
    EstimateTracker tracker;
    EXPECT_THROW((void)tracker.predictGrid(0.0, 60.0, 3), Error);
    (void)tracker.update({kEpoch, leo(), diagonal(100.0, 1e-4)});
    EXPECT_NO_THROW((void)tracker.predictGrid(0.0, 60.0, 3));
    EXPECT_THROW((void)tracker.predictGrid(0.0, 0.0, 3), Error);
}

TEST(EstimateTracker, DeliveriesInTheOtherCelestialFramesGiveTheSameCheck) {
    // As TrueOfDateDeliveries for the frames it does not cover: J2000, MOD, TEME and CIRS. The same
    // truth, delivered in frame f of each epoch (converted with the full IAU 2006/2000A model), must
    // give the residuals and sigmas of the GCRS deliveries, and predictions that are the GCRS ones in
    // f of the output time.
    const Time        t0 = Time::fromIso8601("2026-11-02T00:00:00Z");
    const StateVector sg = leo();
    TrackerOptions    tight;
    tight.integrator.relTol = 1e-13;
    tight.integrator.absTol = 1e-10;
    const StateVector truth = trackerPropagator({t0, sg, Mat6{}}, tight).propagate(1800.0);
    const Mat6        P     = diagonal(100.0, 1e-4);
    EstimateTracker   gcrs;
    (void)gcrs.update({t0, sg, P});
    const auto cg = gcrs.update({t0 + 1800.0, truth, P});
    ASSERT_TRUE(cg);
    const auto pg = gcrs.predictGrid(0.0, 1800.0, 3);
    for (Frame f : {Frame::J2000, Frame::MOD, Frame::TEME, Frame::CIRS}) {
        SCOPED_TRACE(testing::Message() << f);
        const auto into = [&](const StateVector& s, double dt) {
            const Mat3 B = FrameTransformer(t0 + dt, ConstantEop{}).gcrsTo(f);
            return StateVector{B * s.r, B * s.v};
        };
        TrackerOptions options;
        options.frame = f;
        EstimateTracker tracker(options);
        (void)tracker.update({t0, into(sg, 0.0), P});
        const auto c = tracker.update({t0 + 1800.0, into(truth, 1800.0), P});
        ASSERT_TRUE(c);
        // The tracker converts with IAU 2000B, the deliveries here with 2000A (up to 1 mas apart, 3 cm at
        // 7000 km), but in and out alike, so most of it cancels: measured 1.6e-4 m, and 2e-12 in the
        // sigma. Tolerances 1e-3 m and 1e-10.
        test::expectVecNear(c->positionResidualRic, cg->positionResidualRic, 1e-3, "residual");
        EXPECT_NEAR(c->positionSigmaRic.y / cg->positionSigmaRic.y, 1.0, 1e-10);
        // The free function gives the same check.
        const ConsistencyCheck free =
            checkConsistency({t0, into(sg, 0.0), P}, {t0 + 1800.0, into(truth, 1800.0), P}, options);
        test::expectVecNear(free.positionResidualRic, cg->positionResidualRic, 1e-3, "checkConsistency");
        const auto p = tracker.predictGrid(0.0, 1800.0, 3);
        for (std::size_t k = 0; k < p.size(); ++k) {
            test::expectVecNear(p[k].state.r,
                                into(pg[k].state, 1800.0 * static_cast<double>(k + 1)).r,
                                1e-3,
                                "prediction");
        }
    }
}
