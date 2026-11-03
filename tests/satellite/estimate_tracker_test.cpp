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

#include <cmath>
#include <random>
#include <vector>

#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/orbit/elements.h"
#include "astro/satellite/estimate_tracker.h"
#include "astro/satellite/orbital_frames.h"
#include "test_util.h"

using namespace astro;

namespace {

constexpr double km = 1000.0;

const Time kEpoch = Time::fromIso8601("2026-03-01T00:00:00Z");

StateVector leo() {
    return elementsToState({6778 * km, 0.0005, deg2rad(51.6), 0.3, 0.5, 0.7});
}

Mat6 diagonal(double positionVar, double velocityVar) {
    Mat6 P{};
    for (std::size_t i = 0; i < 3; ++i) {
        P[i][i]         = positionVar;
        P[i + 3][i + 3] = velocityVar;
    }
    return P;
}

/// "Truth": the tracker's own force model with a tight tolerance.
StateVector truthAt(const StateVector& s0, double dt) {
    TrackerOptions tight;
    tight.integrator.relTol = 1e-13;
    tight.integrator.absTol = 1e-10;
    return makePropagator({kEpoch, s0, Mat6{}}, tight).propagate(dt);
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

TEST(EstimateTracker, Mahalanobis) {
    const Mat6                  D = diagonal(4.0, 1e-4);
    const std::array<double, 6> d{2.0, 0.0, -4.0, 0.01, 0.0, 0.0};
    EXPECT_NEAR(*mahalanobis2(d, D), 1.0 + 4.0 + 1.0, 1e-12);
    EXPECT_NEAR(*mahalanobis2(d, D, 3), 5.0, 1e-12);
    // Correlated S = A A^T and d = S x give d^T S^-1 d = x^T S x.
    Mat6 A{};
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            A[i][j] = (i == j ? 2.0 : 0.3) * (i < 3 ? 10.0 : 0.01);
        }
    }
    const Mat6                  S = multiply(A, transpose(A));
    const std::array<double, 6> x{0.1, -0.2, 0.05, 30.0, -10.0, 4.0};
    std::array<double, 6>       Sx{};
    double                      xSx = 0.0;
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            Sx[i] += S[i][j] * x[j];
        }
    }
    for (std::size_t i = 0; i < 6; ++i) {
        xSx += x[i] * Sx[i];
    }
    EXPECT_NEAR(*mahalanobis2(Sx, S) / xSx, 1.0, 1e-10);
    EXPECT_FALSE(mahalanobis2(d, Mat6{}).has_value());
    Mat6 bad  = D;
    bad[0][1] = bad[1][0] = 5.0; // |rho| > 1
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
    const StateVector   after = makePropagator({kEpoch + 900.0, burn, Mat6{}}, tight).propagate(900.0);
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
    const StateVector   after = makePropagator({kEpoch + 1780.0, burn, Mat6{}}, tight).propagate(20.0);
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
    const StateVector after = makePropagator({kEpoch + 1500.0, burn, Mat6{}}, tight).propagate(300.0);
    const auto        c =
        checkConsistency({kEpoch, s0, diagonal(100.0, 1e-4)}, {kEpoch + 1800.0, after, diagonal(25.0, 1e-5)});
    EXPECT_NEAR(c.semiMajorAxisResidual / da, 1.0, 0.05);
    EXPECT_GT(c.semiMajorAxisSigma, 0.0);
    EXPECT_GT(c.semiMajorAxisResidual / c.semiMajorAxisSigma, 10.0);
    EXPECT_FALSE(c.consistent);
}

TEST(EstimateTracker, ElementTestIsMirrorInvariant) {
    // Retrograde orbits are mirrored (180 deg about x) before the element conversion. Mirroring the
    // whole problem by hand gives a prograde orbit with the same zonal dynamics (J2 and J4 are even),
    // so the distance must agree.
    const StateVector s0  = elementsToState({7000 * km, 0.002, deg2rad(140.0), 0.5, 1.0, 2.0});
    StateVector       sp  = s0;
    sp.r                 += Vec3{30.0, -10.0, 5.0};
    sp.v                 += Vec3{0.01, 0.02, -0.01};
    Mat6 P                = diagonal(400.0, 4e-4);
    P[0][4] = P[4][0] = 0.05;
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
    const auto c  = checkConsistency(prior, next);
    const auto cm = checkConsistency(mirror(prior), mirror(next));
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
    EXPECT_NEAR(p.covariance[1][1] / grid.back().covariance[1][1], 1.0, 1e-6);
}

TEST(EstimateTracker, InvalidCovarianceIsReported) {
    const StateVector s0 = leo();
    const auto        c  = checkConsistency({kEpoch, s0, Mat6{}}, {kEpoch + 60.0, truthAt(s0, 60.0), Mat6{}});
    EXPECT_FALSE(c.valid);
    EXPECT_FALSE(c.consistent);
    EXPECT_TRUE(std::isnan(c.mahalanobis2));
}
