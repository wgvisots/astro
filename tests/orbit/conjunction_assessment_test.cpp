/**********************************************************************
 * conjunction_assessment_test.cpp
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

// assessConjunction(s): closest approach, covariances propagated to it, and collision probability,
// from states and covariances.

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/orbit/close_approach.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/gravity.h"
#include "astro/orbit/numerical.h"
#include "astro/orbit/two_body.h"

using namespace astro;

namespace {

constexpr double km = 1000.0;

/// Two LEO objects whose two-body trajectories cross at kTca, 1 h after the epoch, 15 m apart
/// along the radius (perpendicular to both velocities, which are 60 deg apart): an exact closest
/// approach.
struct Crossing60 {
    static constexpr double kLead = 3600.0; // s from the epoch to the TCA
    static constexpr double kMiss = 15.0;   // m

    Time        epoch = Time::fromIso8601("2026-11-11T11:00:00Z");
    Time        tca   = epoch + kLead;
    StateVector first;  // at the epoch
    StateVector second; // at the epoch

    Crossing60() {
        first                = elementsToState({7000 * km, 0.0, deg2rad(51.6), 0.3, 0.0, 0.0});
        const StateVector a  = propagateKepler(first, kLead);
        const Vec3        up = normalized(a.r), along = normalized(a.v), across = cross(up, along);
        const Vec3        v2 = norm(a.v) * (std::cos(deg2rad(60.0)) * along + std::sin(deg2rad(60.0)) * across);
        second               = propagateKepler({a.r + kMiss * normalized(cross(a.v, v2)), v2}, -kLead);
    }
};

Mat6 diagonalCovariance(double positionSigma, double velocitySigma) {
    Mat6 P{};
    for (std::size_t i = 0; i < 3; ++i) {
        P[i][i]         = positionSigma * positionSigma;
        P[i + 3][i + 3] = velocitySigma * velocitySigma;
    }
    return P;
}

void expectMatNear(const Mat6& a, const Mat6& b, double relTol, const char* what) {
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            EXPECT_NEAR(a[i][j], b[i][j], relTol * std::sqrt(b[i][i] * b[j][j])) << what << " " << i << " " << j;
        }
    }
}

} // namespace

TEST(ConjunctionAssessment, TwoBodyIsTheChainOfItsParts) {
    const Crossing60 c;
    const Mat6       P1 = diagonalCovariance(10.0, 0.005), P2 = diagonalCovariance(30.0, 0.02);
    const Time       start = c.tca - 600.0, end = c.tca + 600.0;

    const ConjunctionAssessment a =
        assessConjunction({c.epoch, c.first, P1}, {c.epoch, c.second, P2}, start, end, 20.0);

    // The same approach as closestApproach, the covariances of propagateKeplerWithCovariance, and
    // the probability of collisionProbability.
    const CloseApproach ca = closestApproach({c.epoch, c.first}, {c.epoch, c.second}, start, end);
    EXPECT_EQ(a.approach.time - ca.time, 0.0);
    EXPECT_EQ(a.approach.distance, ca.distance);
    EXPECT_NEAR(a.approach.time - c.tca, 0.0, 1e-6);
    EXPECT_NEAR(a.approach.distance, Crossing60::kMiss, 1e-6);
    const StateCovariance s1 = propagateKeplerWithCovariance(c.first, a.approach.time - c.epoch, P1);
    const StateCovariance s2 = propagateKeplerWithCovariance(c.second, a.approach.time - c.epoch, P2);
    EXPECT_EQ(a.firstCovariance, s1.covariance);
    EXPECT_EQ(a.secondCovariance, s2.covariance);
    EXPECT_LT(norm(a.approach.first.r - s1.state.r), 1e-6); // the same time, by two roundings
    EXPECT_EQ(a.collision.probability, collisionProbability(ca, s1.covariance, s2.covariance, 20.0).probability);
    EXPECT_GT(a.collision.probability, 1e-3);
}

TEST(ConjunctionAssessment, MatchesMonteCarloOfInitialStates) {
    // End to end: draw both initial states from their covariances, propagate each draw, find its
    // closest approach, and count the draws that pass within the hard-body radius. Checks the
    // propagated covariances (linear over the hour) and the short-term encounter model together.
    // The radius is large enough for Pc ~ 0.3, so that 100 000 draws resolve it to 0.5 %.
    const Crossing60            c;
    const double                sr = 10.0, sv = 0.005, hbr = 60.0;
    const Mat6                  P = diagonalCovariance(sr, sv);
    const ConjunctionAssessment a =
        assessConjunction({c.epoch, c.first, P}, {c.epoch, c.second, P}, c.tca - 600.0, c.tca + 600.0, hbr);
    const double pc = a.collision.probability;

    std::mt19937_64                  rng(3);
    std::normal_distribution<double> n01;
    const auto                       draw = [&](const StateVector& s) {
        return StateVector{
            s.r + sr * Vec3{n01(rng), n01(rng), n01(rng)},
            s.v + sv * Vec3{n01(rng), n01(rng), n01(rng)}
        };
    };
    const int N    = 100000;
    int       hits = 0;
    for (int k = 0; k < N; ++k) {
        const StateVector a1 = draw(c.first), a2 = draw(c.second);
        hits += closestApproach({c.epoch, a1}, {c.epoch, a2}, c.tca - 60.0, c.tca + 60.0).distance < hbr ? 1 : 0;
    }
    const double mc    = static_cast<double>(hits) / N;
    const double sigma = std::sqrt(pc * (1.0 - pc) / N);
    EXPECT_NEAR(mc, pc, 4.0 * sigma) << "Pc " << pc; // 4 sigma is ~2 % of Pc
    EXPECT_GT(pc, 0.2);
}

TEST(ConjunctionAssessment, SecularModel) {
    const Crossing60                 c;
    const Mat6                       P = diagonalCovariance(10.0, 0.005);
    const UncertainTrajectory        first{c.epoch, c.first, P, TrajectoryModel::SecularJ2J4};
    const UncertainTrajectory        second{c.epoch, c.second, P, TrajectoryModel::SecularJ2J4};
    const SecularPropagatorCorrected prop(c.first);

    // The covariance is SecularPropagatorCorrected's, and its state is the trajectory's.
    const Time            t  = c.epoch + 1234.0;
    const StateCovariance sc = first.stateCovariance(t);
    EXPECT_EQ(sc.state, first.trajectory().state(t));
    expectMatNear(sc.covariance, prop.propagateWithCovariance(1234.0, P).covariance, 1e-8, "covariance");

    // J2 moves the approach; the assessment follows the secular trajectories.
    const ConjunctionAssessment a = assessConjunction(first, second, c.epoch, c.epoch + 7200.0, 20.0);
    const CloseApproach ca        = closestApproach(first.trajectory(), second.trajectory(), c.epoch, c.epoch + 7200.0);
    EXPECT_EQ(a.approach.time - ca.time, 0.0);
    EXPECT_EQ(a.firstCovariance, first.stateCovariance(ca.time).covariance);
    EXPECT_EQ(a.collision.probability,
              collisionProbability(ca, a.firstCovariance, a.secondCovariance, 20.0).probability);
}

TEST(ConjunctionAssessment, AnyTrajectoryWithACovarianceFunction) {
    // A point-mass NumericalPropagator through the function constructor gives the two-body answer.
    const Crossing60 c;
    const Mat6       P = diagonalCovariance(10.0, 0.005);
    AdaptiveOptions  tight;
    tight.relTol                   = 1e-13;
    tight.absTol                   = 1e-10;
    const auto                mass = std::make_shared<GravityField>(CentralBody::pointMass());
    const auto                prop = std::make_shared<NumericalPropagator>(c.epoch,
                                                            c.second,
                                                            std::vector<std::shared_ptr<const ForceModel>>{mass},
                                                            tight);
    const UncertainTrajectory numerical{Trajectory([prop](const Time& t) { return prop->propagate(t); }),
                                        [prop, P, epoch = c.epoch](const Time& t) {
                                            return prop->propagateGridWithCovariance(t - epoch, 1.0, 1, P)[0];
                                        }};

    const UncertainTrajectory   first{c.epoch, c.first, P};
    const ConjunctionAssessment viaFunction = assessConjunction(first, numerical, c.tca - 600.0, c.tca + 600.0, 20.0);
    const ConjunctionAssessment twoBody =
        assessConjunction(first, {c.epoch, c.second, P}, c.tca - 600.0, c.tca + 600.0, 20.0);
    EXPECT_NEAR(viaFunction.approach.time - twoBody.approach.time, 0.0, 1e-6);
    expectMatNear(viaFunction.secondCovariance, twoBody.secondCovariance, 1e-7, "covariance");
    EXPECT_NEAR(viaFunction.collision.probability, twoBody.collision.probability, 1e-6 * twoBody.collision.probability);
}

TEST(ConjunctionAssessment, ListAssessesEveryApproach) {
    const Crossing60          c;
    const Mat6                P = diagonalCovariance(10.0, 0.005);
    const UncertainTrajectory first{c.epoch, c.first, P}, second{c.epoch, c.second, P};
    const Time                end = c.epoch + 86400.0;

    const std::vector<ConjunctionAssessment> all = assessConjunctions(first, second, c.epoch, end, 20.0);
    const std::vector<CloseApproach> list = closeApproaches(first.trajectory(), second.trajectory(), c.epoch, end);
    ASSERT_EQ(all.size(), list.size());
    ASSERT_GT(all.size(), 20u); // two crossings per revolution
    for (std::size_t k = 0; k < all.size(); ++k) {
        EXPECT_EQ(all[k].approach.time - list[k].time, 0.0) << k;
        EXPECT_EQ(all[k].secondCovariance, second.stateCovariance(list[k].time).covariance) << k;
    }
    // With a threshold: the same orbit radius and period bring the objects back together at
    // either node every half revolution, so several pass; the designed one is among them.
    const std::vector<ConjunctionAssessment> close = assessConjunctions(first, second, c.epoch, end, 20.0, 1 * km);
    ASSERT_GE(close.size(), 1u);
    ASSERT_LT(close.size(), all.size());
    bool found = false;
    for (const ConjunctionAssessment& a : close) {
        EXPECT_LE(a.approach.distance, 1 * km);
        found = found || std::abs(a.approach.time - c.tca) < 1e-6;
    }
    EXPECT_TRUE(found);
}

TEST(ConjunctionAssessment, Errors) {
    const Crossing60 c;
    const Mat6       P   = diagonalCovariance(10.0, 0.005);
    Mat6             bad = P;
    bad[2][4]            = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW((void)UncertainTrajectory(c.epoch, c.first, bad), Error);
    EXPECT_THROW(
        (void)UncertainTrajectory(Trajectory{c.epoch, c.first}, UncertainTrajectory::StateCovarianceFunction{}),
        Error);
    const UncertainTrajectory a{c.epoch, c.first, P}, b{c.epoch, c.second, P};
    // The radius is checked even when no approach is found.
    EXPECT_THROW((void)assessConjunctions(a, b, c.epoch, c.epoch + 60.0, 0.0, 1.0), Error);
    EXPECT_THROW((void)assessConjunction(a, b, c.epoch, c.epoch + 60.0, -1.0), Error);
    EXPECT_THROW((void)assessConjunction(a, b, c.epoch + 60.0, c.epoch, 20.0), Error); // window backwards
}

TEST(ConjunctionAssessment, ApproachAtAnEndOfTheWindowIsFlaggedAndAssessedThere) {
    // A window that closes 100 s before the designed approach (or opens 100 s after it): the
    // closest distance in it is at that end, flagged atBoundary. The states, distance and
    // covariances are those at the end, from the two-body propagation written in terms of
    // propagateKepler (checked against an extended-precision reference elsewhere), and the
    // probability is the straight-line encounter through them. That projects the 750 km of relative
    // motion out of the encounter plane, so it stays close to the true probability (within 10 %;
    // the covariances are those 100 s away).
    const Crossing60            c;
    const Mat6                  P1 = diagonalCovariance(10.0, 0.005), P2 = diagonalCovariance(30.0, 0.02);
    const UncertainTrajectory   first{c.epoch, c.first, P1}, second{c.epoch, c.second, P2};
    const ConjunctionAssessment inside = assessConjunction(first, second, c.tca - 600.0, c.tca + 600.0, 20.0);
    EXPECT_FALSE(inside.approach.atBoundary);
    for (const auto& [start, end, edge] : {
             std::tuple{c.tca - 600.0, c.tca - 100.0, c.tca - 100.0},
             std::tuple{c.tca + 100.0, c.tca + 600.0, c.tca + 100.0}
    }) {
        const ConjunctionAssessment a = assessConjunction(first, second, start, end, 20.0);
        EXPECT_TRUE(a.approach.atBoundary);
        EXPECT_EQ(a.approach.time - edge, 0.0);
        const double      dt = edge - c.epoch;
        const StateVector s1 = propagateKepler(c.first, dt), s2 = propagateKepler(c.second, dt);
        EXPECT_LT(norm(a.approach.first.r - s1.r), 1e-6);
        EXPECT_LT(norm(a.approach.second.r - s2.r), 1e-6);
        EXPECT_NEAR(a.approach.distance, norm(s1.r - s2.r), 1e-6);
        EXPECT_GT(a.approach.distance, 500 * km); // 100 s at the 7.5 km/s relative speed
        expectMatNear(a.firstCovariance, propagateKeplerWithCovariance(c.first, dt, P1).covariance, 1e-12, "first");
        expectMatNear(a.secondCovariance, propagateKeplerWithCovariance(c.second, dt, P2).covariance, 1e-12, "second");
        EXPECT_EQ(a.collision.probability,
                  collisionProbability(a.approach, a.firstCovariance, a.secondCovariance, 20.0).probability);
        EXPECT_NEAR(a.collision.probability / inside.collision.probability, 1.0, 0.1);
    }
}

TEST(ConjunctionAssessment, CovarianceFunctionThatIsNotFiniteIsAnError) {
    // A covariance from a user function is not checked when the trajectory is built, but the
    // probability refuses a position covariance that is not finite; the velocity block is not used.
    // An exception thrown by the function passes through.
    const Crossing60          c;
    const Mat6                P = diagonalCovariance(10.0, 0.005);
    const UncertainTrajectory first{c.epoch, c.first, P};
    const auto                withCovariance = [&](double value, std::size_t i, std::size_t j) {
        return UncertainTrajectory(Trajectory{c.epoch, c.second}, [=, epoch = c.epoch, s = c.second](const Time& t) {
            StateCovariance out  = propagateKeplerWithCovariance(s, t - epoch, P);
            out.covariance[i][j] = value;
            out.covariance[j][i] = value;
            return out;
        });
    };
    const double nan = std::numeric_limits<double>::quiet_NaN(), inf = std::numeric_limits<double>::infinity();
    const Time   start = c.tca - 600.0, end = c.tca + 600.0;
    EXPECT_THROW((void)assessConjunction(first, withCovariance(nan, 0, 1), start, end, 20.0), Error);
    EXPECT_THROW((void)assessConjunction(first, withCovariance(inf, 2, 2), start, end, 20.0), Error);
    EXPECT_THROW((void)assessConjunctions(first, withCovariance(nan, 1, 1), start, end, 20.0), Error);
    const ConjunctionAssessment clean    = assessConjunction(first, {c.epoch, c.second, P}, start, end, 20.0);
    const ConjunctionAssessment velocity = assessConjunction(first, withCovariance(nan, 4, 4), start, end, 20.0);
    EXPECT_EQ(velocity.collision.probability, clean.collision.probability);
    const UncertainTrajectory throwing(Trajectory{c.epoch, c.second}, [](const Time&) -> StateCovariance {
        throw std::runtime_error("no covariance");
    });
    EXPECT_THROW((void)assessConjunction(first, throwing, start, end, 20.0), std::runtime_error);
}
