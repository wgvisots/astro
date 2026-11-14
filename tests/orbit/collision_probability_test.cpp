/**********************************************************************
 * collision_probability_test.cpp
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
#include <limits>
#include <random>

#include "astro/core/constants.h"
#include "astro/core/error.h"
#include "astro/orbit/close_approach.h"
#include "orbit/conjunction_reference.h"
#include "test_util.h"

using namespace astro;
using namespace astro::test;

namespace {

/// Isotropic Gaussian (sigma) at distance d from the centre: the Rice distribution's CDF at R,
/// integral of (r / s^2) exp(-(r^2 + d^2) / (2 s^2)) I0(r d / s^2) dr, by Simpson.
double riceReference(double d, double sigma, double R) {
    const double s2 = sigma * sigma;
    const auto   f  = [&](double r) {
        return r / s2 * std::exp(-(r * r + d * d) / (2.0 * s2)) * std::cyl_bessel_i(0.0, r * d / s2);
    };
    const int    n   = 20000;
    const double h   = R / n;
    double       sum = f(0.0) + f(R);
    for (int i = 1; i < n; ++i) {
        sum += (i % 2 ? 4.0 : 2.0) * f(i * h);
    }
    return sum * h / 3.0;
}

/// A closest approach at the origin of relative position `dr` and velocity `dv` (object 1 at rest).
CloseApproach approachOf(const Vec3& dr, const Vec3& dv) {
    CloseApproach c;
    c.first = {
        Vec3{7000e3, 0.0,    0.0},
        Vec3{0.0,    7500.0, 0.0}
    };
    c.second        = {c.first.r + dr, c.first.v + dv};
    c.distance      = norm(dr);
    c.relativeSpeed = norm(dv);
    return c;
}

} // namespace

TEST(ProbabilityInCircle, IsotropicClosedForms) {
    // Centred: 1 - exp(-R^2 / (2 sigma^2)), from the tiny to the near-certain.
    for (double ratio : {1e-3, 0.1, 1.0, 3.0, 8.0}) {
        const double sigma = 100.0, R = ratio * sigma;
        const double exact = -std::expm1(-0.5 * ratio * ratio);
        EXPECT_NEAR(probabilityInCircle(0.0, 0.0, sigma * sigma, 0.0, sigma * sigma, R), exact, 1e-11 * exact) << ratio;
    }
    // Offset: the Rice distribution, in any direction of the offset.
    for (double d : {0.5, 2.0, 5.0}) {
        for (double ratio : {0.2, 1.0, 4.0}) {
            const double sigma = 10.0, R = ratio * sigma;
            const double ref = riceReference(d * sigma, sigma, R);
            const double c = std::cos(0.7), s = std::sin(0.7);
            EXPECT_NEAR(probabilityInCircle(d * sigma * c, d * sigma * s, sigma * sigma, 0.0, sigma * sigma, R),
                        ref,
                        1e-10 * ref)
                << d << " " << ratio;
        }
    }
}

TEST(ProbabilityInCircle, AnisotropicMatchesBruteForce) {
    struct Case {
        double mx, my, sx, sy, rho, R;
    };
    const Case cases[] = {
        {100.0,  50.0,   400.0,  100.0, 0.3,  20.0}, // typical: miss ~ sigma, R << sigma
        {30.0,   -5.0,   1000.0, 1.0,   0.0,  10.0}, // aspect ratio 1000
        {0.0,    40.0,   5.0,    3.0,   -0.5, 50.0}, // R >> sigma, mean inside near the edge
        {800.0,  -300.0, 100.0,  60.0,  0.2,  15.0}, // far: probability 6e-26
        {-200.0, 150.0,  300.0,  250.0, 0.95, 5.0 }, // strongly correlated
    };
    for (const Case& c : cases) {
        const double xx = c.sx * c.sx, yy = c.sy * c.sy, xy = c.rho * c.sx * c.sy;
        const double p   = probabilityInCircle(c.mx, c.my, xx, xy, yy, c.R);
        const double ref = polarReference(c.mx, c.my, xx, xy, yy, c.R);
        EXPECT_GT(p, 0.0);
        EXPECT_NEAR(p, ref, 1e-8 * ref) << c.mx << " " << c.my << " " << c.sx << " " << c.sy;
    }
}

TEST(ProbabilityInCircle, LimitsSymmetriesAndSingularCovariances) {
    // Small circle: pi R^2 times the density at the centre, to O(R^2 / sigma^2).
    const double xx = 400.0 * 400.0, yy = 150.0 * 150.0, xy = 0.4 * 400.0 * 150.0;
    const double mx = 300.0, my = -100.0, R = 0.05;
    const double det  = xx * yy - xy * xy;
    const double q    = (yy * mx * mx - 2.0 * xy * mx * my + xx * my * my) / det;
    const double dens = std::exp(-0.5 * q) / (constants::twoPi * std::sqrt(det));
    EXPECT_NEAR(probabilityInCircle(mx, my, xx, xy, yy, R),
                constants::pi * R * R * dens,
                1e-6 * constants::pi * R * R * dens);

    // Rotating mean and covariance together changes nothing.
    const double p = probabilityInCircle(mx, my, xx, xy, yy, 30.0);
    for (double a : {0.3, 1.9, -2.5}) {
        const double c = std::cos(a), s = std::sin(a);
        const double rxx = c * c * xx - 2 * c * s * xy + s * s * yy;
        const double ryy = s * s * xx + 2 * c * s * xy + c * c * yy;
        const double rxy = c * s * (xx - yy) + (c * c - s * s) * xy;
        EXPECT_NEAR(probabilityInCircle(c * mx - s * my, s * mx + c * my, rxx, rxy, ryy, 30.0), p, 1e-11 * p) << a;
    }

    // Singular covariances: no uncertainty, and all of it along one line (x, offset 6 in y).
    EXPECT_EQ(probabilityInCircle(3.0, 4.0, 0.0, 0.0, 0.0, 5.1), 1.0);
    EXPECT_EQ(probabilityInCircle(3.0, 4.0, 0.0, 0.0, 0.0, 4.9), 0.0);
    const double sigma = 20.0, half = std::sqrt(10.0 * 10.0 - 6.0 * 6.0); // chord |x| < 8
    const double line =
        0.5 * (std::erf((half - 5.0) / (std::sqrt(2.0) * sigma)) - std::erf((-half - 5.0) / (std::sqrt(2.0) * sigma)));
    EXPECT_NEAR(probabilityInCircle(5.0, 6.0, sigma * sigma, 0.0, 0.0, 10.0), line, 1e-14);
    EXPECT_EQ(probabilityInCircle(5.0, 11.0, sigma * sigma, 0.0, 0.0, 10.0), 0.0);

    // Far beyond 15 sigma the probability is 0.
    EXPECT_EQ(probabilityInCircle(1e5, 0.0, 100.0, 0.0, 100.0, 10.0), 0.0);

    EXPECT_THROW((void)probabilityInCircle(0, 0, 1, 0, 1, 0.0), Error);
    EXPECT_THROW((void)probabilityInCircle(0, 0, 1, 0, 1, std::numeric_limits<double>::infinity()), Error);
    EXPECT_THROW((void)probabilityInCircle(0, 0, 1, 2, 1, 1.0), Error); // not positive semi-definite
    EXPECT_THROW((void)probabilityInCircle(0, 0, -1, 0, 1, 1.0), Error);
    EXPECT_THROW((void)probabilityInCircle(std::nan(""), 0, 1, 0, 1, 1.0), Error);
}

TEST(CollisionProbability, ProjectsOntoTheEncounterPlane) {
    // Relative velocity along a skewed direction; miss perpendicular to it.
    const Vec3          dv = Vec3{3000.0, -9000.0, 4000.0};
    const Vec3          z  = normalized(dv);
    const Vec3          m  = 120.0 * normalized(cross(z, Vec3{0.0, 0.0, 1.0}));
    const CloseApproach ca = approachOf(m, dv);
    // Combined covariance: 200 m isotropic, plus 5 km along the relative velocity, which cannot
    // matter: an error along the line of motion only shifts the time of the encounter.
    Mat3 c1  = Mat3::identity();
    c1      *= 150.0 * 150.0;
    Mat3 c2  = Mat3::identity();
    c2      *= 200.0 * 200.0 - 150.0 * 150.0;
    Mat3 along;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            along(i, j) = 5000.0 * 5000.0 * z[i] * z[j];
        }
    }
    Mat3 c2Along  = c2;
    c2Along      += along;

    const CollisionProbability pc = collisionProbability(ca, c1, c2, 10.0);
    EXPECT_NEAR(pc.missDistance, 120.0, 1e-9);
    EXPECT_NEAR(pc.xx, 200.0 * 200.0, 1e-6);
    EXPECT_NEAR(pc.yy, 200.0 * 200.0, 1e-6);
    EXPECT_NEAR(pc.xy, 0.0, 1e-6);
    EXPECT_NEAR(pc.sigmaMajor, 200.0, 1e-9);
    EXPECT_NEAR(pc.sigmaMinor, 200.0, 1e-6);
    EXPECT_NEAR(pc.mahalanobis, 0.6, 1e-10); // positions near 7000 km: dr to ~1e-9 m
    const double ref = riceReference(120.0, 200.0, 10.0);
    EXPECT_NEAR(pc.probability, ref, 1e-10 * ref);
    EXPECT_NEAR(collisionProbability(ca, c1, c2Along, 10.0).probability, pc.probability, 1e-9 * pc.probability);

    // The frame: x along the miss, z along the relative velocity, right-handed.
    test::expectVecNear(pc.xAxis, normalized(m), 1e-10);
    test::expectVecNear(pc.zAxis, z, 1e-12);
    test::expectVecNear(cross(pc.xAxis, pc.yAxis), pc.zAxis, 1e-15);

    // The 6x6 overload reads the position blocks.
    Mat6 p1{}, p2{};
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            p1(i, j) = c1(i, j);
            p2(i, j) = c2Along(i, j);
        }
        p1(i + 3, i + 3) = 4.0; // velocity variances are ignored
        p2(i + 3, i + 3) = 9.0;
    }
    EXPECT_EQ(collisionProbability(ca, p1, p2, 10.0).probability,
              collisionProbability(ca, c1, c2Along, 10.0).probability);

    // A zero miss picks some direction in the plane.
    const CollisionProbability head = collisionProbability(approachOf(Vec3{}, dv), c1, c2, 10.0);
    EXPECT_NEAR(dot(head.xAxis, z), 0.0, 1e-15);
    EXPECT_NEAR(head.probability, -std::expm1(-0.5 * 0.05 * 0.05), 1e-12);
    EXPECT_EQ(head.mahalanobis, 0.0);
}

TEST(CollisionProbability, MatchesMonteCarloOfStraightLineEncounters) {
    // Draw relative positions from the 3-D combined covariance, move each along the relative
    // velocity, and count closest distances under R: the short-term model by brute force.
    const Vec3 dv{-4000.0, 10000.0, 1500.0};
    const Vec3 miss = 60.0 * normalized(cross(dv, Vec3{1.0, 0.0, 0.0}));
    // A correlated 3-D covariance through a Cholesky factor L (C = L L^T).
    // clang-format off
    const Mat3   L{ 80.0,  0.0,   0.0,
                    30.0, 50.0,   0.0,
                   -20.0, 10.0, 120.0};
    // clang-format on
    const Mat3   C = L * L.transpose();
    const Mat3   zero{};
    const double R  = 25.0;
    const double pc = collisionProbability(approachOf(miss, dv), C, zero, R).probability;

    std::mt19937                     rng(1992);
    std::normal_distribution<double> n01;
    const Vec3                       u    = normalized(dv);
    const int                        N    = 4000000;
    int                              hits = 0;
    for (int k = 0; k < N; ++k) {
        const Vec3 p   = miss + L * Vec3{n01(rng), n01(rng), n01(rng)};
        const Vec3 pp  = p - dot(p, u) * u; // closest point of the straight line
        hits          += norm2(pp) < R * R ? 1 : 0;
    }
    const double mc    = static_cast<double>(hits) / N;
    const double sigma = std::sqrt(pc * (1.0 - pc) / N);
    EXPECT_NEAR(mc, pc, 4.0 * sigma) << "Pc " << pc; // 4 sigma is 1.2 % of Pc
    EXPECT_GT(pc, 0.01);
}

TEST(CollisionProbability, NasaCaraReferenceCases) {
    // Omitron cases (CARA asserts these to RelTol 1e-10 in its unit test).
    EXPECT_NEAR(caraPc(kCara[0], kOmitronHbr), kOmitronPc[0], 1e-12 * kOmitronPc[0]);
    EXPECT_NEAR(caraPc(kCara[1], kOmitronHbr), kOmitronPc[1], 1e-12 * kOmitronPc[1]);

    // Omitron 2's states are 2.1 deg from TCA: projecting the miss (the default) gives a slightly
    // larger Pc than CARA's convention (0.017 %).
    CloseApproach raw;
    raw.first               = {kCara[1].r1, kCara[1].v1};
    raw.second              = {kCara[1].r2, kCara[1].v2};
    Mat3 c1                 = kCara[1].c1;
    c1                     *= kCara[1].c1Scale;
    const double projected  = collisionProbability(raw, c1, kCara[1].c2, kOmitronHbr).probability;
    EXPECT_NEAR(projected / kOmitronPc[1] - 1.0, 1.65e-4, 0.05e-4);

    // One geometry over a range of hard-body radii. Its plane covariance has sigmas of 123 km and
    // 18.6 m, with the miss 7 sigma along the narrow axis, so CARA's values are not exact:
    //  - at 10 m and 100 m they are 2.0e-8 and 1.9e-9 from an evaluation with exactly computed
    //    principal axes: the rounding of the minor variance formed as a difference,
    //    epsilon * (123 km / 18.6 m)^2 (collisionProbability forms it as determinant / major);
    //  - from 1 km up the default 64-point quadrature is off by 1-3e-4.
    // astro_study_conjunctions prints all of these against a brute-force grid.
    const double tolerance[] = {1e-7, 1e-8, 1e-3, 1e-3, 1e-3, 1e-3};
    for (int k = 0; k < 6; ++k) {
        EXPECT_NEAR(caraPc(kCara[2], kMultiRadiusHbr[k]), kMultiRadiusPc[k], tolerance[k] * kMultiRadiusPc[k])
            << kMultiRadiusHbr[k];
    }

    // Alfano (2009) cases 5-7 against Alfano's linear Pc, which CARA tests to 1e-3. These are slow
    // encounters (0.2-0.5 m/s) whose relative motion is nevertheless close to linear.
    for (int k = 0; k < 3; ++k) {
        EXPECT_NEAR(caraPc(kCara[3 + k], kAlfanoHbr), kAlfanoLinearPc[k], 1e-3 * kAlfanoLinearPc[k])
            << "Alfano case " << 5 + k;
    }
}

TEST(ProbabilityInCircle, NearSingularCovariances) {
    // A covariance far narrower than the circle, with the narrow axis's line at the circle's edge:
    // the integrand along the major axis is then a peak or a step much narrower than the range.
    // (Before the range was cut at those features, the first case returned 0 and the second took
    // minutes.) Reference: the integral taken the other way round.
    struct Case {
        double mx, my, sx, sy, R;
    };
    const Case cases[] = {
        {3.0,       10.00003,  0.5,  1e-5, 10.0}, // 3 sigma outside the edge: Pc 2.2e-13
        {3.0,       10.0003,   0.5,  1e-4, 10.0},
        {3.0,       10.0003,   10.0, 1e-4, 10.0},
        {1.0,       10.003,    1.0,  1e-3, 10.0},
        {3.0,       10.008,    0.5,  2e-3, 10.0},
        {3.0,       10.03,     0.5,  5e-3, 10.0}, // Pc 4.6e-18
        {3.0,       9.99997,   0.5,  1e-5, 10.0}, // 3 sigma inside the edge
        {3.0,       9.5,       0.5,  1e-5, 10.0}, // well inside: the chord is cut by the line
        {7.0710678, 7.0710678, 1e-3, 1e-5, 10.0}, // a point-like cloud on the edge, at 45 deg
    };
    for (const Case& c : cases) {
        const double ref = crossAxisReference(c.mx, c.my, c.sx, c.sy, c.R);
        ASSERT_GT(ref, 0.0);
        EXPECT_NEAR(probabilityInCircle(c.mx, c.my, c.sx * c.sx, 0.0, c.sy * c.sy, c.R), ref, 1e-9 * ref)
            << c.mx << " " << c.my << " " << c.sx << " " << c.sy;
    }
}

TEST(ProbabilityInCircle, RotatedElongatedCovariance) {
    // Standard deviations 6000 and 0.8 at an angle. The entries are integers, so the determinant,
    // and with it the minor variance, is exact in double arithmetic here. Formed as
    // mean - sqrt(...) it is off by 7e-10, which is 1e-9 to 1e-8 in a probability this far out along
    // the minor axis. The reference takes the two standard deviations directly.
    const double xx = 12960001.0, yy = 23040000.0, xy = 17280000.0;
    const double major = 0.5 * (xx + yy) + std::hypot(0.5 * (xx - yy), xy);
    const double minor = (xx * yy - xy * xy) / major; // both products are below 2^53: exact
    EXPECT_NEAR(minor, 0.64, 1e-6);
    const double angle = 0.5 * std::atan2(2.0 * xy, xx - yy);
    const double c = std::cos(angle), s = std::sin(angle);
    for (double v : {4.0, 4.8, 6.4, 8.0}) {                           // offset along the minor axis: 5 to 10 sigma
        const double mx = -s * v + c * 100.0, my = c * v + s * 100.0; // and 100 along the major axis
        const double ref = crossAxisReference(100.0, v, std::sqrt(major), std::sqrt(minor), 3.0);
        ASSERT_GT(ref, 0.0);
        EXPECT_NEAR(probabilityInCircle(mx, my, xx, xy, yy, 3.0), ref, 1e-10 * ref) << v;
        // The same along its own axes.
        EXPECT_NEAR(probabilityInCircle(100.0, v, major, 0.0, minor, 3.0), ref, 1e-10 * ref) << v;
    }
}

TEST(CollisionProbability, Errors) {
    const CloseApproach ca = approachOf(Vec3{100.0, 0.0, 0.0}, Vec3{0.0, 1000.0, 0.0});
    const Mat3          C  = Mat3::identity();
    EXPECT_THROW((void)collisionProbability(ca, C, C, 0.0), Error);
    EXPECT_THROW((void)collisionProbability(ca, C, C, -1.0), Error);
    EXPECT_THROW((void)collisionProbability(approachOf(Vec3{100.0, 0.0, 0.0}, Vec3{}), C, C, 1.0), Error);
    Mat3 bad  = Mat3::identity();
    bad(0, 0) = -5.0;
    EXPECT_THROW((void)collisionProbability(ca, bad, C, 1.0), Error);
    Mat3 nan  = Mat3::identity();
    nan(2, 2) = std::nan("");
    EXPECT_THROW((void)collisionProbability(ca, nan, C, 1.0), Error);
}
