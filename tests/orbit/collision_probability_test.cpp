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
#include "test_util.h"

using namespace astro;

namespace {

/// Brute-force reference: the 2-D Gaussian integrated over the disc in polar coordinates, Simpson
/// in r and the (spectrally accurate for periodic integrands) midpoint rule in theta.
double polarReference(double mx, double my, double xx, double xy, double yy, double R, int nr = 2000, int nt = 2000) {
    const double det = xx * yy - xy * xy;
    const double ixx = yy / det, ixy = -xy / det, iyy = xx / det;
    const double k    = 1.0 / (constants::twoPi * std::sqrt(det));
    const auto   ring = [&](double r) {
        double s = 0.0;
        for (int j = 0; j < nt; ++j) {
            const double t  = constants::twoPi * (j + 0.5) / nt;
            const double dx = r * std::cos(t) - mx, dy = r * std::sin(t) - my;
            s += std::exp(-0.5 * (ixx * dx * dx + 2.0 * ixy * dx * dy + iyy * dy * dy));
        }
        return k * r * s * constants::twoPi / nt;
    };
    const double h   = R / nr;
    double       sum = ring(0.0) + ring(R);
    for (int i = 1; i < nr; ++i) {
        sum += (i % 2 ? 4.0 : 2.0) * ring(i * h);
    }
    return sum * h / 3.0;
}

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
        {800.0,  -300.0, 100.0,  60.0,  0.2,  15.0}, // far: probability ~ 1e-14
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
            p1[i][j] = c1(i, j);
            p2[i][j] = c2Along(i, j);
        }
        p1[i + 3][i + 3] = 4.0; // velocity variances are ignored
        p2[i + 3][i + 3] = 9.0;
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
    const int                        N    = 400000;
    int                              hits = 0;
    for (int k = 0; k < N; ++k) {
        const Vec3 p   = miss + L * Vec3{n01(rng), n01(rng), n01(rng)};
        const Vec3 pp  = p - dot(p, u) * u; // closest point of the straight line
        hits          += norm2(pp) < R * R ? 1 : 0;
    }
    const double mc    = static_cast<double>(hits) / N;
    const double sigma = std::sqrt(pc * (1.0 - pc) / N);
    EXPECT_NEAR(mc, pc, 5.0 * sigma) << "Pc " << pc;
    EXPECT_GT(pc, 0.01);
}

namespace {

/// A conjunction from NASA CARA's test data: both states and position covariances at TCA.
struct CaraGeometry {
    Vec3   r1, v1, r2, v2;
    Mat3   c1, c2;
    double c1Scale; // factor applied to c1, as in the CARA test
};

// Reference inputs from NASA CARA Analysis Tools (github.com/nasa/CARA_Analysis_Tools, commit
// 43f0160): PcCircle_UnitTest.m and the Alfano (2009) test cases. See THIRD_PARTY_NOTICES.md.
// clang-format off
const CaraGeometry kCara[] = {
    // Omitron 1 (km)
    {{   378.39559,  4305.721887,  5752.767554}, { 2.360800244,  5.580331936, -4.322349039},
     { 374.5180598,  4307.560983,  5751.130418}, {-5.388125081, -3.946827739,  3.322820358},
     Mat3{ 44.5757544811362,  81.6751751052616, -67.8687662707124,
           81.6751751052616,  158.453402956163, -128.616921644857,
          -67.8687662707124, -128.616921644858,  105.490542562701},
     Mat3{ 2.31067077720423,  1.69905293875632,  -1.4170164577661,
           1.69905293875632,  1.24957388457206, -1.04174164279599,
           -1.4170164577661, -1.04174164279599, 0.869260558223714},
     1.0},
    // Omitron 2 (km; C1 scaled by 1e-3 as in the CARA test)
    {{-3239.128337196251,  2404.575152356222,  5703.228541709001}, {-3.745768373154199,  5.012339015927846, -4.231864565717194},
     {-3239.138264917246,  2404.568320465936,  5703.235605231182}, { 6.110192790100711,  -1.76732140789483,  4.140369261741708},
     Mat3{ 0.342072996423899, -0.412677096778269,  0.371500417511149,
          -0.412677096778269,  0.609905946319294, -0.540401385544286,
           0.371500417511149, -0.540401385544286,  0.521238634755377},
     Mat3{ 0.028351300975134, -0.008204103437377,  0.019253747960155,
          -0.008204103437377,  0.002404377774847, -0.005586512197914,
           0.019253747960155, -0.005586512197914,  0.013289250260317},
     0.001},
    // Multiple radii (m)
    {{ -6227006.55254931,  -3000839.63874701, -193256.9099062626}, { 2983.983865443772, -5973.188938878403, -3622.262843716804},
     {-6231009.268962734, -2993196.148130401, -176140.1651747054}, {-325.5968208692156,  994.6514541875586, -7507.340539210305},
     Mat3{ 23292.48149898317, -46240.97512511387, -28103.52472705708,
          -46240.97512511387,  91893.73058443272,   55838.5828667566,
          -28103.52472705708,   55838.5828667566,  33965.06000245927},
     Mat3{ 36891498.93147051, -120066974.8208798,   900253524.718121,
          -120066974.8208798,  390773572.2225543, -2929990982.074453,
            900253524.718121, -2929990982.074453,  21968889281.97153},
     1.0},
    // Alfano (2009) case 5 (m)
    {{ 6878090.1622937,   -17948.6785967,   -17948.6785967}, {  28.09377718079,  5382.8902061517,  5382.8902061517},
     { 6878089.1619754, -17946.678929416, -17947.678286012}, {  28.39378131578,  5383.1902163515,  5382.5902081838},
     Mat3{0.06420452636324, -18.990691938367, -18.990691938367,
          -18.990691938367,  7904.0447002725,  7903.9662413044,
          -18.990691938367,  7903.9662413044,  7904.0447002725},
     Mat3{0.062112348981822,  -18.552069769383,  -18.550001998859,
           -18.552069769383,   7905.3548684174,   7904.3953083026,
           -18.550001998859,   7904.3953083026,   7903.5927814582},
     1.0},
    // Alfano (2009) case 6 (m)
    {{ 6877715.634247, 53834.214039436, 53834.214039436}, {-84.26282780536, 5382.5970952289, 5382.5970952289},
     {6877716.6344015, 53835.214140512, 53836.213772982}, {-84.16281513631, 5382.6970892007, 5382.4970647006},
     Mat3{ 430.02827939777, -18393.891381385, -18393.891381385,
          -18393.891381385,  790190.51867865,   790186.9676547,
          -18393.891381385,   790186.9676547,  790190.51867865},
     Mat3{ 429.35284572627, -18379.895353487, -18379.212284332,
          -18379.895353487,  790234.93895821,  790202.01968289,
          -18379.212284332,  790202.01968289,  790176.20377729},
     1.0},
    // Alfano (2009) case 7 (m)
    {{ -6877469.169824,  -67773.18379991,  -67773.18379991}, { 106.08062887158, -5382.4047050208, -5382.4047050208},
     {-6877468.6856734, -67775.155096794,  -67775.63603513}, { 106.21663006962, -5382.5006947728, -5382.3006891309},
     Mat3{ 47752.931002152, -1772101.1180686, -1772101.1180686,
          -1772101.1180686,  65779308.882914,  65778953.792515,
          -1772101.1180686,  65778953.792515,  65779308.882914},
     Mat3{ 47842.556250334,  -1773779.529826, -1773713.6133993,
           -1773779.529826,  65780452.931669,  65777653.286419,
          -1773713.6133993,  65777653.286419,  65775564.028398},
     1.0},
};
// clang-format on

/// CARA's PcCircle places the miss distance |r2 - r1| in the encounter plane (the states are
/// taken to be exactly at TCA); collisionProbability uses the projection of r2 - r1 onto the
/// plane, the closest distance of the straight-line encounter. To compare the integrals, give it
/// the full distance along the projected direction.
double caraPc(const CaraGeometry& g, double hbr) {
    const Vec3    dr = g.r2 - g.r1, z = normalized(g.v2 - g.v1);
    CloseApproach c;
    c.first   = {Vec3{}, g.v1};
    c.second  = {norm(dr) * normalized(dr - dot(dr, z) * z), g.v2};
    Mat3 c1   = g.c1;
    c1       *= g.c1Scale;
    return collisionProbability(c, c1, g.c2, hbr).probability;
}

} // namespace

TEST(CollisionProbability, NasaCaraReferenceCases) {
    // Omitron cases: CARA's highest-accuracy values (RelTol 1e-10 in its test). HBR 20 m.
    EXPECT_NEAR(caraPc(kCara[0], 0.020), 2.706023476569787e-05, 1e-12 * 2.706023476569787e-05);
    EXPECT_NEAR(caraPc(kCara[1], 0.020), 1.807363058494765e-01, 1e-12 * 1.807363058494765e-01);

    // Omitron 2's states are 2.1 deg from TCA: projecting the miss (the default) gives a slightly
    // larger Pc than CARA's convention.
    CloseApproach raw;
    raw.first               = {kCara[1].r1, kCara[1].v1};
    raw.second              = {kCara[1].r2, kCara[1].v2};
    Mat3 c1                 = kCara[1].c1;
    c1                     *= kCara[1].c1Scale;
    const double projected  = collisionProbability(raw, c1, kCara[1].c2, 0.020).probability;
    EXPECT_GT(projected, 1.807363058494765e-01);
    EXPECT_LT(projected, 1.807363058494765e-01 * 1.001);

    // One geometry over a range of hard-body radii (m). CARA's values come from its default
    // 64-point Gauss-Chebyshev quadrature. The plane covariance here has sigmas of 123 km and 19 m,
    // and from 1 km up that quadrature is off by 1-3e-4 (a brute-force grid agrees with
    // collisionProbability to 1e-13): compare to CARA's own stated accuracy, 1e-3, there.
    // clang-format off
    const double hbr[]      = {1e1,                   1e2,                   1e3,                   3e3,                   1e4,                   3e4};
    const double expected[] = {2.514939701718358e-15, 1.232312816404174e-05, 6.343318773183237e-03, 1.917884657463699e-02, 6.389821634646267e-02, 1.900663422510684e-01};
    // clang-format on
    for (int k = 0; k < 6; ++k) {
        const double tol = (hbr[k] < 1e3 ? 1e-12 : 1e-3) * expected[k];
        EXPECT_NEAR(caraPc(kCara[2], hbr[k]), expected[k], tol) << hbr[k];
    }

    // Alfano (2009) cases 5-7, HBR 10 m, against Alfano's linear Pc, which CARA tests to 1e-3.
    // These are slow encounters (0.2-0.5 m/s) whose relative motion is nevertheless close to linear:
    // Alfano's 1e8-sample Monte Carlo gives 0.044499, 0.0043005 and 0.000161462.
    const double alfano[] = {0.044487386, 0.004335455, 0.000158147};
    for (int k = 0; k < 3; ++k) {
        EXPECT_NEAR(caraPc(kCara[3 + k], 10.0), alfano[k], 1e-3 * alfano[k]) << "Alfano case " << 5 + k;
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
