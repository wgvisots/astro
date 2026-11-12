/**********************************************************************
 * conjunction_reference.h
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

#pragma once

// Reference conjunctions and a brute-force Pc, shared by collision_probability_test.cpp and the
// conjunction study (studies/conjunctions.cpp).

#include <cmath>

#include "astro/core/constants.h"
#include "astro/orbit/close_approach.h"

namespace astro::test {

/// Brute-force reference: the 2-D Gaussian integrated over the disc in polar coordinates, Simpson
/// in r and the (spectrally accurate for periodic integrands) midpoint rule in theta.
inline double
polarReference(double mx, double my, double xx, double xy, double yy, double R, int nr = 2000, int nt = 2000) {
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

/// A conjunction from NASA CARA's test data: both states and position covariances at TCA.
struct CaraGeometry {
    Vec3   r1, v1, r2, v2;
    Mat3   c1, c2;
    double c1Scale; // factor applied to c1, as in the CARA test
};

// Reference inputs from NASA CARA Analysis Tools (github.com/nasa/CARA_Analysis_Tools, commit
// 43f0160): PcCircle_UnitTest.m and the Alfano (2009) test cases. See THIRD_PARTY_NOTICES.md.
// clang-format off
inline const CaraGeometry kCara[] = {
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
inline double caraPc(const CaraGeometry& g, double hbr) {
    const Vec3    dr = g.r2 - g.r1, z = normalized(g.v2 - g.v1);
    CloseApproach c;
    c.first   = {Vec3{}, g.v1};
    c.second  = {norm(dr) * normalized(dr - dot(dr, z) * z), g.v2};
    Mat3 c1   = g.c1;
    c1       *= g.c1Scale;
    return collisionProbability(c, c1, g.c2, hbr).probability;
}

// Reference probabilities for kCara. Omitron 1 and 2 (HBR 0.020 km = 20 m, in their km units):
// CARA's highest-accuracy values (RelTol 1e-10 in its test). Multiple radii (m): CARA's default
// 64-point Gauss-Chebyshev values. Alfano (2009) cases 5-7 (HBR 10 m): Alfano's linear Pc and his
// 1e8-sample Monte Carlo.
// clang-format off
inline constexpr double kOmitronPc[]          = {2.706023476569787e-05, 1.807363058494765e-01};
inline constexpr double kOmitronHbr           = 0.020;
inline constexpr double kMultiRadiusHbr[]     = {1e1,                   1e2,                   1e3,                   3e3,                   1e4,                   3e4};
inline constexpr double kMultiRadiusPc[]      = {2.514939701718358e-15, 1.232312816404174e-05, 6.343318773183237e-03, 1.917884657463699e-02, 6.389821634646267e-02, 1.900663422510684e-01};
inline constexpr double kAlfanoHbr            = 10.0;
inline constexpr double kAlfanoLinearPc[]     = {0.044487386, 0.004335455, 0.000158147};
inline constexpr double kAlfanoMonteCarloPc[] = {0.044498913, 0.0043005,   0.000161462};
// clang-format on

} // namespace astro::test
