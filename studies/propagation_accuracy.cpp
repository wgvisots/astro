/**********************************************************************
 * propagation_accuracy.cpp
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

// Propagation accuracy study (docs/accuracy.md):
//  1. numerical J2/J4: integrator tolerance vs steps and worst position error over 4 h and 1 day;
//  2. grid interpolation adds nothing visible: grid error vs the integrator's own error at its steps;
//  3. analytic secular J2/J4 from a state vector, and with least-squares-fitted mean elements,
//     with the original and the corrected (d position / dt) velocity;
//  4. the J2 axis: GCRS z instead of the true pole, and TOD/TEME frames frozen at the epoch;
//  5. state transition matrices and covariances of the analytic propagators: two-body (closed
//     form) and secular (differenced) against finite differences, the integrated variational
//     equations, and the numerical J2/J4 covariance.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

#include "astro/core/angle.h"
#include "astro/frames/eop.h"
#include "astro/frames/transform.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/gravity.h"
#include "astro/orbit/numerical.h"
#include "astro/orbit/secular.h"
#include "astro/orbit/two_body.h"
#include "astro/satellite/estimate_tracker.h"
#include "astro/satellite/orbital_frames.h"

using namespace astro;

namespace {

const Time kEpoch = Time::fromIso8601("2026-03-01T00:00:00Z");

struct Orbit {
    const char*       name;
    KeplerianElements el;
};

const Orbit kOrbits[] = {
    {"LEO 400 km i=51.6",     {6778e3, 0.0005, deg2rad(51.6), 0.3, 0.5, 0.7}},
    {"SSO 800 km e=0.01",     {7178e3, 0.01, deg2rad(98.6), 1.3, 2.5, 0.2}  },
    {"HEO 1000 x 20000 km",   {16878e3, 0.56, deg2rad(63.4), 0.1, 4.7, 0.0} },
    {"near-equatorial i=0.5", {7000e3, 0.002, deg2rad(0.5), 0.0, 0.0, 1.0}  },
    {"polar i=90",            {7000e3, 0.001, deg2rad(90.0), 2.0, 0.3, 0.1} },
    {"retrograde i=140",      {7000e3, 0.001, deg2rad(140.0), 0.5, 1.0, 2.0}},
};

/// Reference ("truth") propagator: J2 + J4 with a tight tolerance.
TrackerOptions tightOptions() {
    TrackerOptions o;
    o.integrator.relTol  = 1e-14;
    o.integrator.absTol  = 1e-11;
    o.integrator.maxStep = 20.0;
    return o;
}

[[maybe_unused]] Mat6 diagonal(double positionSigma, double velocitySigma) {
    Mat6 P{};
    for (std::size_t i = 0; i < 3; ++i) {
        P[i][i]         = positionSigma * positionSigma;
        P[i + 3][i + 3] = velocitySigma * velocitySigma;
    }
    return P;
}

} // namespace

namespace {

double maxDistance(const std::vector<StateVector>& a, const std::vector<StateVector>& b) {
    double m = 0.0;
    for (std::size_t k = 0; k < a.size(); ++k) {
        m = std::max(m, norm(a[k].r - b[k].r));
    }
    return m;
}

double maxVelocityDifference(const std::vector<StateVector>& a, const std::vector<StateVector>& b) {
    double m = 0.0;
    for (std::size_t k = 0; k < a.size(); ++k) {
        m = std::max(m, norm(a[k].v - b[k].v));
    }
    return m;
}

/// Worst |v - dr/dt| over a grid, with dr/dt from central differences of the grid's own positions.
double maxRateInconsistency(const std::vector<StateVector>& g, double step) {
    double m = 0.0;
    for (std::size_t k = 1; k + 1 < g.size(); ++k) {
        m = std::max(m, norm(g[k].v - (g[k + 1].r - g[k - 1].r) / (2.0 * step)));
    }
    return m;
}

std::vector<StateVector> truthGrid(const StateVector& s0, double step, std::size_t n) {
    return makePropagator({kEpoch, s0, Mat6{}}, tightOptions()).propagateGrid(0.0, step, n);
}

void toleranceStudy() {
    std::printf("1. Numerical J2/J4: tolerance vs worst position error (first three orbits), output every 10 s\n\n");
    std::printf("%-8s %-8s %12s %16s %16s\n", "relTol", "absTol", "steps / 4 h", "max err, 4 h", "max err, 1 day");
    std::vector<std::vector<StateVector>> truth4, truthDay;
    for (int i = 0; i < 3; ++i) {
        const StateVector s0 = elementsToState(kOrbits[i].el);
        truth4.push_back(truthGrid(s0, 10.0, 1441));
        truthDay.push_back(truthGrid(s0, 60.0, 1441));
    }
    for (double rt : {1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7}) {
        TrackerOptions o;
        o.integrator.relTol = rt;
        o.integrator.absTol = rt * 1e3;
        double        e4 = 0.0, eDay = 0.0;
        AdaptiveStats st;
        for (int i = 0; i < 3; ++i) {
            const auto    prop = makePropagator({kEpoch, elementsToState(kOrbits[i].el), Mat6{}}, o);
            AdaptiveStats s;
            e4   = std::max(e4, maxDistance(prop.propagateGrid(0.0, 10.0, 1441, &s), truth4[i]));
            eDay = std::max(eDay, maxDistance(prop.propagateGrid(0.0, 60.0, 1441), truthDay[i]));
            if (i == 0) {
                st = s;
            }
        }
        std::printf("%-8.0e %-8.0e %12d %14.3g m %14.3g m\n", rt, rt * 1e3, st.accepted, e4, eDay);
    }
}

void interpolationStudy() {
    std::printf("\n2. Grid interpolation vs the integrator's own error (LEO, default tolerance, 1 day)\n\n");
    const StateVector   s0    = elementsToState(kOrbits[0].el);
    const auto          prop  = makePropagator({kEpoch, s0, Mat6{}}, [] {
        TrackerOptions o;
        o.integrator = AdaptiveOptions{};
        return o;
    }());
    const auto          truth = truthGrid(s0, 60.0, 1441);
    std::vector<double> dts(1441);
    for (std::size_t k = 0; k < dts.size(); ++k) {
        dts[k] = 60.0 * static_cast<double>(k);
    }
    const double gridErr = maxDistance(prop.propagateGrid(0.0, 60.0, 1441), truth);
    const double stopErr = maxDistance(prop.propagate(dts), truth); // stopping on every output
    const double endErr  = norm(prop.propagate(86400.0).r - truth.back().r);
    std::printf("  grid (natural ~96 s steps, interpolated)  %.2e m\n", gridErr);
    std::printf("  integrator alone at the end of the day    %.2e m\n", endErr);
    std::printf("  stopping at every 60 s output (smaller steps) %.2e m\n", stopErr);
}

using P6 = std::array<double, 6>;

KeplerianElements fromEquinoctialParams(const P6& x) {
    const double e = std::hypot(x[1], x[2]), lp = std::atan2(x[1], x[2]);
    const double t = std::hypot(x[3], x[4]), O = std::atan2(x[3], x[4]);
    return elementsFromMeanAnomaly(x[0], e, 2 * std::atan(t), wrapTwoPi(O), wrapTwoPi(lp - O), wrapTwoPi(x[5] - lp));
}

/// Gauss-Newton fit of SecularPropagator mean elements (non-singular parameters) to positions.
P6 fitMeanElements(const std::vector<StateVector>& obs, double step, P6 x, double& rms) {
    const double h[6] = {1.0, 1e-7, 1e-7, 1e-7, 1e-7, 1e-7};
    for (int it = 0; it < 15; ++it) {
        const auto m0 = SecularPropagator(fromEquinoctialParams(x)).propagateGrid(0.0, step, obs.size());
        std::vector<std::array<Vec3, 6>> J(obs.size());
        for (int j = 0; j < 6; ++j) {
            P6 xp          = x;
            xp[j]         += h[j];
            const auto mp  = SecularPropagator(fromEquinoctialParams(xp)).propagateGrid(0.0, step, obs.size());
            for (std::size_t k = 0; k < obs.size(); ++k) {
                J[k][j] = (mp[k].r - m0[k].r) / h[j];
            }
        }
        double N[6][7] = {}, ss = 0.0;
        for (std::size_t k = 0; k < obs.size(); ++k) {
            const Vec3 res  = obs[k].r - m0[k].r;
            ss             += norm2(res);
            for (int a = 0; a < 6; ++a) {
                for (int b = 0; b < 6; ++b) {
                    N[a][b] += dot(J[k][a], J[k][b]);
                }
                N[a][6] += dot(J[k][a], res);
            }
        }
        rms = std::sqrt(ss / static_cast<double>(obs.size()));
        for (int c = 0; c < 6; ++c) {
            int p = c;
            for (int r = c + 1; r < 6; ++r) {
                if (std::abs(N[r][c]) > std::abs(N[p][c])) {
                    p = r;
                }
            }
            for (int j = 0; j < 7; ++j) {
                std::swap(N[c][j], N[p][j]);
            }
            for (int r = 0; r < 6; ++r) {
                if (r != c) {
                    const double f = N[r][c] / N[c][c];
                    for (int j = c; j < 7; ++j) {
                        N[r][j] -= f * N[c][j];
                    }
                }
            }
        }
        for (int j = 0; j < 6; ++j) {
            x[j] += N[j][6] / N[j][j];
        }
    }
    return x;
}

void secularStudy() {
    std::printf("\n3. Analytic secular J2/J4 vs numerical truth (TOD-like frame, z = pole), 4 h\n\n");
    std::printf("%-24s %16s %22s %18s\n", "orbit", "from state", "fitted (2 orbits)", "fit RMS");
    struct VelocityRow {
        double directBase, directCorrected, fittedBase, fittedCorrected, rateBase, rateCorrected;
    };
    std::vector<VelocityRow> velocity;
    for (int i = 0; i < 3; ++i) {
        const StateVector s0           = elementsToState(kOrbits[i].el);
        TrackerOptions    tod          = tightOptions();
        tod.frame                      = Frame::TOD; // the secular theory needs z = pole
        const auto              truth  = makePropagator({kEpoch, s0, Mat6{}}, tod).propagateGrid(0.0, 30.0, 481);
        const auto              direct = SecularPropagator(s0).propagateGrid(0.0, 30.0, 481);
        const double            period = orbitalPeriod(stateToElements(s0).a);
        const std::size_t       nFit   = std::min<std::size_t>(481, static_cast<std::size_t>(2.0 * period / 30.0) + 1);
        const KeplerianElements k      = stateToElements(s0);
        const double            lp = k.raan + k.argp, t = std::tan(k.i / 2);
        double                  rms = 0.0;
        const P6                x = fitMeanElements({truth.begin(), truth.begin() + static_cast<std::ptrdiff_t>(nFit)},
                                     30.0,
                                                    {k.a,
                                                     k.e * std::sin(lp),
                                                     k.e * std::cos(lp),
                                                     t * std::sin(k.raan),
                                                     t * std::cos(k.raan),
                                                     lp + k.meanAnomaly()},
                                     rms);
        const auto              fitted = SecularPropagator(fromEquinoctialParams(x)).propagateGrid(0.0, 30.0, 481);
        std::printf("%-24s %13.1f km %19.1f km %15.1f km\n",
                    kOrbits[i].name,
                    maxDistance(direct, truth) / 1e3,
                    maxDistance(fitted, truth) / 1e3,
                    rms / 1e3);
        // The same positions with the corrected velocity (d position / dt). Consistency on a 1 s grid.
        const auto directC = SecularPropagatorCorrected(s0).propagateGrid(0.0, 30.0, 481);
        const auto fittedC = SecularPropagatorCorrected(fromEquinoctialParams(x)).propagateGrid(0.0, 30.0, 481);
        velocity.push_back(
            {maxVelocityDifference(direct, truth),
             maxVelocityDifference(directC, truth),
             maxVelocityDifference(fitted, truth),
             maxVelocityDifference(fittedC, truth),
             maxRateInconsistency(SecularPropagator(fromEquinoctialParams(x)).propagateGrid(0.0, 1.0, 14401), 1.0),
             maxRateInconsistency(SecularPropagatorCorrected(fromEquinoctialParams(x)).propagateGrid(0.0, 1.0, 14401),
                                  1.0)});
    }
    std::printf("\n   Worst velocity error vs truth over 4 h (m/s), and worst |v - dr/dt| of the secular output\n\n");
    std::printf("%-24s %22s %22s %22s\n", "", "from state", "fitted (2 orbits)", "|v - dr/dt|, fitted");
    std::printf("%-24s %10s %11s %10s %11s %10s %11s\n",
                "orbit",
                "original",
                "corrected",
                "original",
                "corrected",
                "original",
                "corrected");
    for (int i = 0; i < 3; ++i) {
        const VelocityRow& v = velocity[static_cast<std::size_t>(i)];
        std::printf("%-24s %10.2f %11.2f %10.2f %11.2f %10.3f %11.3f\n",
                    kOrbits[i].name,
                    v.directBase,
                    v.directCorrected,
                    v.fittedBase,
                    v.fittedCorrected,
                    v.rateBase,
                    v.rateCorrected);
    }
}

/// Mean elements whose corrected (d position / dt) velocity at the epoch equals the state's velocity:
/// v_mean = (v - raanDot z x r - argpDot w x r) / (Mdot / n), iterated with the rates.
KeplerianElements matchCorrectedVelocity(const StateVector& s0) {
    KeplerianElements el = stateToElements(s0);
    for (int it = 0; it < 8; ++it) {
        const SecularRates rt = secularRates(el.a, el.e, el.i);
        const double       n  = std::sqrt(constants::gmEarth / (el.a * el.a * el.a));
        const Vec3         w  = normalized(cross(s0.r, s0.v));
        const Vec3         vk =
            (s0.v - rt.raanDot * cross(Vec3{0, 0, 1}, s0.r) - rt.argpDot * cross(w, s0.r)) / (rt.meanAnomalyDot / n);
        el = stateToElements({s0.r, vk});
    }
    return el;
}

void secularStartStudy() {
    std::printf("\n3b. Secular from a state: osculating elements as mean, or mean elements whose corrected\n");
    std::printf("    velocity reproduces the state's velocity at the epoch. Worst position error over 4 h\n\n");
    std::printf("%-24s %16s %22s\n", "orbit", "osculating", "velocity-matched");
    for (const Orbit& o : kOrbits) {
        const StateVector s0  = elementsToState(o.el);
        TrackerOptions    tod = tightOptions();
        tod.frame             = Frame::TOD;
        const auto truth      = makePropagator({kEpoch, s0, Mat6{}}, tod).propagateGrid(0.0, 30.0, 481);
        const auto osc        = SecularPropagator(s0).propagateGrid(0.0, 30.0, 481);
        const auto matched    = SecularPropagatorCorrected(matchCorrectedVelocity(s0)).propagateGrid(0.0, 30.0, 481);
        std::printf("%-24s %13.1f km %19.1f km\n",
                    o.name,
                    maxDistance(osc, truth) / 1e3,
                    maxDistance(matched, truth) / 1e3);
    }
}

void poleStudy() {
    std::printf("\n4. The J2 axis (LEO): J2/J4 about GCRS z instead of the true pole, and of-date frames frozen\n\n");
    const Time        t0            = Time::fromIso8601("2026-11-02T00:00:00Z");
    const StateVector s0            = elementsToState(kOrbits[0].el);
    TrackerOptions    o             = tightOptions();
    o.integrator.maxStep            = 1e30;
    const auto                right = makePropagator({t0, s0, Mat6{}}, o).propagateGrid(0.0, 60.0, 241); // pole aligned
    const std::vector<double> J{0, 0, constants::J2, 0, constants::J4};
    const auto                wrong =
        NumericalPropagator(t0, s0, {std::make_shared<GravityField>(J)}, o.integrator).propagateGrid(0.0, 60.0, 241);
    const Vec3 pole = poleOfDateMatrix(Frame::GCRS, t0).transpose() * Vec3{0, 0, 1};
    std::printf("  true pole vs GCRS z on 2026-11-02: %.4f deg\n", rad2deg(std::acos(pole.z)));
    std::printf("  J2/J4 about GCRS z: %.0f m after 1 h, %.0f m after 4 h\n",
                maxDistance({wrong.begin(), wrong.begin() + 61}, {right.begin(), right.begin() + 61}),
                maxDistance(wrong, right));
    const ConstantEop      eop;
    const FrameTransformer f0(t0, eop);
    for (double dt : {3600.0, 14400.0, 86400.0}) {
        const FrameTransformer f1(t0 + dt, eop);
        const Vec3             r{6778e3, 0, 0};
        std::printf("  TOD of the epoch vs TOD of date after %5.0f h: %.2f m (TEME %.2f m)\n",
                    dt / 3600,
                    norm(f1.transform(f0.transform(r, Frame::TOD, Frame::GCRS), Frame::GCRS, Frame::TOD) - r),
                    norm(f1.transform(f0.transform(r, Frame::TEME, Frame::GCRS), Frame::GCRS, Frame::TEME) - r));
    }
}

/// Largest difference relative to the largest element of the same 3x3 block of b, over all blocks.
double maxBlockRelDiff(const Mat6& a, const Mat6& b) {
    double worst = 0.0;
    for (std::size_t br = 0; br < 2; ++br) {
        for (std::size_t bc = 0; bc < 2; ++bc) {
            double num = 0.0, den = 0.0;
            for (std::size_t i = 0; i < 3; ++i) {
                for (std::size_t j = 0; j < 3; ++j) {
                    num = std::max(num, std::abs(a[3 * br + i][3 * bc + j] - b[3 * br + i][3 * bc + j]));
                    den = std::max(den, std::abs(b[3 * br + i][3 * bc + j]));
                }
            }
            worst = std::max(worst, num / den);
        }
    }
    return worst;
}

/// d(propagateKepler(s0, dt)) / d(s0) by fourth-order central differences (steps 1e-6 |r0|, |v0|).
Mat6 keplerDifferences(const StateVector& s0, double dt) {
    constexpr double offset[4] = {1.0, -1.0, 2.0, -2.0};
    constexpr double weight[4] = {8.0, -8.0, -1.0, 1.0};
    Mat6             m{};
    for (std::size_t j = 0; j < 6; ++j) {
        const double h = 1e-6 * (j < 3 ? norm(s0.r) : norm(s0.v));
        for (std::size_t p = 0; p < 4; ++p) {
            StateVector s                     = s0;
            (j < 3 ? &s.r.x : &s.v.x)[j % 3] += offset[p] * h;
            const StateVector o               = propagateKepler(s, dt);
            const double      c[6]            = {o.r.x, o.r.y, o.r.z, o.v.x, o.v.y, o.v.z};
            for (std::size_t i = 0; i < 6; ++i) {
                m[i][j] += weight[p] * c[i] / (12.0 * h);
            }
        }
    }
    return m;
}

Vec3 ricSigmas(const StateCovariance& s) {
    const Mat6 ric = rotateCovariance(s.covariance, ricMatrix(s.state));
    return {std::sqrt(ric[0][0]), std::sqrt(ric[1][1]), std::sqrt(ric[2][2])};
}

void covarianceStudy() {
    std::printf("\n5. Analytic state transition matrices (block-relative differences) and covariances\n\n");
    std::printf("5a. Two-body, closed form (propagateKeplerStm)\n\n");
    std::printf("| Orbit | dt | vs differences | vs point-mass variational equations | max abs(Phi^T J Phi - J) |\n");
    std::printf("|---|---|---|---|---|\n");
    struct Case {
        const char* name;
        StateVector s0;
        double      dt;
    };
    const Case twoBody[] = {
        {"LEO 400 km",          elementsToState(kOrbits[0].el),                                  86400.0       },
        {"LEO 400 km",          elementsToState(kOrbits[0].el),                                  7 * 86400.0   },
        {"HEO 1000 x 20000 km", elementsToState(kOrbits[2].el),                                  86400.0       },
        {"e = 0.9",             elementsToState({60000e3, 0.9, deg2rad(30.0), 1.0, 2.0, 0.1}),   1.5 * 146000.0},
        {"hyperbola e = 1.5",   elementsToState({-20000e3, 1.5, deg2rad(30.0), 1.0, 2.0, -0.5}), 20000.0       },
    };
    AdaptiveOptions tight;
    tight.relTol = 1e-14;
    tight.absTol = 1e-11;
    Mat6 J{};
    for (std::size_t i = 0; i < 3; ++i) {
        J[i][i + 3] = 1.0;
        J[i + 3][i] = -1.0;
    }
    for (const Case& c : twoBody) {
        const Mat6                phi = propagateKeplerStm(c.s0, c.dt).stm;
        const NumericalPropagator mass(Time(), c.s0, {std::make_shared<GravityField>(std::vector<double>{})}, tight);
        std::vector<StateVector>  st(1);
        std::vector<Mat6>         var(1);
        mass.propagateGridStm(c.dt, 1.0, st, var);
        const Mat6 sym   = multiply(transpose(phi), multiply(J, phi));
        double     worst = 0.0;
        for (std::size_t i = 0; i < 6; ++i) {
            for (std::size_t j = 0; j < 6; ++j) {
                worst = std::max(worst, std::abs(sym[i][j] - J[i][j]));
            }
        }
        std::printf("| %s | %.1f d | %.1e | %.1e | %.1e |\n",
                    c.name,
                    c.dt / 86400.0,
                    maxBlockRelDiff(phi, keplerDifferences(c.s0, c.dt)),
                    maxBlockRelDiff(phi, var[0]),
                    worst);
    }

    std::printf("\n5b. Secular, differenced (SecularPropagatorCorrected::propagateStm) with J2 = J4 = 0, against\n"
                "    the closed-form two-body Phi\n\n");
    std::printf("| Orbit | 4 h | 1 day | 1 week |\n|---|---|---|---|\n");
    ZonalParameters none;
    none.J2 = 0.0;
    none.J4 = 0.0;
    for (const Orbit& o : kOrbits) {
        const StateVector s0 = elementsToState(o.el);
        std::printf("| %s |", o.name);
        for (double dt : {14400.0, 86400.0, 7 * 86400.0}) {
            std::printf(" %.1e |",
                        maxBlockRelDiff(SecularPropagatorCorrected(s0, none).propagateStm(dt).stm,
                                        propagateKeplerStm(s0, dt).stm));
        }
        std::printf("\n");
    }
    const StateVector gto = elementsToState({26560e3, 0.74, deg2rad(63.4), 0.3, 4.7, 0.2});
    std::printf("| e = 0.74 |");
    for (double dt : {14400.0, 86400.0, 7 * 86400.0}) {
        std::printf(" %.1e |",
                    maxBlockRelDiff(SecularPropagatorCorrected(gto, none).propagateStm(dt).stm,
                                    propagateKeplerStm(gto, dt).stm));
    }
    std::printf("\n");

    std::printf("\n5c. Covariance (100 m and 0.1 m/s per axis at the epoch, no process noise): RIC sigmas of the\n"
                "    secular J2/J4 and two-body covariances over those of the numerical J2/J4 covariance\n\n");
    std::printf("| Orbit | dt | numerical sigma R / I / C (m) | secular / numerical | two-body / numerical |\n");
    std::printf("|---|---|---|---|---|\n");
    Mat6 P0{};
    for (std::size_t i = 0; i < 3; ++i) {
        P0[i][i]         = 100.0 * 100.0;
        P0[i + 3][i + 3] = 0.1 * 0.1;
    }
    for (const Orbit& o : kOrbits) {
        const StateVector s0                       = elementsToState(o.el); // TOD
        TrackerOptions    tod                      = tightOptions();
        tod.frame                                  = Frame::TOD;
        const auto                       numerical = makePropagator({kEpoch, s0, Mat6{}}, tod);
        const SecularPropagatorCorrected secular(s0);
        for (double dt : {14400.0, 86400.0}) {
            const Vec3 n = ricSigmas(numerical.propagateGridWithCovariance(dt, 1.0, 1, P0)[0]);
            const Vec3 s = ricSigmas(secular.propagateWithCovariance(dt, P0));
            const Vec3 k = ricSigmas(propagateKeplerWithCovariance(s0, P0, dt));
            std::printf("| %s | %.0f h | %.0f / %.0f / %.0f | %.3f / %.3f / %.3f | %.3f / %.3f / %.3f |\n",
                        o.name,
                        dt / 3600.0,
                        n.x,
                        n.y,
                        n.z,
                        s.x / n.x,
                        s.y / n.y,
                        s.z / n.z,
                        k.x / n.x,
                        k.y / n.y,
                        k.z / n.z);
        }
    }
}

} // namespace

int main() {
    toleranceStudy();
    interpolationStudy();
    secularStudy();
    secularStartStudy();
    poleStudy();
    covarianceStudy();
    return 0;
}
