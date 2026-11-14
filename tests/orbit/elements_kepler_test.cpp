/**********************************************************************
 * elements_kepler_test.cpp
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

// Accuracy and domain tests for Kepler's equation (elliptic, hyperbolic, parabolic), the element
// conversions near their singular angles, the equinoctial elements and the secular design helpers.
//
// The truths used here are independent of the library: solutions computed with mpmath at 60 digits,
// residuals evaluated in extended precision (long double) in a form that does not cancel, and
// states built from elements with the perifocal formulas written out in extended precision.

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>

#include "astro/core/angle.h"
#include "astro/core/constants.h"
#include "astro/core/error.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/kepler.h"
#include "astro/orbit/secular.h"
#include "test_util.h"

using namespace astro;

namespace {

constexpr double      kEps    = std::numeric_limits<double>::epsilon();
constexpr double      kInf    = std::numeric_limits<double>::infinity();
constexpr double      kNaN    = std::numeric_limits<double>::quiet_NaN();
constexpr double      kMax    = std::numeric_limits<double>::max();
constexpr double      kPi     = constants::pi;
constexpr double      kTwoPi  = constants::twoPi;
constexpr long double kTwoPiL = 6.283185307179586476925286766559005768L;

} // namespace

// ------------------------------------------------------------------------- Kepler's equation

namespace {

/// sinh(x) - x (sign = +1) or x - sin(x) (sign = -1) in extended precision. Below 0.5 the
/// difference is summed as a Taylor series, because the direct subtraction cancels.
long double cubicTail(long double x, int sign) {
    if (std::abs(x) >= 0.5L) {
        return sign > 0 ? std::sinh(x) - x : x - std::sin(x);
    }
    long double term = x * x * x / 6.0L;
    long double sum  = 0.0L;
    for (int k = 1; k < 30; ++k) {
        sum  += term;
        term *= static_cast<long double>(sign) * x * x / static_cast<long double>((2 * k + 2) * (2 * k + 3));
    }
    return sum;
}

/// Residual of Kepler's equation E - e sin E = M at a returned E, in units of the rounding of the
/// result and of the input: |residual| / (eps (|M| + |E| dM/dE)). A value of order one means E
/// solves the equation to about one unit in the last place of E or of M.
///
/// Evaluated in extended precision as (1 - e) E + e (E - sin E) - M, after taking the same whole
/// number of revolutions off E and M (the equation is periodic), so that nothing cancels for e
/// near 1 close to a perigee passage.
double ellipticResidualInUlps(double E, double e, double M) {
    const long double revs = std::nearbyint(static_cast<long double>(M) / kTwoPiL);
    const long double Er   = static_cast<long double>(E) - revs * kTwoPiL;
    const long double Mr   = static_cast<long double>(M) - revs * kTwoPiL;
    const long double el   = e;
    const long double res  = (1.0L - el) * Er + el * cubicTail(Er, -1) - Mr;
    if (res == 0.0L) {
        return 0.0;
    }
    const long double sh   = std::sin(0.5L * Er);
    const long double dMdE = (1.0L - el) + 2.0L * el * sh * sh;
    const long double unit = static_cast<long double>(kEps)
                           * (std::abs(static_cast<long double>(M)) + std::abs(static_cast<long double>(E)) * dMdE);
    return static_cast<double>(std::abs(res) / unit);
}

/// The same for the hyperbolic equation e sinh F - F = M, as (e - 1) F + e (sinh F - F) - M.
double hyperbolicResidualInUlps(double F, double e, double M) {
    const long double Fl = F, el = e, Ml = M;
    const long double res = (el - 1.0L) * Fl + el * cubicTail(Fl, +1) - Ml;
    if (res == 0.0L) {
        return 0.0;
    }
    const long double sh   = std::sinh(0.5L * Fl);
    const long double dMdF = (el - 1.0L) + 2.0L * el * sh * sh;
    const long double unit = static_cast<long double>(kEps) * (std::abs(Ml) + std::abs(Fl) * dMdF);
    return static_cast<double>(std::abs(res) / unit);
}

struct AnomalyReference {
    double e;
    double M;
    double anomaly;
};

} // namespace

TEST(Kepler, EccentricAnomalyMatchesHighPrecisionReference) {
    // Solutions of E - e sin E = M by bisection in mpmath at 60 digits, for the exact double e and M
    // shown, rounded to double. They include e = 1 - 1e-12 from perigee to apogee, and mean
    // anomalies 16 and 1000 revolutions away from zero.
    // clang-format off
    const AnomalyReference refs[] = {
        {0.5,            1.0,                1.4987011335178484    },
        {0.9,            0.02,               0.18976764205337365   },
        {0.999,          0.001,              0.17085095632357902   },
        {0.999999,       1e-08,              0.003407264597719929  },
        {0.999999999,    1e-06,              0.018171195869132243  },
        {0.999999999999, 1e-15,              1.8061145475683216e-05},
        {0.999999999999, 1e-10,              0.0008434303040921716 },
        {0.999999999999, 0.0001,             0.08435326955645084   },
        {0.999999999999, 0.5,                1.497300389094816     },
        {0.999999999999, -3.0,               -3.0707667271420047   },
        {0.999999999999, 3.141592653589793,  3.141592653589793     },
        {0.7,            100.0,              99.35343692253775     },
        {0.999999999999, -6283.185307079586, -6283.17687282926     },
    };
    // clang-format on
    for (const AnomalyReference& ref : refs) {
        const double E = meanToEccentricAnomaly(ref.M, ref.e);
        // Four units in the last place of E, plus what four units in the last place of M do to E
        // (dE/dM = 1 / (1 - e cos E) is large near the perigee of a very eccentric orbit; it
        // matters only in the last row, where M is known to 9e-13 and dE/dM is 2.8e4).
        const long double dMdE =
            1.0L - static_cast<long double>(ref.e) * std::cos(static_cast<long double>(ref.anomaly));
        const double tol = 4.0 * kEps * (std::abs(ref.anomaly) + std::abs(ref.M) / static_cast<double>(dMdE));
        EXPECT_NEAR(E, ref.anomaly, tol) << "e=1-" << 1.0 - ref.e << " M=" << ref.M;
    }
}

TEST(Kepler, EccentricAnomalyResidualIsAtRoundingLevelUpToNearParabolic) {
    const double eccentricities[] =
        {0.0, 1e-300, 1e-8, 0.1, 0.5, 0.9, 0.99, 0.999999, 1.0 - 1e-9, 1.0 - 1e-12, std::nextafter(1.0, 0.0)};
    const double meanAnomalies[] =
        {1e-300, 1e-15, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 0.5, 1.0, 2.0, 3.0, kPi - 1e-9, kPi};
    double worst = 0.0, worstE = 0.0, worstM = 0.0;
    for (const double e : eccentricities) {
        for (const double m : meanAnomalies) {
            for (const double M : {m, -m}) {
                const double E = meanToEccentricAnomaly(M, e);
                // E has the sign of M and stays within the half revolution.
                EXPECT_LE(std::abs(E), kPi) << "e=1-" << 1.0 - e << " M=" << M;
                EXPECT_EQ(std::signbit(E), std::signbit(M)) << "e=1-" << 1.0 - e << " M=" << M;
                const double res = ellipticResidualInUlps(E, e, M);
                if (!(res <= worst)) {
                    worst  = res;
                    worstE = e;
                    worstM = M;
                }
            }
        }
    }
    // Measured 0.9 over six million random (e, M) with glibc; 4 leaves room for another libm.
    EXPECT_LE(worst, 4.0) << "worst residual at e=1-" << 1.0 - worstE << " M=" << worstM;
}

TEST(Kepler, EccentricAnomalyKeepsTheRevolutionAndItsAccuracy) {
    double worst = 0.0, worstE = 0.0, worstM = 0.0;
    for (const double e : {0.3, 0.9, 1.0 - 1e-9, 1.0 - 1e-12}) {
        for (const double revs : {1.0, 7.0, 1000.0, -1000.0, 123456.0}) {
            for (const double offset : {-1.0, -1e-3, -1e-7, 0.0, 1e-7, 1e-3, 1.0, 3.0}) {
                const double M = kTwoPi * revs + offset;
                const double E = meanToEccentricAnomaly(M, e);
                // E - M = e sin E, so E is within e of M: the revolution count is kept.
                EXPECT_LE(std::abs(E - M), e * (1.0 + 4.0 * kEps) + 4.0 * kEps * std::abs(M))
                    << "e=1-" << 1.0 - e << " M=" << M;
                const double res = ellipticResidualInUlps(E, e, M);
                if (!(res <= worst)) {
                    worst  = res;
                    worstE = e;
                    worstM = M;
                }
            }
        }
    }
    // The reduction of M by whole revolutions uses 2 pi rounded to double, which is the same as an
    // error of 0.2 units in the last place of M; the rest is as in the test above.
    EXPECT_LE(worst, 4.0) << "worst residual at e=1-" << 1.0 - worstE << " M=" << worstM;
}

TEST(Kepler, HyperbolicAnomalyMatchesHighPrecisionReference) {
    // Solutions of e sinh F - F = M by bisection in mpmath at 60 digits, for the exact double e and
    // M shown, rounded to double. e = 1.005, M = 0.5 is the audit's reproducer (it threw).
    // clang-format off
    const AnomalyReference refs[] = {
        {1.5,            1.0,       1.1616354445046073    },
        {1.005,          0.5,       1.3879821326276007    },
        {1.001,          0.9,       1.6731177001539475    },
        {1.000001,       1e-08,     0.0034072615353025817 },
        {1.000000000001, 1e-15,     1.8061133256342496e-05},
        {1.000000000001, 1e-10,     0.0008434302838285713 },
        {1.000000000001, 0.0001,    0.08433326955400465   },
        {1.000000000001, 0.5,       1.396250871729208     },
        {1.000000000001, -3.0,      -2.3853380234835835   },
        {1.000000000001, 1000000.0, 14.508672247090466    },
        {2.0,            -50.0,     -3.989125544758968    },
        {1.2,            1e+19,     44.25994239065286     },
        {1.2,            2e+19,     44.95308957121281     },
        {10.0,           1e+300,    689.1660899857796     },
    };
    // clang-format on
    for (const AnomalyReference& ref : refs) {
        double F = kNaN;
        EXPECT_NO_THROW(F = meanToHyperbolicAnomaly(ref.M, ref.e)) << "e-1=" << ref.e - 1.0 << " M=" << ref.M;
        // Four units in the last place of F (one unit in the last place of M moves F by less than
        // one in its own: (M / F) dF/dM <= 1).
        EXPECT_NEAR(F, ref.anomaly, 4.0 * kEps * std::abs(ref.anomaly)) << "e-1=" << ref.e - 1.0 << " M=" << ref.M;
    }
}

TEST(Kepler, HyperbolicAnomalyConvergesForNearParabolicOrbits) {
    // The audit's reproducer: at e = 1.005, 52 % of M in (0, 1) threw "did not converge".
    int thrown = 0;
    for (int k = 1; k < 1000; ++k) {
        const double M = k / 1000.0;
        try {
            const double F = meanToHyperbolicAnomaly(M, 1.005);
            EXPECT_LE(hyperbolicResidualInUlps(F, 1.005, M), 4.0) << "M=" << M;
        } catch (const Error&) {
            ++thrown;
        }
    }
    EXPECT_EQ(thrown, 0) << "of 999 mean anomalies in (0, 1) at e = 1.005";

    // All eccentricities down to e - 1 = 1e-12 and beyond, and mean anomalies of any size.
    const double excesses[] = {kEps, 1e-12, 1e-9, 1e-6, 1e-3, 5e-3, 1e-2, 0.1, 1.0, 49.0, 1e6};
    const double meanAnomalies[] =
        {1e-300, 1e-15, 1e-12, 1e-9, 1e-6, 1e-3, 0.1, 0.5, 0.9, 1.0, 1.5, 10.0, 1e4, 1e8, 1e19, 1e300, kMax};
    double worst = 0.0, worstE = 0.0, worstM = 0.0;
    for (const double excess : excesses) {
        const double e = 1.0 + excess;
        for (const double m : meanAnomalies) {
            double F = kNaN;
            ASSERT_NO_THROW(F = meanToHyperbolicAnomaly(m, e)) << "e-1=" << excess << " M=" << m;
            EXPECT_TRUE(std::isfinite(F) && F > 0.0) << "e-1=" << excess << " M=" << m << " F=" << F;
            // Odd in M, exactly.
            EXPECT_EQ(meanToHyperbolicAnomaly(-m, e), -F) << "e-1=" << excess << " M=" << m;
            const double res = hyperbolicResidualInUlps(F, e, m);
            if (!(res <= worst)) {
                worst  = res;
                worstE = e;
                worstM = m;
            }
        }
    }
    // Measured 1.0 over six million random (e, M) with glibc; 4 leaves room for another libm.
    EXPECT_LE(worst, 4.0) << "worst residual at e-1=" << worstE - 1.0 << " M=" << worstM;
}

TEST(Kepler, SolversRejectEccentricitiesOutsideTheirDomainAndPassNonFiniteMeanAnomaliesThrough) {
    EXPECT_THROW((void)meanToEccentricAnomaly(1.0, 1.0), Error);
    EXPECT_THROW((void)meanToEccentricAnomaly(1.0, -0.1), Error);
    EXPECT_THROW((void)meanToEccentricAnomaly(1.0, kNaN), Error);
    EXPECT_TRUE(std::isnan(meanToEccentricAnomaly(kNaN, 0.3)));
    EXPECT_TRUE(std::isnan(meanToEccentricAnomaly(kInf, 0.3))); // an infinite angle has no revolution

    EXPECT_THROW((void)meanToHyperbolicAnomaly(1.0, 1.0), Error);
    EXPECT_THROW((void)meanToHyperbolicAnomaly(1.0, 0.5), Error);
    EXPECT_THROW((void)meanToHyperbolicAnomaly(1.0, kNaN), Error);
    EXPECT_THROW((void)meanToHyperbolicAnomaly(1.0, kInf), Error);
    EXPECT_TRUE(std::isnan(meanToHyperbolicAnomaly(kNaN, 1.5)));
    EXPECT_EQ(meanToHyperbolicAnomaly(kInf, 1.5), kInf);
    EXPECT_EQ(meanToHyperbolicAnomaly(-kInf, 1.5), -kInf);
}

namespace {

/// Residual of Barker's equation D + D^3/3 = M at a returned D, in units of rounding as above
/// (dM/dD = 1 + D^2), evaluated in extended precision.
double parabolicResidualInUlps(double D, double M) {
    const long double Dl = D, Ml = M;
    const long double res = Dl + Dl * Dl * Dl / 3.0L - Ml;
    if (res == 0.0L) {
        return 0.0;
    }
    const long double unit = static_cast<long double>(kEps) * (std::abs(Ml) + std::abs(Dl) * (1.0L + Dl * Dl));
    return static_cast<double>(std::abs(res) / unit);
}

} // namespace

TEST(Kepler, ParabolicAnomalySolvesBarkersEquationForBothSigns) {
    // Solutions of D + D^3/3 = M by bisection in mpmath at 60 digits, rounded to double. The old
    // closed form returned 0.99964e-12 for M = 1e-12 and -infinity for M = -1e8.
    // clang-format off
    const double refs[][2] = {
        {1e-12,  1e-12                },
        {-0.001, -0.000999999666667   },
        {0.5,    0.46622052391077345  },
        {-3.0,   -1.6096954940166688  },
        {1000.0, 14.353160112373454   },
        {-1e8,   -669.4314562805873   },
        {1e30,   14422495703.074083   },
        {-1e200, -6.694329500821695e66},
    };
    // clang-format on
    for (const auto& ref : refs) {
        EXPECT_NEAR(meanToParabolicAnomaly(ref[0]), ref[1], 4.0 * kEps * std::abs(ref[1])) << "M=" << ref[0];
    }

    double worst = 0.0, worstM = 0.0;
    for (const double m : {5e-324, 1e-300, 1e-12, 1e-9, 1e-6, 1e-3, 0.1,  0.5,   1.0,   3.0,
                           10.0,   1e3,    1e6,   1e8,  1e9,  1e15, 1e30, 1e100, 1e300, kMax}) {
        const double D = meanToParabolicAnomaly(m);
        EXPECT_TRUE(std::isfinite(D) && D > 0.0) << "M=" << m << " D=" << D;
        // Odd in M, exactly: a negative M is as accurate as a positive one.
        EXPECT_EQ(meanToParabolicAnomaly(-m), -D) << "M=" << m;
        const double res = parabolicResidualInUlps(D, m);
        if (!(res <= worst)) {
            worst  = res;
            worstM = m;
        }
    }
    // Measured 1.8 over four million random M of both signs from 1e-300 to 1.8e308 with glibc.
    EXPECT_LE(worst, 4.0) << "worst residual of Barker's equation at |M|=" << worstM;

    EXPECT_EQ(meanToParabolicAnomaly(0.0), 0.0);
    EXPECT_TRUE(std::isnan(meanToParabolicAnomaly(kNaN)));
    EXPECT_EQ(meanToParabolicAnomaly(kInf), kInf);
    EXPECT_EQ(meanToParabolicAnomaly(-kInf), -kInf);
}

// ------------------------------------------------------------------------- Classical elements

namespace {

constexpr double km = 1000.0;

// Tolerances of the element tests: a few units of rounding (eps = 2.2e-16), about four times what
// was measured. The acos formulas missed each of them by seven orders of magnitude or more.
constexpr double kInclinationTolerance = 8.0 * kEps;  // relative; measured 2 eps
constexpr double kAngleTolerance       = 16.0 * kEps; // rad; measured 4 eps
// Relative error of a state rebuilt from its elements: a, e and four angles each carry a rounding
// error. Measured 2e-15; 1e-14 is 0.07 micrometres at 7000 km (acos: 4e-8, 0.3 m).
constexpr double kRoundTripTolerance = 1e-14;
// SecularPropagator at its epoch: its Kepler solver stops at a residual of 2e-15 rad, 1.4e-8 m at
// 7000 km, on top of the element conversion. Measured 6e-9 m and 5e-12 m/s (acos: 5e-3 m).
constexpr double kEpochPositionTolerance = 5e-8;  // m
constexpr double kEpochVelocityTolerance = 5e-11; // m/s

/// Difference of two angles, wrapped to [-pi, pi].
double angleDifference(double a, double b) {
    return std::remainder(a - b, kTwoPi);
}

/// Position and velocity from classical elements by the perifocal construction (unit vectors P
/// towards perigee and Q a quarter turn ahead of it in the orbit plane), evaluated in extended
/// precision and rounded to double. Independent of elementsToState.
StateVector
stateFromElements(double a, double e, double i, double raan, double argp, double nu, double mu = constants::gmEarth) {
    using L        = long double;
    const L     cO = std::cos(static_cast<L>(raan)), sO = std::sin(static_cast<L>(raan));
    const L     cw = std::cos(static_cast<L>(argp)), sw = std::sin(static_cast<L>(argp));
    const L     ci = std::cos(static_cast<L>(i)), si = std::sin(static_cast<L>(i));
    const L     cn = std::cos(static_cast<L>(nu)), sn = std::sin(static_cast<L>(nu));
    const L     P[3] = {cO * cw - sO * sw * ci, sO * cw + cO * sw * ci, sw * si};
    const L     Q[3] = {-cO * sw - sO * cw * ci, -sO * sw + cO * cw * ci, cw * si};
    const L     el   = e;
    const L     p    = static_cast<L>(a) * (1.0L - el * el);
    const L     rm   = p / (1.0L + el * cn);
    const L     k    = std::sqrt(static_cast<L>(mu) / p);
    StateVector s;
    for (std::size_t j = 0; j < 3; ++j) {
        s.r[j] = static_cast<double>(rm * (cn * P[j] + sn * Q[j]));
        s.v[j] = static_cast<double>(k * (-sn * P[j] + (el + cn) * Q[j]));
    }
    return s;
}

std::string describe(double a, double e, double i, double raan, double argp, double nu) {
    std::ostringstream os;
    os.precision(17);
    os << "a=" << a / km << " km, e=" << e << ", i=" << i << " (pi-i=" << kPi - i << "), raan=" << raan
       << ", argp=" << argp << ", nu=" << nu;
    return os.str();
}

} // namespace

TEST(Elements, AnglesNearZeroAndPiAreRecoveredToRounding) {
    // States built from known elements; stateToElements must return them. Formulas with acos
    // cannot: near 0 and pi they resolve an angle only to ~1.5e-8 rad, and an inclination of 1e-9
    // came back as exactly 0.
    const double inclinations[] = {1e-9, 1e-6, 0.4, 0.5 * kPi, 2.5, kPi - 1e-6, kPi - 1e-9};
    const double perigees[]     = {1e-9, 1.0, kPi - 1e-9, kPi + 1e-9, 5.0, kTwoPi - 1e-9};
    const double anomalies[]    = {1e-9, 1e-6, 2.0, kPi - 1e-9, kPi + 1e-9, 4.0, kTwoPi - 1e-6, kTwoPi - 1e-9};
    const double nodes[]        = {0.3, 3.0, 5.5};
    // clang-format off
    const double shapes[][2]    = {{7000.0 * km,  0.1},
                                   {26600.0 * km, 0.7}};
    // clang-format on
    double      worstI = 0.0, worstNode = 0.0, worstPerigee = 0.0, worstAnomaly = 0.0, worstLatitude = 0.0;
    std::string atI, atNode, atPerigee, atAnomaly, atLatitude;
    std::size_t count = 0;
    for (const auto& shape : shapes) {
        const double a = shape[0], e = shape[1];
        for (const double i : inclinations) {
            for (const double argp : perigees) {
                for (const double nu : anomalies) {
                    const double            raan = nodes[count++ % 3];
                    const KeplerianElements el   = stateToElements(stateFromElements(a, e, i, raan, argp, nu));
                    const auto              note = [&](double error, double& worst, std::string& at) {
                        if (!(error <= worst)) {
                            worst = error;
                            at    = describe(a, e, i, raan, argp, nu);
                        }
                    };
                    EXPECT_TRUE(el.raan >= 0.0 && el.raan < kTwoPi && el.argp >= 0.0 && el.argp < kTwoPi && el.nu >= 0.0
                                && el.nu < kTwoPi)
                        << describe(a, e, i, raan, argp, nu);
                    // The inclination relative to itself: 1e-9 rad is recovered as 1e-9 rad.
                    note(std::abs(el.i - i) / i, worstI, atI);
                    note(std::abs(angleDifference(el.raan, raan)), worstNode, atNode);
                    // Perigee direction and true anomaly carry the rounding of the eccentricity
                    // vector, eps / e; their sum, the argument of latitude, does not.
                    note(std::abs(angleDifference(el.argp, argp)) * e, worstPerigee, atPerigee);
                    note(std::abs(angleDifference(el.nu, nu)) * e, worstAnomaly, atAnomaly);
                    note(std::abs(angleDifference(el.argp + el.nu, argp + nu)), worstLatitude, atLatitude);
                }
            }
        }
    }
    EXPECT_LE(worstI, kInclinationTolerance) << "relative inclination error, at " << atI;
    EXPECT_LE(worstNode, kAngleTolerance) << "RAAN error (rad), at " << atNode;
    EXPECT_LE(worstPerigee, kAngleTolerance) << "argument of perigee error times e (rad), at " << atPerigee;
    EXPECT_LE(worstAnomaly, kAngleTolerance) << "true anomaly error times e (rad), at " << atAnomaly;
    EXPECT_LE(worstLatitude, kAngleTolerance) << "argument of latitude error (rad), at " << atLatitude;
}

TEST(Elements, ArgumentOfLatitudeOfCircularOrbitsIsRecoveredNearZeroAndPi) {
    for (const double i : {1e-6, 0.4, 2.5, kPi - 1e-6}) {
        for (const double u : {1e-9, 1e-6, 1.0, kPi - 1e-9, kPi + 1e-9, 4.0, kTwoPi - 1e-9}) {
            const KeplerianElements el = stateToElements(stateFromElements(7000.0 * km, 0.0, i, 2.0, 0.0, u));
            EXPECT_LT(el.e, 1e-11);
            EXPECT_EQ(el.argp, 0.0);
            EXPECT_NEAR(angleDifference(el.nu, u), 0.0, kAngleTolerance) << "i=" << i << " u=" << u;
        }
    }
}

TEST(Elements, RoundTripThroughElementsIsExactToRounding) {
    // state -> elements -> state for ordinary orbits with every angle at and near 0 and pi. The
    // truth is the starting state.
    const double inclinations[] = {0.0, 1e-9, 1e-6, 0.4, 0.5 * kPi, 2.5, kPi - 1e-6, kPi - 1e-9, kPi};
    const double perigees[]     = {0.0, 1e-9, 1.0, kPi - 1e-9, kPi + 1e-9, 5.0, kTwoPi - 1e-9};
    const double anomalies[] = {0.0, 1e-9, 1e-6, 2.0, kPi - 1e-9, kPi, kPi + 1e-9, 4.0, kTwoPi - 1e-6, kTwoPi - 1e-9};
    const double nodes[]     = {0.0, 0.3, 3.0, 5.5};
    // clang-format off
    const double shapes[][2]    = {{7000.0 * km,  0.0  },
                                   {7000.0 * km,  0.001},
                                   {7000.0 * km,  0.1  },
                                   {26600.0 * km, 0.7  },
                                   {42164.0 * km, 0.0  }};
    // clang-format on
    double      worstR = 0.0, worstV = 0.0;
    std::string atR, atV;
    std::size_t count = 0;
    for (const auto& shape : shapes) {
        const double a = shape[0], e = shape[1];
        for (const double i : inclinations) {
            for (const double argp : perigees) {
                for (const double nu : anomalies) {
                    const double      raan = nodes[count++ % 4];
                    const StateVector s    = stateFromElements(a, e, i, raan, argp, nu);
                    const StateVector back = elementsToState(stateToElements(s));
                    const double      dr = norm(back.r - s.r) / norm(s.r), dv = norm(back.v - s.v) / norm(s.v);
                    if (!(dr <= worstR)) {
                        worstR = dr;
                        atR    = describe(a, e, i, raan, argp, nu);
                    }
                    if (!(dv <= worstV)) {
                        worstV = dv;
                        atV    = describe(a, e, i, raan, argp, nu);
                    }
                }
            }
        }
    }
    EXPECT_LE(worstR, kRoundTripTolerance) << "relative position error, at " << atR;
    EXPECT_LE(worstV, kRoundTripTolerance) << "relative velocity error, at " << atV;

    // Hyperbolic orbits: true anomaly in [-pi, pi), accurate through perigee.
    for (const double nu : {-2.0, -1e-9, 0.0, 1e-9, 1.0}) {
        const StateVector       s  = stateFromElements(-20000.0 * km, 1.5, 0.5, 1.0, 2.0, nu);
        const KeplerianElements el = stateToElements(s);
        EXPECT_NEAR(el.nu, nu, kAngleTolerance) << "nu=" << nu;
        const StateVector back = elementsToState(el);
        EXPECT_LE(norm(back.r - s.r), kRoundTripTolerance * norm(s.r)) << "nu=" << nu;
    }
}

TEST(Elements, SpecialCaseConventionsAndThresholdsAreKept) {
    const double a = 8000.0 * km;
    // Circular inclined: argument of perigee 0, nu is the argument of latitude.
    {
        const KeplerianElements el = stateToElements(stateFromElements(a, 0.0, 0.5, 1.0, 0.0, 0.7));
        EXPECT_EQ(el.argp, 0.0);
        EXPECT_NEAR(el.raan, 1.0, 1e-15);
        EXPECT_NEAR(el.nu, 0.7, 1e-15);
    }
    // Elliptic equatorial: RAAN 0, argp is the longitude of periapsis.
    {
        const KeplerianElements el = stateToElements(stateFromElements(a, 0.2, 0.0, 0.0, 1.1, 0.4));
        EXPECT_EQ(el.i, 0.0);
        EXPECT_EQ(el.raan, 0.0);
        EXPECT_NEAR(el.argp, 1.1, 1e-15);
        EXPECT_NEAR(el.nu, 0.4, 1e-15);
    }
    // Circular equatorial: RAAN and argp 0, nu is the true longitude.
    {
        const KeplerianElements el = stateToElements(stateFromElements(a, 0.0, 0.0, 0.0, 0.0, 2.0));
        EXPECT_EQ(el.raan, 0.0);
        EXPECT_EQ(el.argp, 0.0);
        EXPECT_NEAR(el.nu, 2.0, 1e-15);
    }
    // Retrograde equatorial (i = pi): the same, with longitudes measured in the direction of motion
    // (clockwise seen from +z). Each state is the mirror image in the x axis of a prograde one, so
    // its perigee is at longitude -1.1 (or the satellite at -2.0) and it moves clockwise.
    {
        const StateVector       p  = stateFromElements(a, 0.2, 0.0, 0.0, 1.1, 0.4);
        const KeplerianElements el = stateToElements({
            {p.r.x, -p.r.y, 0.0},
            {p.v.x, -p.v.y, 0.0}
        });
        EXPECT_EQ(el.i, kPi);
        EXPECT_EQ(el.raan, 0.0);
        EXPECT_NEAR(el.argp, 1.1, 4e-15); // wrapped twice, each rounding to 4.4e-16
        EXPECT_NEAR(el.nu, 0.4, 1e-15);
        const StateVector       c  = stateFromElements(a, 0.0, 0.0, 0.0, 0.0, 2.0);
        const KeplerianElements ec = stateToElements({
            {c.r.x, -c.r.y, 0.0},
            {c.v.x, -c.v.y, 0.0}
        });
        EXPECT_EQ(ec.argp, 0.0);
        EXPECT_NEAR(ec.nu, 2.0, 1e-15);
    }
    // The documented thresholds: e and i below 1e-11 are the special cases, above it they are not.
    {
        const KeplerianElements below = stateToElements(stateFromElements(a, 5e-12, 0.5, 1.0, 2.0, 0.7));
        EXPECT_EQ(below.argp, 0.0);
        EXPECT_NEAR(below.nu, 2.7, 4e-15);
        const KeplerianElements above = stateToElements(stateFromElements(a, 1e-9, 0.5, 1.0, 2.0, 0.7));
        EXPECT_NEAR(above.argp, 2.0, 1e-5); // perigee direction of e = 1e-9: good to eps / e
        EXPECT_NEAR(above.argp + above.nu, 2.7, 4e-15);
    }
    {
        const KeplerianElements below = stateToElements(stateFromElements(a, 0.2, 5e-12, 1.0, 2.0, 0.7));
        EXPECT_EQ(below.raan, 0.0);
        EXPECT_NEAR(below.argp, 3.0, 1e-10); // longitude of periapsis = RAAN + argp
        // An inclination of 2e-11 used to come back as 0, and so as equatorial.
        const KeplerianElements above = stateToElements(stateFromElements(a, 0.2, 2e-11, 1.0, 2.0, 0.7));
        EXPECT_NEAR(above.i, 2e-11, 1e-20);
        EXPECT_NEAR(above.raan, 1.0, 1e-9);
        EXPECT_NEAR(above.argp, 2.0, 1e-9);
    }
}

TEST(Secular, PropagatorReproducesItsInputStateAtTheEpochNearSingularAngles) {
    // "The input state is reproduced at the epoch" (secular.h), also where acos could not resolve
    // the angles: true anomaly and argument of perigee within 1e-7 of 0 and pi (it gave 5 mm).
    for (const double nu : {1e-9, 1e-7, kPi - 1e-7, kTwoPi - 1e-8}) {
        for (const double argp : {1e-8, kPi - 1e-9, 4.0}) {
            const StateVector s = stateFromElements(7000.0 * km, 0.01, 0.9, 0.3, argp, nu);
            const StateVector b = SecularPropagator(s).propagate(0.0);
            test::expectVecNear(b.r, s.r, kEpochPositionTolerance, "position");
            test::expectVecNear(b.v, s.v, kEpochVelocityTolerance, "velocity");
        }
    }
}

// ------------------------------------------------------------------------- Equinoctial elements

namespace {

// Measured 1.2e-14 in a, h, k, p, q (worst in p and q at i = 2.8 rad, where 1 / (1 + cos i) = 17
// multiplies the rounding), 1.8e-15 rad in the mean longitude, and 1.8e-14 of |r| for the state
// rebuilt from the elements (worst at e = 0.9).
constexpr double kEquinoctialElementTolerance   = 5e-14;
constexpr double kMeanLongitudeTolerance        = 1e-14; // rad
constexpr double kEquinoctialRoundTripTolerance = 1e-13;

/// Mean anomaly of an elliptic orbit from the true anomaly, in extended precision.
long double meanFromTrueAnomaly(double nu, double e) {
    using L    = long double;
    const L el = e;
    const L E  = 2.0L * std::atan(std::sqrt((1.0L - el) / (1.0L + el)) * std::tan(0.5L * static_cast<L>(nu)));
    return E - el * std::sin(E);
}

} // namespace

TEST(Elements, EquinoctialElementsMatchTheirDefinitionForEllipses) {
    // h = e sin(raan + argp), k = e cos(raan + argp), p = tan(i/2) sin(raan), q = tan(i/2) cos(raan),
    // lambda = raan + argp + M (Broucke & Cefola 1972), from the elements the state was built with.
    double worstElement = 0.0, worstLambda = 0.0, worstR = 0.0;
    for (const double e : {0.0, 0.01, 0.3, 0.7, 0.9}) {
        for (const double i : {0.0, 1e-9, 0.4, 1.5, 2.8}) {
            for (const double nu : {0.0, 1e-9, 1.0, 3.0, kPi, 4.0, 6.0}) {
                const double              a = 12000.0 * km / (1.0 - e), raan = 1.0, argp = 2.0;
                const StateVector         s = stateFromElements(a, e, i, raan, argp, nu);
                const EquinoctialElements q = stateToEquinoctial(s);
                const double              t = std::tan(0.5 * i);
                worstElement                = std::max({worstElement,
                                                        std::abs(q.a / a - 1.0),
                                                        std::abs(q.h - e * std::sin(raan + argp)),
                                                        std::abs(q.k - e * std::cos(raan + argp)),
                                                        std::abs(q.p - t * std::sin(raan)),
                                                        std::abs(q.q - t * std::cos(raan))});
                const double lambda =
                    static_cast<double>(static_cast<long double>(raan) + argp + meanFromTrueAnomaly(nu, e));
                EXPECT_TRUE(q.lambda >= 0.0 && q.lambda < kTwoPi);
                worstLambda = std::max(worstLambda, std::abs(angleDifference(q.lambda, lambda)));
                // ... and back to the state.
                const StateVector back = equinoctialToState(q);
                worstR                 = std::max(worstR, norm(back.r - s.r) / norm(s.r));
            }
        }
    }
    EXPECT_LE(worstElement, kEquinoctialElementTolerance);
    EXPECT_LE(worstLambda, kMeanLongitudeTolerance);
    EXPECT_LE(worstR, kEquinoctialRoundTripTolerance);
}

TEST(Elements, EquinoctialElementsAreDefinedForEllipsesOnly) {
    // The mean longitude is an angle only when the mean anomaly is one. For a hyperbola its wrapped
    // value was returned without warning and the round trip was off by 1.6e8 m (audit reproducer:
    // a = -20000 km, e = 1.5, i = 0.5, raan = 1, argp = 2, nu = 1.5).
    const KeplerianElements hyperbola{-20000.0 * km, 1.5, 0.5, 1.0, 2.0, 1.5};
    const StateVector       s =
        stateFromElements(hyperbola.a, hyperbola.e, hyperbola.i, hyperbola.raan, hyperbola.argp, hyperbola.nu);
    EXPECT_THROW((void)stateToEquinoctial(s), Error);
    EXPECT_THROW((void)stateToEquinoctialJacobian(s), Error);
    // Escape speed exactly: parabolic.
    const double      r = 7000.0 * km;
    const StateVector parabola{
        {r,   0.0,                                     0.0},
        {0.0, std::sqrt(2.0 * constants::gmEarth / r), 0.0}
    };
    EXPECT_THROW((void)stateToEquinoctial(parabola), Error);
    // h^2 + k^2 >= 1, or a <= 0, is not an ellipse.
    EXPECT_THROW((void)equinoctialToKeplerian({-20000.0 * km, 1.2, 0.9, 0.1, 0.2, 1.0}), Error);
    EXPECT_THROW((void)equinoctialToKeplerian({20000.0 * km, 0.6, 0.8, 0.1, 0.2, 1.0}), Error);
    EXPECT_THROW((void)equinoctialToKeplerian({-20000.0 * km, 0.1, 0.2, 0.1, 0.2, 1.0}), Error);
    EXPECT_THROW((void)equinoctialToKeplerian({kNaN, 0.1, 0.2, 0.1, 0.2, 1.0}), Error);
    EXPECT_THROW((void)equinoctialToState({-20000.0 * km, 1.2, 0.9, 0.1, 0.2, 1.0}), Error);
    // keplerianToEquinoctial throws as equinoctialToKeplerian does: for a hyperbola, a parabola, a <= 0 and NaN.
    EXPECT_THROW((void)keplerianToEquinoctial(hyperbola), Error);
    EXPECT_THROW((void)keplerianToEquinoctial({20000.0 * km, 1.0, 0.5, 1.0, 2.0, 1.5}), Error);  // parabola
    EXPECT_THROW((void)keplerianToEquinoctial({0.0, 0.1, 0.5, 1.0, 2.0, 1.5}), Error);           // a = 0
    EXPECT_THROW((void)keplerianToEquinoctial({-20000.0 * km, 0.1, 0.5, 1.0, 2.0, 1.5}), Error); // a < 0
    EXPECT_THROW((void)keplerianToEquinoctial({kNaN, 0.1, 0.5, 1.0, 2.0, 1.5}), Error);
    EXPECT_THROW((void)keplerianToEquinoctial({20000.0 * km, kNaN, 0.5, 1.0, 2.0, 1.5}), Error);
    EXPECT_NO_THROW((void)keplerianToEquinoctial({20000.0 * km, 0.999, 0.5, 1.0, 2.0, 1.5}));
    // An ellipse close to the parabola is still converted.
    const StateVector         nearParabolic = stateFromElements(7000.0 * km / (1.0 - 0.999), 0.999, 0.5, 1.0, 2.0, 1.5);
    const EquinoctialElements qe            = stateToEquinoctial(nearParabolic);
    EXPECT_NEAR(std::hypot(qe.h, qe.k), 0.999, 1e-14);
    EXPECT_NO_THROW((void)equinoctialToKeplerian(qe));
}

// ------------------------------------------------------------------------- Secular helpers

namespace {

/// The first-order J2 secular rates (Brouwer 1959) written out in their textbook form, for checking
/// the design helpers without the library's secularRates.
struct J2Rates {
    double raanDot, argpDot, meanAnomalyDot;
};

J2Rates j2Rates(double a, double e, double i, const CentralBody& body) {
    const double n  = std::sqrt(body.mu / (a * a * a));
    const double p  = a * (1.0 - e * e);
    const double k  = 1.5 * n * body.zonal(2) * (body.re / p) * (body.re / p);
    const double ci = std::cos(i);
    return {-k * ci, 0.5 * k * (5.0 * ci * ci - 1.0), n + 0.5 * k * std::sqrt(1.0 - e * e) * (3.0 * ci * ci - 1.0)};
}

} // namespace

TEST(Secular, RepeatGroundTrackSatisfiesTheRepeatCondition) {
    // revs nodal periods must equal days nodal days.
    const CentralBody body;
    struct Case {
        int    revs, days;
        double e, i;
    };
    // clang-format off
    const Case cases[] = {{233, 16, 0.0,   degToRad(98.2)},
                          {15,  1,  0.001, degToRad(51.6)},
                          {2,   1,  0.72,  degToRad(63.4)},
                          {1,   1,  0.0,   0.0          }};
    // clang-format on
    for (const Case& c : cases) {
        const double  a = repeatGroundTrackSemiMajorAxis(c.revs, c.days, c.e, c.i, body);
        const J2Rates r = j2Rates(a, c.e, c.i, body);
        // The iteration stops when a changes by less than 1e-13 of itself.
        EXPECT_NEAR((r.meanAnomalyDot + r.argpDot) * c.days / ((constants::earthRotationRate - r.raanDot) * c.revs),
                    1.0,
                    1e-12)
            << c.revs << " revolutions in " << c.days << " days";
    }
}

TEST(Secular, RepeatGroundTrackThrowsInsteadOfReturningAnUnconvergedValue) {
    // The audit's reproducer returned 0.000 km: with e = 0.99 the perigee of a 12 h orbit is deep
    // inside the Earth, the J2 node rate is so large that no semi-major axis satisfies the repeat
    // condition, and the iteration ran away towards a = 0.
    EXPECT_THROW((void)repeatGroundTrackSemiMajorAxis(2, 1, 0.99, 1.1), Error);
    // Not an ellipse.
    EXPECT_THROW((void)repeatGroundTrackSemiMajorAxis(2, 1, 1.0, 1.1), Error);
    EXPECT_THROW((void)repeatGroundTrackSemiMajorAxis(2, 1, 1.5, 1.1), Error);
    EXPECT_THROW((void)repeatGroundTrackSemiMajorAxis(2, 1, -0.1, 1.1), Error);
    EXPECT_THROW((void)repeatGroundTrackSemiMajorAxis(2, 1, kNaN, 1.1), Error);
    EXPECT_THROW((void)repeatGroundTrackSemiMajorAxis(2, 1, 0.1, kNaN), Error);
    EXPECT_THROW((void)repeatGroundTrackSemiMajorAxis(0, 1, 0.1, 1.1), Error);
    EXPECT_THROW((void)repeatGroundTrackSemiMajorAxis(2, -1, 0.1, 1.1), Error);
}

TEST(Secular, SunSynchronousInclinationRejectsOrbitsThatAreNotEllipses) {
    EXPECT_THROW((void)sunSynchronousInclination(0.0), Error); // returned pi/2
    EXPECT_THROW((void)sunSynchronousInclination(-7200.0 * km), Error);
    EXPECT_THROW((void)sunSynchronousInclination(kNaN), Error);
    EXPECT_THROW((void)sunSynchronousInclination(7200.0 * km, 1.5), Error); // returned an angle
    EXPECT_THROW((void)sunSynchronousInclination(7200.0 * km, 1.0), Error);
    EXPECT_THROW((void)sunSynchronousInclination(7200.0 * km, -0.1), Error);
    EXPECT_THROW((void)sunSynchronousInclination(7200.0 * km, kNaN), Error);
    // A valid orbit still has its solution: the node advances 2 pi in a tropical year.
    const CentralBody body;
    const double      a = 7078.0 * km, e = 0.001;
    const double      i = sunSynchronousInclination(a, e, body);
    EXPECT_NEAR(j2Rates(a, e, i, body).raanDot * 365.242190 * 86400.0 / kTwoPi, 1.0, 1e-14);
}

TEST(Secular, PropagateSecularRejectsTheSameElementsAsThePropagatorClass) {
    // a < 0 with e < 1 returned NaN angles from the free function while the class threw.
    // clang-format off
    const KeplerianElements invalid[] = {{-7000.0 * km, 0.1,  1.0, 0.0, 0.0, 0.0},
                                         {0.0,          0.1,  1.0, 0.0, 0.0, 0.0},
                                         {kNaN,         0.1,  1.0, 0.0, 0.0, 0.0},
                                         {7000.0 * km,  -0.1, 1.0, 0.0, 0.0, 0.0},
                                         {7000.0 * km,  1.0,  1.0, 0.0, 0.0, 0.0},
                                         {-7000.0 * km, 1.2,  1.0, 0.0, 0.0, 0.0},
                                         {7000.0 * km,  kNaN, 1.0, 0.0, 0.0, 0.0}};
    // clang-format on
    for (const KeplerianElements& el : invalid) {
        EXPECT_THROW((void)propagateSecular(el, 100.0), Error) << "a=" << el.a << " e=" << el.e;
        EXPECT_THROW(SecularPropagator{el}, Error) << "a=" << el.a << " e=" << el.e;
    }
    // A valid ellipse: the node moves at the J2 rate.
    const CentralBody       body;
    const KeplerianElements el{7000.0 * km, 0.01, 1.0, 0.5, 1.0, 2.0};
    const double            raan = propagateSecular(el, 1000.0, body, SecularOrder::J2).raan;
    EXPECT_NEAR(raan, 0.5 + j2Rates(el.a, el.e, el.i, body).raanDot * 1000.0, 1e-15);
}

// ------------------------------------------------------------------------- Smaller helpers

namespace {

/// Eccentric (e < 1) or hyperbolic (e > 1) anomaly from the true anomaly, in extended precision,
/// from tan(E/2) = sqrt((1 - e) / (1 + e)) tan(nu/2) and tanh(F/2) = sqrt((e - 1) / (e + 1)) tan(nu/2).
long double anomalyFromTrue(double nu, double e) {
    using L    = long double;
    const L el = e, t = std::tan(0.5L * static_cast<L>(nu));
    return el < 1.0L ? 2.0L * std::atan(std::sqrt((1.0L - el) / (1.0L + el)) * t)
                     : 2.0L * std::atanh(std::sqrt((el - 1.0L) / (el + 1.0L)) * t);
}

} // namespace

TEST(Elements, EccentricAndMeanAnomalyOfTheElementsFollowTheirDefinitions) {
    // E and M = E - e sin E (ellipse), F and M = e sinh F - F (hyperbola), from the half-angle forms
    // in extended precision. Angles are compared modulo 2 pi. Allowed: 8 rounding errors of E and of
    // nu carried into E (dE/dnu = sqrt(|1 - e^2|) / (1 + e cos nu), 14 near the apoapsis at e = 0.99),
    // times dM/dE for M (measured: up to 0.15 of these).
    for (const double e : {0.0, 0.1, 0.7, 0.99, 1.2, 3.0}) {
        const double nuMax = e < 1.0 ? kPi : 0.95 * std::acos(-1.0 / e);
        for (const double f : {-0.9, -0.3, 0.0, 0.2, 0.6, 0.99}) {
            const double            nu = f * nuMax;
            const KeplerianElements el{(e < 1.0 ? 1.0 : -1.0) * 9000.0 * km,
                                       e,
                                       0.5,
                                       1.0,
                                       2.0,
                                       e < 1.0 ? wrapTwoPi(nu) : nu};
            const long double       E    = anomalyFromTrue(nu, e);
            const long double       M    = e < 1.0 ? E - e * std::sin(E) : e * std::sinh(E) - E;
            const double            dEdv = std::sqrt(std::abs(1.0 - e * e)) / (1.0 + e * std::cos(nu));
            const double            dMdE =
                e < 1.0 ? 1.0 - e * std::cos(static_cast<double>(E)) : e * std::cosh(static_cast<double>(E)) - 1.0;
            const double tolE = 8.0 * kEps * (std::abs(static_cast<double>(E)) + kTwoPi * dEdv);
            const double tolM = tolE * dMdE + 8.0 * kEps * std::abs(static_cast<double>(M));
            EXPECT_NEAR(angleDifference(el.eccentricAnomaly(), static_cast<double>(E)), 0.0, tolE)
                << "e=" << e << " nu=" << nu;
            EXPECT_NEAR(angleDifference(el.meanAnomaly(), static_cast<double>(M)), 0.0, tolM)
                << "e=" << e << " nu=" << nu;
        }
    }
}

TEST(Elements, ElementsFromMeanAnomalySolveKeplersEquation) {
    // The true anomaly returned must put the mean anomaly back, evaluated in extended precision from
    // the half-angle form (independent of the library's anomaly functions), and the other elements
    // are passed through. Allowed: 8 rounding errors of M and e (measured: up to 0.2 of these).
    for (const double e : {0.0, 0.3, 0.9, 1.5}) {
        for (const double M : {-2.0, -1e-6, 0.0, 0.4, 3.0, 5.5}) {
            const KeplerianElements el = elementsFromMeanAnomaly(7000.0 * km, e, 0.5, 1.0, 2.0, M);
            EXPECT_EQ(el.a, 7000.0 * km);
            EXPECT_EQ(el.e, e);
            EXPECT_EQ(el.i, 0.5);
            EXPECT_EQ(el.raan, 1.0);
            EXPECT_EQ(el.argp, 2.0);
            const long double E    = anomalyFromTrue(el.nu, e);
            const long double back = e < 1.0 ? E - e * std::sin(E) : e * std::sinh(E) - E;
            const double      tol  = 8.0 * kEps * (std::abs(M) + 1.0) * (1.0 + e);
            EXPECT_NEAR(e < 1.0 ? angleDifference(static_cast<double>(back), M) : static_cast<double>(back) - M,
                        0.0,
                        tol)
                << "e=" << e << " M=" << M;
        }
    }
    EXPECT_THROW((void)elementsFromMeanAnomaly(7000.0 * km, -0.1, 0.5, 1.0, 2.0, 1.0), Error);
    EXPECT_THROW((void)elementsFromMeanAnomaly(7000.0 * km, kNaN, 0.5, 1.0, 2.0, 1.0), Error);
}

TEST(Elements, FlightPathAngleAndApsidesMatchTheConic) {
    // tan(gamma) = e sin nu / (1 + e cos nu) for any conic; |r| at nu = 0 and pi is a (1 - e) and
    // a (1 + e). States from the extended-precision perifocal construction.
    for (const double e : {0.0, 0.2, 0.8, 1.5}) {
        const double a = (e < 1.0 ? 1.0 : -1.0) * 12000.0 * km;
        for (const double nu : {-1.5, -0.4, 0.0, 0.7, 1.9}) {
            if (e > 1.0 && std::abs(nu) > 0.95 * std::acos(-1.0 / e)) {
                continue;
            }
            const StateVector s     = stateFromElements(a, e, 0.9, 1.0, 2.0, nu);
            const double      gamma = std::atan2(e * std::sin(nu), 1.0 + e * std::cos(nu));
            EXPECT_NEAR(flightPathAngle(s), gamma, 8.0 * kEps) << "e=" << e << " nu=" << nu;
        }
        const double rp = norm(stateFromElements(a, e, 0.9, 1.0, 2.0, 0.0).r);
        EXPECT_NEAR(periapsisRadius(a, e), rp, 4.0 * kEps * rp) << "e=" << e;
        if (e < 1.0) {
            const double ra = norm(stateFromElements(a, e, 0.9, 1.0, 2.0, kPi).r);
            EXPECT_NEAR(apoapsisRadius(a, e), ra, 4.0 * kEps * ra) << "e=" << e;
        }
    }
}

TEST(Kepler, ParabolicAnomalyIsTheTangentOfHalfTheTrueAnomaly) {
    // D = tan(nu / 2) and back, and the time law: on a parabola of perifocal distance q the radius is
    // q (1 + D^2) = 2 q / (1 + cos nu).
    for (const double nu : {-3.0, -1.0, -1e-9, 0.0, 0.5, 2.0, 3.1}) {
        const double D = trueToParabolicAnomaly(nu);
        EXPECT_NEAR(D, std::tan(0.5 * nu), 4.0 * kEps * std::abs(D)) << "nu=" << nu;
        EXPECT_NEAR(parabolicToTrueAnomaly(D), nu, 4.0 * kEps * std::max(1.0, std::abs(nu))) << "nu=" << nu;
        const long double radius = 2.0L / (1.0L + std::cos(static_cast<long double>(nu))); // r / q
        EXPECT_NEAR(1.0 + D * D, static_cast<double>(radius), 8.0 * kEps * (1.0 + D * D)) << "nu=" << nu;
    }
    EXPECT_NEAR(parabolicToTrueAnomaly(1e300), kPi, 4.0 * kEps);
    EXPECT_NEAR(parabolicToTrueAnomaly(-1e300), -kPi, 4.0 * kEps);
}

TEST(Elements, ConversionsRefuseStatesAndElementsWithoutAnOrbit) {
    const double      r = 7000.0 * km, v = std::sqrt(constants::gmEarth / r);
    const StateVector circular{
        {r,   0.0, 0.0},
        {0.0, v,   0.0}
    };
    // Rectilinear: r = 0, or r parallel to v.
    EXPECT_THROW((void)stateToElements({Vec3{}, circular.v}), Error);
    EXPECT_THROW((void)stateToElements({
                     circular.r,
                     Vec3{2.0 * v, 0.0, 0.0}
    }),
                 Error);
    // Parabolic: escape speed exactly.
    EXPECT_THROW((void)stateToElements({circular.r, circular.v * std::sqrt(2.0)}), Error);
    // A semi-latus rectum that is not positive.
    EXPECT_THROW((void)elementsToState({-7000.0 * km, 0.5, 0.5, 1.0, 2.0, 0.3}), Error);
    EXPECT_THROW((void)elementsToState({7000.0 * km, 1.5, 0.5, 1.0, 2.0, 0.3}), Error);
    const StateVector retrograde{
        {r,   0.0, 0.0},
        {0.0, -v,  0.0}
    };
    // Retrograde equatorial: the equinoctial elements are singular there, and the message says so.
    const auto refusedAsRetrograde = [](auto&& convert) {
        try {
            convert();
            ADD_FAILURE() << "no exception for a retrograde equatorial orbit";
        } catch (const Error& e) {
            EXPECT_NE(std::string(e.what()).find("retrograde"), std::string::npos) << e.what();
        }
    };
    refusedAsRetrograde([&] { (void)stateToEquinoctial(retrograde); });
    refusedAsRetrograde([&] { (void)stateToEquinoctialJacobian(retrograde); });
    // ... while a slightly inclined retrograde orbit still converts.
    EXPECT_NO_THROW((void)stateToEquinoctial({
        circular.r,
        Vec3{0.0, -v, 1e-3 * v}
    }));
}
