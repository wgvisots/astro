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
