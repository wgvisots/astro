/**********************************************************************
 * kepler.cpp
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

#include "astro/orbit/kepler.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>

#include "astro/core/constants.h"
#include "astro/core/error.h"

namespace astro {

namespace {

// Newton's method ends when a step is below this fraction of the iterate. It converges
// quadratically, so what is left after such a step is of the order of the tolerance squared, below
// rounding (the solvers give the bound).
constexpr double kStepTolerance = 1e-9;

// Newton's method needs at most 4 steps here (measured over 12 million eccentricities and mean
// anomalies, including e within 1e-16 of 1); the limit only guards against a defect.
constexpr int kMaxIterations = 50;

// Below this argument x - sin x, 1 - cos x, sinh x - x and cosh x - 1 are summed as Taylor series;
// from it on the differences are formed directly (they lose about one bit there).
constexpr double kSeriesLimit = 2.0;

// Ratios of successive terms of those series, without the power of x: 1 / (2k (2k + 1)) from the
// term x^3/3! on, and 1 / ((2k - 1) 2k) from the term x^2/2! on. They end where the first term left
// out is below 2e-18 of the sum for x < 2 (x^25/25! and x^26/26!).
// clang-format off
constexpr double kOddRatios[]  = {1.0 /  20.0, 1.0 /  42.0, 1.0 /  72.0, 1.0 / 110.0, 1.0 / 156.0, 1.0 / 210.0,
                                  1.0 / 272.0, 1.0 / 342.0, 1.0 / 420.0, 1.0 / 506.0};
constexpr double kEvenRatios[] = {1.0 /  12.0, 1.0 /  30.0, 1.0 /  56.0, 1.0 /  90.0, 1.0 / 132.0, 1.0 / 182.0,
                                  1.0 / 240.0, 1.0 / 306.0, 1.0 / 380.0, 1.0 / 462.0, 1.0 / 552.0};
// clang-format on

/// 1 + z r0 (1 + z r1 (1 + ...)): a series whose successive terms have the ratios z r0, z r1, ...
template <std::size_t N>
inline double nestedSum(double z, const double (&ratios)[N]) noexcept {
    double sum = 1.0;
    for (std::size_t k = N; k-- > 0;) {
        sum = 1.0 + z * ratios[k] * sum;
    }
    return sum;
}

/// x - sin x and 1 - cos x for x >= 0, accurate to rounding. Formed directly they have absolute
/// errors of eps x and eps, which is all there is of them for small x.
inline void sinRemainders(double x, double& xMinusSin, double& oneMinusCos) noexcept {
    if (x < kSeriesLimit) {
        const double z = -x * x;
        xMinusSin      = x * x * x / 6.0 * nestedSum(z, kOddRatios);
        oneMinusCos    = 0.5 * x * x * nestedSum(z, kEvenRatios);
    } else {
        xMinusSin   = x - std::sin(x);
        oneMinusCos = 1.0 - std::cos(x);
    }
}

/// sinh x - x and cosh x - 1 for x >= 0, accurate to rounding (see sinRemainders).
inline void sinhRemainders(double x, double& sinhMinusX, double& coshMinusOne) noexcept {
    if (x < kSeriesLimit) {
        const double z = x * x;
        sinhMinusX     = x * x * x / 6.0 * nestedSum(z, kOddRatios);
        coshMinusOne   = 0.5 * x * x * nestedSum(z, kEvenRatios);
    } else {
        sinhMinusX   = std::sinh(x) - x;
        coshMinusOne = std::cosh(x) - 1.0;
    }
}

/// The positive root of c3 x^3 + c1 x = m, for c1 > 0, c3 >= 0 and m >= 0: Kepler's equation with
/// the sine, or the hyperbolic sine, cut off after its cubic term. Used for starting values.
inline double cubicRoot(double c1, double c3, double m) noexcept {
    // With x = (m / c1) y the equation is kappa y^3 + y = 1, kappa = c3 m^2 / c1^3.
    const double x1    = m / c1;
    const double kappa = c3 * x1 * x1 / c1;
    if (!(kappa > 1e-100)) {
        return x1; // the cubic term does not count (and 1 / kappa^3 below would overflow)
    }
    // Cardano: y = A - B with A B = P / 3 and A^3 - B^3 = P, where P = 1 / kappa. Written as
    // P / (A^2 + A B + B^2), a sum of positive terms, because A - B cancels for small kappa.
    const double P = 1.0 / kappa;
    const double A = std::cbrt(0.5 * P + std::sqrt(0.25 * P * P + P * P * P / 27.0));
    const double B = P / (3.0 * A);
    return x1 * P / (A * A + P / 3.0 + B * B);
}

} // namespace

double trueToEccentricAnomaly(double nu, double e) noexcept {
    return std::atan2(std::sqrt(1.0 - e * e) * std::sin(nu), e + std::cos(nu));
}

double eccentricToTrueAnomaly(double E, double e) noexcept {
    return std::atan2(std::sqrt(1.0 - e * e) * std::sin(E), std::cos(E) - e);
}

double eccentricToMeanAnomaly(double E, double e) noexcept {
    return E - e * std::sin(E);
}

double meanToEccentricAnomaly(double M, double e) {
    if (!(e >= 0.0 && e < 1.0)) {
        throw Error("meanToEccentricAnomaly: requires 0 <= e < 1");
    }
    if (!std::isfinite(M)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    // Solve within half a revolution of perigee, |Mw| <= pi, for |Mw| (the equation is odd), then
    // restore the sign and the revolutions. std::remainder is exact and leaves an M that is already
    // in range untouched; adding and subtracting pi instead costs 4e-16 of absolute error in Mw,
    // which 1 / (1 - e) then multiplies in E.
    const double Mw = std::remainder(M, constants::twoPi);
    const double m  = std::abs(Mw);

    // f(E) = E - e sin E - m is evaluated as (1 - e) E + e (E - sin E) - m, and its derivative
    // 1 - e cos E as (1 - e) + e (1 - cos E): sums of positive terms, good to rounding where the
    // plain forms cancel (e near 1, E small).
    //
    // f increases and is convex on [0, pi] (f'' = e sin E >= 0). A Newton step from any point of
    // [0, pi] therefore lands on or above the solution E*, and from there the iterates come down
    // to E* without crossing it: the iteration cannot fail. A step d from above leaves an error of
    // at most d^2 f'' / (2 f'(E*)) <= 2.5 d^2 / E* (by sin x <= x and 1 - cos x >= 2 x^2 / pi^2), so
    // after a step below kStepTolerance the iterate is E* to rounding. A step that is zero or
    // negative means that rounding has taken over, which is the same.
    //
    // Starting value, a lower bound of E*: m itself (E = m + e sin E >= m), or the root of
    // (1 - e) E + e E^3/6 = m (because sin E >= E - E^3/6), whichever is larger. The cubic's root is
    // close where a start is hard to find, near perigee for e near 1.
    const double q = 1.0 - e;
    double       E = m * m > 6.0 ? m : cubicRoot(q, e / 6.0, m);
    for (int it = 0; it < kMaxIterations; ++it) {
        double xMinusSin = 0.0, oneMinusCos = 0.0;
        sinRemainders(E, xMinusSin, oneMinusCos);
        const double step = (q * E + e * xMinusSin - m) / (q + e * oneMinusCos);
        E                 = std::min(E - step, constants::pi); // f(pi) >= 0: E* is not beyond pi
        // The first step starts below E*, where a small step would prove nothing: always take it.
        if (it > 0 && !(step > kStepTolerance * E)) {
            return std::copysign(E, Mw) + (M - Mw);
        }
    }
    throw Error("meanToEccentricAnomaly: did not converge");
}

double trueToHyperbolicAnomaly(double nu, double e) noexcept {
    return 2.0 * std::atanh(std::sqrt((e - 1.0) / (e + 1.0)) * std::tan(0.5 * nu));
}

double hyperbolicToTrueAnomaly(double F, double e) noexcept {
    return 2.0 * std::atan(std::sqrt((e + 1.0) / (e - 1.0)) * std::tanh(0.5 * F));
}

double hyperbolicToMeanAnomaly(double F, double e) noexcept {
    return e * std::sinh(F) - F;
}

double meanToHyperbolicAnomaly(double M, double e) {
    if (!(e > 1.0) || std::isinf(e)) {
        throw Error("meanToHyperbolicAnomaly: requires a finite e > 1");
    }
    if (std::isnan(M)) {
        return M;
    }
    // Solve for |M| (the equation is odd) and restore the sign.
    const double m = std::abs(M);
    if (m >= 0x1p64) {
        // F is at most 710 (the largest argument with a finite sinh), less than half a unit in the
        // last place of M here: e sinh F = M + F is M to rounding. Infinite M gives infinite F.
        return std::copysign(std::asinh(m / e), M);
    }

    // g(F) = e sinh F - F - m is evaluated as (e - 1) F + e (sinh F - F) - m, and its derivative
    // e cosh F - 1 as (e - 1) + e (cosh F - 1): sums of positive terms, good to rounding where the
    // plain forms cancel (e near 1, F small).
    //
    // g increases and is convex for F >= 0, so Newton's method comes down to the solution F* from
    // above without crossing it, as in meanToEccentricAnomaly: it cannot fail. Here a step d leaves
    // an error of at most d^2 g'' / (2 g'(F*)) <= (F*/2) coth(F*/2) d^2 / F*, and F* < 46, so again
    // the iterate after a step below kStepTolerance is F* to rounding.
    //
    // Starting value, an upper bound of F*: the root Fc of (e - 1) F + e F^3/6 = m is one (because
    // sinh F >= F + F^3/6), and then so is asinh((m + Fc) / e), because
    // sinh F* = (m + F*) / e <= (m + Fc) / e. The usual start m / (e - 1) is unbounded as e -> 1,
    // and Newton's method then comes back by only one unit of F per step.
    const double q = e - 1.0;
    double       F = std::asinh((m + cubicRoot(q, e / 6.0, m)) / e);
    for (int it = 0; it < kMaxIterations; ++it) {
        double sinhMinusF = 0.0, coshMinusOne = 0.0;
        sinhRemainders(F, sinhMinusF, coshMinusOne);
        const double step  = (q * F + e * sinhMinusF - m) / (q + e * coshMinusOne);
        F                 -= step;
        // Always take the first step: rounding may have put the start just below F*.
        if (it > 0 && !(step > kStepTolerance * F)) {
            return std::copysign(F, M);
        }
    }
    throw Error("meanToHyperbolicAnomaly: did not converge");
}

double trueToParabolicAnomaly(double nu) noexcept {
    return std::tan(0.5 * nu);
}
double parabolicToTrueAnomaly(double D) noexcept {
    return 2.0 * std::atan(D);
}
double parabolicToMeanAnomaly(double D) noexcept {
    return D + D * D * D / 3.0;
}

double meanToParabolicAnomaly(double M) noexcept {
    // Barker's equation D + D^3/3 = M is odd: solve for |M| and restore the sign. With
    // w^3 = 3|M|/2 + sqrt(1 + 9 M^2/4) the solution is D = w - 1/w. For M < 0 the cube itself
    // cancelled (to zero from M = -1e8 on, giving -infinity), and w - 1/w cancels for small |M|.
    const double m = std::abs(M);
    if (m >= 0x1p100) {
        // D = cbrt(3 |M|) to rounding (1/w is below 1e-20 of w), scaled so that 3 |M| cannot
        // overflow. Infinite M gives infinite D.
        return std::copysign(2.0 * std::cbrt(0.375 * m), M);
    }
    // From w^3 - 1/w^3 = 3|M|: D = 3|M| / (w^2 + 1 + 1/w^2), a quotient of positive terms.
    const double a  = 1.5 * m;
    const double w  = std::cbrt(a + std::sqrt(1.0 + a * a));
    const double w2 = w * w;
    double       D  = 3.0 * m / (w2 + 1.0 + 1.0 / w2);
    // One Newton step takes out the rounding of the cube root (a few units in the last place).
    D -= (D + D * D * D / 3.0 - m) / (1.0 + D * D);
    return std::copysign(D, M);
}

double meanToTrueAnomaly(double M, double e) {
    if (e < 1.0) {
        return eccentricToTrueAnomaly(meanToEccentricAnomaly(M, e), e);
    }
    if (e == 1.0) {
        return parabolicToTrueAnomaly(meanToParabolicAnomaly(M));
    }
    return hyperbolicToTrueAnomaly(meanToHyperbolicAnomaly(M, e), e);
}

double trueToMeanAnomaly(double nu, double e) noexcept {
    if (e < 1.0) {
        return eccentricToMeanAnomaly(trueToEccentricAnomaly(nu, e), e);
    }
    if (e == 1.0) {
        return parabolicToMeanAnomaly(trueToParabolicAnomaly(nu));
    }
    return hyperbolicToMeanAnomaly(trueToHyperbolicAnomaly(nu, e), e);
}

} // namespace astro
