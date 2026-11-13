/**********************************************************************
 * rootfind.h
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

#include <cmath>
#include <limits>
#include <utility>

#include "astro/core/error.h"

/// \file
/// Scalar root finding and minimization (header-only templates).
namespace astro {

/// Find a root of f in [a, b] by bisection. f(a) and f(b) must differ in sign.
template <class F>
[[nodiscard]] double bisect(F&& f, double a, double b, double tol = 1e-12, int maxIter = 200) {
    double       fa = f(a);
    const double fb = f(b);
    if (fa == 0.0) {
        return a;
    }
    if (fb == 0.0) {
        return b;
    }
    if ((fa > 0.0) == (fb > 0.0)) {
        throw Error("bisect: root not bracketed");
    }
    for (int it = 0; it < maxIter; ++it) {
        const double m  = 0.5 * (a + b);
        const double fm = f(m);
        if (fm == 0.0 || 0.5 * std::abs(b - a) < tol) {
            return m;
        }
        if ((fm > 0.0) == (fa > 0.0)) {
            a  = m;
            fa = fm;
        } else {
            b = m;
        }
    }
    return 0.5 * (a + b);
}

/// Brent's method root finder on [a, b]. f(a) and f(b) must differ in sign.
/// Converges to |x - root| <= tol.
template <class F>
[[nodiscard]] double brentRoot(F&& f, double a, double b, double tol = 1e-12, int maxIter = 200) {
    double fa = f(a), fb = f(b);
    if (fa == 0.0) {
        return a;
    }
    if (fb == 0.0) {
        return b;
    }
    if ((fa > 0.0) == (fb > 0.0)) {
        throw Error("brentRoot: root not bracketed");
    }
    double           c = a, fc = fa, d = b - a, e = d;
    constexpr double eps = std::numeric_limits<double>::epsilon();
    for (int it = 0; it < maxIter; ++it) {
        if ((fb > 0.0) == (fc > 0.0)) {
            c  = a;
            fc = fa;
            d = e = b - a;
        }
        if (std::abs(fc) < std::abs(fb)) {
            a  = b;
            b  = c;
            c  = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        const double tol1 = 2.0 * eps * std::abs(b) + 0.5 * tol;
        const double xm   = 0.5 * (c - b);
        if (std::abs(xm) <= tol1 || fb == 0.0) {
            return b;
        }
        if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb)) {
            double       p, q, r;
            const double s = fb / fa;
            if (a == c) {
                p = 2.0 * xm * s;
                q = 1.0 - s;
            } else {
                q = fa / fc;
                r = fb / fc;
                p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0.0) {
                q = -q;
            }
            p                 = std::abs(p);
            const double min1 = 3.0 * xm * q - std::abs(tol1 * q);
            const double min2 = std::abs(e * q);
            if (2.0 * p < (min1 < min2 ? min1 : min2)) {
                e = d;
                d = p / q;
            } else {
                d = xm;
                e = d;
            }
        } else {
            d = xm;
            e = d;
        }
        a   = b;
        fa  = fb;
        b  += (std::abs(d) > tol1) ? d : (xm > 0.0 ? tol1 : -tol1);
        fb  = f(b);
    }
    throw Error("brentRoot: did not converge");
}

/// Result of a scalar minimization.
struct MinimizeResult {
    double x;  ///< abscissa of the minimum found
    double fx; ///< f(x)
};

/// Brent's method minimizer of f on the interval between a and b, which may be given in either
/// order (no derivatives). Finds a local minimum.
///
/// `tol` is an absolute tolerance on x, in the units of x. The minimum is located to within about
/// tol + 3 sqrt(eps) |x|, where eps is the machine epsilon (sqrt(eps) = 1.5e-8). The second term is a
/// floor that no tolerance lowers: closer to a minimum than that, f differs from its least value by
/// less than its own rounding error. A `tol` below the floor, such as the default, is therefore safe
/// and costs no extra evaluations.
template <class F>
[[nodiscard]] MinimizeResult brentMinimize(F&& f, double a, double b, double tol = 1e-10, int maxIter = 200) {
    constexpr double golden  = 0.3819660112501051;    // (3 - sqrt(5)) / 2
    constexpr double sqrtEps = 1.4901161193847656e-8; // sqrt(machine epsilon) = 2^-26
    // Everything below takes a as the lower end. With a > b the first convergence test, which
    // subtracts half of (b - a), was met at once and the first trial point was returned.
    if (a > b) {
        std::swap(a, b);
    }
    double x = a + golden * (b - a), w = x, v = x;
    double fx = f(x), fw = fx, fv = fx;
    double d = 0.0, e = 0.0;
    for (int it = 0; it < maxIter; ++it) {
        const double xm = 0.5 * (a + b);
        // Brent's tolerance (Brent 1973, "Algorithms for Minimization without Derivatives",
        // chapter 5): a relative part of sqrt(eps), below which f cannot tell two points apart, plus
        // the caller's absolute part. No two evaluations are closer than tol1 and the search ends
        // with x within 2 tol1 of both ends, so the minimum is within 3 tol1 = 3 sqrt(eps) |x| + tol
        // of x; hence the tol / 3. (With `tol` as the relative part and no floor, a tolerance of
        // 1e-10 spent three to six times the evaluations moving about inside the rounding error of f.)
        const double tol1 = sqrtEps * std::abs(x) + tol / 3.0;
        const double tol2 = 2.0 * tol1;
        if (std::abs(x - xm) <= tol2 - 0.5 * (b - a)) {
            return {x, fx};
        }
        bool golden_step = true;
        if (std::abs(e) > tol1) {
            const double r = (x - w) * (fx - fv);
            double       q = (x - v) * (fx - fw);
            double       p = (x - v) * q - (x - w) * r;
            q              = 2.0 * (q - r);
            if (q > 0.0) {
                p = -p;
            }
            q                  = std::abs(q);
            const double etemp = e;
            e                  = d;
            if (!(std::abs(p) >= std::abs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x))) {
                d              = p / q;
                const double u = x + d;
                if (u - a < tol2 || b - u < tol2) {
                    d = (xm - x >= 0.0) ? tol1 : -tol1;
                }
                golden_step = false;
            }
        }
        if (golden_step) {
            e = (x >= xm) ? a - x : b - x;
            d = golden * e;
        }
        const double u  = (std::abs(d) >= tol1) ? x + d : x + (d >= 0.0 ? tol1 : -tol1);
        const double fu = f(u);
        if (fu <= fx) {
            if (u >= x) {
                a = x;
            } else {
                b = x;
            }
            v  = w;
            fv = fw;
            w  = x;
            fw = fx;
            x  = u;
            fx = fu;
        } else {
            if (u < x) {
                a = u;
            } else {
                b = u;
            }
            if (fu <= fw || w == x) {
                v  = w;
                fv = fw;
                w  = u;
                fw = fu;
            } else if (fu <= fv || v == x || v == w) {
                v  = u;
                fv = fu;
            }
        }
    }
    return {x, fx};
}

/// Newton-Raphson with a derivative. `fdf(x)` returns std::pair{f(x), f'(x)}.
template <class FdF>
[[nodiscard]] double newtonRoot(FdF&& fdf, double x0, double tol = 1e-14, int maxIter = 50) {
    double x = x0;
    for (int it = 0; it < maxIter; ++it) {
        const auto [fx, dfx] = fdf(x);
        if (dfx == 0.0) {
            throw Error("newtonRoot: zero derivative");
        }
        const double dx  = fx / dfx;
        x               -= dx;
        if (std::abs(dx) <= tol * (1.0 + std::abs(x))) {
            return x;
        }
    }
    throw Error("newtonRoot: did not converge");
}

} // namespace astro
