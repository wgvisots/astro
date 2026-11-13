/**********************************************************************
 * lambert.cpp
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

// Lambert's problem by the method of Izzo (2015), "Revisiting Lambert's problem", Celest. Mech.
// Dyn. Astr. 121, 1-15 (arXiv:1403.2705). Written from the paper; equation numbers below are the
// paper's. Where the paper leaves a choice open, the choice made here is explained at the place:
//  - T(x) is evaluated from eq. 18 with psi from an arctangent or an inverse hyperbolic sine of
//    eq. 17, and near the parabola from Battin's series (eq. 20);
//  - near the parabola the derivatives come from the series too (derived here), because eq. 22
//    divides a vanishing numerator by 1 - x^2 there;
//  - y -+ lambda x and x -+ lambda y are formed without cancellation, and y uses 1 - lambda^2 = c / s
//    from the geometry, as the paper suggests in section 5;
//  - the iteration runs to convergence, every root is checked against the time of flight and
//    against the side of the curve it belongs to, and a bracketed search stands behind it;
//  - the minimum of a multi-revolution curve, which decides whether its solutions exist, is
//    verified to be one (Halley's iteration from x = 0 leaves the domain for lambda near -1).

#include "astro/orbit/two_body.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <vector>

#include "astro/core/constants.h"
#include "astro/core/error.h"
#include "astro/core/rootfind.h"

namespace astro {

namespace {

using constants::pi;

/// Gauss's hypergeometric function 2F1(a, b; c; z) by its series. Used for |z| <= kSeriesLimit,
/// where at most ~25 terms reach rounding.
[[nodiscard]] double hypergeometric(double a, double b, double c, double z) noexcept {
    constexpr double eps  = std::numeric_limits<double>::epsilon();
    double           term = 1.0, sum = 1.0;
    for (int n = 0; n < 200; ++n) {
        term *= (a + n) * (b + n) / ((c + n) * (n + 1.0)) * z;
        sum  += term;
        if (std::abs(term) <= eps * std::abs(sum)) {
            break;
        }
    }
    return sum;
}

/// |S1| below which the single-revolution time of flight is taken from Battin's series (eq. 20).
/// S1 = sin^2(psi / 2) (minus sinh^2 for a hyperbola), so the series covers psi up to ~0.64 rad.
/// Eq. 18 loses about epsilon / S1 there, and the series needs more terms as |S1| grows; at 0.1
/// both are at rounding level (measured against extended precision).
constexpr double kSeriesLimit = 0.1;

/// The non-dimensional time of flight and its first three derivatives with respect to x.
struct Flight {
    double T{0.0}, d1{0.0}, d2{0.0}, d3{0.0};
};

/// The time-of-flight curves T(x; M) of one geometry (eq. 18): x is the Lancaster-Blanchard
/// variable (x < 1 ellipse, x = 1 parabola, x > 1 hyperbola), M the number of revolutions.
class FlightCurve {
public:
    /// lambda (eq. after 13: negative for a transfer angle above 180 deg) and c / s = 1 - lambda^2.
    FlightCurve(double lambda, double chordOverS) noexcept
        : l_(lambda)
        , q_(chordOverS) {}

    [[nodiscard]] double lambda() const noexcept { return l_; }

    /// y = sqrt(1 - lambda^2 (1 - x^2)), written with c / s so that it keeps its precision when
    /// lambda^2 is close to 1.
    [[nodiscard]] double y(double x) const noexcept { return std::sqrt(q_ + l_ * l_ * x * x); }

    /// @{
    /// y -+ lambda x and x -+ lambda y, which enter the time of flight and the velocities. Whichever
    /// of each pair is a difference of like-signed terms cancels when lambda^2 is close to 1 or x is
    /// large (a fast hyperbola); it is formed from the other and their product instead:
    ///   (y - lambda x)(y + lambda x) = c / s,
    ///   (x - lambda y)(x + lambda y) = (c / s)(x^2 (1 + lambda^2) - lambda^2).
    [[nodiscard]] double yMinusLambdaX(double x, double yv) const noexcept {
        const double lx = l_ * x;
        return lx > 0.0 ? q_ / (yv + lx) : yv - lx;
    }
    [[nodiscard]] double yPlusLambdaX(double x, double yv) const noexcept {
        const double lx = l_ * x;
        return lx < 0.0 ? q_ / (yv - lx) : yv + lx;
    }
    [[nodiscard]] double xMinusLambdaY(double x, double yv) const noexcept {
        const double ly = l_ * yv;
        return (x > 0.0 && ly > 0.0) || (x < 0.0 && ly < 0.0) ? crossProduct(x) / (x + ly) : x - ly;
    }
    [[nodiscard]] double xPlusLambdaY(double x, double yv) const noexcept {
        const double ly = l_ * yv;
        return (x > 0.0 && ly < 0.0) || (x < 0.0 && ly > 0.0) ? crossProduct(x) / (x - ly) : x + ly;
    }
    /// @}

    /// T(0; 0) (eq. 19) and T(1; 0) (eq. 21).
    [[nodiscard]] double T00() const noexcept { return std::acos(l_) + l_ * std::sqrt(q_); }
    [[nodiscard]] double T1() const noexcept { return 2.0 / 3.0 * (1.0 - l_ * l_ * l_); }

    /// T and its derivatives at x (x > -1; x < 1 for M > 0).
    [[nodiscard]] Flight at(double x, int M) const noexcept {
        const double yv  = y(x);
        const double eta = yMinusLambdaX(x, yv);
        if (M == 0) {
            const double S1 = 0.5 * (1.0 - l_ - x * eta);
            if (std::abs(S1) < kSeriesLimit) {
                return nearParabola(x, yv, eta, S1);
            }
        }
        // 1 - x^2 = a_min / a: positive on an ellipse, negative on a hyperbola. As a product it keeps
        // its relative precision near x = +-1, where 1 + x or 1 - x is exact.
        const double u = (1.0 - x) * (1.0 + x);
        Flight       f;
        if (u > 0.0) { // eq. 17: sin psi and cos psi, psi in [0, pi]
            const double root = std::sqrt(u);
            const double psi  = std::atan2(eta * root, x * yv + l_ * u);
            f.T               = ((psi + M * pi) / root - xMinusLambdaY(x, yv)) / u;
        } else { // eq. 17: sinh psi
            const double root = std::sqrt(-u);
            const double psi  = std::asinh(eta * root);
            f.T               = (xMinusLambdaY(x, yv) - psi / root) / (-u);
        }
        // Eq. 22.
        const double l2 = l_ * l_, l3 = l2 * l_, y2 = yv * yv;
        f.d1 = (3.0 * f.T * x - 2.0 + 2.0 * l3 * x / yv) / u;
        f.d2 = (3.0 * f.T + 5.0 * x * f.d1 + 2.0 * q_ * l3 / (y2 * yv)) / u;
        f.d3 = (7.0 * x * f.d2 + 8.0 * f.d1 - 6.0 * q_ * l2 * l3 * x / (y2 * y2 * yv)) / u;
        return f;
    }

private:
    /// Eq. 20: 2 T = eta^3 Q + 4 lambda eta, Q = 4/3 F(S1), F(z) = 2F1(3, 1; 5/2; z), with
    /// eta = y - lambda x and S1 = (1 - lambda - x eta) / 2. Its derivatives follow from the chain
    /// rule with
    ///   y' = lambda^2 x / y,  y'' = lambda^2 (1 - lambda^2) / y^3,  y''' = -3 lambda^4 (1 - lambda^2) x / y^5,
    ///   eta' = y' - lambda,   eta'' = y'',   eta''' = y''',
    ///   S1' = -(eta + x eta') / 2,  S1'' = -(2 eta' + x eta'') / 2,  S1''' = -(3 eta'' + x eta''') / 2,
    /// and the derivative of a hypergeometric function, d/dz 2F1(a, b; c; z) = (a b / c) 2F1(a+1, b+1; c+1; z).
    [[nodiscard]] Flight nearParabola(double x, double yv, double eta, double S1) const noexcept {
        const double l2 = l_ * l_, y2 = yv * yv;
        const double yy = l2 * q_ / (y2 * yv);
        const double y3 = -3.0 * l2 * l2 * q_ * x / (y2 * y2 * yv);
        const double e1 = -l_ * eta / yv, e2 = yy, e3 = y3; // y' - lambda = -lambda (y - lambda x) / y
        const double s1 = -0.5 * (eta + x * e1), s2 = -0.5 * (2.0 * e1 + x * e2), s3 = -0.5 * (3.0 * e2 + x * e3);

        const double F0 = hypergeometric(3.0, 1.0, 2.5, S1);
        const double F1 = 6.0 / 5.0 * hypergeometric(4.0, 2.0, 3.5, S1);                           // F'
        const double F2 = 6.0 / 5.0 * 16.0 / 7.0 * hypergeometric(5.0, 3.0, 4.5, S1);              // F''
        const double F3 = 6.0 / 5.0 * 16.0 / 7.0 * 10.0 / 3.0 * hypergeometric(6.0, 4.0, 5.5, S1); // F'''
        const double Q  = 4.0 / 3.0 * F0;
        const double Q1 = 4.0 / 3.0 * F1 * s1;
        const double Q2 = 4.0 / 3.0 * (F2 * s1 * s1 + F1 * s2);
        const double Q3 = 4.0 / 3.0 * (F3 * s1 * s1 * s1 + 3.0 * F2 * s1 * s2 + F1 * s3);

        const double n2 = eta * eta, n3 = n2 * eta;
        Flight       f;
        f.T  = 0.5 * (n3 * Q + 4.0 * l_ * eta);
        f.d1 = 0.5 * (3.0 * n2 * e1 * Q + n3 * Q1 + 4.0 * l_ * e1);
        f.d2 = 0.5 * (6.0 * eta * e1 * e1 * Q + 3.0 * n2 * e2 * Q + 6.0 * n2 * e1 * Q1 + n3 * Q2 + 4.0 * l_ * e2);
        f.d3 = 0.5
             * (6.0 * e1 * e1 * e1 * Q + 18.0 * eta * e1 * e2 * Q + 18.0 * eta * e1 * e1 * Q1 + 3.0 * n2 * e3 * Q
                + 9.0 * n2 * e2 * Q1 + 9.0 * n2 * e1 * Q2 + n3 * Q3 + 4.0 * l_ * e3);
        return f;
    }

    /// x^2 - lambda^2 y^2.
    [[nodiscard]] double crossProduct(double x) const noexcept { return q_ * (x * x * (1.0 + l_ * l_) - l_ * l_); }

    double l_;
    double q_;
};

/// Distance kept from the ends of the domain, x = -1 (and x = 1 for M > 0), where T is infinite.
constexpr double kEdge = 1e-13;

/// x moved into the open domain of the curve of M revolutions.
[[nodiscard]] double intoDomain(double x, int M) noexcept {
    x = std::max(x, -1.0 + kEdge);
    return M > 0 ? std::min(x, 1.0 - kEdge) : x;
}

/// The minimum of the curve of M >= 1 revolutions: x and T there. Halley's iteration on dT/dx = 0
/// from x = 0, as in the paper's Algorithm 2. For lambda near -1 that iteration runs out of the
/// domain (the curve bends sharply near x = 0 there), so its result is used only if it is a
/// stationary point with positive curvature inside the domain; otherwise the minimum is bracketed
/// by a scan and refined with Brent's minimizer.
[[nodiscard]] std::pair<double, double> minimumOf(const FlightCurve& curve, int M) {
    double x = 0.0;
    for (int it = 0; it < 30; ++it) {
        const Flight f = curve.at(x, M);
        if (f.d2 > 0.0 && std::abs(f.d1) <= 1e-10 * f.d2) { // the Newton step to the minimum is negligible
            return {x, f.T};
        }
        const double denom = 2.0 * f.d2 * f.d2 - f.d1 * f.d3;
        const double next  = x - 2.0 * f.d1 * f.d2 / denom;
        if (!(next > -1.0 + 1e-6 && next < 1.0 - 1e-6)) { // also catches a zero or non-finite denominator
            break;
        }
        x = next;
    }
    constexpr int kScan = 64;
    const auto    T     = [&](double xx) { return curve.at(xx, M).T; };
    double        best = 0.0, bestT = T(0.0);
    for (int k = 1; k < kScan; ++k) {
        const double xx = -1.0 + 2.0 * k / kScan;
        if (const double t = T(xx); t < bestT) {
            best  = xx;
            bestT = t;
        }
    }
    const MinimizeResult m = brentMinimize(T, best - 2.0 / kScan, best + 2.0 / kScan, 1.5e-8);
    return {m.x, m.fx};
}

/// Which root of T(x; M) = target is wanted.
enum class Branch {
    Single, ///< M = 0: the curve falls from infinity at x = -1 to zero
    Left,   ///< M > 0, left of the minimum: the curve is falling
    Right,  ///< M > 0, right of the minimum: the curve is rising
};

/// True if the curve value f at an x inside the domain solves T = target on the wanted branch.
[[nodiscard]] bool accepted(const Flight& f, double target, Branch branch) noexcept {
    if (!(std::abs(f.T - target) <= 1e-10 * target)) {
        return false;
    }
    // Near the minimum the slope is close to zero and may carry either sign.
    const double slack = 1e-6 * target;
    return branch == Branch::Single || (branch == Branch::Left ? f.d1 <= slack : f.d1 >= -slack);
}

/// Householder's iteration (section 4.1) on T(x; M) - target from x0. It stops when the residual
/// is at rounding level, or when the step is negligible against the distance to the nearer end of
/// the domain (the iteration is of fourth order, so a step of 1e-9 of that distance leaves an error
/// far below rounding); the value at the final x is then checked by `accepted`. A step that would
/// leave the domain goes half-way to the end instead.
[[nodiscard]] std::optional<double>
householder(const FlightCurve& curve, double x0, int M, double target, Branch branch) noexcept {
    double x    = intoDomain(x0, M);
    double step = std::numeric_limits<double>::infinity();
    for (int it = 0; it < 60; ++it) {
        const Flight f    = curve.at(x, M);
        const double g    = f.T - target;
        const double left = 1.0 + x, right = 1.0 - x; // distances to the ends
        const double scale = M > 0 ? std::min(left, right) : std::min(left, 1.0 + std::abs(x));
        if (std::abs(g) <= 4e-15 * target || step <= 1e-9 * scale) {
            return accepted(f, target, branch) ? std::optional<double>(x) : std::nullopt;
        }
        const double d1sq  = f.d1 * f.d1;
        const double denom = f.d1 * (d1sq - g * f.d2) + f.d3 * g * g / 6.0;
        double       dx    = -g * (d1sq - 0.5 * g * f.d2) / denom;
        if (!std::isfinite(dx)) {
            return std::nullopt;
        }
        if (dx <= -left) {
            dx = -0.5 * left;
        } else if (M > 0 && dx >= right) {
            dx = 0.5 * right;
        }
        step = std::abs(dx);
        x    = intoDomain(x + dx, M);
    }
    return std::nullopt;
}

/// The root on `branch` by bracketing and Brent's method, for when Householder's iteration does not
/// deliver an accepted root. In the accuracy study that happens only for transfer angles within
/// ~1e-2 rad of 0 or 360 deg (about 3 solutions in 10 000 there), where lambda^2 is close to 1 and
/// the curves bend sharply near x = 0. nullopt if the branch has no root.
[[nodiscard]] std::optional<double> bracketed(const FlightCurve& curve, int M, double target, Branch branch) {
    const auto g  = [&](double x) { return curve.at(x, M).T - target; };
    double     lo = -1.0 + kEdge, hi = 1.0 - kEdge;
    if (branch == Branch::Single) {
        hi = 1.0;
        while (g(hi) > 0.0) { // T falls to zero as x grows
            hi *= 2.0;
            if (!(hi < 1e300)) {
                return std::nullopt;
            }
        }
    } else {
        const double xMin                  = minimumOf(curve, M).first;
        (branch == Branch::Left ? hi : lo) = xMin;
    }
    const double glo = g(lo), ghi = g(hi);
    if (!std::isfinite(glo) || !std::isfinite(ghi) || (glo > 0.0) == (ghi > 0.0)) {
        return glo == 0.0 ? std::optional<double>(lo) : ghi == 0.0 ? std::optional<double>(hi) : std::nullopt;
    }
    return brentRoot(g, lo, hi, 1e-15);
}

/// The root of T(x; M) = target on `branch`, starting Householder's iteration from x0.
[[nodiscard]] std::optional<double> solve(const FlightCurve& curve, double x0, int M, double target, Branch branch) {
    if (const std::optional<double> x = householder(curve, x0, M, target, branch)) {
        return x;
    }
    return bracketed(curve, M, target, branch);
}

/// The starter of the single-revolution case (eq. 30, with the improved middle line of section
/// 4.2.1): the time-of-flight curve is close to straight lines in the plane of xi = log(1 + x) and
/// tau = log T, through (x = 0, T00) and (x = 1, T1) with the slopes of its asymptotes.
[[nodiscard]] double singleRevolutionStarter(const FlightCurve& curve, double T) noexcept {
    const double T0 = curve.T00(), T1 = curve.T1(), l = curve.lambda();
    if (T >= T0) {
        return std::pow(T0 / T, 2.0 / 3.0) - 1.0;
    }
    if (T < T1) {
        const double l5 = l * l * l * l * l;
        return 2.5 * T1 * (T1 - T) / (T * (1.0 - l5)) + 1.0;
    }
    // Between the two: xi0 = log 2 (tau - tau0) / (tau1 - tau0), so 1 + x0 = 2^((tau - tau0) / (tau1 - tau0)).
    return std::exp2(std::log(T / T0) / std::log(T1 / T0)) - 1.0;
}

} // namespace

std::vector<LambertSolution>
lambertIzzo(const Vec3& r1, const Vec3& r2, double tof, double mu, bool prograde, int maxRevolutions) {
    if (!(tof > 0.0) || !std::isfinite(tof)) {
        throw Error("lambertIzzo: time of flight must be positive and finite");
    }
    if (!(mu > 0.0) || !std::isfinite(mu)) {
        throw Error("lambertIzzo: gravitational parameter must be positive and finite");
    }
    const double R1 = norm(r1), R2 = norm(r2);
    if (!(R1 > 0.0) || !(R2 > 0.0) || !std::isfinite(R1) || !std::isfinite(R2)) {
        throw Error("lambertIzzo: r1 and r2 must be finite and non-zero");
    }

    // Geometry (Algorithm 1). The unit vector along r1 x r2 is the orbit normal of the transfer
    // through less than 180 deg; the other way round the normal is its opposite and lambda is
    // negative.
    const Vec3   ir1 = r1 / R1, ir2 = r2 / R2;
    const Vec3   h      = cross(ir1, ir2);
    const double sinAng = norm(h);
    if (!(sinAng > 1e-12)) {
        throw Error("lambertIzzo: r1 and r2 are collinear; transfer plane undefined");
    }
    const double c = norm(r2 - r1);
    const double s = 0.5 * (R1 + R2 + c);
    // Prograde: counter-clockwise seen from +z. That is the long way round when r1 x r2 points down.
    // If the plane contains the z axis the z component of r1 x r2 is zero, or just its own rounding
    // error (a few units in the last place of its two products): its sign then says nothing, and the
    // shorter way round is taken whichever direction was asked for.
    const double zNoise =
        4.0 * std::numeric_limits<double>::epsilon() * (std::abs(ir1.x * ir2.y) + std::abs(ir1.y * ir2.x));
    const bool   polar   = std::abs(h.z) <= zNoise;
    const bool   longWay = !polar && (prograde ? h.z < 0.0 : h.z > 0.0);
    const Vec3   normal  = (longWay ? -1.0 : 1.0) / sinAng * h;
    const Vec3   it1 = cross(normal, ir1), it2 = cross(normal, ir2);
    const double chordOverS = c / s; // 1 - lambda^2
    const double lambda     = (longWay ? -1.0 : 1.0) * std::sqrt(std::max(0.0, 1.0 - chordOverS));
    const double T          = std::sqrt(2.0 * mu / (s * s * s)) * tof;

    const FlightCurve curve(lambda, chordOverS);

    // The largest number of revolutions with a solution (Algorithm 2): T(0; M) = T00 + M pi, and the
    // minimum of the curve of M revolutions is below that.
    int          maxM  = 0;
    const double turns = std::floor(T / pi);
    if (maxRevolutions > 0 && turns >= 1.0) {
        if (turns > maxRevolutions) {
            maxM = maxRevolutions; // T >= (maxM + 1) pi >= T00 + maxM pi: a solution exists
        } else {
            maxM = static_cast<int>(turns);
            if (T < curve.T00() + maxM * pi && minimumOf(curve, maxM).second > T) {
                --maxM;
            }
        }
    }

    // The roots: x and the solution they belong to.
    struct Root {
        double x;
        int    revolutions;
        bool   right;
    };
    std::vector<Root> roots;
    const auto        keep = [&](std::optional<double> x, int M, bool right) {
        if (!x) {
            throw Error("lambertIzzo: no solution found for " + std::to_string(M) + " revolutions");
        }
        roots.push_back({*x, M, right});
    };
    keep(solve(curve, singleRevolutionStarter(curve, T), 0, T, Branch::Single), 0, false);
    for (int M = 1; M <= maxM; ++M) {
        // Eq. 31: the starters from the two asymptotes of the curve in the (xi, tau) plane.
        const double left  = std::pow((M * pi + pi) / (8.0 * T), 2.0 / 3.0);
        const double right = std::pow(8.0 * T / (M * pi), 2.0 / 3.0);
        keep(solve(curve, (left - 1.0) / (left + 1.0), M, T, Branch::Left), M, false);
        keep(solve(curve, (right - 1.0) / (right + 1.0), M, T, Branch::Right), M, true);
    }

    // Terminal velocities (Algorithm 1): radial and transverse components at both ends.
    const double                 gamma = std::sqrt(0.5 * mu * s);
    const double                 rho   = (R1 - R2) / c;
    const double                 sigma = std::sqrt(std::max(0.0, (1.0 - rho) * (1.0 + rho)));
    std::vector<LambertSolution> out;
    out.reserve(roots.size());
    for (const Root& root : roots) {
        const double    x = root.x, y = curve.y(x);
        const double    a = -curve.xMinusLambdaY(x, y), b = curve.xPlusLambdaY(x, y); // lambda y -+ x
        const double    vr1 = gamma * (a - rho * b) / R1;
        const double    vr2 = -gamma * (a + rho * b) / R2;
        const double    vt  = gamma * sigma * curve.yPlusLambdaX(x, y);
        LambertSolution sol;
        sol.v1          = vr1 * ir1 + (vt / R1) * it1;
        sol.v2          = vr2 * ir2 + (vt / R2) * it2;
        sol.revolutions = root.revolutions;
        sol.rightBranch = root.right;
        out.push_back(sol);
    }
    return out;
}

LambertSolution lambert(const Vec3& r1, const Vec3& r2, double tof, double mu, bool prograde) {
    return lambertIzzo(r1, r2, tof, mu, prograde, 0).front();
}

} // namespace astro
