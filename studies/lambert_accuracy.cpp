/**********************************************************************
 * lambert_accuracy.cpp
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

// Lambert accuracy study (docs/accuracy.md): lambertIzzo on random problems in ordinary and extreme
// regimes. Every solution is judged by propagating the departure state over the time of flight in
// extended precision (universal variables, long double) and comparing with the arrival point: the
// "closure" |r(t) - r2| / |r2|.
//
// Closure measures the problem as well as the solver: where the arrival point is very sensitive to
// the departure velocity, a velocity that is right to the last bit still closes poorly. So each
// closure is also compared with what a change of the velocity by one unit in the last place does:
// a solution is "beyond conditioning" if it closes more than 20 times worse than that (plus 1e-12).
//
// A fast hyperbola that swings past the centre (the long way round in a very short time) cannot be
// checked this way at all: errors grow along it by the exponential of its hyperbolic anomaly, e^50
// and more, which no arithmetic here survives. Those solutions are counted as "not checkable".
//
//   astro_study_lambert [problems per regime]   (default 10000)

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <optional>
#include <random>
#include <vector>

#include "astro/core/constants.h"
#include "astro/orbit/two_body.h"

using namespace astro;

namespace {

using Real = long double;

struct Vec3L {
    Real x, y, z;
};
Vec3L operator+(const Vec3L& a, const Vec3L& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
Vec3L operator-(const Vec3L& a, const Vec3L& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Vec3L operator*(Real s, const Vec3L& a) {
    return {s * a.x, s * a.y, s * a.z};
}
Real dot(const Vec3L& a, const Vec3L& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
Real norm(const Vec3L& a) {
    return std::sqrt(dot(a, a));
}
Vec3L lift(const Vec3& a) {
    return {a.x, a.y, a.z};
}

/// The universal functions U0..U3 at (chi, alpha): U_k = chi^k c_k(alpha chi^2).
void universal(Real chi, Real alpha, Real& U0, Real& U1, Real& U2, Real& U3) {
    const Real z = alpha * chi * chi;
    Real       c0, c1, c2 = 0, c3 = 0;
    if (std::abs(z) < 0.5L) { // series for c2 and c3; c0 and c1 from the recurrence
        Real t2 = 0.5L, t3 = 1.0L / 6.0L;
        for (int j = 0; j < 40; ++j) {
            c2 += t2;
            c3 += t3;
            t2 *= -z / ((3 + 2 * j) * (4.0L + 2 * j));
            t3 *= -z / ((4 + 2 * j) * (5.0L + 2 * j));
        }
        c1 = 1 - z * c3;
        c0 = 1 - z * c2;
    } else if (z > 0) {
        const Real q = std::sqrt(z);
        c0           = std::cos(q);
        c1           = std::sin(q) / q;
        c2           = (1 - c0) / z;
        c3           = (1 - c1) / z;
    } else {
        const Real q = std::sqrt(-z);
        c0           = std::cosh(q);
        c1           = std::sinh(q) / q;
        c2           = (1 - c0) / z;
        c3           = (1 - c1) / z;
    }
    U0 = c0;
    U1 = chi * c1;
    U2 = chi * chi * c2;
    U3 = chi * chi * chi * c3;
}

/// Hyperbolic anomaly span beyond which the reference propagation loses more than ~1e-9: it works to
/// ~1e-19 and errors grow by e^span.
constexpr long double kMaxSpan = 23.0L;

/// |r(tof) - r2| / |r2| for the two-body orbit leaving r1 with velocity v1; nullopt if the orbit
/// cannot be propagated accurately enough to say (see kMaxSpan).
std::optional<double> closure(const Vec3& r1, const Vec3& r2, const Vec3& v1, double tof, double mu) {
    const Vec3L r0 = lift(r1), v0 = lift(v1), target = lift(r2);
    const Real  r0m = norm(r0), rootMu = std::sqrt(static_cast<Real>(mu));
    const Real  sigma = dot(r0, v0) / rootMu, alpha = 2 / r0m - dot(v0, v0) / mu, goal = rootMu * tof;
    Real        U0, U1, U2, U3;
    // Kepler's equation in the universal anomaly: f increases with chi (df/dchi = r > 0).
    const auto f = [&](Real chi) {
        universal(chi, alpha, U0, U1, U2, U3);
        return r0m * U1 + sigma * U2 + U3 - goal;
    };
    Real lo = 0, hi = 1e-6L * goal / r0m + 1e-300L;
    for (int guard = 0; f(hi) < 0 && guard < 2000; ++guard) {
        lo  = hi;
        hi *= 2;
    }
    Real chi = 0.5L * (lo + hi);
    for (int it = 0; it < 300; ++it) { // Newton's method inside the bracket
        const Real value      = f(chi);
        (value > 0 ? hi : lo) = chi;
        Real next             = chi - value / (r0m * U0 + sigma * U1 + U2);
        if (!(next > lo && next < hi)) {
            next = 0.5L * (lo + hi);
        }
        const bool done = std::abs(next - chi) <= 1e-19L * std::abs(chi);
        chi             = next;
        if (done) {
            break;
        }
    }
    (void)f(chi);
    if (alpha < 0 && std::sqrt(-alpha) * chi > kMaxSpan) {
        return std::nullopt;
    }
    const Vec3L  r = (1 - U2 / r0m) * r0 + ((r0m * U1 + sigma * U2) / rootMu) * v0;
    const double e = static_cast<double>(norm(r - target) / norm(target));
    return std::isfinite(e) ? e : 1.0;
}

struct Tally {
    long   problems{0}, solutions{0}, nonFinite{0}, unchecked{0}, beyond{0}, above10{0}, above6{0};
    double worst{0.0};
};

void judge(Tally& t, const LambertSolution& s, const Vec3& r1, const Vec3& r2, double tof, double mu) {
    ++t.solutions;
    if (!std::isfinite(s.v1.x + s.v1.y + s.v1.z + s.v2.x + s.v2.y + s.v2.z)) {
        ++t.nonFinite;
        return;
    }
    const std::optional<double> checked = closure(r1, r2, s.v1, tof, mu);
    if (!checked) {
        ++t.unchecked;
        return;
    }
    const double e  = *checked;
    t.worst         = std::max(t.worst, e);
    t.above10      += e > 1e-10;
    t.above6       += e > 1e-6;
    double floor    = 0.0; // closure after moving one component of v1 by one unit in the last place
    for (int i = 0; i < 3; ++i) {
        for (double direction : {-1.0, 1.0}) {
            Vec3    v = s.v1;
            double& c = i == 0 ? v.x : i == 1 ? v.y : v.z;
            c         = std::nextafter(c, direction * 1e300);
            floor     = std::max(floor, closure(r1, r2, v, tof, mu).value_or(1.0));
        }
    }
    t.beyond += e > 20.0 * floor + 1e-12;
}

/// One regime: `generate(r1, r2, tof, mu, prograde)` fills a problem and returns false to skip it.
template <class Generate>
void run(const char* name, int problems, int maxRevolutions, Generate&& generate) {
    Tally t;
    for (int k = 0; k < problems; ++k) {
        Vec3   r1, r2;
        double tof = 0.0, mu = 1.0;
        bool   prograde = true;
        if (!generate(r1, r2, tof, mu, prograde)) {
            continue;
        }
        ++t.problems;
        for (const LambertSolution& s : lambertIzzo(r1, r2, tof, mu, prograde, maxRevolutions)) {
            judge(t, s, r1, r2, tof, mu);
        }
    }
    std::printf("| %s | %d | %ld | %ld | %.1e | %ld | %ld | %ld | %ld |\n",
                name,
                maxRevolutions,
                t.solutions,
                t.unchecked,
                t.worst,
                t.above10,
                t.above6,
                t.beyond,
                t.nonFinite);
}

} // namespace

int main(int argc, char** argv) {
    const int                              n = argc > 1 ? std::atoi(argv[1]) : 10000;
    std::mt19937_64                        rng(77);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    const auto                             unit = [&] { // a random direction
        const double z = 2.0 * U(rng) - 1.0, a = constants::twoPi * U(rng), s = std::sqrt(1.0 - z * z);
        return Vec3{s * std::cos(a), s * std::sin(a), z};
    };
    const auto turned = [](const Vec3& v, const Vec3& axis, double angle) { // Rodrigues' rotation
        return std::cos(angle) * v + std::sin(angle) * cross(axis, v) + (1.0 - std::cos(angle)) * dot(axis, v) * axis;
    };
    const auto notCollinear = [](const Vec3& a, const Vec3& b) {
        return ::astro::norm(cross(a, b)) > 1e-6 * ::astro::norm(a) * ::astro::norm(b);
    };
    const auto radius = [&] { return 0.5 + 3.0 * U(rng); };

    std::printf("Lambert solutions against extended-precision propagation (closure = |r(t) - r2| / |r2|)\n\n");
    std::printf(
        "| Regime (mu = 1 unless stated) | Max revolutions | Solutions | Not checkable | Worst closure | > 1e-10 "
        "| > 1e-6 | Beyond conditioning | Non-finite |\n|---|---|---|---|---|---|---|---|---|\n");

    run("Izzo's test box: components in [-4, 4], t in [0.1, 100]",
        2 * n,
        5,
        [&](Vec3& r1, Vec3& r2, double& tof, double&, bool& prograde) {
            r1       = {8 * U(rng) - 4, 8 * U(rng) - 4, 8 * U(rng) - 4};
            r2       = {8 * U(rng) - 4, 8 * U(rng) - 4, 8 * U(rng) - 4};
            tof      = 0.1 + 99.9 * U(rng);
            prograde = U(rng) < 0.5;
            return ::astro::norm(r1) > 0.05 && ::astro::norm(r2) > 0.05 && notCollinear(r1, r2);
        });
    run("Earth orbits, SI units: radii 6600-66600 km, 5 min to 3.5 d",
        n,
        10,
        [&](Vec3& r1, Vec3& r2, double& tof, double& mu, bool& prograde) {
            r1       = (6.6e6 + 6e7 * U(rng)) * unit();
            r2       = (6.6e6 + 6e7 * U(rng)) * unit();
            tof      = 300.0 + 3e5 * U(rng);
            mu       = constants::gmEarth;
            prograde = U(rng) < 0.5;
            return notCollinear(r1, r2);
        });
    run("Many revolutions: t in [10, 1500]", n / 20, 50, [&](Vec3& r1, Vec3& r2, double& tof, double&, bool& prograde) {
        r1       = radius() * unit();
        r2       = radius() * unit();
        tof      = 10.0 + 1500.0 * U(rng);
        prograde = U(rng) < 0.5;
        return notCollinear(r1, r2);
    });
    run("Transfer angle within 1e-2 to 1e-6 rad of 0 or 360 deg",
        n,
        5,
        [&](Vec3& r1, Vec3& r2, double& tof, double&, bool& prograde) {
            r1  = radius() * unit();
            r2  = (0.9 + 0.2 * U(rng)) * turned(r1, normalized(cross(r1, unit())), std::pow(10.0, -2.0 - 4.0 * U(rng)));
            tof = 0.1 + 99.9 * U(rng);
            prograde = U(rng) < 0.5;
            return true;
        });
    run("Transfer angle within 1e-2 to 1e-6 rad of 180 deg",
        n,
        5,
        [&](Vec3& r1, Vec3& r2, double& tof, double&, bool& prograde) {
            r1 = radius() * unit();
            r2 = radius() / ::astro::norm(r1)
               * turned(r1, normalized(cross(r1, unit())), constants::pi - std::pow(10.0, -2.0 - 4.0 * U(rng)));
            tof      = 0.1 + 99.9 * U(rng);
            prograde = U(rng) < 0.5;
            return true;
        });
    run("Within 1e-2 to 1e-9 of the parabolic time",
        2 * n,
        0,
        [&](Vec3& r1, Vec3& r2, double& tof, double&, bool& prograde) {
            r1       = radius() * unit();
            r2       = radius() * unit();
            prograde = U(rng) < 0.5;
            if (::astro::norm(cross(r1, r2)) < 1e-3 * ::astro::norm(r1) * ::astro::norm(r2)) {
                return false;
            }
            // Euler's equation: + the long way round.
            const double R1 = ::astro::norm(r1), R2 = ::astro::norm(r2), c = ::astro::norm(r2 - r1);
            const bool   longWay = prograde ? cross(r1, r2).z < 0.0 : cross(r1, r2).z > 0.0;
            const double parabolic =
                (std::pow(R1 + R2 + c, 1.5) + (longWay ? 1.0 : -1.0) * std::pow(R1 + R2 - c, 1.5)) / 6.0;
            tof = parabolic * (1.0 + (2.0 * U(rng) - 1.0) * std::pow(10.0, -2.0 - 7.0 * U(rng)));
            return true;
        });
    for (bool shortWay : {true, false}) {
        run(shortWay ? "Very short flights, t in [1e-6, 0.1], the short way round"
                     : "Very short flights, t in [1e-6, 0.1], the long way round (past the centre)",
            n,
            0,
            [&](Vec3& r1, Vec3& r2, double& tof, double&, bool& prograde) {
                r1       = radius() * unit();
                r2       = radius() * unit();
                tof      = std::pow(10.0, -6.0 + 5.0 * U(rng));
                prograde = (cross(r1, r2).z > 0.0) == shortWay;
                return notCollinear(r1, r2);
            });
    }
    run("Very long flights, t in [1e2, 1e8]", n, 0, [&](Vec3& r1, Vec3& r2, double& tof, double&, bool& prograde) {
        r1       = radius() * unit();
        r2       = radius() * unit();
        tof      = std::pow(10.0, 2.0 + 6.0 * U(rng));
        prograde = U(rng) < 0.5;
        return notCollinear(r1, r2);
    });
    return 0;
}
