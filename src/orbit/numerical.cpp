/**********************************************************************
 * numerical.cpp
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

#include "astro/orbit/numerical.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>

#include "astro/core/error.h"
#include "astro/orbit/gravity.h"
#include "gravity_detail.h"

namespace astro {

namespace {

/// Central difference of accel(t, r, v) along (1, v, a): the time derivative along the motion.
template <class F>
Vec3 derivativeAlongMotion(const F& accel, double t, const Vec3& r, const Vec3& v, const Vec3& a) {
    const double vm = norm(v);
    const double h  = vm > 0.0 ? 1e-5 * norm(r) / vm : 1.0;
    return (accel(t + h, r + v * h, v + a * h) - accel(t - h, r - v * h, v - a * h)) / (2.0 * h);
}

} // namespace

Vec3 ForceModel::partials(double t, const Vec3& r, const Vec3& v, Mat3& dadr, Mat3& dadv) const {
    // Central differences, one column per position and velocity component.
    const double hr = 1e-6 * std::max(norm(r), 1.0);
    const double hv = 1e-6 * std::max(norm(v), 1.0);
    for (std::size_t j = 0; j < 3; ++j) {
        Vec3 dr, dv;
        (j == 0 ? dr.x : j == 1 ? dr.y : dr.z) = hr;
        (j == 0 ? dv.x : j == 1 ? dv.y : dv.z) = hv;
        const Vec3 ar                          = (acceleration(t, r + dr, v) - acceleration(t, r - dr, v)) / (2.0 * hr);
        const Vec3 av                          = (acceleration(t, r, v + dv) - acceleration(t, r, v - dv)) / (2.0 * hv);
        dadr(0, j)                             = ar.x;
        dadr(1, j)                             = ar.y;
        dadr(2, j)                             = ar.z;
        dadv(0, j)                             = av.x;
        dadv(1, j)                             = av.y;
        dadv(2, j)                             = av.z;
    }
    return acceleration(t, r, v);
}

Vec3 ForceModel::jerk(double t, const Vec3& r, const Vec3& v, const Vec3& totalAcceleration) const {
    return derivativeAlongMotion([this](double tt, const Vec3& rr, const Vec3& vv) { return acceleration(tt, rr, vv); },
                                 t,
                                 r,
                                 v,
                                 totalAcceleration);
}

GravityField::GravityField(std::vector<double> J, double mu, double re, const Mat3& toBodyFrame)
    : J_(std::move(J))
    , mu_(mu)
    , re_(re)
    , R_(toBodyFrame)
    , identity_(toBodyFrame == Mat3::identity()) {}

Vec3 GravityField::acceleration(double /*t*/, const Vec3& r, const Vec3& /*v*/) const {
    Vec3 a;
    if (identity_) {
        detail::gravity<false>(r, J_.data(), J_.size(), mu_, re_, true, &a, nullptr);
        return a;
    }
    detail::gravity<false>(R_ * r, J_.data(), J_.size(), mu_, re_, true, &a, nullptr);
    return R_.transpose() * a;
}

Vec3 GravityField::jerk(double /*t*/, const Vec3& r, const Vec3& v, const Vec3& /*totalAcceleration*/) const {
    detail::GravityPartials p;
    if (identity_) {
        detail::gravity<true>(r, J_.data(), J_.size(), mu_, re_, true, nullptr, &p);
        return p.jerk(v);
    }
    detail::gravity<true>(R_ * r, J_.data(), J_.size(), mu_, re_, true, nullptr, &p);
    return R_.transpose() * p.jerk(R_ * v);
}

Vec3 GravityField::partials(double /*t*/, const Vec3& r, const Vec3& /*v*/, Mat3& dadr, Mat3& dadv) const {
    detail::GravityPartials p;
    dadv = Mat3{};
    if (identity_) {
        detail::gravity<true>(r, J_.data(), J_.size(), mu_, re_, true, nullptr, &p);
        dadr = p.gradient();
        return p.acceleration;
    }
    detail::gravity<true>(R_ * r, J_.data(), J_.size(), mu_, re_, true, nullptr, &p);
    dadr = R_.transpose() * p.gradient() * R_;
    return R_.transpose() * p.acceleration;
}

Vec3 UnmodeledAcceleration::acceleration(double t, const Vec3& /*r*/, const Vec3& /*v*/) const {
    return a0_ + j0_ * t;
}

Vec3 UnmodeledAcceleration::jerk(double /*t*/, const Vec3& /*r*/, const Vec3& /*v*/, const Vec3& /*a*/) const {
    return j0_;
}

Vec3 UnmodeledAcceleration::partials(double t, const Vec3& /*r*/, const Vec3& /*v*/, Mat3& dadr, Mat3& dadv) const {
    dadr = Mat3{};
    dadv = Mat3{};
    return a0_ + j0_ * t;
}

ThirdBodyForce::ThirdBodyForce(Body body, double muBody, std::shared_ptr<const Ephemeris> ephemeris, const Time& epoch)
    : body_(body)
    , mu_(muBody)
    , eph_(std::move(ephemeris))
    , epoch_(epoch) {}

Vec3 ThirdBodyForce::acceleration(double t, const Vec3& r, const Vec3& /*v*/) const {
    return thirdBodyAcceleration(r, eph_->position(body_, Body::Earth, epoch_ + t), mu_);
}

Vec3 ThirdBodyForce::partials(double t, const Vec3& r, const Vec3& /*v*/, Mat3& dadr, Mat3& dadv) const {
    // a = mu (d / |d|^3 - s / |s|^3) with d = s - r, so da/dr = mu / |d|^3 (3 d_hat d_hat^T - I).
    const Vec3 s = eph_->position(body_, Body::Earth, epoch_ + t);
    dadr         = pointMassGradient(s - r, mu_);
    dadv         = Mat3{};
    return thirdBodyAcceleration(r, s, mu_);
}

NumericalPropagator::NumericalPropagator(const Time&                                    epoch,
                                         const StateVector&                             initial,
                                         std::vector<std::shared_ptr<const ForceModel>> forces,
                                         const AdaptiveOptions&                         options)
    : epoch_(epoch)
    , s0_(initial)
    , forces_(std::move(forces))
    , opt_(options) {
    if (forces_.empty()) {
        throw Error("NumericalPropagator: at least one force model is required");
    }
}

NumericalPropagator::NumericalPropagator(const Time&                                    epoch,
                                         const StateVector12&                           initial,
                                         std::vector<std::shared_ptr<const ForceModel>> forces,
                                         HigherOrderMode                                mode,
                                         const AdaptiveOptions&                         options)
    : NumericalPropagator(epoch, initial.rv(), std::move(forces), options) {
    Vec3 a;
    for (const auto& f : forces_) {
        a += f->acceleration(0.0, initial.r, initial.v);
    }
    da0_ = initial.a - a;
    // In UnmodeledForce mode the total acceleration at the epoch is the input one; velocity-dependent
    // forces (drag) see it through the jerk's dependence on the acceleration.
    const Vec3& aTotal = mode == HigherOrderMode::UnmodeledForce ? initial.a : a;
    Vec3        j;
    for (const auto& f : forces_) {
        j += f->jerk(0.0, initial.r, initial.v, aTotal);
    }
    dj0_ = initial.j - j;
    if (mode == HigherOrderMode::UnmodeledForce) {
        forces_.push_back(std::make_shared<UnmodeledAcceleration>(da0_, dj0_));
    }
}

namespace {

std::shared_ptr<const ForceModel> j2j4Field(Frame frame, const Time& epoch) {
    return std::make_shared<GravityField>(std::vector<double>{0.0, 0.0, constants::J2, 0.0, constants::J4},
                                          constants::gmEarth,
                                          constants::earthRadiusEgm2008,
                                          poleOfDateMatrix(frame, epoch));
}

} // namespace

NumericalPropagator
NumericalPropagator::j2j4(const Time& epoch, const StateVector& initial, Frame frame, const AdaptiveOptions& options) {
    return NumericalPropagator(epoch, initial, {j2j4Field(frame, epoch)}, options);
}

NumericalPropagator NumericalPropagator::j2j4(const Time&            epoch,
                                              const StateVector12&   initial,
                                              Frame                  frame,
                                              HigherOrderMode        mode,
                                              const AdaptiveOptions& options) {
    return NumericalPropagator(epoch, initial, {j2j4Field(frame, epoch)}, mode, options);
}

NumericalPropagator NumericalPropagator::j2j4(const Time&          epoch,
                                              const StateVector12& initial,
                                              HigherOrderMode      mode,
                                              double               mu,
                                              double               re,
                                              double               J2,
                                              double               J4) {
    return NumericalPropagator(epoch,
                               initial,
                               {std::make_shared<GravityField>(std::vector<double>{0.0, 0.0, J2, 0.0, J4}, mu, re)},
                               mode);
}

NumericalPropagator
NumericalPropagator::j2j4(const Time& epoch, const StateVector& initial, double mu, double re, double J2, double J4) {
    return NumericalPropagator(epoch,
                               initial,
                               {std::make_shared<GravityField>(std::vector<double>{0.0, 0.0, J2, 0.0, J4}, mu, re)});
}

OdeState<6> NumericalPropagator::derivative(double t, const OdeState<6>& y) const {
    const Vec3 r{y[0], y[1], y[2]}, v{y[3], y[4], y[5]};
    Vec3       a;
    for (const auto& f : forces_) {
        a += f->acceleration(t, r, v);
    }
    return {v.x, v.y, v.z, a.x, a.y, a.z};
}

namespace {

StateVector toStateVector(const OdeState<6>& y) {
    return {
        {y[0], y[1], y[2]},
        {y[3], y[4], y[5]}
    };
}

/// Degree-7 Hermite interpolant of position over one step [ta, ta + h], matching position, velocity,
/// acceleration and jerk at both ends: p(tau) = sum c_k tau^k with tau = (t - ta) / h.
struct SepticHermite {
    Vec3   c[8];
    double h{0.0};

    SepticHermite(const StateVector12& A, const StateVector12& B, double step) noexcept
        : h(step) {
        const double h2 = h * h, h3 = h2 * h;
        c[0] = A.r;
        c[1] = A.v * h;
        c[2] = A.a * (0.5 * h2);
        c[3] = A.j * (h3 / 6.0);
        // Remaining conditions at tau = 1, then the exact inverse of the 4x4 system for c4..c7.
        const Vec3 D0 = B.r - (c[0] + c[1] + c[2] + c[3]);
        const Vec3 D1 = B.v * h - (c[1] + 2.0 * c[2] + 3.0 * c[3]);
        const Vec3 D2 = B.a * (0.5 * h2) - (c[2] + 3.0 * c[3]);
        const Vec3 D3 = B.j * (h3 / 6.0) - c[3];
        c[4]          = 35.0 * D0 - 15.0 * D1 + 5.0 * D2 - D3;
        c[5]          = -84.0 * D0 + 39.0 * D1 - 14.0 * D2 + 3.0 * D3;
        c[6]          = 70.0 * D0 - 34.0 * D1 + 13.0 * D2 - 3.0 * D3;
        c[7]          = -20.0 * D0 + 10.0 * D1 - 4.0 * D2 + D3;
    }

    [[nodiscard]] StateVector at(double tau) const noexcept {
        Vec3 p = c[7], dp = 7.0 * c[7];
        for (int k = 6; k >= 1; --k) {
            p  = p * tau + c[k];
            dp = dp * tau + static_cast<double>(k) * c[k];
        }
        p = p * tau + c[0];
        return {p, dp / h};
    }
};

} // namespace

template <class Emit>
void NumericalPropagator::grid(double start, double step, std::size_t count, AdaptiveStats* stats, Emit&& emit) const {
    if (count == 0) {
        return;
    }
    if (!(std::abs(step) > 0.0) || !std::isfinite(step) || !std::isfinite(start)) {
        throw Error("NumericalPropagator::propagateGrid: step must be finite and non-zero");
    }
    auto f = [this](double t, const OdeState<6>& y) { return derivative(t, y); };

    AdaptiveStats total, st;
    auto          add = [&] {
        total.accepted    += st.accepted;
        total.rejected    += st.rejected;
        total.evaluations += st.evaluations;
        st                 = {};
    };

    // Integrate (exactly) to the first output time.
    OdeState<6> y{s0_.r.x, s0_.r.y, s0_.r.z, s0_.v.x, s0_.v.y, s0_.v.z};
    if (start != 0.0) {
        y = integrateRkf78<6>(f, 0.0, y, start, opt_, &st);
        add();
    }
    emit(std::size_t{0}, start, toStateVector(y));

    if (count > 1) {
        const double  dir    = step > 0.0 ? 1.0 : -1.0;
        const double  tEnd   = start + step * static_cast<double>(count - 1);
        std::size_t   k      = 1;
        StateVector12 A      = withDerivatives(start, toStateVector(y));
        auto          onStep = [&](double ta, const OdeState<6>& /*ya*/, double tb, const OdeState<6>& yb) {
            const StateVector12 B   = withDerivatives(tb, toStateVector(yb));
            const double        h   = tb - ta;
            const double        tk0 = start + step * static_cast<double>(k);
            if (k < count && dir * (tk0 - tb) <= 0.0) {
                const SepticHermite poly(A, B, h);
                for (; k < count; ++k) {
                    const double tk = start + step * static_cast<double>(k);
                    if (dir * (tk - tb) > 0.0) {
                        break;
                    }
                    emit(k, tk, tk == tb ? B.rv() : poly.at((tk - ta) / h));
                }
            }
            A = B;
        };
        (void)integrateRkf78Steps<6>(f, start, y, tEnd, onStep, opt_, &st);
        add();
    }
    if (stats) {
        *stats = total;
    }
}

// ------------------------------------------------------------ variational equations
//
// The extended state y holds r, v (0-5), Phi row-major (6-41) and, with process noise, Q row-major
// (42-77). With A = [[0, I], [G, D]], G = da/dr and D = da/dv:
//     Phi' = A Phi,     Q' = A Q + Q A^T + diag(0, 0, 0, q, q, q).

namespace {

constexpr std::size_t kPhi = 6;  // offset of Phi in the extended state
constexpr std::size_t kQ   = 42; // offset of Q

/// Writes A X into d for the 6x6 row-major block X at y (A = [[0, I], [G, D]]).
void applyA(const double* y, double* d, const Mat3& G, const Mat3& D) noexcept {
    for (std::size_t j = 0; j < 6; ++j) {
        for (std::size_t i = 0; i < 3; ++i) {
            d[6 * i + j] = y[6 * (i + 3) + j];
            double s     = 0.0;
            for (std::size_t m = 0; m < 3; ++m) {
                s += G(i, m) * y[6 * m + j] + D(i, m) * y[6 * (m + 3) + j];
            }
            d[6 * (i + 3) + j] = s;
        }
    }
}

/// Adds Adot X to d for the 6x6 row-major block X at y, with Adot = [[0, 0], [Gd, Dd]].
void addAdot(const double* y, double* d, const Mat3& Gd, const Mat3& Dd) noexcept {
    for (std::size_t j = 0; j < 6; ++j) {
        for (std::size_t i = 0; i < 3; ++i) {
            double s = 0.0;
            for (std::size_t m = 0; m < 3; ++m) {
                s += Gd(i, m) * y[6 * m + j] + Dd(i, m) * y[6 * (m + 3) + j];
            }
            d[6 * (i + 3) + j] += s;
        }
    }
}

/// The 6x6 row-major block at y.
Mat6 unpack(const double* y) noexcept {
    Mat6 m;
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            m[i][j] = y[6 * i + j];
        }
    }
    return m;
}

} // namespace

template <std::size_t N, class Emit>
void NumericalPropagator::variationalGrid(double         start,
                                          double         step,
                                          std::size_t    count,
                                          double         q,
                                          AdaptiveStats* stats,
                                          Emit&&         emit) const {
    static_assert(N == 42 || N == 78);
    constexpr bool withNoise = N == 78;
    if (count == 0) {
        return;
    }
    if (!(std::abs(step) > 0.0) || !std::isfinite(step) || !std::isfinite(start)) {
        throw Error("NumericalPropagator: grid step must be finite and non-zero");
    }
    if (!(q >= 0.0) || !std::isfinite(q)) {
        throw Error("NumericalPropagator: acceleration noise PSD must be finite and non-negative");
    }
    AdaptiveOptions opt = opt_;
    opt.errorComponents = 6; // same steps as the state alone

    // Sum of the force models' acceleration and partials.
    auto partialsAt = [this](double t, const Vec3& r, const Vec3& v, Mat3& G, Mat3& D) {
        Vec3 a;
        G = Mat3{};
        D = Mat3{};
        for (const auto& force : forces_) {
            Mat3 gr, gv;
            a += force->partials(t, r, v, gr, gv);
            G += gr;
            D += gv;
        }
        return a;
    };
    auto f = [&](double t, const OdeState<N>& y) {
        const Vec3  r{y[0], y[1], y[2]}, v{y[3], y[4], y[5]};
        Mat3        G, D;
        const Vec3  a = partialsAt(t, r, v, G, D);
        OdeState<N> d{};
        d[0] = v.x;
        d[1] = v.y;
        d[2] = v.z;
        d[3] = a.x;
        d[4] = a.y;
        d[5] = a.z;
        applyA(y.data() + kPhi, d.data() + kPhi, G, D);
        if constexpr (withNoise) {
            double m[36]; // A Q
            applyA(y.data() + kQ, m, G, D);
            for (std::size_t i = 0; i < 6; ++i) {
                for (std::size_t j = 0; j < 6; ++j) {
                    d[kQ + 6 * i + j] = m[6 * i + j] + m[6 * j + i] + (i == j && i >= 3 ? q : 0.0);
                }
            }
        }
        return d;
    };
    // Second derivatives of Phi and Q for the quintic Hermite interpolation:
    //     Phi'' = A Phi' + Adot Phi,   Q'' = A Q' + Q' A^T + Adot Q + Q Adot^T,
    // with Adot from a central difference of G and D along the motion.
    auto second = [&partialsAt](double t, const OdeState<N>& y, const OdeState<N>& dy) {
        const Vec3   r{y[0], y[1], y[2]}, v{y[3], y[4], y[5]}, a{dy[3], dy[4], dy[5]};
        const double vm = norm(v);
        const double h  = vm > 0.0 ? 1e-5 * norm(r) / vm : 1.0;
        Mat3         G, D, Gp, Dp, Gm, Dm;
        (void)partialsAt(t, r, v, G, D);
        (void)partialsAt(t + h, r + v * h, v + a * h, Gp, Dp);
        (void)partialsAt(t - h, r - v * h, v - a * h, Gm, Dm);
        const Mat3  Gd = (Gp - Gm) * (0.5 / h);
        const Mat3  Dd = (Dp - Dm) * (0.5 / h);
        OdeState<N> dd{};
        applyA(dy.data() + kPhi, dd.data() + kPhi, G, D);
        addAdot(y.data() + kPhi, dd.data() + kPhi, Gd, Dd);
        if constexpr (withNoise) {
            double m[36]; // A Q' + Adot Q
            applyA(dy.data() + kQ, m, G, D);
            addAdot(y.data() + kQ, m, Gd, Dd);
            for (std::size_t i = 0; i < 6; ++i) {
                for (std::size_t j = 0; j < 6; ++j) {
                    dd[kQ + 6 * i + j] = m[6 * i + j] + m[6 * j + i];
                }
            }
        }
        return dd;
    };

    AdaptiveStats total, st;
    auto          add = [&] {
        total.accepted    += st.accepted;
        total.rejected    += st.rejected;
        total.evaluations += st.evaluations;
        st                 = {};
    };
    auto noiseOf = [](const OdeState<N>& y) {
        if constexpr (withNoise) {
            return unpack(y.data() + kQ);
        } else {
            (void)y;
            return Mat6{};
        }
    };
    auto stateOf = [](const OdeState<N>& y) {
        return StateVector{
            {y[0], y[1], y[2]},
            {y[3], y[4], y[5]}
        };
    };

    OdeState<N> y{};
    y[0] = s0_.r.x;
    y[1] = s0_.r.y;
    y[2] = s0_.r.z;
    y[3] = s0_.v.x;
    y[4] = s0_.v.y;
    y[5] = s0_.v.z;
    for (std::size_t i = 0; i < 6; ++i) {
        y[kPhi + 7 * i] = 1.0;
    }
    if (start != 0.0) {
        y = integrateRkf78<N>(f, 0.0, y, start, opt, &st);
        add();
    }
    emit(std::size_t{0}, stateOf(y), unpack(y.data() + kPhi), noiseOf(y));

    if (count > 1) {
        const double  dir    = step > 0.0 ? 1.0 : -1.0;
        const double  tEnd   = start + step * static_cast<double>(count - 1);
        std::size_t   k      = 1;
        StateVector12 A      = withDerivatives(start, stateOf(y));
        OdeState<N>   dA     = f(start, y);
        OdeState<N>   ddA    = second(start, y, dA);
        auto          onStep = [&](double ta, const OdeState<N>& ya, double tb, const OdeState<N>& yb) {
            const StateVector12 B   = withDerivatives(tb, stateOf(yb));
            const OdeState<N>   dB  = f(tb, yb);
            const OdeState<N>   ddB = second(tb, yb, dB);
            const double        h   = tb - ta;
            const double        tk0 = start + step * static_cast<double>(k);
            if (k < count && dir * (tk0 - tb) <= 0.0) {
                const SepticHermite poly(A, B, h);
                for (; k < count; ++k) {
                    const double tk = start + step * static_cast<double>(k);
                    if (dir * (tk - tb) > 0.0) {
                        break;
                    }
                    if (tk == tb) {
                        emit(k, B.rv(), unpack(yb.data() + kPhi), noiseOf(yb));
                        continue;
                    }
                    // Quintic Hermite for Phi and Q from values, first and second derivatives.
                    const double tau = (tk - ta) / h, t2 = tau * tau, t3 = t2 * tau, t4 = t3 * tau, t5 = t4 * tau;
                    const double H0 = 1.0 - 10.0 * t3 + 15.0 * t4 - 6.0 * t5;
                    const double H1 = (tau - 6.0 * t3 + 8.0 * t4 - 3.0 * t5) * h;
                    const double H2 = 0.5 * (t2 - 3.0 * t3 + 3.0 * t4 - t5) * h * h;
                    const double H3 = 10.0 * t3 - 15.0 * t4 + 6.0 * t5;
                    const double H4 = (-4.0 * t3 + 7.0 * t4 - 3.0 * t5) * h;
                    const double H5 = 0.5 * (t3 - 2.0 * t4 + t5) * h * h;
                    OdeState<N>  yi{};
                    for (std::size_t i = kPhi; i < N; ++i) {
                        yi[i] = H0 * ya[i] + H1 * dA[i] + H2 * ddA[i] + H3 * yb[i] + H4 * dB[i] + H5 * ddB[i];
                    }
                    emit(k, poly.at(tau), unpack(yi.data() + kPhi), noiseOf(yi));
                }
            }
            A   = B;
            dA  = dB;
            ddA = ddB;
        };
        (void)integrateRkf78Steps<N>(f, start, y, tEnd, onStep, opt, &st);
        add();
    }
    if (stats) {
        *stats = total;
    }
}

void NumericalPropagator::propagateGridWithCovariance(double                     start,
                                                      double                     step,
                                                      const Mat6&                P0,
                                                      std::span<StateCovariance> out,
                                                      double                     accelerationNoisePsd,
                                                      AdaptiveStats*             stats) const {
    auto emit = [&](std::size_t k, const StateVector& s, const Mat6& phi, const Mat6& Q) {
        out[k].state      = s;
        out[k].covariance = add(congruence(phi, P0), Q);
    };
    if (accelerationNoisePsd > 0.0) {
        variationalGrid<78>(start, step, out.size(), accelerationNoisePsd, stats, emit);
    } else {
        variationalGrid<42>(start, step, out.size(), accelerationNoisePsd, stats, emit);
    }
}

std::vector<StateCovariance> NumericalPropagator::propagateGridWithCovariance(double         start,
                                                                              double         step,
                                                                              std::size_t    count,
                                                                              const Mat6&    P0,
                                                                              double         accelerationNoisePsd,
                                                                              AdaptiveStats* stats) const {
    std::vector<StateCovariance> out(count);
    propagateGridWithCovariance(start, step, P0, out, accelerationNoisePsd, stats);
    return out;
}

void NumericalPropagator::propagateGridStm(double                 start,
                                           double                 step,
                                           std::span<StateVector> states,
                                           std::span<Mat6>        stm,
                                           AdaptiveStats*         stats) const {
    if (states.size() != stm.size()) {
        throw Error("NumericalPropagator::propagateGridStm: states and stm must have the same size");
    }
    variationalGrid<42>(start,
                        step,
                        states.size(),
                        0.0,
                        stats,
                        [&](std::size_t k, const StateVector& s, const Mat6& phi, const Mat6&) {
                            states[k] = s;
                            stm[k]    = phi;
                        });
}

void NumericalPropagator::propagateGrid(double                 start,
                                        double                 step,
                                        std::span<StateVector> out,
                                        AdaptiveStats*         stats) const {
    grid(start, step, out.size(), stats, [&](std::size_t k, double, const StateVector& s) { out[k] = s; });
}

void NumericalPropagator::propagateGrid12(double                   start,
                                          double                   step,
                                          std::span<StateVector12> out,
                                          AdaptiveStats*           stats) const {
    grid(start, step, out.size(), stats, [&](std::size_t k, double t, const StateVector& s) {
        out[k] = withDerivatives(t, s);
    });
}

std::vector<StateVector>
NumericalPropagator::propagateGrid(double start, double step, std::size_t count, AdaptiveStats* stats) const {
    std::vector<StateVector> out(count);
    propagateGrid(start, step, out, stats);
    return out;
}

std::vector<StateVector12>
NumericalPropagator::propagateGrid12(double start, double step, std::size_t count, AdaptiveStats* stats) const {
    std::vector<StateVector12> out(count);
    propagateGrid12(start, step, out, stats);
    return out;
}

StateVector12 NumericalPropagator::withDerivatives(double t, const StateVector& s) const {
    Vec3 a;
    for (const auto& f : forces_) {
        a += f->acceleration(t, s.r, s.v);
    }
    Vec3 j;
    for (const auto& f : forces_) {
        j += f->jerk(t, s.r, s.v, a);
    }
    return {s.r, s.v, a, j};
}

StateVector12 NumericalPropagator::propagate12(double dt, AdaptiveStats* stats) const {
    return withDerivatives(dt, propagate(dt, stats));
}

std::vector<StateVector12> NumericalPropagator::propagate12(const std::vector<double>& dts,
                                                            AdaptiveStats*             stats) const {
    const std::vector<StateVector> states = propagate(dts, stats);
    std::vector<StateVector12>     out;
    out.reserve(states.size());
    for (std::size_t i = 0; i < states.size(); ++i) {
        out.push_back(withDerivatives(dts[i], states[i]));
    }
    return out;
}

StateVector NumericalPropagator::propagate(double dt, AdaptiveStats* stats) const {
    return propagate(std::vector<double>{dt}, stats).front();
}

std::vector<StateVector> NumericalPropagator::propagate(const std::vector<double>& dts, AdaptiveStats* stats) const {
    std::vector<StateVector> out(dts.size());
    std::vector<std::size_t> order(dts.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return dts[a] < dts[b]; });
    auto              f = [this](double t, const OdeState<6>& y) { return derivative(t, y); };
    const OdeState<6> y0{s0_.r.x, s0_.r.y, s0_.r.z, s0_.v.x, s0_.v.y, s0_.v.z};
    auto              toState = [](const OdeState<6>& y) {
        return StateVector{
                         {y[0], y[1], y[2]},
                         {y[3], y[4], y[5]}
        };
    };
    AdaptiveStats total;
    auto          run = [&](auto begin, auto end) {
        OdeState<6> y = y0;
        double      t = 0.0;
        for (auto it = begin; it != end; ++it) {
            AdaptiveStats st;
            y                  = integrateRkf78<6>(f, t, y, dts[*it], opt_, &st);
            total.accepted    += st.accepted;
            total.rejected    += st.rejected;
            total.evaluations += st.evaluations;
            t                  = dts[*it];
            out[*it]           = toState(y);
        }
    };
    // Forward for dt >= 0 in ascending order, backward for dt < 0 in descending order.
    const auto split = std::find_if(order.begin(), order.end(), [&](std::size_t i) { return dts[i] >= 0.0; });
    run(split, order.end());
    run(std::make_reverse_iterator(split), std::make_reverse_iterator(order.begin()));
    if (stats) {
        *stats = total;
    }
    return out;
}

} // namespace astro
