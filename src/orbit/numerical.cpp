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
#include <numeric>
#include <utility>

#include "astro/orbit/gravity.h"

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
    Vec3 a = pointMassAcceleration(r, mu_);
    if (J_.size() > 2) {
        if (identity_) {
            a += zonalAcceleration(r, J_, mu_, re_);
        } else {
            a += R_.transpose() * zonalAcceleration(R_ * r, J_, mu_, re_);
        }
    }
    return a;
}

Vec3 GravityField::jerk(double t, const Vec3& r, const Vec3& v, const Vec3& totalAcceleration) const {
    Vec3 j = pointMassJerk(r, v, mu_);
    if (J_.size() > 2) {
        auto zonal = [this](double, const Vec3& rr, const Vec3&) {
            return identity_ ? zonalAcceleration(rr, J_, mu_, re_)
                             : R_.transpose() * zonalAcceleration(R_ * rr, J_, mu_, re_);
        };
        j += derivativeAlongMotion(zonal, t, r, v, totalAcceleration);
    }
    return j;
}

Vec3 UnmodeledAcceleration::acceleration(double t, const Vec3& /*r*/, const Vec3& /*v*/) const {
    return a0_ + j0_ * t;
}

Vec3 UnmodeledAcceleration::jerk(double /*t*/, const Vec3& /*r*/, const Vec3& /*v*/, const Vec3& /*a*/) const {
    return j0_;
}

ThirdBodyForce::ThirdBodyForce(Body body, double muBody, std::shared_ptr<const Ephemeris> ephemeris, const Time& epoch)
    : body_(body)
    , mu_(muBody)
    , eph_(std::move(ephemeris))
    , epoch_(epoch) {}

Vec3 ThirdBodyForce::acceleration(double t, const Vec3& r, const Vec3& /*v*/) const {
    return thirdBodyAcceleration(r, eph_->position(body_, Body::Earth, epoch_ + t), mu_);
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

StateVector12 NumericalPropagator::propagate12(double dt) const {
    return withDerivatives(dt, propagate(dt));
}

std::vector<StateVector12> NumericalPropagator::propagate12(const std::vector<double>& dts) const {
    const std::vector<StateVector> states = propagate(dts);
    std::vector<StateVector12>     out;
    out.reserve(states.size());
    for (std::size_t i = 0; i < states.size(); ++i) {
        out.push_back(withDerivatives(dts[i], states[i]));
    }
    return out;
}

StateVector NumericalPropagator::propagate(double dt) const {
    return propagate(std::vector<double>{dt}).front();
}

std::vector<StateVector> NumericalPropagator::propagate(const std::vector<double>& dts) const {
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
    stats_ = total;
    return out;
}

} // namespace astro
