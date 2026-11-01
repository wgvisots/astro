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

#include "astro/orbit/gravity.h"

namespace astro {

GravityField::GravityField(std::vector<double> J, double mu, double re, const Mat3& toBodyFrame)
    : J_(std::move(J)), mu_(mu), re_(re), R_(toBodyFrame), identity_(toBodyFrame == Mat3::identity()) {}

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

ThirdBodyForce::ThirdBodyForce(Body body, double muBody, std::shared_ptr<const Ephemeris> ephemeris, const Time& epoch)
    : body_(body), mu_(muBody), eph_(std::move(ephemeris)), epoch_(epoch) {}

Vec3 ThirdBodyForce::acceleration(double t, const Vec3& r, const Vec3& /*v*/) const {
    return thirdBodyAcceleration(r, eph_->position(body_, Body::Earth, epoch_ + t), mu_);
}

NumericalPropagator::NumericalPropagator(const Time& epoch, const StateVector& initial,
                                         std::vector<std::shared_ptr<const ForceModel>> forces,
                                         const AdaptiveOptions& options)
    : epoch_(epoch), s0_(initial), forces_(std::move(forces)), opt_(options) {
    if (forces_.empty()) throw Error("NumericalPropagator: at least one force model is required");
}

NumericalPropagator NumericalPropagator::j2j4(const Time& epoch, const StateVector& initial, double mu, double re,
                                              double J2, double J4) {
    return NumericalPropagator(epoch, initial,
                               {std::make_shared<GravityField>(std::vector<double>{0.0, 0.0, J2, 0.0, J4}, mu, re)});
}

OdeState<6> NumericalPropagator::derivative(double t, const OdeState<6>& y) const {
    const Vec3 r{y[0], y[1], y[2]}, v{y[3], y[4], y[5]};
    Vec3 a;
    for (const auto& f : forces_) a += f->acceleration(t, r, v);
    return {v.x, v.y, v.z, a.x, a.y, a.z};
}

StateVector NumericalPropagator::propagate(double dt) const {
    return propagate(std::vector<double>{dt}).front();
}

std::vector<StateVector> NumericalPropagator::propagate(const std::vector<double>& dts) const {
    std::vector<StateVector> out(dts.size());
    std::vector<std::size_t> order(dts.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return dts[a] < dts[b]; });
    auto f = [this](double t, const OdeState<6>& y) { return derivative(t, y); };
    const OdeState<6> y0{s0_.r.x, s0_.r.y, s0_.r.z, s0_.v.x, s0_.v.y, s0_.v.z};
    auto toState = [](const OdeState<6>& y) { return StateVector{{y[0], y[1], y[2]}, {y[3], y[4], y[5]}}; };
    AdaptiveStats total;
    auto run = [&](auto begin, auto end) {
        OdeState<6> y = y0;
        double t = 0.0;
        for (auto it = begin; it != end; ++it) {
            AdaptiveStats st;
            y = integrateRkf78<6>(f, t, y, dts[*it], opt_, &st);
            total.accepted += st.accepted;
            total.rejected += st.rejected;
            total.evaluations += st.evaluations;
            t = dts[*it];
            out[*it] = toState(y);
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
