/**********************************************************************
 * close_approach.cpp
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

#include "astro/orbit/close_approach.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <optional>
#include <utility>

#include "../core/enum_io.h"
#include "astro/core/error.h"
#include "astro/core/rootfind.h"
#include "astro/orbit/two_body.h"

namespace astro {

std::ostream& operator<<(std::ostream& os, TrajectoryModel value) {
    return detail::writeEnum(os, "TrajectoryModel", value);
}

namespace {

/// Shortest sampling step (s), so that a trajectory through the centre (|r| = 0) still advances.
constexpr double kMinStep = 1e-3;

/// The SecularOrder of a secular model; throws astro::Error for any other value.
SecularOrder secularOrder(TrajectoryModel model) {
    switch (model) {
    case TrajectoryModel::SecularJ2:   return SecularOrder::J2;
    case TrajectoryModel::SecularJ2J4: return SecularOrder::J2SquaredJ4;
    case TrajectoryModel::TwoBody:     break;
    }
    throw Error("Trajectory: invalid TrajectoryModel");
}

/// Two-body propagation from s0: seconds after the epoch -> state.
std::function<StateVector(double)> twoBody(const StateVector& s0, double mu) {
    return [s0, mu](double dt) { return propagateKepler(s0, dt, mu); };
}

/// Secular propagation with velocity = d(position)/dt, which the range rate needs (see
/// SecularPropagatorCorrected). Shared, so copies of the function stay cheap.
std::function<StateVector(double)> secular(std::shared_ptr<const SecularPropagatorCorrected> p) {
    return [p = std::move(p)](double dt) { return p->propagate(dt); };
}

[[nodiscard]] bool isFinite(const StateVector& s) noexcept {
    return std::isfinite(s.r.x) && std::isfinite(s.r.y) && std::isfinite(s.r.z) && std::isfinite(s.v.x)
        && std::isfinite(s.v.y) && std::isfinite(s.v.z);
}

/// Both states at s seconds after the window start, and the range rate there.
struct Sample {
    double      s{0.0}; // seconds after the window start
    StateVector a;      // first object
    StateVector b;      // second object
    double      g{0.0}; // (rb - ra) . (vb - va): negative while closing (m^2/s)
};

/// The sampled range-rate search shared by closestApproach and closeApproaches.
class Search {
public:
    Search(const Trajectory& first, const Trajectory& second, const Time& start, const CloseApproachOptions& opt)
        : fa_(first.secondsAfter(start))
        , fb_(second.secondsAfter(start))
        , opt_(opt) {}

    [[nodiscard]] Sample at(double s) const {
        Sample x{s, fa_(s), fb_(s), 0.0};
        if (!isFinite(x.a) || !isFinite(x.b)) {
            throw Error("closestApproach: a trajectory returned a non-finite state");
        }
        x.g = dot(x.b.r - x.a.r, x.b.v - x.a.v);
        return x;
    }

    /// Samples [0, T] and calls onMinimum with each refined local minimum, in time order. Returns
    /// the samples at 0 and T.
    template <class OnMinimum>
    std::pair<Sample, Sample> run(double T, OnMinimum&& onMinimum) const {
        const Sample first = at(0.0);
        Sample       prev  = first;
        while (prev.s < T) {
            Sample cur = at(std::min(prev.s + step(prev), T));
            if (prev.g < 0.0 && cur.g >= 0.0) { // closing, then separating: a minimum in (prev, cur]
                const double tca = brentRoot([this](double s) { return at(s).g; }, prev.s, cur.s, opt_.timeTolerance);
                onMinimum(at(tca));
            }
            prev = std::move(cur);
        }
        return {first, prev};
    }

private:
    /// |r| / |v|: the time scale on which the object's position turns (infinite at rest).
    [[nodiscard]] static double timeScale(const StateVector& s) noexcept {
        const double v = norm(s.v);
        return v > 0.0 ? norm(s.r) / v : std::numeric_limits<double>::infinity();
    }

    [[nodiscard]] double step(const Sample& x) const noexcept {
        return std::max(opt_.stepFraction * std::min(timeScale(x.a), timeScale(x.b)), kMinStep);
    }

    std::function<StateVector(double)> fa_;
    std::function<StateVector(double)> fb_;
    CloseApproachOptions               opt_;
};

[[nodiscard]] CloseApproach approach(const Sample& x, const Time& start, bool atBoundary) {
    CloseApproach c;
    c.time          = start + x.s;
    c.distance      = norm(x.b.r - x.a.r);
    c.relativeSpeed = norm(x.b.v - x.a.v);
    c.first         = x.a;
    c.second        = x.b;
    c.atBoundary    = atBoundary;
    return c;
}

/// Window length (s), after checking the window and the options.
[[nodiscard]] double windowLength(const Time& start, const Time& end, const CloseApproachOptions& opt) {
    if (!(opt.stepFraction > 0.0) || !std::isfinite(opt.stepFraction)) {
        throw Error("closestApproach: stepFraction must be positive and finite");
    }
    if (!(opt.timeTolerance > 0.0) || !std::isfinite(opt.timeTolerance)) {
        throw Error("closestApproach: timeTolerance must be positive and finite");
    }
    const double T = end - start;
    if (!(T >= 0.0) || !std::isfinite(T)) {
        throw Error("closestApproach: end must not be before start");
    }
    return T;
}

} // namespace

Trajectory::Trajectory(const Time& epoch, const StateVector& state, TrajectoryModel model, const ZonalParameters& zp)
    : epoch_(epoch) {
    if (model == TrajectoryModel::TwoBody) {
        fromEpoch_ = twoBody(state, zp.mu);
    } else {
        fromEpoch_ = secular(std::make_shared<const SecularPropagatorCorrected>(state, zp, secularOrder(model)));
    }
}

Trajectory::Trajectory(const Time&              epoch,
                       const KeplerianElements& elements,
                       TrajectoryModel          model,
                       const ZonalParameters&   zp)
    : epoch_(epoch) {
    if (model == TrajectoryModel::TwoBody) {
        fromEpoch_ = twoBody(elementsToState(elements, zp.mu), zp.mu);
    } else {
        fromEpoch_ = secular(std::make_shared<const SecularPropagatorCorrected>(elements, zp, secularOrder(model)));
    }
}

Trajectory::Trajectory(const Time&                epoch,
                       const EquinoctialElements& elements,
                       TrajectoryModel            model,
                       const ZonalParameters&     zp)
    : Trajectory(epoch, equinoctialToKeplerian(elements), model, zp) {}

Trajectory::Trajectory(StateFunction stateAt)
    : stateAt_(std::move(stateAt)) {
    if (!stateAt_) {
        throw Error("Trajectory: empty state function");
    }
}

StateVector Trajectory::state(const Time& t) const {
    return fromEpoch_ ? fromEpoch_(t - epoch_) : stateAt_(t);
}

std::function<StateVector(double)> Trajectory::secondsAfter(const Time& origin) const {
    if (fromEpoch_) {
        return [f = fromEpoch_, offset = origin - epoch_](double s) { return f(offset + s); };
    }
    return [f = stateAt_, origin](double s) { return f(origin + s); };
}

CloseApproach closestApproach(const Trajectory&           first,
                              const Trajectory&           second,
                              const Time&                 start,
                              const Time&                 end,
                              const CloseApproachOptions& options) {
    const double                 T = windowLength(start, end, options);
    const Search                 search(first, second, start, options);
    std::optional<CloseApproach> best;
    const auto                   consider = [&best](CloseApproach c) {
        if (!best || c.distance < best->distance) {
            best = std::move(c);
        }
    };
    // Minima first, so a true closest approach wins a tie with a window end.
    const auto [s0, sT] = search.run(T, [&](const Sample& x) { consider(approach(x, start, false)); });
    consider(approach(s0, start, true));
    consider(approach(sT, start, true));
    return *best;
}

std::vector<CloseApproach> closeApproaches(const Trajectory&           first,
                                           const Trajectory&           second,
                                           const Time&                 start,
                                           const Time&                 end,
                                           double                      maxDistance,
                                           const CloseApproachOptions& options) {
    const double T = windowLength(start, end, options);
    if (std::isnan(maxDistance)) {
        throw Error("closeApproaches: maxDistance is NaN");
    }
    const Search               search(first, second, start, options);
    std::vector<CloseApproach> out;
    (void)search.run(T, [&](const Sample& x) {
        CloseApproach c = approach(x, start, false);
        if (c.distance <= maxDistance) {
            out.push_back(std::move(c));
        }
    });
    return out;
}

} // namespace astro
