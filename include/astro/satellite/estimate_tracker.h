/**********************************************************************
 * estimate_tracker.h
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

#include <iosfwd>
#include <optional>
#include <span>
#include <vector>

#include <boost/describe/enum.hpp>

#include "astro/core/constants.h"
#include "astro/core/mat6.h"
#include "astro/core/vec.h"
#include "astro/orbit/integrators.h"
#include "astro/orbit/numerical.h"
#include "astro/time/time.h"

/// Following a stream of state estimates from an upstream filter.
///
/// Each delivery from a filter already contains everything before it, so successive deliveries are
/// correlated and must not be fused as independent measurements. EstimateTracker keeps the latest
/// one, propagates it (with its covariance) on request, and checks each new delivery against the
/// previous one propagated to the new epoch, to flag manoeuvres, bad data or unrealistic
/// covariances.
namespace astro {

/// A state estimate: inertial position and velocity (z axis along the Earth's pole, e.g. GCRS) and
/// their 6x6 covariance (m^2, m^2/s, m^2/s^2) at an epoch.
struct StateEstimate {
    Time        epoch;
    StateVector state;
    Mat6        covariance;
};

/// Which residual decides ConsistencyCheck::consistent. Both statistics are always reported.
enum class ConsistencyTest {
    /// Position residual against the predicted position covariance (3 degrees of freedom). Holds its
    /// nominal false-alarm rate even 4 h out with 100 m / 0.1 m/s priors (measured, LEO). Misses a
    /// velocity change until it has moved the position, e.g. a burn just before the delivery.
    Position,
    /// Position and velocity residual against the full predicted covariance (6 degrees of freedom).
    /// Also sees velocity changes, but Cartesian linear covariance propagation stops describing the
    /// error distribution once the propagated uncertainty grows: measured false-alarm rates at a
    /// nominal 1e-3 were 1e-3 up to 2 h and 0.06 at 4 h for 20 m / 2 cm/s priors, and 0.011 at
    /// 2 h for 100 m / 0.1 m/s priors.
    PositionVelocity,
};

BOOST_DESCRIBE_ENUM(ConsistencyTest, Position, PositionVelocity)

/// Writes "ConsistencyTest::<enumerator>", or "ConsistencyTest::INVALID_VALUE" for any other value.
std::ostream& operator<<(std::ostream& os, ConsistencyTest value);

/// Settings for propagating estimates.
struct TrackerOptions {
    /// Zonal coefficients by degree (J[2] = J2, ...) for the gravity model; point mass included.
    std::vector<double> J{0.0, 0.0, constants::J2, 0.0, constants::J4};
    double              mu{constants::gmEarth};
    double              re{constants::earthRadiusEgm2008};
    /// Integrator settings. The default (relTol 1e-9) keeps LEO positions within ~0.2 m over 4 h.
    AdaptiveOptions integrator{[] {
        AdaptiveOptions o;
        o.relTol = 1e-9;
        o.absTol = 1e-6;
        return o;
    }()};
    /// White acceleration process noise for the covariance (m^2/s^3); 0 for none. Tune it so that
    /// consistent deliveries pass the check (see ConsistencyCheck).
    double accelerationNoisePsd{0.0};
    /// A delivery is flagged inconsistent when the selected test's p-value falls below this.
    double          alarmProbability{1e-3};
    ConsistencyTest test{ConsistencyTest::Position};
};

/// How a new estimate compares with the previous one propagated to its epoch.
///
/// The distances use the predicted covariance only: d^2 = delta^T P_pred^-1 delta with
/// delta = new - predicted. If both came from the same consistent filter, delta has covariance at
/// most P_pred (the new estimate is better than the prediction), so the chi-square p-values are
/// conservative. A large d^2 means the new state is not where the previous estimate said it would
/// be. See ConsistencyTest for which statistic sets `consistent`.
struct ConsistencyCheck {
    double      dt{0.0};   ///< new epoch minus previous epoch (s)
    StateVector predicted; ///< previous estimate propagated to the new epoch
    Mat6        predictedCovariance{};
    Vec3        positionResidualRic;       ///< new - predicted position along radial, in-track, cross-track (m)
    Vec3        velocityResidualRic;       ///< new - predicted inertial velocity in the same axes (m/s)
    Vec3        positionSigmaRic;          ///< predicted 1-sigma position along R, I, C (m)
    Vec3        velocitySigmaRic;          ///< predicted 1-sigma velocity along R, I, C (m/s)
    double      mahalanobis2{0.0};         ///< 6-D squared distance
    double      positionMahalanobis2{0.0}; ///< position-only (3-D) squared distance
    double      pValue{1.0};               ///< P(chi^2 with 6 dof >= mahalanobis2)
    double      positionPValue{1.0};       ///< P(chi^2 with 3 dof >= positionMahalanobis2)
    bool        valid{false};              ///< false if the predicted covariance is not positive definite
    bool        consistent{false};         ///< valid and the selected test's p-value >= alarmProbability
    bool        replaced{false};           ///< whether the tracker took the new estimate (it is not older)
};

/// Survival function of the chi-square distribution with 3 and 6 degrees of freedom, P(X >= x).
[[nodiscard]] double chiSquareSurvival3(double x) noexcept;
[[nodiscard]] double chiSquareSurvival6(double x) noexcept;

/// The numerical propagator TrackerOptions describes, started from `e`.
[[nodiscard]] NumericalPropagator makePropagator(const StateEstimate& e, const TrackerOptions& options = {});

/// Propagate `prior` to next.epoch and compare it with `next` (see ConsistencyCheck).
[[nodiscard]] ConsistencyCheck
checkConsistency(const StateEstimate& prior, const StateEstimate& next, const TrackerOptions& options = {});

/// Latest-wins tracking of an upstream filter's estimates with a consistency check per delivery.
/// Not safe to use from several threads at once; use one tracker per thread.
class EstimateTracker {
public:
    explicit EstimateTracker(TrackerOptions options = {});

    /// Take a new estimate. If there is a current one, `next` is first checked against it; the
    /// result is returned (nullopt for the first estimate). `next` then replaces the current
    /// estimate, consistent or not, unless it is older than the current one (check.replaced).
    std::optional<ConsistencyCheck> update(const StateEstimate& next);

    [[nodiscard]] bool hasEstimate() const noexcept { return current_.has_value(); }
    /// The current estimate. Throws astro::Error if there is none yet.
    [[nodiscard]] const StateEstimate& current() const;

    /// The current estimate propagated to t, with its covariance.
    [[nodiscard]] StateCovariance predict(const Time& t) const;
    /// States and covariances at current().epoch + start + k * step for k = 0 .. out.size() - 1.
    void predictGrid(double start, double step, std::span<StateCovariance> out) const;
    [[nodiscard]] std::vector<StateCovariance> predictGrid(double start, double step, std::size_t count) const;

    [[nodiscard]] const TrackerOptions& options() const noexcept { return opt_; }

private:
    [[nodiscard]] const NumericalPropagator& propagator() const;

    TrackerOptions                     opt_;
    std::optional<StateEstimate>       current_;
    std::optional<NumericalPropagator> prop_;
};

} // namespace astro
