/**********************************************************************
 * impact.cpp
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

#include "astro/orbit/impact.h"

#include <algorithm>
#include <cmath>

#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/geodesy/local.h"

namespace astro {

double ellipseScale(double probability) {
    if (!(probability > 0.0 && probability < 1.0)) {
        throw Error("ellipseScale: probability must be in (0, 1)");
    }
    return std::sqrt(-2.0 * std::log1p(-probability));
}

ErrorEllipse ErrorEllipse::scaled(double probability) const {
    const double k = ellipseScale(probability);
    return {k * semiMajor, k * semiMinor, rotation};
}

ErrorEllipse errorEllipse(double eastVariance, double northVariance, double eastNorthCovariance) noexcept {
    // Eigenvalues of [[ee, en], [en, nn]]: mean +- radius.
    const double mean   = 0.5 * (eastVariance + northVariance);
    const double half   = 0.5 * (eastVariance - northVariance);
    const double radius = std::hypot(half, eastNorthCovariance);
    ErrorEllipse e;
    e.semiMajor = std::sqrt(std::max(mean + radius, 0.0));
    e.semiMinor = std::sqrt(std::max(mean - radius, 0.0));
    if (radius > 0.0) {
        // Major axis counterclockwise from east, in [-pi/2, pi/2]; as an azimuth, clockwise from
        // north, in [0, pi].
        const double fromEast = 0.5 * std::atan2(2.0 * eastNorthCovariance, eastVariance - northVariance);
        const double azimuth  = constants::halfPi - fromEast;
        e.rotation            = azimuth >= constants::pi ? azimuth - constants::pi : azimuth;
    }
    return e;
}

namespace {

void checkOptions(const ImpactOptions& o, Frame frame) {
    if (isEarthFixed(frame)) {
        throw Error("predictImpact: the state must be in a celestial frame (not TIRS, PEF or ITRS)");
    }
    if (!(o.maxFlightTime > 0.0) || !std::isfinite(o.maxFlightTime)) {
        throw Error("predictImpact: maxFlightTime must be finite and positive");
    }
    if (!std::isfinite(o.surfaceHeight)) {
        throw Error("predictImpact: surfaceHeight must be finite");
    }
    if (!(o.accelerationNoisePsd >= 0.0) || !std::isfinite(o.accelerationNoisePsd)) {
        throw Error("predictImpact: accelerationNoisePsd must be finite and non-negative");
    }
}

std::optional<Impact> impact(const Time&          epoch,
                             const StateVector&   state,
                             const Mat6*          covariance,
                             Frame                frame,
                             const EopProvider&   eop,
                             const ImpactOptions& o) {
    checkOptions(o, frame);

    // Propagate in GCRS, with the zonal field about the true pole at the epoch.
    const Mat3        toFrame0 = gcrsToCelestialMatrix(frame, epoch); // GCRS -> frame at the epoch
    const Mat3        toGcrs0  = toFrame0.transpose();
    const StateVector s0{toGcrs0 * state.r, toGcrs0 * state.v};
    std::vector<std::shared_ptr<const ForceModel>> forces{
        std::make_shared<GravityField>(o.body, poleOfDateMatrix(Frame::GCRS, epoch))};
    forces.insert(forces.end(), o.additionalForces.begin(), o.additionalForces.end());
    const NumericalPropagator prop(epoch, s0, forces, o.integrator);

    // GCRS -> ITRS over the flight: the caller's series, or one for this call.
    std::optional<ItrsRotationSeries> own;
    const ItrsRotationSeries*         rot    = o.rotations;
    double                            offset = 0.0; // series time of the epoch
    if (rot != nullptr) {
        offset = epoch - rot->start();
        if (offset < 0.0 || offset + o.maxFlightTime > rot->duration()) {
            throw Error("predictImpact: ImpactOptions::rotations must cover [epoch, epoch + maxFlightTime]");
        }
    } else {
        own.emplace(epoch, o.maxFlightTime, eop, o.frameOptions);
        rot = &*own;
    }

    const auto height = [&](double dt, const StateVector& s) {
        return ecefToGeodetic(rot->gcrsToItrs(s.r, offset + dt), o.ellipsoid).h - o.surfaceHeight;
    };
    const auto event = prop.propagateToEvent(height, o.maxFlightTime, Crossing::Falling, o.timeTolerance);
    if (!event) {
        return std::nullopt;
    }

    Impact out;
    out.flightTime     = event->dt;
    out.time           = epoch + event->dt;
    out.itrs           = rot->gcrsToItrs(event->state, offset + event->dt);
    out.point          = ecefToGeodetic(out.itrs.r, o.ellipsoid);
    const Mat3 enu     = ecefToEnuMatrix(out.point.lat, out.point.lon); // rows: east, north, up
    const Vec3 v       = enu * out.itrs.v;
    out.speed          = norm(out.itrs.v);
    out.angle          = std::atan2(-v.z, std::hypot(v.x, v.y));
    out.heading        = wrapTwoPi(std::atan2(v.x, v.y));
    const Mat3 toFrame = gcrsToCelestialMatrix(frame, out.time);
    out.state          = {toFrame * event->state.r, toFrame * event->state.v};

    if (covariance != nullptr) {
        // Covariance at the nominal impact time, in GCRS.
        const StateCovariance pc =
            prop.propagateWithCovariance(event->dt, rotateCovariance(*covariance, toGcrs0), o.accelerationNoisePsd);
        const Mat6& P   = pc.covariance;
        const Mat3  Prr = block(P, 0, 0);

        // A position error dr (GCRS) at the nominal time moves the impact point (ITRS) by
        //     dp = (I - v u^T / (u . v)) R dr,
        // R = GCRS -> ITRS, u the local up (the gradient of the height), v the ground-relative
        // velocity: the trajectory slides along v until its height is zero again, which takes
        //     dt = -(u . R dr) / (u . v).
        // dp lies in the horizontal plane; velocity errors only enter at second order.
        const Mat3   R     = rot->gcrsToItrsMatrix(offset + event->dt);
        const Vec3   up    = Vec3{enu(2, 0), enu(2, 1), enu(2, 2)};
        const double vUp   = dot(up, out.itrs.v); // negative: falling
        const Mat3   slide = Mat3::identity() - outer(out.itrs.v, up) * (1.0 / vUp);
        const Mat3   A     = enu * slide * R; // GCRS position error -> east, north, (up = 0)
        const Mat3   C     = A * Prr * A.transpose();
        const Vec3   gradT = R.transpose() * up / vUp; // dt = -gradT . dr

        ImpactUncertainty u;
        u.covariance          = rotateCovariance(P, toFrame);
        u.eastVariance        = C(0, 0);
        u.northVariance       = C(1, 1);
        u.eastNorthCovariance = 0.5 * (C(0, 1) + C(1, 0));
        u.ellipse             = errorEllipse(u.eastVariance, u.northVariance, u.eastNorthCovariance);
        u.timeSigma           = std::sqrt(std::max(dot(gradT, Prr * gradT), 0.0));
        out.uncertainty       = u;
    }
    return out;
}

} // namespace

std::optional<Impact> predictImpact(const Time&          epoch,
                                    const StateVector&   state,
                                    Frame                frame,
                                    const EopProvider&   eop,
                                    const ImpactOptions& options) {
    return impact(epoch, state, nullptr, frame, eop, options);
}

std::optional<Impact> predictImpact(const Time&          epoch,
                                    const StateVector&   state,
                                    const Mat6&          covariance,
                                    Frame                frame,
                                    const EopProvider&   eop,
                                    const ImpactOptions& options) {
    return impact(epoch, state, &covariance, frame, eop, options);
}

} // namespace astro
