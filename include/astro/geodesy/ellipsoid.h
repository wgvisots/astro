/**********************************************************************
 * ellipsoid.h
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

#include "astro/core/constants.h"

namespace astro {

/// Oblate ellipsoid of revolution, defined by equatorial radius a (m) and flattening f.
struct Ellipsoid {
    double a{constants::wgs84A};
    double f{constants::wgs84F};

    [[nodiscard]] static constexpr Ellipsoid wgs84() noexcept { return {constants::wgs84A, constants::wgs84F}; }
    [[nodiscard]] static constexpr Ellipsoid grs80() noexcept { return {constants::grs80A, constants::grs80F}; }
    [[nodiscard]] static constexpr Ellipsoid iers2010() noexcept {
        return {constants::iers2010A, constants::iers2010F};
    }
    /// Sphere of radius r.
    [[nodiscard]] static constexpr Ellipsoid sphere(double r) noexcept { return {r, 0.0}; }

    /// Polar radius b = a (1 - f).
    [[nodiscard]] constexpr double b() const noexcept { return a * (1.0 - f); }
    /// First eccentricity squared e^2 = f (2 - f).
    [[nodiscard]] constexpr double e2() const noexcept { return f * (2.0 - f); }
    /// Second eccentricity squared e'^2 = e^2 / (1 - e^2).
    [[nodiscard]] constexpr double ep2() const noexcept { return e2() / (1.0 - e2()); }
    /// Third flattening n = f / (2 - f).
    [[nodiscard]] constexpr double n() const noexcept { return f / (2.0 - f); }

    /// Prime-vertical radius of curvature N(lat).
    [[nodiscard]] double primeVerticalRadius(double lat) const noexcept {
        const double s = std::sin(lat);
        return a / std::sqrt(1.0 - e2() * s * s);
    }
    /// Meridional radius of curvature M(lat).
    [[nodiscard]] double meridionalRadius(double lat) const noexcept {
        const double s = std::sin(lat);
        const double w2 = 1.0 - e2() * s * s;
        return a * (1.0 - e2()) / (w2 * std::sqrt(w2));
    }
    /// Gaussian mean radius sqrt(M N) at latitude.
    [[nodiscard]] double gaussianRadius(double lat) const noexcept {
        return std::sqrt(meridionalRadius(lat) * primeVerticalRadius(lat));
    }
    /// Geocentric radius of the ellipsoid surface at geodetic latitude.
    [[nodiscard]] double geocentricRadius(double lat) const noexcept {
        const double c = std::cos(lat), s = std::sin(lat);
        const double bb = b();
        const double num = (a * a * c) * (a * a * c) + (bb * bb * s) * (bb * bb * s);
        const double den = (a * c) * (a * c) + (bb * s) * (bb * s);
        return std::sqrt(num / den);
    }
    /// Mean radius R1 = (2a + b) / 3.
    [[nodiscard]] constexpr double meanRadius() const noexcept { return (2.0 * a + b()) / 3.0; }
};

} // namespace astro
