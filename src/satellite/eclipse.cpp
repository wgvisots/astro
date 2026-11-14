/**********************************************************************
 * eclipse.cpp
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

#include "astro/satellite/eclipse.h"

#include <algorithm>
#include <cmath>

#include "../core/enum_io.h"
#include "astro/core/error.h"
#include "astro/core/rootfind.h"

namespace astro {

namespace {

// Apparent radii of Sun (a) and Earth (b), and their angular separation (c), seen from the satellite.
struct DiskGeometry {
    double a, b, c;
};

DiskGeometry disks(const Vec3& sat, const Vec3& sun, double re, double rs) {
    const Vec3   toSun = sun - sat;
    const double ds = norm(toSun), dr = norm(sat);
    return {std::asin(std::min(1.0, rs / ds)),
            std::asin(std::min(1.0, re / dr)),
            std::acos(std::clamp(dot(-sat, toSun) / (dr * ds), -1.0, 1.0))};
}

// Signed functions whose zeros are the shadow boundaries (negative = inside).
double penumbraFunction(const Vec3& sat, const Vec3& sun, ShadowModel model, double re, double rs) {
    if (model == ShadowModel::Cylindrical) {
        const Vec3   s     = normalized(sun);
        const double along = dot(sat, s);
        const double perp  = norm(sat - along * s);
        // Outside if on the sunward side; otherwise compare the distance from the shadow axis.
        return along > 0.0 ? perp + along : perp - re;
    }
    const DiskGeometry g = disks(sat, sun, re, rs);
    return g.c - (g.a + g.b);
}

double umbraFunction(const Vec3& sat, const Vec3& sun, ShadowModel model, double re, double rs) {
    if (model == ShadowModel::Cylindrical) {
        return penumbraFunction(sat, sun, model, re, rs);
    }
    const DiskGeometry g = disks(sat, sun, re, rs);
    return g.c - (g.b - g.a);
}

} // namespace

ShadowResult earthShadow(const Vec3& sat, const Vec3& sun, ShadowModel model, double re, double rs) noexcept {
    if (model == ShadowModel::Cylindrical) {
        return penumbraFunction(sat, sun, model, re, rs) < 0.0 ? ShadowResult{ShadowState::Umbra, 0.0}
                                                               : ShadowResult{ShadowState::Sunlit, 1.0};
    }
    const DiskGeometry g = disks(sat, sun, re, rs);
    const double       a = g.a, b = g.b, c = g.c;
    if (c >= a + b) {
        return {ShadowState::Sunlit, 1.0};
    }
    if (c < b - a) {
        return {ShadowState::Umbra, 0.0};
    }
    if (c < a - b) {
        return {ShadowState::Annular, 1.0 - (b * b) / (a * a)};
    }
    // Partial overlap of the two disks.
    const double x = (c * c + a * a - b * b) / (2.0 * c);
    const double y = std::sqrt(std::max(0.0, a * a - x * x));
    const double area =
        a * a * std::acos(std::clamp(x / a, -1.0, 1.0)) + b * b * std::acos(std::clamp((c - x) / b, -1.0, 1.0)) - c * y;
    return {ShadowState::Penumbra, std::clamp(1.0 - area / (constants::pi * a * a), 0.0, 1.0)};
}

double betaAngle(const StateVector& sat, const Vec3& sun) noexcept {
    return std::asin(std::clamp(dot(normalized(cross(sat.r, sat.v)), normalized(sun)), -1.0, 1.0));
}

std::vector<EclipseInterval> findEclipses(const std::function<Vec3(const Time&)>& satPos,
                                          const std::function<Vec3(const Time&)>& sunPos,
                                          const Time&                             start,
                                          const Time&                             end,
                                          double                                  step,
                                          ShadowModel                             model,
                                          double                                  re) {
    const double rs   = constants::sunRadius;
    const double span = end - start;
    if (!std::isfinite(span) || span < 0.0) {
        throw Error("findEclipses: start and end must be finite, with end not before start");
    }
    if (!(step > 0.0) || !std::isfinite(step)) {
        throw Error("findEclipses: step must be positive and finite");
    }
    if (!(span / step <= 2e9)) { // the count of steps is kept in an int
        throw Error("findEclipses: the window holds more than 2e9 steps");
    }
    std::vector<EclipseInterval> out;
    if (span == 0.0) {
        return out; // an empty window holds no eclipse
    }
    auto fp = [&](double s) {
        const Time t = start + s;
        return penumbraFunction(satPos(t), sunPos(t), model, re, rs);
    };
    auto fu = [&](double s) {
        const Time t = start + s;
        return umbraFunction(satPos(t), sunPos(t), model, re, rs);
    };
    const int    n = std::max(2, static_cast<int>(std::ceil(span / step)));
    const double h = span / n;

    // Crossing times of a boundary function.
    auto crossings = [&](auto&& f) {
        std::vector<std::pair<double, bool>> c; // time, entering
        double                               prev = f(0.0);
        for (int i = 1; i <= n; ++i) {
            const double cur = f(i * h);
            if ((prev < 0.0) != (cur < 0.0)) {
                c.push_back({brentRoot(f, (i - 1) * h, i * h, 1e-3), cur < 0.0});
            }
            prev = cur;
        }
        return c;
    };
    const auto pen = crossings(fp);
    const auto umb = crossings(fu);

    // Assemble intervals from penumbra entries/exits; attach the umbra crossings inside each.
    double entry = fp(0.0) < 0.0 ? 0.0 : -1.0;
    auto   close = [&](double exitS) {
        EclipseInterval e;
        e.penumbraEntry = start + entry;
        e.penumbraExit  = start + exitS;
        double     ue = -1.0, ux = -1.0;
        const bool inUmbraAtEntry = fu(entry) < 0.0;
        if (inUmbraAtEntry) {
            ue = entry;
        }
        for (const auto& [t, entering] : umb) {
            if (t < entry || t > exitS) {
                continue;
            }
            if (entering && ue < 0.0) {
                ue = t;
            }
            if (!entering) {
                ux = t;
            }
        }
        if (ue >= 0.0 && ux < 0.0) {
            ux = exitS;
        }
        e.umbral     = ue >= 0.0;
        e.umbraEntry = start + (e.umbral ? ue : entry);
        e.umbraExit  = start + (e.umbral ? ux : exitS);
        out.push_back(e);
    };
    for (const auto& [t, entering] : pen) {
        if (entering) {
            entry = t;
        } else if (entry >= 0.0) {
            close(t);
            entry = -1.0;
        }
    }
    if (entry >= 0.0) {
        close(span);
    }
    return out;
}

std::ostream& operator<<(std::ostream& os, ShadowModel value) {
    return detail::writeEnum(os, "ShadowModel", value);
}

std::ostream& operator<<(std::ostream& os, ShadowState value) {
    return detail::writeEnum(os, "ShadowState", value);
}

} // namespace astro
