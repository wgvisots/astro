/**********************************************************************
 * integrators.h
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

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

#include "astro/core/error.h"

/// Explicit Runge-Kutta ODE integrators (header-only) for y' = f(t, y),
/// y an std::array<double, N>.
namespace astro {

template <std::size_t N>
using OdeState = std::array<double, N>;

namespace detail {
template <std::size_t N>
OdeState<N> axpy(const OdeState<N>& y, double h, std::initializer_list<std::pair<double, const OdeState<N>*>> terms) {
    OdeState<N> out = y;
    for (const auto& [c, k] : terms) {
        if (c != 0.0) {
            for (std::size_t i = 0; i < N; ++i) {
                out[i] += h * c * (*k)[i];
            }
        }
    }
    return out;
}
} // namespace detail

/// One classical 4th-order Runge-Kutta step.
template <std::size_t N, class F>
[[nodiscard]] OdeState<N> rk4Step(F&& f, double t, const OdeState<N>& y, double h) {
    const OdeState<N> k1  = f(t, y);
    const OdeState<N> k2  = f(t + 0.5 * h,
                             detail::axpy<N>(y,
                                             h,
                                              {
                                                 {0.5, &k1}
    }));
    const OdeState<N> k3  = f(t + 0.5 * h,
                             detail::axpy<N>(y,
                                             h,
                                              {
                                                 {0.5, &k2}
    }));
    const OdeState<N> k4  = f(t + h,
                             detail::axpy<N>(y,
                                             h,
                                              {
                                                 {1.0, &k3}
    }));
    OdeState<N>       out = y;
    for (std::size_t i = 0; i < N; ++i) {
        out[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
    return out;
}

/// Integrate with fixed-step RK4 from t0 to t1 (the last step is shortened to land on t1).
template <std::size_t N, class F>
[[nodiscard]] OdeState<N> integrateRk4(F&& f, double t0, OdeState<N> y, double t1, double step) {
    if (step <= 0.0) {
        throw Error("integrateRk4: step must be positive");
    }
    const double dir = t1 >= t0 ? 1.0 : -1.0;
    double       t   = t0;
    while (dir * (t1 - t) > 0.0) {
        const double h  = dir * std::min(step, dir * (t1 - t));
        y               = rk4Step<N>(f, t, y, h);
        t              += h;
    }
    return y;
}

struct AdaptiveOptions {
    double relTol{1e-12};
    double absTol{1e-9};
    double initialStep{60.0};
    double minStep{1e-6};
    double maxStep{1e30};
    int    maxSteps{10000000};
};

struct AdaptiveStats {
    int accepted{0};
    int rejected{0};
    int evaluations{0};
};

/// Runge-Kutta-Fehlberg 7(8) adaptive integration from t0 to t1, landing exactly on t1.
/// Error control per component: |err_i| <= absTol + relTol * max(|y_i|, |y_new_i|).
template <std::size_t N, class F>
[[nodiscard]] OdeState<N> integrateRkf78(F&&                    f,
                                         double                 t0,
                                         OdeState<N>            y,
                                         double                 t1,
                                         const AdaptiveOptions& opt   = {},
                                         AdaptiveStats*         stats = nullptr) {
    // Fehlberg (1968) 7(8) coefficients.
    // clang-format off
    static constexpr double c[13] = {0.0, 2.0 / 27.0, 1.0 / 9.0, 1.0 / 6.0, 5.0 / 12.0, 0.5, 5.0 / 6.0,
                                     1.0 / 6.0, 2.0 / 3.0, 1.0 / 3.0, 1.0, 0.0, 1.0};
    static constexpr double a[13][12] = {
        {},
        {2.0 / 27.0},
        {1.0 / 36.0, 1.0 / 12.0},
        {1.0 / 24.0, 0.0, 1.0 / 8.0},
        {5.0 / 12.0, 0.0, -25.0 / 16.0, 25.0 / 16.0},
        {1.0 / 20.0, 0.0, 0.0, 1.0 / 4.0, 1.0 / 5.0},
        {-25.0 / 108.0, 0.0, 0.0, 125.0 / 108.0, -65.0 / 27.0, 125.0 / 54.0},
        {31.0 / 300.0, 0.0, 0.0, 0.0, 61.0 / 225.0, -2.0 / 9.0, 13.0 / 900.0},
        {2.0, 0.0, 0.0, -53.0 / 6.0, 704.0 / 45.0, -107.0 / 9.0, 67.0 / 90.0, 3.0},
        {-91.0 / 108.0, 0.0, 0.0, 23.0 / 108.0, -976.0 / 135.0, 311.0 / 54.0, -19.0 / 60.0, 17.0 / 6.0, -1.0 / 12.0},
        {2383.0 / 4100.0, 0.0, 0.0, -341.0 / 164.0, 4496.0 / 1025.0, -301.0 / 82.0, 2133.0 / 4100.0, 45.0 / 82.0,
         45.0 / 164.0, 18.0 / 41.0},
        {3.0 / 205.0, 0.0, 0.0, 0.0, 0.0, -6.0 / 41.0, -3.0 / 205.0, -3.0 / 41.0, 3.0 / 41.0, 6.0 / 41.0, 0.0},
        {-1777.0 / 4100.0, 0.0, 0.0, -341.0 / 164.0, 4496.0 / 1025.0, -289.0 / 82.0, 2193.0 / 4100.0, 51.0 / 82.0,
         33.0 / 164.0, 12.0 / 41.0, 0.0, 1.0},
    };
    // 8th-order weights (propagated solution); the error estimate is 41/840 (k1 + k11 - k12 - k13).
    static constexpr double b8[13] = {0.0, 0.0, 0.0, 0.0, 0.0,
                                      34.0 / 105.0, 9.0 / 35.0, 9.0 / 35.0, 9.0 / 280.0, 9.0 / 280.0,
                                      0.0, 41.0 / 840.0, 41.0 / 840.0};
    // clang-format on

    if (t1 == t0) {
        return y;
    }
    const double                dir = t1 > t0 ? 1.0 : -1.0;
    double                      h   = dir * std::min(std::abs(opt.initialStep), std::abs(t1 - t0));
    double                      t   = t0;
    std::array<OdeState<N>, 13> k;
    AdaptiveStats               st;
    for (int step = 0; step < opt.maxSteps; ++step) {
        if (dir * (t + h - t1) > 0.0) {
            h = t1 - t;
        }
        for (int s = 0; s < 13; ++s) {
            OdeState<N> ys = y;
            for (int j = 0; j < s; ++j) {
                if (a[s][j] != 0.0) {
                    for (std::size_t i = 0; i < N; ++i) {
                        ys[i] += h * a[s][j] * k[static_cast<std::size_t>(j)][i];
                    }
                }
            }
            k[static_cast<std::size_t>(s)] = f(t + c[s] * h, ys);
        }
        st.evaluations      += 13;
        OdeState<N> yn       = y;
        double      errNorm  = 0.0;
        for (std::size_t i = 0; i < N; ++i) {
            double acc = 0.0;
            for (int s = 0; s < 13; ++s) {
                acc += b8[s] * k[static_cast<std::size_t>(s)][i];
            }
            yn[i]            += h * acc;
            const double err  = h * 41.0 / 840.0 * (k[0][i] + k[10][i] - k[11][i] - k[12][i]);
            const double sc   = opt.absTol + opt.relTol * std::max(std::abs(y[i]), std::abs(yn[i]));
            errNorm           = std::max(errNorm, std::abs(err) / sc);
        }
        if (errNorm <= 1.0) {
            t += h;
            y  = yn;
            ++st.accepted;
            if (dir * (t1 - t) <= 0.0) {
                if (stats) {
                    *stats = st;
                }
                return y;
            }
        } else {
            ++st.rejected;
        }
        const double factor = errNorm == 0.0 ? 4.0 : std::clamp(0.9 * std::pow(errNorm, -1.0 / 8.0), 0.2, 4.0);
        h                   = dir * std::clamp(std::abs(h) * factor, opt.minStep, opt.maxStep);
        if (std::abs(h) <= opt.minStep && errNorm > 1.0) {
            throw Error("integrateRkf78: step size underflow");
        }
    }
    throw Error("integrateRkf78: maximum number of steps exceeded");
}

} // namespace astro
