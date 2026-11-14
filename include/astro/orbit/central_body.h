/**********************************************************************
 * central_body.h
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

#include <cstddef>
#include <vector>

#include "astro/core/constants.h"

/// \file
/// The gravity of a central body (gravitational parameter and zonal harmonics), shared by the
/// secular and numerical propagators, the gravity functions, the tracker and impact prediction.
namespace astro {

/// A central body's gravity: its gravitational parameter and its zonal harmonics.
///
/// The defaults are the Earth's, as described in core/constants.h: mu is constants::gmEarth (the
/// TCG-compatible value of IERS Conventions 2010, 7.5e-10 larger than EGM2008's own), re is
/// constants::earthRadiusEgm2008, and J2, J3, J4 are the EGM2008 tide-free values constants::J2,
/// constants::J3 and constants::J4. A default-constructed CentralBody equals CentralBody::earth().
///
/// Each user takes what its model needs: the numerical gravity field and the gravity functions use
/// every entry of J; the secular theory uses mu, re, zonal(2) and zonal(4) only.
struct CentralBody {
    double mu{constants::gmEarth};            ///< gravitational parameter (m^3/s^2)
    double re{constants::earthRadiusEgm2008}; ///< reference radius of the harmonics J (m)
    /// Unnormalized zonal coefficients indexed by degree: J[2] = J2, J[3] = J3, ... Entries 0 and 1
    /// are ignored (the point mass is always included through mu); an empty list is a point mass.
    std::vector<double> J{0.0, 0.0, constants::J2, constants::J3, constants::J4};

    /// The zonal coefficient Jn: J[n], or 0 for a degree beyond the list.
    [[nodiscard]] double zonal(std::size_t n) const noexcept { return n < J.size() ? J[n] : 0.0; }

    /// The Earth with EGM2008 J2, J3 and J4 (the default).
    [[nodiscard]] static CentralBody earth() { return {}; }
    /// The Earth with EGM2008 J2 and J4, without J3.
    [[nodiscard]] static CentralBody earthJ2J4() {
        return {
            constants::gmEarth,
            constants::earthRadiusEgm2008,
            {0.0, 0.0, constants::J2, 0.0, constants::J4}
        };
    }
    /// A point mass of gravitational parameter `mu` (m^3/s^2): no zonal harmonics (re keeps the
    /// Earth's value, which then plays no part).
    [[nodiscard]] static CentralBody pointMass(double mu = constants::gmEarth) {
        return {mu, constants::earthRadiusEgm2008, {}};
    }
};

} // namespace astro
