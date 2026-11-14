/**********************************************************************
 * conventions.cpp
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

// Library conventions: units, enums with Boost.Describe and operator<<, and errors.

#include <boost/describe/enumerators.hpp>
#include <boost/mp11/algorithm.hpp>

#include <iostream>

#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/frames/transform.h"
#include "astro/time/time.h"

using namespace astro;

int main() {
    //! [enums]
    // Every public enum prints as Enum::Value and is reflected with Boost.Describe.
    std::cout << TimeScale::TT << '\n';          // TimeScale::TT
    std::cout << static_cast<Frame>(42) << '\n'; // Frame::INVALID_VALUE
    boost::mp11::mp_for_each<boost::describe::describe_enumerators<Frame>>(
        [](auto d) { std::cout << d.name << (isEarthFixed(d.value) ? " (Earth-fixed)\n" : "\n"); });
    //! [enums]

    //! [errors]
    // Bad input and missing data throw astro::Error (a std::runtime_error).
    try {
        (void)Time::fromIso8601("not a time");
    } catch (const Error& e) {
        std::cout << "astro::Error: " << e.what() << '\n';
    }
    //! [errors]

    //! [units]
    // SI units and radians everywhere; convert at the edges.
    const double lat = degToRad(38.8895);
    std::cout << formatDms(lat) << '\n';
    //! [units]
    return 0;
}
