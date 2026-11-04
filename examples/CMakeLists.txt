###########################################################################
# examples/CMakeLists.txt
#
# Copyright 2024 W. Greg Visotski
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
###########################################################################

# Each example is a small program that the user guide quotes (between //! [tag] markers). They are
# built and run as tests so the documented code keeps working.
set(ASTRO_EXAMPLES time_scales geodesy frames sgp4_passes propagation covariance tracker sun_moon conventions)
foreach(name IN LISTS ASTRO_EXAMPLES)
    add_executable(astro_example_${name} ${name}.cpp)
    target_link_libraries(astro_example_${name} PRIVATE astro::astro)
    target_compile_options(astro_example_${name} PRIVATE ${ASTRO_WARNING_FLAGS})
    add_test(NAME example.${name} COMMAND astro_example_${name})
endforeach()

# Run the ephemeris example on the test kernel too, if it is in the source tree.
set(ASTRO_EXAMPLE_KERNEL ${PROJECT_SOURCE_DIR}/tests/data/de440s_2024q1_le.bsp)
if(EXISTS ${ASTRO_EXAMPLE_KERNEL})
    add_test(NAME example.sun_moon_spk COMMAND astro_example_sun_moon ${ASTRO_EXAMPLE_KERNEL})
endif()
