###########################################################################
# studies/CMakeLists.txt
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

# Accuracy studies behind docs/accuracy.md. They take seconds to minutes, so they are not tests.
add_executable(astro_study_propagation propagation_accuracy.cpp)
add_executable(astro_study_consistency consistency_monte_carlo.cpp)
add_executable(astro_study_impact impact_monte_carlo.cpp)
foreach(t astro_study_propagation astro_study_consistency astro_study_impact)
    target_link_libraries(${t} PRIVATE astro::astro)
    target_compile_options(${t} PRIVATE ${ASTRO_WARNING_FLAGS})
endforeach()
