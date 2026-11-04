###########################################################################
# tests/CMakeLists.txt
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

find_package(GTest QUIET)
if(NOT GTest_FOUND)
    include(FetchContent)
    FetchContent_Declare(googletest
            URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.tar.gz
            DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
endif()

include(GoogleTest)

file(GLOB_RECURSE ASTRO_TEST_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/*_test.cpp)

add_executable(astro_tests ${ASTRO_TEST_SOURCES})
find_package(Threads REQUIRED)
target_link_libraries(astro_tests PRIVATE astro::astro GTest::gtest_main Threads::Threads)
target_include_directories(astro_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(astro_tests PRIVATE ASTRO_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
target_compile_options(astro_tests PRIVATE ${ASTRO_WARNING_FLAGS})

gtest_discover_tests(astro_tests DISCOVERY_TIMEOUT 60)
