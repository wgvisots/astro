###########################################################################
# CMakeLists.txt
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

cmake_minimum_required(VERSION 3.20)

project(astro
        VERSION 0.1.0
        DESCRIPTION "Time, geodesy, reference frame, ephemeris and orbit utilities"
        LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

option(ASTRO_BUILD_TESTS "Build the astro unit tests" ${PROJECT_IS_TOP_LEVEL})
option(ASTRO_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(ASTRO_BUILD_BENCHMARKS "Build the astro benchmarks" OFF)

include(GNUInstallDirs)

file(GLOB_RECURSE ASTRO_SOURCES CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB_RECURSE ASTRO_HEADERS CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/include/astro/*.h)

add_library(astro ${ASTRO_SOURCES} ${ASTRO_HEADERS})
add_library(astro::astro ALIAS astro)

target_compile_features(astro PUBLIC cxx_std_20)
set_target_properties(astro PROPERTIES
        CXX_EXTENSIONS OFF
        POSITION_INDEPENDENT_CODE ON
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR})

target_include_directories(astro PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# Boost.Describe (header-only, Boost >= 1.77) provides enum reflection; it is the only dependency.
find_package(Boost 1.77 CONFIG REQUIRED)
target_link_libraries(astro PUBLIC Boost::headers)

# Warning flags shared by the library and the tests (-Wall -Wextra -Wpedantic).
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(ASTRO_WARNING_FLAGS -Wall -Wextra -Wpedantic)
    set(ASTRO_WERROR_FLAG -Werror)
elseif(MSVC)
    set(ASTRO_WARNING_FLAGS /W4 /permissive-)
    set(ASTRO_WERROR_FLAG /WX)
endif()
if(ASTRO_WARNINGS_AS_ERRORS)
    list(APPEND ASTRO_WARNING_FLAGS ${ASTRO_WERROR_FLAG})
endif()

target_compile_options(astro PRIVATE ${ASTRO_WARNING_FLAGS})
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # Extra checks for the library sources only.
    target_compile_options(astro PRIVATE -Wshadow -Wconversion -Wno-sign-conversion)
endif()

# ---------------------------------------------------------------- install
install(TARGETS astro EXPORT astroTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/astro DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT astroTargets
        NAMESPACE astro::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/astro)

include(CMakePackageConfigHelpers)
file(WRITE ${PROJECT_BINARY_DIR}/astroConfig.cmake.in
        "@PACKAGE_INIT@\ninclude(CMakeFindDependencyMacro)\nfind_dependency(Boost 1.77 CONFIG)\ninclude(\"\${CMAKE_CURRENT_LIST_DIR}/astroTargets.cmake\")\n")
configure_package_config_file(${PROJECT_BINARY_DIR}/astroConfig.cmake.in
        ${PROJECT_BINARY_DIR}/astroConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/astro)
write_basic_package_version_file(${PROJECT_BINARY_DIR}/astroConfigVersion.cmake
        COMPATIBILITY SameMajorVersion)
install(FILES ${PROJECT_BINARY_DIR}/astroConfig.cmake ${PROJECT_BINARY_DIR}/astroConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/astro)

# ---------------------------------------------------------------- tests
if(ASTRO_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(ASTRO_BUILD_BENCHMARKS)
    add_subdirectory(bench)
endif()
