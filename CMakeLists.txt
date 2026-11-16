###########################################################################
# CMakeLists.txt
#
# Copyright 2024-2026 W. Greg Visotski
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
        VERSION 0.4.0
        DESCRIPTION "Time, geodesy, reference frame, ephemeris and orbit utilities"
        LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

option(ASTRO_BUILD_TESTS "Build the astro unit tests" ${PROJECT_IS_TOP_LEVEL})
option(ASTRO_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(ASTRO_BUILD_BENCHMARKS "Build the astro benchmarks" OFF)
option(ASTRO_BUILD_STUDIES "Build the accuracy studies behind docs/accuracy.md" OFF)
option(ASTRO_DOCS_WARNINGS_AS_ERRORS "Fail the docs target on Doxygen warnings" OFF)
option(ASTRO_BUILD_EXAMPLES "Build (and test) the astro examples" ${PROJECT_IS_TOP_LEVEL})
option(ASTRO_BUILD_FUZZERS "Build libFuzzer targets for the parsers (needs Clang)" OFF)
option(ASTRO_COVERAGE "Instrument the build for line coverage (GCC; see tools/coverage_report.py)" OFF)

include(GNUInstallDirs)

# Windows finds a DLL next to the program that loads it, not through a run path. So that the tests and
# examples run against a shared astro (and gtest_discover_tests can list the tests at build time),
# every program and DLL of the build goes to one directory.
if(WIN32 AND PROJECT_IS_TOP_LEVEL AND NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
endif()

if(ASTRO_COVERAGE)
    # Counts how often each line runs (gcov). Use an unoptimized build so that lines map to code.
    if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        message(FATAL_ERROR "ASTRO_COVERAGE needs GCC (gcov)")
    endif()
    add_compile_options(--coverage -O0)
    add_link_options(--coverage)
endif()

if(ASTRO_BUILD_FUZZERS)
    # Fuzzing instruments everything: coverage feedback for libFuzzer, and the address and
    # undefined-behaviour sanitizers so that a memory error or an overflow stops the run.
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        message(FATAL_ERROR "ASTRO_BUILD_FUZZERS needs Clang (libFuzzer); set CMAKE_CXX_COMPILER=clang++")
    endif()
    add_compile_options(-fsanitize=fuzzer-no-link,address,undefined -fno-sanitize-recover=undefined
                        -fno-omit-frame-pointer)
    add_link_options(-fsanitize=address,undefined)
endif()

file(GLOB_RECURSE ASTRO_SOURCES CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB_RECURSE ASTRO_HEADERS CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/include/astro/*.h)

add_library(astro ${ASTRO_SOURCES} ${ASTRO_HEADERS})
add_library(astro::astro ALIAS astro)

target_compile_features(astro PUBLIC cxx_std_20)
# The soname carries the major and minor version (libastro.so.0.3). astro promises source compatibility
# within a major version, not binary compatibility: a minor release may change the application binary
# interface (ABI), so it gets a new soname and a program built against 0.3 never loads 0.4 by mistake.
# A patch release keeps the ABI, and keeps the soname.
set_target_properties(astro PROPERTIES
        CXX_EXTENSIONS OFF
        POSITION_INDEPENDENT_CODE ON
        CXX_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN ON
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})

# astro/version.h is generated from include/astro/version.h.in.
configure_file(include/astro/version.h.in ${PROJECT_BINARY_DIR}/include/astro/version.h @ONLY)

# astro/export.h defines ASTRO_EXPORT, which marks the public API. The library is built with hidden
# visibility, so a shared build exports only what is marked. In a static build ASTRO_STATIC_DEFINE
# makes ASTRO_EXPORT empty, for the library and for every program that links it.
#
# A class exported from a DLL makes MSVC warn about each member or base that is not exported itself
# (C4251, C4275): std::vector, std::string, std::runtime_error, and astro's own Vec3 and Mat3. Both
# warnings are harmless when the DLL and the program use the same compiler and runtime library, which
# a C++ DLL requires anyway. Each header that exports a class silences them between
# ASTRO_DLL_WARNINGS_PUSH and ASTRO_DLL_WARNINGS_POP, so the warnings stay on in the user's own code.
include(GenerateExportHeader)
set(ASTRO_EXPORT_EXTRA [=[
/* MSVC: no C4251/C4275 for astro's exported classes (see CMakeLists.txt). */
#if defined(_MSC_VER) && !defined(ASTRO_STATIC_DEFINE)
#  define ASTRO_DLL_WARNINGS_PUSH __pragma(warning(push)) __pragma(warning(disable : 4251 4275))
#  define ASTRO_DLL_WARNINGS_POP __pragma(warning(pop))
#else
#  define ASTRO_DLL_WARNINGS_PUSH
#  define ASTRO_DLL_WARNINGS_POP
#endif]=])
generate_export_header(astro
        BASE_NAME astro
        EXPORT_FILE_NAME ${PROJECT_BINARY_DIR}/include/astro/export.h
        CUSTOM_CONTENT_FROM_VARIABLE ASTRO_EXPORT_EXTRA)
get_target_property(ASTRO_LIBRARY_TYPE astro TYPE)
if(ASTRO_LIBRARY_TYPE STREQUAL "STATIC_LIBRARY")
    target_compile_definitions(astro PUBLIC ASTRO_STATIC_DEFINE)
endif()
target_include_directories(astro PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
# src/ holds the internal headers: not installed, so not part of the public API.
target_include_directories(astro PRIVATE ${PROJECT_SOURCE_DIR}/src)

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
install(DIRECTORY include/astro DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} PATTERN "*.in" EXCLUDE)
install(FILES ${PROJECT_BINARY_DIR}/include/astro/version.h ${PROJECT_BINARY_DIR}/include/astro/export.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/astro)
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
install(FILES LICENSE THIRD_PARTY_NOTICES.md README.md CHANGELOG.md DESTINATION ${CMAKE_INSTALL_DOCDIR})

# ---------------------------------------------------------------- packages
# `cpack -G DEB` or `cpack -G RPM` in the build directory (see cmake/AstroPackaging.cmake). Only when
# astro is the top-level project: a project that includes it packages itself.
set(ASTRO_PACKAGE_VENDOR "W. Greg Visotski" CACHE STRING "Vendor named in the packages")
set(ASTRO_PACKAGE_CONTACT "W. Greg Visotski" CACHE STRING "Maintainer named in the packages")
set(ASTRO_PACKAGE_HOMEPAGE "https://github.com/wgvisots/astro" CACHE STRING "Home page named in the packages")
if(PROJECT_IS_TOP_LEVEL)
    include(cmake/AstroPackaging.cmake)
endif()

# ---------------------------------------------------------------- tests
if(ASTRO_BUILD_TESTS OR ASTRO_BUILD_EXAMPLES)
    enable_testing()
endif()
if(ASTRO_BUILD_TESTS)
    add_subdirectory(tests)
endif()
if(ASTRO_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(ASTRO_BUILD_BENCHMARKS)
    add_subdirectory(bench)
endif()
if(ASTRO_BUILD_STUDIES)
    add_subdirectory(studies)
endif()
if(ASTRO_BUILD_FUZZERS)
    enable_testing()
    add_subdirectory(fuzz)
endif()

# ---------------------------------------------------------------- docs
# `cmake --build build --target docs` builds the API reference and the Markdown pages into
# build/docs/html. The code in the Markdown pages must match the examples (test docs.snippets).
find_package(Doxygen 1.9 OPTIONAL_COMPONENTS dot)
if(DOXYGEN_FOUND)
    set(ASTRO_DOCS_OUTPUT ${PROJECT_BINARY_DIR}/docs)
    if(DOXYGEN_DOT_FOUND)
        set(ASTRO_DOCS_HAVE_DOT YES)
    else()
        set(ASTRO_DOCS_HAVE_DOT NO)
    endif()
    if(ASTRO_DOCS_WARNINGS_AS_ERRORS)
        set(ASTRO_DOCS_WARN_AS_ERROR FAIL_ON_WARNINGS)
    else()
        set(ASTRO_DOCS_WARN_AS_ERROR NO)
    endif()
    configure_file(docs/Doxyfile.in ${ASTRO_DOCS_OUTPUT}/Doxyfile @ONLY)
    add_custom_target(docs
            COMMAND Doxygen::doxygen ${ASTRO_DOCS_OUTPUT}/Doxyfile
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            COMMENT "Building HTML documentation in ${ASTRO_DOCS_OUTPUT}/html"
            VERBATIM)
endif()
find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND AND (ASTRO_BUILD_TESTS OR ASTRO_BUILD_EXAMPLES))
    add_test(NAME docs.snippets COMMAND Python3::Interpreter ${PROJECT_SOURCE_DIR}/tools/check_doc_snippets.py ${PROJECT_SOURCE_DIR})
endif()
if(Python3_FOUND AND ASTRO_BUILD_TESTS)
    # astro/astro.h includes every public header.
    add_test(NAME headers.astro_h COMMAND Python3::Interpreter ${PROJECT_SOURCE_DIR}/tools/check_astro_header.py ${PROJECT_SOURCE_DIR})
    # A shared build exports no internal symbol (detail namespaces, anonymous namespaces).
    if(ASTRO_LIBRARY_TYPE STREQUAL "SHARED_LIBRARY" AND CMAKE_NM)
        add_test(NAME exports.internal
                 COMMAND Python3::Interpreter ${PROJECT_SOURCE_DIR}/tools/check_exports.py $<TARGET_FILE:astro> --nm ${CMAKE_NM})
    endif()
endif()
