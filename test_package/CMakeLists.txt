# Conan's test of the astro package (conan create runs it): the consumer program CI also builds
# against the Debian package, here found through the Conan-generated CMake files.
cmake_minimum_required(VERSION 3.20)
project(test_package LANGUAGES CXX)

find_package(astro REQUIRED CONFIG)

add_executable(test_package ${CMAKE_CURRENT_SOURCE_DIR}/../tests/package/main.cpp)
target_link_libraries(test_package PRIVATE astro::astro)
