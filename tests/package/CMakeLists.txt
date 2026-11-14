# A program outside the astro build that uses an installed astro through find_package. CI builds it
# against the installed Debian package; to try it on any installation:
#
#   cmake -S tests/package -B build-consumer -DCMAKE_PREFIX_PATH=<install prefix>
#   cmake --build build-consumer && build-consumer/astro_consumer
cmake_minimum_required(VERSION 3.20)
project(astro_consumer LANGUAGES CXX)

find_package(astro 0.1 REQUIRED)

add_executable(astro_consumer main.cpp)
target_link_libraries(astro_consumer PRIVATE astro::astro)
