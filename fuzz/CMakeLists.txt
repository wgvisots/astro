# Fuzz targets for the parsers that SECURITY.md names as the library's attack surface, and for the
# two text parsers beside them. They need Clang (libFuzzer). The top-level CMakeLists.txt has
# already instrumented the whole build for coverage feedback and for the address and
# undefined-behaviour sanitizers.
#
#   cmake -B build-fuzz -DCMAKE_CXX_COMPILER=clang++ -DASTRO_BUILD_FUZZERS=ON -DASTRO_BUILD_TESTS=OFF \
#         -DASTRO_BUILD_EXAMPLES=OFF
#   cmake --build build-fuzz -j
#   build-fuzz/fuzz/astro_fuzz_tle fuzz/corpus/tle -max_total_time=60
#
# Each target stops at the first crash, sanitizer report or failed property and saves the input.
# `ctest -L fuzz` runs every target for a few seconds over its seed corpus, as a smoke test.

set(ASTRO_FUZZ_TARGETS angle eop iso8601 leap_seconds spk tle)

foreach(name IN LISTS ASTRO_FUZZ_TARGETS)
    add_executable(astro_fuzz_${name} ${name}.cpp)
    target_link_libraries(astro_fuzz_${name} PRIVATE astro::astro)
    target_compile_options(astro_fuzz_${name} PRIVATE ${ASTRO_WARNING_FLAGS})
    target_link_options(astro_fuzz_${name} PRIVATE -fsanitize=fuzzer)
    # libFuzzer writes what it finds into the first directory: keep that out of the source tree.
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/corpus/${name})
    add_test(NAME fuzz.${name}
             COMMAND astro_fuzz_${name} ${CMAKE_CURRENT_BINARY_DIR}/corpus/${name}
                     ${CMAKE_CURRENT_SOURCE_DIR}/corpus/${name} -max_total_time=5 -max_len=65536 -timeout=10
                     -rss_limit_mb=2048 -artifact_prefix=${CMAKE_CURRENT_BINARY_DIR}/${name}-)
    set_tests_properties(fuzz.${name} PROPERTIES LABELS fuzz TIMEOUT 120)
endforeach()
