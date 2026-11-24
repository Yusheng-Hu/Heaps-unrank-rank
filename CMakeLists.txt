# CMake build for heaps_unrank (rank / unrank for the Heap permutation order)
#
# Two implementations are built into one library:
#   heaps_unrank.cpp       main                   Theta(n^2)
#   heaps_unrank_o_n3.cpp  baseline (dense)       Theta(n^3)
#
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
#   cmake --build build --config Release
#   ctest --test-dir build -C Release --output-on-failure

cmake_minimum_required(VERSION 3.14)

project(heaps_unrank LANGUAGES CXX VERSION 2.2.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if(MSVC)
  set(HEAPS_WARNINGS /W4)
else()
  set(HEAPS_WARNINGS -Wall -Wextra -Wpedantic)
endif()

add_library(heaps_unrank STATIC heaps_unrank.cpp heaps_unrank_o_n3.cpp heaps_generator.cpp
                              heaps_unrank_variant.cpp)
target_include_directories(heaps_unrank PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(heaps_unrank PRIVATE ${HEAPS_WARNINGS})

enable_testing()

add_executable(test_heaps_unrank test_heaps_unrank.cpp)
target_link_libraries(test_heaps_unrank PRIVATE heaps_unrank)
target_compile_options(test_heaps_unrank PRIVATE ${HEAPS_WARNINGS})

# The full suite includes the exhaustive n = 10 check (3,628,800 codes) and
# exercises both implementations.
add_test(NAME heaps_unrank_full COMMAND test_heaps_unrank --full)

# Resumable generation (paper, Section 5): restore the classical counter array,
# or drive the loop by the factorial counter.  Both are compared with
# heaps::unrank() for n <= 8.
add_executable(test_generators test_generators.cpp)
target_link_libraries(test_generators PRIVATE heaps_unrank)
target_compile_options(test_generators PRIVATE ${HEAPS_WARNINGS})
add_test(NAME generators COMMAND test_generators)

# Swap-after-call variant of the recursion (heaps_unrank_variant.hpp): a
# different order from Heap's, still Theta(n^2) in both directions.  Exhaustive
# to n = 9 by default; `test_variant --full` goes to n = 11 (43,954,710
# permutations, ~47 s).
add_executable(test_variant test_variant.cpp)
target_link_libraries(test_variant PRIVATE heaps_unrank)
target_compile_options(test_variant PRIVATE ${HEAPS_WARNINGS})
add_test(NAME variant COMMAND test_variant)

# Wall-clock benchmark (not a test); the paper's operation counts come from an
# instrumented build.
add_executable(bench_heaps bench_heaps.cpp)
target_link_libraries(bench_heaps PRIVATE heaps_unrank)
target_compile_options(bench_heaps PRIVATE ${HEAPS_WARNINGS})
