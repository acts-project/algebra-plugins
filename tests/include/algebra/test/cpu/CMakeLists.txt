# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message(STATUS "Building 'algebra::test_cpu' component")

# Set the CPU build flags.
include(algebra-plugins-compiler-options-cpp)

# Set up a test library, which the "new style" benchmarks and tests could use.
add_library(
    algebra_test_cpu
    INTERFACE
    "matrix_test.hpp"
    "vector_test.hpp"
    "transform_test.hpp"
)

add_library(algebra::test_cpu ALIAS algebra_test_cpu)

target_link_libraries(
    algebra_test_cpu
    INTERFACE algebra::common algebra::test_framework algebra::utils
)
