# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message(STATUS "Building the algebra plugins test suite")

# Set the default C++ compiler flags.
include(algebra-plugins-compiler-options-cpp)

# Algebra-plugins unit test framework (depends on google test)
add_library(
    algebra_test_framework
    INTERFACE
    "include/algebra/test/framework/test_base.hpp"
    "include/algebra/test/framework/types.hpp"
)
add_library(algebra::test_framework ALIAS algebra_test_framework)

target_include_directories(
    algebra_test_framework
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_link_libraries(
    algebra_test_framework
    INTERFACE GTest::gtest algebra::utils
)

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")
    target_compile_definitions(
        algebra_test_framework
        INTERFACE -D_USE_MATH_DEFINES
    )
endif()

# Build the local test libraries
add_subdirectory(include/algebra/test)

# Build the unittests
add_subdirectory(unit_tests)
