# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message(STATUS "Building 'algebra::test_sycl' component")

# Enable SYCL as a language.
enable_language(SYCL)

# Set the default SYCL compiler flags.
include(algebra-plugins-compiler-options-sycl)

if(${CMAKE_SYCL_STANDARD} LESS 20)
    message(
        SEND_ERROR
        "CMAKE_SYCL_STANDARD=${CMAKE_SYCL_STANDARD}, but algebra-plugins requires C++>=20"
    )
endif()

# Helper library for the SYCL tests.
add_library(
    algebra_test_sycl
    INTERFACE
    "algebra_test_suite.hpp"
    "execute_sycl_test.hpp"
    "sycl_test_fixture.hpp"
)

target_include_directories(
    algebra_test_sycl
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(
    algebra_test_sycl
    INTERFACE
        GTest::gtest
        algebra::test_framework
        algebra::test_device
        vecmem::core
        vecmem::sycl
)

add_library(algebra::test_sycl ALIAS algebra_test_sycl)
