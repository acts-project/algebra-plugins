# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message(STATUS "Building 'algebra::test_device' component")

# Set up a "common library" for the device tests.
add_library(
    algebra_test_device
    INTERFACE
    "execute_host_test.hpp"
    "device_fixture.hpp"
    "matrix_fixture.hpp"
    "transform_fixture.hpp"
    "vector_fixture.hpp"
)
add_library(algebra::test_device ALIAS algebra_test_device)

target_include_directories(
    algebra_test_device
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(
    algebra_test_device
    INTERFACE GTest::gtest vecmem::core algebra::common algebra::test_framework
)

if(ALGEBRA_BUILD_CUDA)
    add_subdirectory(cuda)
endif()

if(ALGEBRA_BUILD_SYCL)
    add_subdirectory(sycl)
endif()
