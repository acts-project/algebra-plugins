# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up a "common library" for the device tests.
add_library( algebra_tests_accel_common INTERFACE )
target_include_directories( algebra_tests_accel_common INTERFACE
   "${CMAKE_CURRENT_SOURCE_DIR}/common" )
target_link_libraries( algebra_tests_accel_common
   INTERFACE GTest::gtest algebra::common algebra::tests_common )
add_library( algebra::tests_accel_common ALIAS algebra_tests_accel_common )

# Run the CUDA test(s) by default, if CUDA is available.
option( ALGEBRA_PLUGINS_TEST_CUDA "Test the code in CUDA device code" FALSE )
if( ALGEBRA_PLUGINS_TEST_CUDA )
   add_subdirectory( cuda )
endif()

option( ALGEBRA_PLUGINS_TEST_SYCL "Test the code in SYCL device code" FALSE )
if( ALGEBRA_PLUGINS_TEST_SYCL )
   add_subdirectory( sycl )
endif()
