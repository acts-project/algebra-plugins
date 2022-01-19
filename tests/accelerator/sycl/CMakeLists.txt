# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Enable SYCL as a language.
enable_language( SYCL )

# Helper library for the SYCL tests.
add_library( algebra_tests_sycl_common INTERFACE )
target_include_directories( algebra_tests_sycl_common INTERFACE
   "${CMAKE_CURRENT_SOURCE_DIR}/common" )
target_link_libraries( algebra_tests_sycl_common
   INTERFACE GTest::gtest algebra::tests_common algebra::tests_accel_common
             vecmem::core vecmem::sycl )
add_library( algebra::tests_sycl_common ALIAS algebra_tests_sycl_common )

# Set up all of the (available) CUDA tests.
algebra_add_test( array_sycl "array_cmath.sycl"
   LINK_LIBRARIES algebra::array_cmath algebra::tests_sycl_common
                  GTest::gtest_main )

if( ALGEBRA_PLUGINS_INCLUDE_VECMEM )
   algebra_add_test( vecmem_sycl "vecmem_cmath.sycl"
      LINK_LIBRARIES algebra::vecmem_cmath algebra::tests_sycl_common
                     GTest::gtest_main )
endif()
