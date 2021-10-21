# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up an INTERFACE library for the common header(s).
add_library( algebra_tests_common INTERFACE )
target_include_directories( algebra_tests_common INTERFACE
   "${CMAKE_CURRENT_SOURCE_DIR}/common" )
target_link_libraries( algebra_tests_common INTERFACE
   algebra::common )
add_library( algebra::tests_common ALIAS algebra_tests_common )

# Set up all of the (available) tests.
algebra_add_test( array_cmath "array/array_cmath.cpp"
   LINK_LIBRARIES GTest::gtest_main algebra::tests_common algebra::array_cmath )

if( ALGEBRA_PLUGIN_INCLUDE_EIGEN )
   algebra_add_test( eigen_cmath "eigen/eigen_cmath.cpp"
      LINK_LIBRARIES GTest::gtest_main algebra::tests_common
                     algebra::eigen_cmath )
endif()

if( ALGEBRA_PLUGIN_INCLUDE_VC )
   algebra_add_test( vc_cmath "vc/vc_cmath.cpp"
      LINK_LIBRARIES GTest::gtest_main algebra::tests_common algebra::vc_cmath )
endif()

if( ALGEBRA_PLUGIN_INCLUDE_VECMEM )
   algebra_add_test( vecmem_cmath "vecmem/vecmem_cmath.cpp"
      LINK_LIBRARIES GTest::gtest_main algebra::tests_common
                     algebra::vecmem_cmath )
endif()
