# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Helper library for the SYCL tests.
add_library( algebra_tests_sycl_common INTERFACE )
target_include_directories( algebra_tests_sycl_common INTERFACE
   "${CMAKE_CURRENT_SOURCE_DIR}/common" )
target_link_libraries( algebra_tests_sycl_common
   INTERFACE GTest::gtest algebra::tests_common algebra::tests_accel_common )
add_library( algebra::tests_sycl_common ALIAS algebra_tests_sycl_common )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsycl")

# Set up all of the (available) SYCL tests.
algebra_add_test( array_sycl "array_cmath.cpp"
   LINK_LIBRARIES algebra::array_cmath algebra::tests_sycl_common
                  GTest::gtest_main )

target_link_libraries(algebra_test_array_sycl PRIVATE sycl)

if( ALGEBRA_PLUGINS_INCLUDE_EIGEN )
   algebra_add_test( eigen_sycl "eigen_eigen.cpp"
      LINK_LIBRARIES algebra::eigen_eigen
                     algebra::tests_sycl_common GTest::gtest_main )

   target_link_libraries(algebra_test_eigen_sycl PRIVATE sycl)
   # Prevent Eigen from considering the SYCL device compilation as either:
   #  - host compilation;
   #  - CUDA compilation;
   #  - HIP compilation.
   target_compile_definitions( algebra_test_eigen_sycl
      PRIVATE EIGEN_NO_CUDA EIGEN_NO_HIP EIGEN_DONT_VECTORIZE )
endif()

