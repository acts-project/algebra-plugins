/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "algebra/smatrix_cmath.hpp"

// Test include(s).
#include "test_host_basics.hpp"

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <string>

/// Struct providing a readable name for the test
struct test_specialisation_name {
  template <typename T>
  static std::string GetName(int i) {
    switch (i) {
      case 0:
        return "smatrix_cmath<float>";
      case 1:
        return "smatrix_cmath<double>";
      default:
        return "unknown";
    }
  }
};

// Register the tests
REGISTER_TYPED_TEST_SUITE_P(test_host_basics_vector, local_vectors, vector3,
                            getter);
REGISTER_TYPED_TEST_SUITE_P(test_host_basics_matrix, matrix3, matrix64,
                            matrix22);
REGISTER_TYPED_TEST_SUITE_P(test_host_basics_transform, transform3,
                            global_transformations);

// Instantiate the test(s).
typedef testing::Types<
    test_types<
        float, algebra::smatrix::point2<float>, algebra::smatrix::point3<float>,
        algebra::smatrix::vector2<float>, algebra::smatrix::vector3<float>,
        algebra::smatrix::transform3<float>, unsigned int,
        algebra::smatrix::matrix_type,
        algebra::matrix::actor<float,
                               algebra::matrix::determinant::preset0<float>,
                               algebra::matrix::inverse::preset0<float>>>,
    test_types<
        double, algebra::smatrix::point2<double>,
        algebra::smatrix::point3<double>, algebra::smatrix::vector2<double>,
        algebra::smatrix::vector3<double>, algebra::smatrix::transform3<double>,
        unsigned int, algebra::smatrix::matrix_type,
        algebra::matrix::actor<double,
                               algebra::matrix::determinant::preset0<double>,
                               algebra::matrix::inverse::preset0<double>>>>
    smatrix_cmath_types;
INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, test_host_basics_vector,
                               smatrix_cmath_types, test_specialisation_name);
INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, test_host_basics_matrix,
                               smatrix_cmath_types, test_specialisation_name);
INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, test_host_basics_transform,
                               smatrix_cmath_types, test_specialisation_name);
