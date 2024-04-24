/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "algebra/vc_cmath.hpp"

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
        return "vc_cmath<float>";
      case 1:
        return "vc_cmath<double>";
      default:
        return "unknown";
    }
  }
};

// Register the tests
REGISTER_TYPED_TEST_SUITE_P(test_host_basics_vector, local_vectors, vector3,
                            getter);
/*REGISTER_TYPED_TEST_SUITE_P(test_host_basics_matrix, matrix3, matrix64,
matrix22); REGISTER_TYPED_TEST_SUITE_P(test_host_basics_transform, transform3,
                            global_transformations);*/

// Instantiate the test(s).
typedef testing::Types<
    test_types<float, algebra::vc_aos::point2<float>,
               algebra::vc_aos::point3<float>, algebra::vc_aos::vector2<float>,
               algebra::vc_aos::vector3<float>,
               algebra::vc_aos::transform3<float>, std::size_t,
               algebra::vc_aos::matrix_type, float>,
    test_types<
        double, algebra::vc_aos::point2<double>,
        algebra::vc_aos::point3<double>, algebra::vc_aos::vector2<double>,
        algebra::vc_aos::vector3<double>, algebra::vc_aos::transform3<double>,
        std::size_t, algebra::vc_aos::matrix_type, double>>
    vc_cmath_types;
INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, test_host_basics_vector,
                               vc_cmath_types, test_specialisation_name);
/*INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, test_host_basics_matrix,
                               vc_cmath_types, test_specialisation_name);
INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, test_host_basics_transform,
                               vc_cmath_types, test_specialisation_name);*/
