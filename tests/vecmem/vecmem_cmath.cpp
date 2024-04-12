/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "algebra/vecmem_cmath.hpp"

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
        return "vecmem_cmath<float>";
      case 1:
        return "vecmem_cmath<double>";
      default:
        return "unknown";
    }
  }
};

// Register the tests
REGISTER_TYPED_TEST_SUITE_P(test_host_basics_vector, local_vectors, vector3,
                            getter);
TEST_HOST_BASICS_MATRIX_TESTS();
REGISTER_TYPED_TEST_SUITE_P(test_host_basics_transform, transform3,
                            global_transformations);

// Instantiate the test(s).
typedef testing::Types<
    test_types<float, algebra::vecmem::point2<float>,
               algebra::vecmem::point3<float>, algebra::vecmem::vector2<float>,
               algebra::vecmem::vector3<float>,
               algebra::vecmem::transform3<float>, std::size_t,
               algebra::vecmem::matrix_type>,
    test_types<
        double, algebra::vecmem::point2<double>,
        algebra::vecmem::point3<double>, algebra::vecmem::vector2<double>,
        algebra::vecmem::vector3<double>, algebra::vecmem::transform3<double>,
        std::size_t, algebra::vecmem::matrix_type>>
    vecmem_cmath_types;
INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, test_host_basics_vector,
                               vecmem_cmath_types, test_specialisation_name);
INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, test_host_basics_matrix,
                               vecmem_cmath_types, test_specialisation_name);
INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, test_host_basics_transform,
                               vecmem_cmath_types, test_specialisation_name);
