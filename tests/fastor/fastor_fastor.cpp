/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "algebra/fastor_fastor.hpp"

// Test include(s).
#include "test_host_basics.hpp"

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <cstddef>  // for the std::size_t data type
#include <string>

/// Struct providing a readable name for the test
struct test_specialisation_name {
  template <typename T>
  static std::string GetName(int i) {
    switch (i) {
      case 0:
        return "fastor_fastor<float>";
      case 1:
        return "fastor_fastor<double>";
      default:
        return "unknown";
    }
  }
};


// This test checks to see if the `dot` function can handle when one of its operands is a sum or difference of two vectors.
// We also test to see if the `dot` function plays nicely with scalar multiplication (just to be safe).
TYPED_TEST_P(test_host_basics_vector, dot_product_with_ops) {
  typename TypeParam::vector3 v1{1.f, 2.f, 3.f};
  typename TypeParam::vector3 v2{3.f, 4.f, 5.f};

  ASSERT_NEAR(algebra::vector::dot(v1 + v2, v2), 76.f, this->m_epsilon);
  ASSERT_NEAR(algebra::vector::dot(v1, v2 - v1), 12.f, this->m_epsilon);
  ASSERT_NEAR(algebra::vector::dot(v1 + v2, v1 - v2), -36.f, this->m_epsilon);
  ASSERT_NEAR(algebra::vector::dot(v1 + v2, 2 * v2), 152.f, this->m_epsilon);  
}

// This test checks to see if the `cross` function can handle when one of its operands is a sum or difference of two vectors.
// We also test to see if the `cross` function plays nicely with scalar multiplication (just to be safe).
TYPED_TEST_P(test_host_basics_vector, cross_product_add_sub) {
  typename TypeParam::vector3 v1{1.f, 2.f, 3.f};
  typename TypeParam::vector3 v2{3.f, 4.f, 5.f};
  typename TypeParam::vector3 v3{-6.f, 7.f, -9.f};

  typename TypeParam::vector3 v = algebra::vector::cross(v1 + v2, v3);
  typename TypeParam::vector3 ans{-110.f, -12.f, 64.f};

  ASSERT_NEAR(v[0], ans[0], this->m_epsilon);
  ASSERT_NEAR(v[1], ans[1], this->m_epsilon);
  ASSERT_NEAR(v[2], ans[2], this->m_epsilon);

  v = algebra::vector::cross(v3 - v1, v1 + v2);
  ans = {112.f, 8.f, -62.f};

  ASSERT_NEAR(v[0], ans[0], this->m_epsilon);
  ASSERT_NEAR(v[1], ans[1], this->m_epsilon);
  ASSERT_NEAR(v[2], ans[2], this->m_epsilon);
}

// Register the tests
REGISTER_TYPED_TEST_SUITE_P(test_host_basics_vector, local_vectors, vector3,
                            getter, dot_product_with_ops, cross_product_add_sub);
TEST_HOST_BASICS_MATRIX_TESTS();
REGISTER_TYPED_TEST_SUITE_P(test_host_basics_transform, transform3,
                            global_transformations);

// Instantiate the test(s).
typedef testing::Types<
    test_types<float, algebra::fastor::point2<float>,
               algebra::fastor::point3<float>, algebra::fastor::vector2<float>,
               algebra::fastor::vector3<float>,
               algebra::fastor::transform3<float>, std::size_t,
               algebra::fastor::matrix_type>,
    test_types<
        double, algebra::fastor::point2<double>,
        algebra::fastor::point3<double>, algebra::fastor::vector2<double>,
        algebra::fastor::vector3<double>, algebra::fastor::transform3<double>,
        std::size_t, algebra::fastor::matrix_type>>
    fastor_fastor_types;
INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, test_host_basics_vector,
                               fastor_fastor_types, test_specialisation_name);
INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, test_host_basics_matrix,
                               fastor_fastor_types, test_specialisation_name);
INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, test_host_basics_transform,
                               fastor_fastor_types, test_specialisation_name);
