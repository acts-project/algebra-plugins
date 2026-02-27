/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/concepts.hpp"
#include "algebra/type_traits.hpp"
#include "algebra/utils/approximately_equal.hpp"
#include "algebra/utils/casts.hpp"
#include "algebra/utils/print.hpp"

// Test include(s).
#include "algebra/test/framework/test_base.hpp"
#include "algebra/test/framework/types.hpp"

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>

namespace algebra::test {

/// Test case class, to be specialised for the different plugins - vectors
template <algebra::concepts::algebra A>
class vector_test : public testing::Test, public test_base<A> {};

TYPED_TEST_SUITE_P(vector_test);

// This defines the local frame test suite
TYPED_TEST_P(vector_test, vector2D) {
  // Print the linear algebra types of this backend
  using algebra::operator<<;

  // Construction
  algebra::get_point2D_t<TypeParam> vA{0.f, 1.f};
  ASSERT_EQ(vA[0], 0.f);
  ASSERT_EQ(vA[1], 1.f);

  // Test printing
  std::cout << vA << std::endl;

  // Test comparison
  constexpr auto epsilon{
      std::numeric_limits<algebra::get_scalar_t<TypeParam>>::epsilon()};
  ASSERT_TRUE(algebra::approx_equal(vA, vA));
  ASSERT_TRUE(algebra::approx_equal(vA, vA, epsilon));

  algebra::get_scalar_t<TypeParam> rel_err{1.f + 10.f * epsilon};
  algebra::get_point2D_t<TypeParam> vA_err = rel_err * vA;
  ASSERT_TRUE(algebra::approx_equal(vA, vA_err, 11.f * epsilon));
  ASSERT_FALSE(algebra::approx_equal(vA, vA_err, 9.f * epsilon));

  rel_err = 1.f + 17.f * epsilon;
  vA_err = rel_err * vA;
  ASSERT_TRUE(algebra::approx_equal(vA, vA_err, 18.f * epsilon));
  ASSERT_FALSE(algebra::approx_equal(vA, vA_err, 16.f * epsilon));
  // Cast to (different) precision
  const auto vA_cast_f = algebra::cast_to<float>(vA);

  for (algebra::get_size_t<TypeParam> i = 0; i < 2; ++i) {
    auto elem_i = vA_cast_f[i];

    static_assert(std::same_as<decltype(elem_i), float>);
    ASSERT_FLOAT_EQ(elem_i, static_cast<float>(vA[i]));
  }

  const auto vA_cast_d = algebra::cast_to<double>(vA);

  for (algebra::get_size_t<TypeParam> i = 0; i < 2; ++i) {
    auto elem_i = vA_cast_d[i];

    static_assert(std::same_as<decltype(elem_i), double>);
    ASSERT_DOUBLE_EQ(elem_i, static_cast<double>(vA[i]));
  }

  const auto vA_cast_i = algebra::cast_to<int>(vA);

  for (algebra::get_size_t<TypeParam> i = 0; i < 2; ++i) {
    auto elem_i = vA_cast_i[i];

    static_assert(std::same_as<decltype(elem_i), int>);
    ASSERT_EQ(elem_i, static_cast<int>(vA[i]));
  }

  // Assignment
  algebra::get_point2D_t<TypeParam> vB = vA;
  ASSERT_EQ(vB[0], 0.f);
  ASSERT_EQ(vB[1], 1.f);

  // Addition
  algebra::get_point2D_t<TypeParam> vC = vA + vB;
  ASSERT_EQ(vC[0], 0.f);
  ASSERT_EQ(vC[1], 2.f);

  // Multiplication by scalar
  algebra::get_point2D_t<TypeParam> vC2 = vC * 2.f;
  ASSERT_EQ(vC2[0], 0.f);
  ASSERT_EQ(vC2[1], 4.f);

  // Cast operations to phi, theta, eta, perp
  algebra::get_vector2D_t<TypeParam> vD{1.f, 1.f};
  algebra::get_scalar_t<TypeParam> phi = algebra::vector::phi(vD);
  ASSERT_NEAR(phi, M_PI_4, this->m_epsilon);

  algebra::get_scalar_t<TypeParam> perp = algebra::vector::perp(vD);
  ASSERT_NEAR(perp, std::sqrt(2.), this->m_epsilon);

  algebra::get_scalar_t<TypeParam> norm = algebra::vector::norm(vD);
  ASSERT_NEAR(norm, std::sqrt(2.), this->m_epsilon);

  algebra::get_vector2D_t<TypeParam> vDnorm = algebra::vector::normalize(vD);
  ASSERT_NEAR(vDnorm[0], 1. / std::sqrt(2.), this->m_epsilon);
  ASSERT_NEAR(vDnorm[1], 1. / std::sqrt(2.), this->m_epsilon);

  // Test comparison
  algebra::get_point2D_t<TypeParam> vE{100.f, 462809.f};

  ASSERT_TRUE(algebra::approx_equal(vE, vE));
  ASSERT_TRUE(algebra::approx_equal(vE, vE, epsilon));

  rel_err = 1.f + 10.f * epsilon;
  algebra::get_point2D_t<TypeParam> vE_err = rel_err * vE;
  ASSERT_TRUE(algebra::approx_equal(vE, vE_err, 11.f * epsilon));
  ASSERT_FALSE(algebra::approx_equal(vE, vE_err, 9.f * epsilon));

  algebra::get_point2D_t<TypeParam> vE_abs_err{100.00001f, 462809.05f};
  ASSERT_TRUE(algebra::approx_equal(vE, vE_abs_err, 0.00001f));
  ASSERT_FALSE(algebra::approx_equal(vE, vE_abs_err, epsilon));

  rel_err = 1.f + 17.f * epsilon;
  vE_err = rel_err * vE;
  ASSERT_TRUE(algebra::approx_equal(vE, vE_err, 18.f * epsilon));
  ASSERT_FALSE(algebra::approx_equal(vE, vE_err, 16.f * epsilon));
}

// This defines the vector3 test suite
TYPED_TEST_P(vector_test, vector3D) {
  // Print the linear algebra types of this backend
  using algebra::operator<<;

  // Test concepts
  static_assert(algebra::concepts::scalar<algebra::get_scalar_t<TypeParam>>);
  static_assert(!algebra::concepts::vector<algebra::get_scalar_t<TypeParam>>);

  static_assert(!algebra::concepts::scalar<algebra::get_vector3D_t<TypeParam>>);
  static_assert(algebra::concepts::vector<algebra::get_vector3D_t<TypeParam>>);
  static_assert(
      algebra::concepts::vector3D<algebra::get_vector3D_t<TypeParam>>);
  static_assert(
      !algebra::concepts::vector2D<algebra::get_vector3D_t<TypeParam>>);

  static_assert(!algebra::concepts::scalar<algebra::get_vector2D_t<TypeParam>>);
  static_assert(algebra::concepts::vector<algebra::get_vector2D_t<TypeParam>>);
  static_assert(
      algebra::concepts::vector2D<algebra::get_vector2D_t<TypeParam>>);
  static_assert(
      !algebra::concepts::vector3D<algebra::get_vector2D_t<TypeParam>>);

  // Construction
  algebra::get_vector3D_t<TypeParam> vA{0.f, 1.f, 2.f};
  ASSERT_EQ(vA[0], 0.f);
  ASSERT_EQ(vA[1], 1.f);
  ASSERT_EQ(vA[2], 2.f);

  // Test printing
  std::cout << vA << std::endl;

  // Cast to (different) precision
  const auto vA_cast_f = algebra::cast_to<float>(vA);

  for (algebra::get_size_t<TypeParam> i = 0; i < 3; ++i) {
    auto elem_i = vA_cast_f[i];

    static_assert(std::same_as<decltype(elem_i), float>);
    ASSERT_FLOAT_EQ(elem_i, static_cast<float>(vA[i]));
  }

  const auto vA_cast_d = algebra::cast_to<double>(vA);

  for (algebra::get_size_t<TypeParam> i = 0; i < 3; ++i) {
    auto elem_i = vA_cast_d[i];

    static_assert(std::same_as<decltype(elem_i), double>);
    ASSERT_DOUBLE_EQ(elem_i, static_cast<double>(vA[i]));
  }

  const auto vA_cast_i = algebra::cast_to<int>(vA);

  for (algebra::get_size_t<TypeParam> i = 0; i < 3; ++i) {
    auto elem_i = vA_cast_i[i];

    static_assert(std::same_as<decltype(elem_i), int>);
    ASSERT_EQ(elem_i, static_cast<int>(vA[i]));
  }

  // Assignment
  algebra::get_vector3D_t<TypeParam> vB = vA;
  ASSERT_EQ(vB[0], 0.f);
  ASSERT_EQ(vB[1], 1.f);
  ASSERT_EQ(vB[2], 2.f);

  // Addition
  algebra::get_vector3D_t<TypeParam> vC = vA + vB;
  ASSERT_EQ(vC[0], 0.f);
  ASSERT_EQ(vC[1], 2.f);
  ASSERT_EQ(vC[2], 4.f);

  // Multiplication by scalar
  algebra::get_vector3D_t<TypeParam> vC2 = vC * 2.0f;
  ASSERT_EQ(vC2[0], 0.f);
  ASSERT_EQ(vC2[1], 4.f);
  ASSERT_EQ(vC2[2], 8.f);

  // Cast operations to phi, theta, eta, perp
  algebra::get_vector3D_t<TypeParam> vD{1.f, 1.f, 1.f};
  algebra::get_scalar_t<TypeParam> phi = algebra::vector::phi(vD);
  ASSERT_NEAR(phi, M_PI_4, this->m_epsilon);

  algebra::get_scalar_t<TypeParam> theta = algebra::vector::theta(vD);
  ASSERT_NEAR(theta, std::atan2(std::sqrt(2.), 1.), this->m_epsilon);

  algebra::get_scalar_t<TypeParam> eta = algebra::vector::eta(vD);
  ASSERT_NEAR(eta, 0.65847891569137573, this->m_isclose);

  algebra::get_scalar_t<TypeParam> perp = algebra::vector::perp(vD);
  ASSERT_NEAR(perp, std::sqrt(2.), this->m_epsilon);

  algebra::get_scalar_t<TypeParam> norm = algebra::vector::norm(vD);
  ASSERT_NEAR(norm, std::sqrt(3.), this->m_epsilon);
}

// This defines the vector operation test suite
TYPED_TEST_P(vector_test, element_getter) {

  algebra::get_vector3D_t<TypeParam> v3{1.f, 1.f, 1.f};

  // Normalization
  algebra::get_vector3D_t<TypeParam> v3n = algebra::vector::normalize(v3);
  ASSERT_NEAR(v3n[0], 1. / std::sqrt(3.), this->m_epsilon);
  ASSERT_NEAR(v3n[1], 1. / std::sqrt(3.), this->m_epsilon);
  ASSERT_NEAR(v3n[2], 1. / std::sqrt(3.), this->m_epsilon);

  // Cross product
  algebra::get_vector3D_t<TypeParam> z = algebra::vector::normalize(
      algebra::get_vector3D_t<TypeParam>{3.f, 2.f, 1.f});
  algebra::get_vector3D_t<TypeParam> x = algebra::vector::normalize(
      algebra::get_vector3D_t<TypeParam>{2.f, -3.f, 0.f});
  algebra::get_vector3D_t<TypeParam> y =
      algebra::vector::cross(z, algebra::vector::normalize(x));

  // Check with dot product
  ASSERT_NEAR(algebra::vector::dot(x, y), 0.f, this->m_epsilon);
  ASSERT_NEAR(algebra::vector::dot(y, z), 0.f, this->m_epsilon);
  ASSERT_NEAR(algebra::vector::dot(z, x), 0.f, this->m_epsilon);
}

// This test checks to see if the `dot` function can handle when one of its
// operands is a sum or difference of two vectors. We also test to see if the
// `dot` function plays nicely with scalar multiplication (just to be safe).
TYPED_TEST_P(vector_test, dot_product_with_ops) {
  get_vector3D_t<TypeParam> v1{1.f, 2.f, 3.f};
  get_vector3D_t<TypeParam> v2{3.f, 4.f, 5.f};

  ASSERT_NEAR(algebra::vector::dot(v1 + v2, v2), 76.f, this->m_epsilon);
  ASSERT_NEAR(algebra::vector::dot(v1, v2 - v1), 12.f, this->m_epsilon);
  ASSERT_NEAR(algebra::vector::dot(v1 + v2, v1 - v2), -36.f, this->m_epsilon);
  ASSERT_NEAR(algebra::vector::dot(v1 + v2, 2 * v2), 152.f, this->m_epsilon);
}

// This test checks to see if the `cross` function can handle when one of its
// operands is a sum or difference of two vectors. We also test to see if the
// `cross` function plays nicely with scalar multiplication (just to be safe).
TYPED_TEST_P(vector_test, cross_product_add_sub) {
  get_vector3D_t<TypeParam> v1{1.f, 2.f, 3.f};
  get_vector3D_t<TypeParam> v2{3.f, 4.f, 5.f};
  get_vector3D_t<TypeParam> v3{-6.f, 7.f, -9.f};

  get_vector3D_t<TypeParam> v = algebra::vector::cross(v1 + v2, v3);
  get_vector3D_t<TypeParam> ans{-110.f, -12.f, 64.f};

  ASSERT_NEAR(v[0], ans[0], this->m_epsilon);
  ASSERT_NEAR(v[1], ans[1], this->m_epsilon);
  ASSERT_NEAR(v[2], ans[2], this->m_epsilon);

  v = algebra::vector::cross(v3 - 2 * v1, 3 * (v1 + v2));
  ans = {342.f, 12.f, -180.f};

  ASSERT_NEAR(v[0], ans[0], this->m_epsilon);
  ASSERT_NEAR(v[1], ans[1], this->m_epsilon);
  ASSERT_NEAR(v[2], ans[2], this->m_epsilon);
}

// clang-format off
#define ALGEBRA_REGISTER_VECTOR_TESTS(...) \
  REGISTER_TYPED_TEST_SUITE_P(vector_test \
    , vector2D \
    , vector3D \
    , element_getter \
    , dot_product_with_ops \
    , cross_product_add_sub \
)
// clang-format on

}  // namespace algebra::test
