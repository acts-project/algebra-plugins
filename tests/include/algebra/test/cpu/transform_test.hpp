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

/// Test case class, to be specialised for the different plugins - transforms
template <algebra::concepts::algebra A>
class transform_test : public testing::Test, public test_base<A> {};

TYPED_TEST_SUITE_P(transform_test);

// This defines the transform3 test suite
TYPED_TEST_P(transform_test, transform3D) {
  // Print the linear algebra types of this backend
  using algebra::operator<<;

  static_assert(
      algebra::concepts::transform3D<algebra::get_transform3D_t<TypeParam>>);

  // Preparation work
  algebra::get_vector3D_t<TypeParam> z = algebra::vector::normalize(
      algebra::get_vector3D_t<TypeParam>{3.f, 2.f, 1.f});
  algebra::get_vector3D_t<TypeParam> x = algebra::vector::normalize(
      algebra::get_vector3D_t<TypeParam>{2.f, -3.f, 0.f});
  algebra::get_vector3D_t<TypeParam> y = algebra::vector::cross(z, x);
  algebra::get_point3D_t<TypeParam> t = {2.f, 3.f, 4.f};

  // Test constructor from t, z, x
  algebra::get_transform3D_t<TypeParam> trf1(t, z, x);
  ASSERT_TRUE(trf1 == trf1);
  algebra::get_transform3D_t<TypeParam> trf2;
  trf2 = trf1;

  // Test printing
  std::cout << trf1 << std::endl;

  // Test comparison
  constexpr auto epsilon{
      std::numeric_limits<algebra::get_scalar_t<TypeParam>>::epsilon()};

  ASSERT_TRUE(algebra::approx_equal(trf1, trf1));
  ASSERT_TRUE(algebra::approx_equal(trf1, trf1, epsilon));

  algebra::get_scalar_t<TypeParam> rel_err{1.f + 10.f * epsilon};
  algebra::get_transform3D_t<TypeParam> trf1_err(rel_err * t, rel_err * z,
                                                 rel_err * x);
  ASSERT_TRUE(algebra::approx_equal(trf1, trf1_err, 200.f * epsilon));
  ASSERT_FALSE(algebra::approx_equal(trf1, trf1_err, 10.f * epsilon));
  // Cast to (different) precision
  const auto trf1_cast_f = algebra::cast_to<float>(trf1);
  const auto& mat_f = trf1_cast_f.matrix();
  const auto& mat_inv_f = trf1_cast_f.matrix_inverse();

  for (algebra::get_size_t<TypeParam> j = 0; j < 3; ++j) {
    for (algebra::get_size_t<TypeParam> i = 0; i < 2; ++i) {

      auto elem_i = algebra::getter::element(mat_f, i, j);
      auto elem_inv_i = algebra::getter::element(mat_inv_f, i, j);

      static_assert(std::same_as<decltype(elem_i), float>);
      static_assert(std::same_as<decltype(elem_inv_i), float>);
      ASSERT_FLOAT_EQ(
          elem_i, static_cast<float>(algebra::getter::element(mat_f, i, j)));
      ASSERT_FLOAT_EQ(elem_inv_i, static_cast<float>(algebra::getter::element(
                                      mat_inv_f, i, j)));
    }
  }

  const auto trf1_cast_d = algebra::cast_to<double>(trf1);
  const auto& mat_d = trf1_cast_d.matrix();
  const auto& mat_inv_d = trf1_cast_d.matrix_inverse();

  for (algebra::get_size_t<TypeParam> j = 0; j < 3; ++j) {
    for (algebra::get_size_t<TypeParam> i = 0; i < 2; ++i) {

      auto elem_i = algebra::getter::element(mat_d, i, j);
      auto elem_inv_i = algebra::getter::element(mat_inv_d, i, j);

      static_assert(std::same_as<decltype(elem_i), double>);
      static_assert(std::same_as<decltype(elem_inv_i), double>);
      ASSERT_DOUBLE_EQ(
          elem_i, static_cast<double>(algebra::getter::element(mat_d, i, j)));
      ASSERT_DOUBLE_EQ(elem_inv_i, static_cast<double>(algebra::getter::element(
                                       mat_inv_d, i, j)));
    }
  }

  const auto trf1_cast_i = algebra::cast_to<int>(trf1);
  const auto& mat_i = trf1_cast_i.matrix();
  const auto& mat_inv_i = trf1_cast_i.matrix_inverse();

  for (algebra::get_size_t<TypeParam> j = 0; j < 3; ++j) {
    for (algebra::get_size_t<TypeParam> i = 0; i < 2; ++i) {

      auto elem_i = algebra::getter::element(mat_i, i, j);
      auto elem_inv_i = algebra::getter::element(mat_inv_i, i, j);

      static_assert(std::same_as<decltype(elem_i), int>);
      static_assert(std::same_as<decltype(elem_inv_i), int>);
      ASSERT_EQ(elem_i,
                static_cast<int>(algebra::getter::element(mat_i, i, j)));
      ASSERT_EQ(elem_inv_i,
                static_cast<int>(algebra::getter::element(mat_inv_i, i, j)));
    }
  }

  const auto rot = trf2.rotation();
  ASSERT_NEAR(algebra::getter::element(rot, 0, 0), x[0], this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(rot, 1, 0), x[1], this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(rot, 2, 0), x[2], this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(rot, 0, 1), y[0], this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(rot, 1, 1), y[1], this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(rot, 2, 1), y[2], this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(rot, 0, 2), z[0], this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(rot, 1, 2), z[1], this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(rot, 2, 2), z[2], this->m_epsilon);

  auto trn = trf2.translation();
  ASSERT_NEAR(trn[0], 2.f, this->m_epsilon);
  ASSERT_NEAR(trn[1], 3.f, this->m_epsilon);
  ASSERT_NEAR(trn[2], 4.f, this->m_epsilon);

  // Test constructor from matrix
  auto m44 = trf2.matrix();
  algebra::get_transform3D_t<TypeParam> trfm(m44);

  // Make sure that algebra::getter:vector can be called.
  [[maybe_unused]] algebra::get_vector3D_t<TypeParam>
      test_vector;  // we need to declare a variable in order to use the
                    // [[maybe_unused]] attribute here

  test_vector = algebra::getter::vector<3>(m44, 0, 0);
  test_vector = algebra::getter::vector<3>(m44, 0, 1);
  test_vector = algebra::getter::vector<3>(m44, 0, 2);

  // Test constructor from inverse matrix
  auto m44_inv = trf2.matrix_inverse();

  // Make sure that algebra::getter:vector can be called.
  test_vector = algebra::getter::vector<3>(m44_inv, 0, 0);
  test_vector = algebra::getter::vector<3>(m44_inv, 0, 1);
  test_vector = algebra::getter::vector<3>(m44_inv, 0, 2);

  // Re-evaluate rot and trn
  auto rotm = trfm.rotation();
  ASSERT_NEAR(algebra::getter::element(rotm, 0, 0), x[0], this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(rotm, 1, 0), x[1], this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(rotm, 2, 0), x[2], this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(rotm, 0, 1), y[0], this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(rotm, 1, 1), y[1], this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(rotm, 2, 1), y[2], this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(rotm, 0, 2), z[0], this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(rotm, 1, 2), z[1], this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(rotm, 2, 2), z[2], this->m_epsilon);

  auto trnm = trfm.translation();
  ASSERT_NEAR(trnm[0], 2.f, this->m_epsilon);
  ASSERT_NEAR(trnm[1], 3.f, this->m_epsilon);
  ASSERT_NEAR(trnm[2], 4.f, this->m_epsilon);

  // Check a contruction from an array[16]
  std::array<algebra::get_scalar_t<TypeParam>, 16> matray_helper = {
      1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f,
      0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 0.f};
  typename algebra::get_transform3D_t<TypeParam>::template array_type<16>
      matray;
  for (unsigned int i = 0u; i < 16u; ++i) {
    matray[i] = matray_helper[i];
  }
  algebra::get_transform3D_t<TypeParam> trfma(matray);

  // Re-evaluate rot and trn
  auto rotma = trfma.rotation();
  ASSERT_NEAR(algebra::getter::element(rotma, 0, 0), 1.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(rotma, 1, 0), 0.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(rotma, 2, 0), 0.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(rotma, 0, 1), 0.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(rotma, 1, 1), 1.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(rotma, 2, 1), 0.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(rotma, 0, 2), 0.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(rotma, 1, 2), 0.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(rotma, 2, 2), 1.f, this->m_epsilon);

  auto trnma = trfma.translation();
  ASSERT_NEAR(trnma[0], 0.f, this->m_epsilon);
  ASSERT_NEAR(trnma[1], 0.f, this->m_epsilon);
  ASSERT_NEAR(trnma[2], 0.f, this->m_epsilon);
}

// This test global coordinate transforms
TYPED_TEST_P(transform_test, global_transformations) {

  // Preparation work
  algebra::get_vector3D_t<TypeParam> z = algebra::vector::normalize(
      algebra::get_vector3D_t<TypeParam>{3.f, 2.f, 1.f});
  algebra::get_vector3D_t<TypeParam> x = algebra::vector::normalize(
      algebra::get_vector3D_t<TypeParam>{2.f, -3.f, 0.f});
  [[maybe_unused]] algebra::get_vector3D_t<TypeParam> y =
      algebra::vector::cross(z, x);
  algebra::get_point3D_t<TypeParam> t = {2.f, 3.f, 4.f};
  algebra::get_transform3D_t<TypeParam> trf(t, z, x);

  // Check that local origin translates into global translation
  algebra::get_point3D_t<TypeParam> lzero = {0.f, 0.f, 0.f};
  algebra::get_point3D_t<TypeParam> gzero = trf.point_to_global(lzero);
  ASSERT_NEAR(gzero[0], t[0], this->m_epsilon);
  ASSERT_NEAR(gzero[1], t[1], this->m_epsilon);
  ASSERT_NEAR(gzero[2], t[2], this->m_epsilon);

  // Check a round trip for point
  algebra::get_point3D_t<TypeParam> lpoint = {3.f, 4.f, 5.f};
  algebra::get_point3D_t<TypeParam> gpoint = trf.point_to_global(lpoint);
  algebra::get_point3D_t<TypeParam> lpoint_r = trf.point_to_local(gpoint);
  ASSERT_NEAR(lpoint[0], lpoint_r[0], this->m_isclose);
  ASSERT_NEAR(lpoint[1], lpoint_r[1], this->m_isclose);
  ASSERT_NEAR(lpoint[2], lpoint_r[2], this->m_isclose);

  // Check a point versus vector transform
  // vector should not change if transformed by a pure translation
  algebra::get_transform3D_t<TypeParam> ttrf(t);

  algebra::get_vector3D_t<TypeParam> gvector = {1.f, 1.f, 1.f};
  algebra::get_vector3D_t<TypeParam> lvector = ttrf.vector_to_local(gvector);
  ASSERT_NEAR(gvector[0], lvector[0], this->m_isclose);
  ASSERT_NEAR(gvector[1], lvector[1], this->m_isclose);
  ASSERT_NEAR(gvector[2], lvector[2], this->m_isclose);

  // Check a round trip for vector
  algebra::get_vector3D_t<TypeParam> lvectorB = {7.f, 8.f, 9.f};
  algebra::get_vector3D_t<TypeParam> gvectorB = trf.vector_to_local(lvectorB);
  algebra::get_vector3D_t<TypeParam> lvectorC = trf.vector_to_global(gvectorB);
  ASSERT_NEAR(lvectorB[0], lvectorC[0], this->m_isclose);
  ASSERT_NEAR(lvectorB[1], lvectorC[1], this->m_isclose);
  ASSERT_NEAR(lvectorB[2], lvectorC[2], this->m_isclose);
}

// clang-format off
#define ALGEBRA_REGISTER_TRANSFORM_TESTS(...) \
  REGISTER_TYPED_TEST_SUITE_P(transform_test \
    , transform3D \
    , global_transformations \
)
// clang-format on

}  // namespace algebra::test
