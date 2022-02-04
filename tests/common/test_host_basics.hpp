/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "test_base.hpp"

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <array>
#include <climits>
#include <cmath>

/// Test case class, to be specialised for the different plugins
template <typename T>
class test_host_basics : public testing::Test, public test_base<T> {};
TYPED_TEST_SUITE_P(test_host_basics);

// This defines the local frame test suite
TYPED_TEST_P(test_host_basics, local_vectors) {

  // Construction
  typename TypeParam::point2 vA{0., 1.};
  ASSERT_EQ(vA[0], 0.);
  ASSERT_EQ(vA[1], 1.);

  // Assignment
  typename TypeParam::point2 vB = vA;
  ASSERT_EQ(vB[0], 0.);
  ASSERT_EQ(vB[1], 1.);

  // Addition
  typename TypeParam::point2 vC = vA + vB;
  ASSERT_EQ(vC[0], 0.);
  ASSERT_EQ(vC[1], 2.);

  // Multiplication by scalar
  typename TypeParam::point2 vC2 = vC * 2.;
  ASSERT_EQ(vC2[0], 0.);
  ASSERT_EQ(vC2[1], 4.);

  // Cast operations to phi, theta, eta, perp
  typename TypeParam::vector2 vD{1., 1.};
  typename TypeParam::scalar phi = algebra::getter::phi(vD);
  ASSERT_NEAR(phi, M_PI_4, this->m_epsilon);

  typename TypeParam::scalar perp = algebra::getter::perp(vD);
  ASSERT_NEAR(perp, std::sqrt(2.), this->m_epsilon);

  typename TypeParam::scalar norm = algebra::getter::norm(vD);
  ASSERT_NEAR(norm, std::sqrt(2.), this->m_epsilon);

  typename TypeParam::vector2 vDnorm = algebra::vector::normalize(vD);
  ASSERT_NEAR(vDnorm[0], 1. / std::sqrt(2.), this->m_epsilon);
  ASSERT_NEAR(vDnorm[1], 1. / std::sqrt(2.), this->m_epsilon);
}

// This defines the vector3 test suite
TYPED_TEST_P(test_host_basics, vector3) {

  // Construction
  typename TypeParam::vector3 vA{0., 1., 2.};
  ASSERT_EQ(vA[0], 0.);
  ASSERT_EQ(vA[1], 1.);
  ASSERT_EQ(vA[2], 2.);

  // Assignment
  typename TypeParam::vector3 vB = vA;
  ASSERT_EQ(vB[0], 0.);
  ASSERT_EQ(vB[1], 1.);
  ASSERT_EQ(vB[2], 2.);

  // Addition
  typename TypeParam::vector3 vC = vA + vB;
  ASSERT_EQ(vC[0], 0.);
  ASSERT_EQ(vC[1], 2.);
  ASSERT_EQ(vC[2], 4.);

  // Multiplication by scalar
  typename TypeParam::vector3 vC2 = vC * 2.0;
  ASSERT_EQ(vC2[0], 0.);
  ASSERT_EQ(vC2[1], 4.);
  ASSERT_EQ(vC2[2], 8.);

  // Cast operations to phi, theta, eta, perp
  typename TypeParam::vector3 vD{1., 1., 1.};
  typename TypeParam::scalar phi = algebra::getter::phi(vD);
  ASSERT_NEAR(phi, M_PI_4, this->m_epsilon);

  typename TypeParam::scalar theta = algebra::getter::theta(vD);
  ASSERT_NEAR(theta, std::atan2(std::sqrt(2.), 1.), this->m_epsilon);

  typename TypeParam::scalar eta = algebra::getter::eta(vD);
  ASSERT_NEAR(eta, 0.65847891569137573, this->m_isclose);

  typename TypeParam::scalar perp = algebra::getter::perp(vD);
  ASSERT_NEAR(perp, std::sqrt(2.), this->m_epsilon);

  typename TypeParam::scalar norm = algebra::getter::norm(vD);
  ASSERT_NEAR(norm, std::sqrt(3.), this->m_epsilon);
}

// This defines the vector operation test suite
TYPED_TEST_P(test_host_basics, getter) {

  typename TypeParam::vector3 v3{1., 1., 1.};

  // Normalization
  typename TypeParam::vector3 v3n = algebra::vector::normalize(v3);
  ASSERT_NEAR(v3n[0], 1. / std::sqrt(3.), this->m_epsilon);
  ASSERT_NEAR(v3n[1], 1. / std::sqrt(3.), this->m_epsilon);
  ASSERT_NEAR(v3n[2], 1. / std::sqrt(3.), this->m_epsilon);

  // Cross product
  typename TypeParam::vector3 z =
      algebra::vector::normalize(typename TypeParam::vector3{3., 2., 1.});
  typename TypeParam::vector3 x =
      algebra::vector::normalize(typename TypeParam::vector3{2., -3., 0.});
  typename TypeParam::vector3 y = algebra::vector::cross(z, x);

  // Check with dot product
  ASSERT_NEAR(algebra::vector::dot(x, y), 0., this->m_epsilon);
  ASSERT_NEAR(algebra::vector::dot(y, z), 0., this->m_epsilon);
  ASSERT_NEAR(algebra::vector::dot(z, x), 0., this->m_epsilon);
}

// This defines the transform3 test suite
TYPED_TEST_P(test_host_basics, transform3) {

  // Preparatioon work
  typename TypeParam::vector3 z =
      algebra::vector::normalize(typename TypeParam::vector3{3., 2., 1.});
  typename TypeParam::vector3 x =
      algebra::vector::normalize(typename TypeParam::vector3{2., -3., 0.});
  typename TypeParam::vector3 y = algebra::vector::cross(z, x);
  typename TypeParam::point3 t = {2., 3., 4.};

  // Test constructor from t, z, x
  typename TypeParam::transform3 trf1(t, z, x);
  ASSERT_TRUE(trf1 == trf1);
  typename TypeParam::transform3 trf2;
  trf2 = trf1;

  // Helper object for the matrix checks.
  auto element_getter = typename TypeParam::transform3::element_getter();

  const auto rot = trf2.rotation();
  ASSERT_NEAR(element_getter(rot, 0, 0), x[0], this->m_epsilon);
  ASSERT_NEAR(element_getter(rot, 1, 0), x[1], this->m_epsilon);
  ASSERT_NEAR(element_getter(rot, 2, 0), x[2], this->m_epsilon);
  ASSERT_NEAR(element_getter(rot, 0, 1), y[0], this->m_epsilon);
  ASSERT_NEAR(element_getter(rot, 1, 1), y[1], this->m_epsilon);
  ASSERT_NEAR(element_getter(rot, 2, 1), y[2], this->m_epsilon);
  ASSERT_NEAR(element_getter(rot, 0, 2), z[0], this->m_epsilon);
  ASSERT_NEAR(element_getter(rot, 1, 2), z[1], this->m_epsilon);
  ASSERT_NEAR(element_getter(rot, 2, 2), z[2], this->m_epsilon);

  auto trn = trf2.translation();
  ASSERT_NEAR(trn[0], 2., this->m_epsilon);
  ASSERT_NEAR(trn[1], 3., this->m_epsilon);
  ASSERT_NEAR(trn[2], 4., this->m_epsilon);

  // Test constructor from matrix
  auto m44 = trf2.matrix();
  typename TypeParam::transform3 trfm(m44);

  // Make sure that algebra::getter:vector can be called.
  (void)algebra::getter::vector<3>(m44, 0, 0);
  (void)algebra::getter::vector<3>(m44, 0, 1);
  (void)algebra::getter::vector<3>(m44, 0, 2);

  // Re-evaluate rot and trn
  auto rotm = trfm.rotation();
  ASSERT_NEAR(element_getter(rotm, 0, 0), x[0], this->m_epsilon);
  ASSERT_NEAR(element_getter(rotm, 1, 0), x[1], this->m_epsilon);
  ASSERT_NEAR(element_getter(rotm, 2, 0), x[2], this->m_epsilon);
  ASSERT_NEAR(element_getter(rotm, 0, 1), y[0], this->m_epsilon);
  ASSERT_NEAR(element_getter(rotm, 1, 1), y[1], this->m_epsilon);
  ASSERT_NEAR(element_getter(rotm, 2, 1), y[2], this->m_epsilon);
  ASSERT_NEAR(element_getter(rotm, 0, 2), z[0], this->m_epsilon);
  ASSERT_NEAR(element_getter(rotm, 1, 2), z[1], this->m_epsilon);
  ASSERT_NEAR(element_getter(rotm, 2, 2), z[2], this->m_epsilon);

  auto trnm = trfm.translation();
  ASSERT_NEAR(trnm[0], 2., this->m_epsilon);
  ASSERT_NEAR(trnm[1], 3., this->m_epsilon);
  ASSERT_NEAR(trnm[2], 4., this->m_epsilon);

  // Check a contruction from an array[16]
  std::array<typename TypeParam::scalar, 16> matray_helper = {
      1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0};
  typename TypeParam::transform3::template array_type<
      typename TypeParam::scalar, 16>
      matray;
  for (int i = 0; i < 16; ++i) {
    matray[i] = matray_helper[i];
  }
  typename TypeParam::transform3 trfma(matray);

  // Re-evaluate rot and trn
  auto rotma = trfma.rotation();
  ASSERT_NEAR(element_getter(rotma, 0, 0), 1., this->m_epsilon);
  ASSERT_NEAR(element_getter(rotma, 1, 0), 0., this->m_epsilon);
  ASSERT_NEAR(element_getter(rotma, 2, 0), 0., this->m_epsilon);
  ASSERT_NEAR(element_getter(rotma, 0, 1), 0., this->m_epsilon);
  ASSERT_NEAR(element_getter(rotma, 1, 1), 1., this->m_epsilon);
  ASSERT_NEAR(element_getter(rotma, 2, 1), 0., this->m_epsilon);
  ASSERT_NEAR(element_getter(rotma, 0, 2), 0., this->m_epsilon);
  ASSERT_NEAR(element_getter(rotma, 1, 2), 0., this->m_epsilon);
  ASSERT_NEAR(element_getter(rotma, 2, 2), 1., this->m_epsilon);

  auto trnma = trfma.translation();
  ASSERT_NEAR(trnma[0], 0., this->m_epsilon);
  ASSERT_NEAR(trnma[1], 0., this->m_epsilon);
  ASSERT_NEAR(trnma[2], 0., this->m_epsilon);
}

// This test global coordinate transforms
TYPED_TEST_P(test_host_basics, global_transformations) {

  // Preparatioon work
  typename TypeParam::vector3 z =
      algebra::vector::normalize(typename TypeParam::vector3{3., 2., 1.});
  typename TypeParam::vector3 x =
      algebra::vector::normalize(typename TypeParam::vector3{2., -3., 0.});
  typename TypeParam::vector3 y = algebra::vector::cross(z, x);
  (void)y;
  typename TypeParam::point3 t = {2., 3., 4.};
  typename TypeParam::transform3 trf(t, z, x);

  // Check that local origin translates into global translation
  typename TypeParam::point3 lzero = {0., 0., 0.};
  typename TypeParam::point3 gzero = trf.point_to_global(lzero);
  ASSERT_NEAR(gzero[0], t[0], this->m_epsilon);
  ASSERT_NEAR(gzero[1], t[1], this->m_epsilon);
  ASSERT_NEAR(gzero[2], t[2], this->m_epsilon);

  // Check a round trip for point
  typename TypeParam::point3 lpoint = {3., 4., 5.};
  typename TypeParam::point3 gpoint = trf.point_to_global(lpoint);
  typename TypeParam::point3 lpoint_r = trf.point_to_local(gpoint);
  ASSERT_NEAR(lpoint[0], lpoint_r[0], this->m_isclose);
  ASSERT_NEAR(lpoint[1], lpoint_r[1], this->m_isclose);
  ASSERT_NEAR(lpoint[2], lpoint_r[2], this->m_isclose);

  // Check a point versus vector transform
  // vector should not change if transformed by a pure translation
  typename TypeParam::transform3 ttrf(t);

  typename TypeParam::vector3 gvector = {1., 1., 1};
  typename TypeParam::vector3 lvector = ttrf.vector_to_local(gvector);
  ASSERT_NEAR(gvector[0], lvector[0], this->m_isclose);
  ASSERT_NEAR(gvector[1], lvector[1], this->m_isclose);
  ASSERT_NEAR(gvector[2], lvector[2], this->m_isclose);

  // Check a round trip for vector
  typename TypeParam::vector3 lvectorB = {7., 8., 9};
  typename TypeParam::vector3 gvectorB = trf.vector_to_local(lvectorB);
  typename TypeParam::vector3 lvectorC = trf.vector_to_global(gvectorB);
  ASSERT_NEAR(lvectorB[0], lvectorC[0], this->m_isclose);
  ASSERT_NEAR(lvectorB[1], lvectorC[1], this->m_isclose);
  ASSERT_NEAR(lvectorB[2], lvectorC[2], this->m_isclose);
}

// This test local coordinate transforms
TYPED_TEST_P(test_host_basics, local_transformations) {

  typename TypeParam::point2 p2 = {3., 3.};
  typename TypeParam::point3 p3 = {3., 3., 5.};

  typename TypeParam::cartesian2 c2;
  auto cart2fromp3 = c2(p3);
  ASSERT_NEAR(p2[0], cart2fromp3[0], this->m_epsilon);
  ASSERT_NEAR(p2[1], cart2fromp3[1], this->m_epsilon);

  typename TypeParam::polar2 pol2;
  auto polfrom2 = pol2(p2);
  auto polfrom3 = pol2(p3);

  // Check r-phi
  ASSERT_NEAR(polfrom2[0], std::sqrt(18.), this->m_isclose);
  ASSERT_NEAR(polfrom2[1], M_PI_4, this->m_isclose);

  // Need to be identical
  ASSERT_NEAR(polfrom2[0], polfrom3[0], this->m_epsilon);
  ASSERT_NEAR(polfrom2[1], polfrom3[1], this->m_epsilon);
}

REGISTER_TYPED_TEST_SUITE_P(test_host_basics, local_vectors, vector3, getter,
                            transform3, global_transformations,
                            local_transformations);
