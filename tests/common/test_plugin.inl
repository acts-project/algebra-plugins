/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/types.hpp"

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <climits>
#include <cmath>

/// @note __plugin has to be defined with a preprocessor command
using namespace algebra;

// Two-dimensional definitions
using point2cart = __plugin::point2;
using point2pol = __plugin::point2;
using point2cyl = __plugin::point2;
__plugin::cartesian2 cartesian2;
__plugin::polar2 polar2;
__plugin::cylindrical2 cylindrical2;

// Three-dimensional definitions
using transform3 = __plugin::transform3;
using vector3 = __plugin::vector3;
using point3 = __plugin::point3;

constexpr scalar epsilon = std::numeric_limits<scalar>::epsilon();
constexpr scalar isclose = 1e-5;

// This defines the local frame test suite
TEST(ALGEBRA_PLUGIN, local_vectors) {
  // Construction
  point2cart vA{0., 1.};
  ASSERT_EQ(vA[0], 0.);
  ASSERT_EQ(vA[1], 1.);

  // Assignment
  point2cart vB = vA;
  ASSERT_EQ(vB[0], 0.);
  ASSERT_EQ(vB[1], 1.);

  // Addition
  point2cart vC = vA + vB;
  ASSERT_EQ(vC[0], 0.);
  ASSERT_EQ(vC[1], 2.);

  // Multiplication by scalar
  point2cart vC2 = vC * 2.;
  ASSERT_EQ(vC2[0], 0.);
  ASSERT_EQ(vC2[1], 4.);

  // Cast operations to phi, theta, eta, perp
  point2cart vD{1., 1.};
  scalar phi = getter::phi(vD);
  ASSERT_NEAR(phi, M_PI_4, epsilon);

  scalar perp = getter::perp(vD);
  ASSERT_NEAR(perp, std::sqrt(2.), epsilon);

  scalar norm = getter::norm(vD);
  ASSERT_NEAR(norm, std::sqrt(2.), epsilon);

  auto vDnorm = vector::normalize(vD);
  ASSERT_NEAR(vDnorm[0], 1. / std::sqrt(2.), epsilon);
  ASSERT_NEAR(vDnorm[1], 1. / std::sqrt(2.), epsilon);
}

// This defines the vector3 test suite
TEST(ALGEBRA_PLUGIN, vector3) {
  // Construction
  vector3 vA{0., 1., 2.};
  ASSERT_EQ(vA[0], 0.);
  ASSERT_EQ(vA[1], 1.);
  ASSERT_EQ(vA[2], 2.);

  // Assignment
  vector3 vB = vA;
  ASSERT_EQ(vB[0], 0.);
  ASSERT_EQ(vB[1], 1.);
  ASSERT_EQ(vB[2], 2.);

  // Addition
  vector3 vC = vA + vB;
  ASSERT_EQ(vC[0], 0.);
  ASSERT_EQ(vC[1], 2.);
  ASSERT_EQ(vC[2], 4.);

  // Multiplication by scalar
  vector3 vC2 = vC * 2.0;
  ASSERT_EQ(vC2[0], 0.);
  ASSERT_EQ(vC2[1], 4.);
  ASSERT_EQ(vC2[2], 8.);

  // Cast operations to phi, theta, eta, perp
  vector3 vD{1., 1., 1.};
  scalar phi = getter::phi(vD);
  ASSERT_NEAR(phi, M_PI_4, epsilon);

  scalar theta = getter::theta(vD);
  ASSERT_NEAR(theta, std::atan2(std::sqrt(2.), 1.), epsilon);

  scalar eta = getter::eta(vD);
  ASSERT_NEAR(eta, 0.65847891569137573, isclose);

  scalar perp = getter::perp(vD);
  ASSERT_NEAR(perp, std::sqrt(2.), epsilon);

  scalar norm = getter::norm(vD);
  ASSERT_NEAR(norm, std::sqrt(3.), epsilon);
}

// This defines the vector operation test suite
TEST(ALGEBRA_PLUGIN, getter) {
  vector3 v3{1., 1., 1.};

  // Normalization
  auto v3n = vector::normalize(v3);
  ASSERT_NEAR(v3n[0], 1. / std::sqrt(3.), epsilon);
  ASSERT_NEAR(v3n[1], 1. / std::sqrt(3.), epsilon);
  ASSERT_NEAR(v3n[2], 1. / std::sqrt(3.), epsilon);

  // Cross product
  vector3 z = vector::normalize(vector3{3., 2., 1.});
  vector3 x = vector::normalize(vector3{2., -3., 0.});
  vector3 y = vector::cross(z, x);

  // Check with dot product
  ASSERT_NEAR(vector::dot(x, y), 0., epsilon);
  ASSERT_NEAR(vector::dot(y, z), 0., epsilon);
  ASSERT_NEAR(vector::dot(z, x), 0., epsilon);
}

// This defines the transform3 test suite
TEST(ALGEBRA_PLUGIN, transform3) {
  // Preparatioon work
  vector3 z = vector::normalize(vector3{3., 2., 1.});
  vector3 x = vector::normalize(vector3{2., -3., 0.});
  vector3 y = vector::cross(z, x);
  point3 t = {2., 3., 4.};

  // Test constructor from t, z, x
  transform3 trf(t, z, x);

  ASSERT_TRUE(trf == trf);

  const auto rot = trf.rotation();
#ifndef __plugin_without_matrix_element_accessor
  ASSERT_NEAR(rot(0, 0), x[0], epsilon);
  ASSERT_NEAR(rot(1, 0), x[1], epsilon);
  ASSERT_NEAR(rot(2, 0), x[2], epsilon);
  ASSERT_NEAR(rot(0, 1), y[0], epsilon);
  ASSERT_NEAR(rot(1, 1), y[1], epsilon);
  ASSERT_NEAR(rot(2, 1), y[2], epsilon);
  ASSERT_NEAR(rot(0, 2), z[0], epsilon);
  ASSERT_NEAR(rot(1, 2), z[1], epsilon);
  ASSERT_NEAR(rot(2, 2), z[2], epsilon);
#endif

  auto trn = trf.translation();
  ASSERT_NEAR(trn[0], 2., epsilon);
  ASSERT_NEAR(trn[1], 3., epsilon);
  ASSERT_NEAR(trn[2], 4., epsilon);

  // Test constructor from matrix
  auto m44 = trf.matrix();
  transform3 trfm(m44);

  // Re-evaluate rot and trn
  auto rotm = trfm.rotation();
#ifndef __plugin_without_matrix_element_accessor
  ASSERT_NEAR(rotm(0, 0), x[0], epsilon);
  ASSERT_NEAR(rotm(1, 0), x[1], epsilon);
  ASSERT_NEAR(rotm(2, 0), x[2], epsilon);
  ASSERT_NEAR(rotm(0, 1), y[0], epsilon);
  ASSERT_NEAR(rotm(1, 1), y[1], epsilon);
  ASSERT_NEAR(rotm(2, 1), y[2], epsilon);
  ASSERT_NEAR(rotm(0, 2), z[0], epsilon);
  ASSERT_NEAR(rotm(1, 2), z[1], epsilon);
  ASSERT_NEAR(rotm(2, 2), z[2], epsilon);
#endif

  auto trnm = trfm.translation();
  ASSERT_NEAR(trnm[0], 2., epsilon);
  ASSERT_NEAR(trnm[1], 3., epsilon);
  ASSERT_NEAR(trnm[2], 4., epsilon);

  // Check a contruction from an array[16]
  array_s<scalar, 16> matray = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0};
  transform3 trfma(matray);

// Re-evaluate rot and trn
#ifndef __plugin_without_matrix_element_accessor
  auto rotma = trfma.rotation();
  ASSERT_NEAR(rotma(0, 0), 1., epsilon);
  ASSERT_NEAR(rotma(1, 0), 0., epsilon);
  ASSERT_NEAR(rotma(2, 0), 0., epsilon);
  ASSERT_NEAR(rotma(0, 1), 0., epsilon);
  ASSERT_NEAR(rotma(1, 1), 1., epsilon);
  ASSERT_NEAR(rotma(2, 1), 0., epsilon);
  ASSERT_NEAR(rotma(0, 2), 0., epsilon);
  ASSERT_NEAR(rotma(1, 2), 0., epsilon);
  ASSERT_NEAR(rotma(2, 2), 1., epsilon);
#endif

  auto trnma = trfma.translation();
  ASSERT_NEAR(trnma[0], 0., epsilon);
  ASSERT_NEAR(trnma[1], 0., epsilon);
  ASSERT_NEAR(trnma[2], 0., epsilon);
}

// This test global coordinate transforms
TEST(ALGEBRA_PLUGIN, global_transformations) {
  // Preparatioon work
  vector3 z = vector::normalize(vector3{3., 2., 1.});
  vector3 x = vector::normalize(vector3{2., -3., 0.});
  vector3 y = vector::cross(z, x);
  point3 t = {2., 3., 4.};
  transform3 trf(t, z, x);

  // Check that local origin translates into global translation
  point3 lzero = {0., 0., 0.};
  auto gzero = trf.point_to_global(lzero);
  ASSERT_NEAR(gzero[0], t[0], epsilon);
  ASSERT_NEAR(gzero[1], t[1], epsilon);
  ASSERT_NEAR(gzero[2], t[2], epsilon);

  // Check a round trip for point
  point3 lpoint = {3., 4., 5.};
  auto gpoint = trf.point_to_global(lpoint);
  auto lpoint_r = trf.point_to_local(gpoint);
  ASSERT_NEAR(lpoint[0], lpoint_r[0], isclose);
  ASSERT_NEAR(lpoint[1], lpoint_r[1], isclose);
  ASSERT_NEAR(lpoint[2], lpoint_r[2], isclose);

  // Check a point versus vector transform
  // vector should not change if transformed by a pure translation
  point3 tt = {2., 3., 4.};
  transform3 ttrf(t);

  vector3 gvector = {1., 1., 1};
  auto lvector = ttrf.vector_to_local(gvector);
  ASSERT_NEAR(gvector[0], lvector[0], isclose);
  ASSERT_NEAR(gvector[1], lvector[1], isclose);
  ASSERT_NEAR(gvector[2], lvector[2], isclose);

  // Check a round trip for vector
  vector3 lvectorB = {7., 8., 9};
  vector3 gvectorB = trf.vector_to_local(lvectorB);
  vector3 lvectorC = trf.vector_to_global(gvectorB);
  ASSERT_NEAR(lvectorB[0], lvectorC[0], isclose);
  ASSERT_NEAR(lvectorB[1], lvectorC[1], isclose);
  ASSERT_NEAR(lvectorB[2], lvectorC[2], isclose);
}

// This test local coordinate transforms
TEST(ALGEBRA_PLUGIN, local_transformations) {
  point2cart p2 = {3., 3.};
  point3 p3 = {3., 3., 5.};

  auto cart2fromp3 = cartesian2(p3);
  ASSERT_NEAR(p2[0], cart2fromp3[0], epsilon);
  ASSERT_NEAR(p2[1], cart2fromp3[1], epsilon);

  auto polfrom2 = polar2(p2);
  auto polfrom3 = polar2(p3);

  // Check r-phi
  ASSERT_NEAR(polfrom2[0], sqrt(18.), isclose);
  ASSERT_NEAR(polfrom2[1], M_PI_4, isclose);

  // Need to be identical
  ASSERT_NEAR(polfrom2[0], polfrom3[0], epsilon);
  ASSERT_NEAR(polfrom2[1], polfrom3[1], epsilon);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}
