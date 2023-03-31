/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.f
 */

// Project include(s).
#include "algebra/vc_soa.hpp"

// GoogleTest include(s).
#include <gtest/gtest.h>

#include <array>
#include <iostream>

using namespace algebra;

using scalar_t = float;

constexpr float tol{1e-6f};

/// This test the vector functions on an SoA (Vc::Vector) based vector
TEST(test_vc_host, vc_soa_vector) {

  using vector3_v = vc_soa::vector3<scalar_t>;
  // Value type is Vc::Vector<float>
  using value_t = typename vector3_v::value_type;

  vector3_v a{value_t::One(), 2.f * value_t::One(), 3.f * value_t::One()};
  vector3_v b{4.f * value_t::One(), 5.f * value_t::One(), 6.f * value_t::One()};

  EXPECT_TRUE((a[0] == value_t(1.f)).isFull());
  EXPECT_TRUE((a[1] == value_t(2.f)).isFull());
  EXPECT_TRUE((a[2] == value_t(3.f)).isFull());

  // Masked comparison
  auto m = a.compare(a);
  EXPECT_TRUE(m[0].isFull());
  EXPECT_TRUE(m[1].isFull());
  EXPECT_TRUE(m[2].isFull());

  m = a.compare(b);
  EXPECT_FALSE(m[0].isFull());
  EXPECT_FALSE(m[1].isFull());
  EXPECT_FALSE(m[2].isFull());

  // Full comparisons
  EXPECT_TRUE(a == a);
  EXPECT_FALSE(a == b);

  // Addition
  auto v_add = a + b;
  EXPECT_TRUE((v_add[0] == value_t(5.f)).isFull());
  EXPECT_TRUE((v_add[1] == value_t(7.f)).isFull());
  EXPECT_TRUE((v_add[2] == value_t(9.f)).isFull());

  // Subration
  auto v_sub = a - b;
  EXPECT_TRUE((v_sub[0] == value_t(-3.f)).isFull());
  EXPECT_TRUE((v_sub[1] == value_t(-3.f)).isFull());
  EXPECT_TRUE((v_sub[2] == value_t(-3.f)).isFull());

  // Multiplication
  auto v_mul = a * b;
  EXPECT_TRUE((v_mul[0] == value_t(4.f)).isFull());
  EXPECT_TRUE((v_mul[1] == value_t(10.f)).isFull());
  EXPECT_TRUE((v_mul[2] == value_t(18.f)).isFull());

  // Division
  auto v_div = a / b;
  EXPECT_TRUE((v_div[0] == value_t(0.25f)).isFull());
  EXPECT_TRUE((v_div[1] == value_t(0.4f)).isFull());
  EXPECT_TRUE((v_div[2] == value_t(0.5f)).isFull());

  // Scalar multiplication
  auto v_smul = 2.f * b;
  EXPECT_TRUE((v_smul[0] == value_t(8.f)).isFull());
  EXPECT_TRUE((v_smul[1] == value_t(10.f)).isFull());
  EXPECT_TRUE((v_smul[2] == value_t(12.f)).isFull());

  // Expression
  auto v_expr = (b / a) - (2.5f * b) + vector3_v{};
  EXPECT_TRUE((v_expr[0] == value_t(-6.f)).isFull());
  EXPECT_TRUE((v_expr[1] == value_t(-10.f)).isFull());
  EXPECT_TRUE((v_expr[2] == value_t(-13.f)).isFull());

  auto d{vector::dot(a, b)};
  EXPECT_TRUE((d == value_t(32.f)).isFull());

  value_t norms_a{getter::norm(vector::normalize(a))};
  value_t norms_b{getter::norm(vector::normalize(b))};
  for (unsigned int i{0u}; i < norms_a.size(); ++i) {
    EXPECT_NEAR(norms_a[i], 1.f, tol);
    EXPECT_NEAR(norms_b[i], 1.f, tol);
  }

  auto cr{vector::cross(a, b)};
  EXPECT_TRUE((cr[0] == value_t(-3.f)).isFull());
  EXPECT_TRUE((cr[1] == value_t(6.f)).isFull());
  EXPECT_TRUE((cr[2] == value_t(-3.f)).isFull());

  static_assert(std::is_convertible_v<decltype(v_expr), vector3_v>,
                "expression type not convertible");
}

/// This test the getter functions on an SoA (Vc::Vector) based vector
TEST(test_vc_host, vc_soa_getter) {

  using vector3_v = vc_soa::vector3<scalar_t>;
  // Value type is Vc::Vector<float>
  using value_t = typename vector3_v::value_type;

  vector3_v a{value_t::One(), 2.f * value_t::One(), 3.f * value_t::One()};

  // All results in the vector are the same, so only check the first one

  // Phi angle
  auto v_phi = getter::phi(a);
  EXPECT_NEAR(v_phi[0], static_cast<scalar_t>(std::atan2(2., 1.)), tol);

  // Perpendicular projection
  auto v_perp = getter::perp(a);
  EXPECT_NEAR(v_perp[0], std::sqrt(5.), tol);

  // Theta angle
  auto v_theta = getter::theta(a);
  EXPECT_NEAR(v_theta[0], static_cast<scalar_t>(std::atan2(std::sqrt(5.), 3.)),
              tol);

  // Norm of the vector
  auto v_norm = getter::norm(a);
  EXPECT_NEAR(v_norm[0], std::sqrt(14.), tol);

  // Eta of the vector
  auto v_eta = getter::eta(a);
  EXPECT_NEAR(v_eta[0],
              static_cast<scalar_t>(std::atanh(1. / std::sqrt(14.) * 3.)), tol);
}