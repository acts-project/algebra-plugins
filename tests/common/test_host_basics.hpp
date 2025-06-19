/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/utils/approximately_equal.hpp"
#include "algebra/utils/casts.hpp"
#include "algebra/utils/print.hpp"

// Local include(s).
#include "test_base.hpp"

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>

using namespace algebra;

/// Test case class, to be specialised for the different plugins - vectors
template <typename T>
class test_host_basics_vector : public testing::Test, public test_base<T> {};
TYPED_TEST_SUITE_P(test_host_basics_vector);

/// Test case class, to be specialised for the different plugins - matrices
template <typename T>
class test_host_basics_matrix : public testing::Test, public test_base<T> {
 protected:
  template <typename A, std::size_t ROWS, std::size_t COLS>
  void test_matrix_ops_any_matrix() {
    // Test the set_product method.
    {
      typename A::template matrix<ROWS, ROWS> m1;
      typename A::template matrix<ROWS, COLS> m2;

      for (std::size_t i = 0; i < ROWS; ++i) {
        for (std::size_t j = 0; j < ROWS; ++j) {
          algebra::getter::element(m1, i, j) =
              static_cast<algebra::traits::scalar_t<decltype(m1)>>(i * ROWS +
                                                                   j);
        }
      }

      for (std::size_t i = 0; i < ROWS; ++i) {
        for (std::size_t j = 0; j < COLS; ++j) {
          algebra::getter::element(m2, i, j) =
              static_cast<algebra::traits::scalar_t<decltype(m1)>>(i * COLS +
                                                                   j);
        }
      }

      {
        typename A::template matrix<ROWS, COLS> r1 = m1 * m2;
        typename A::template matrix<ROWS, COLS> r2;
        algebra::matrix::set_product(r2, m1, m2);

        for (std::size_t i = 0; i < ROWS; ++i) {
          for (std::size_t j = 0; j < COLS; ++j) {
            ASSERT_NEAR(algebra::getter::element(r1, i, j),
                        algebra::getter::element(r2, i, j), this->m_epsilon);
          }
        }
      }
    }

    // Test the set_product_right_transpose method.
    {
      typename A::template matrix<ROWS, ROWS> m1;
      typename A::template matrix<COLS, ROWS> m2;

      for (std::size_t i = 0; i < ROWS; ++i) {
        for (std::size_t j = 0; j < ROWS; ++j) {
          algebra::getter::element(m1, i, j) =
              static_cast<algebra::traits::scalar_t<decltype(m1)>>(i * ROWS +
                                                                   j);
        }
      }

      for (std::size_t i = 0; i < COLS; ++i) {
        for (std::size_t j = 0; j < ROWS; ++j) {
          algebra::getter::element(m2, i, j) =
              static_cast<algebra::traits::scalar_t<decltype(m1)>>(i * COLS +
                                                                   j);
        }
      }

      {
        typename A::template matrix<ROWS, COLS> r1 =
            m1 * algebra::matrix::transpose(m2);
        typename A::template matrix<ROWS, COLS> r2;
        algebra::matrix::set_product_right_transpose(r2, m1, m2);

        for (std::size_t i = 0; i < ROWS; ++i) {
          for (std::size_t j = 0; j < COLS; ++j) {
            ASSERT_NEAR(algebra::getter::element(r1, i, j),
                        algebra::getter::element(r2, i, j), this->m_epsilon);
          }
        }
      }
    }

    // Test the set_product_left_transpose method.
    {
      typename A::template matrix<ROWS, ROWS> m1;
      typename A::template matrix<ROWS, COLS> m2;

      for (std::size_t i = 0; i < ROWS; ++i) {
        for (std::size_t j = 0; j < ROWS; ++j) {
          algebra::getter::element(m1, i, j) =
              static_cast<algebra::traits::scalar_t<decltype(m1)>>(i * ROWS +
                                                                   j);
        }
      }

      for (std::size_t i = 0; i < ROWS; ++i) {
        for (std::size_t j = 0; j < COLS; ++j) {
          algebra::getter::element(m2, i, j) =
              static_cast<algebra::traits::scalar_t<decltype(m1)>>(i * COLS +
                                                                   j);
        }
      }

      {
        typename A::template matrix<ROWS, COLS> r1 =
            algebra::matrix::transpose(m1) * m2;
        typename A::template matrix<ROWS, COLS> r2;
        algebra::matrix::set_product_left_transpose(r2, m1, m2);

        for (std::size_t i = 0; i < ROWS; ++i) {
          for (std::size_t j = 0; j < COLS; ++j) {
            ASSERT_NEAR(algebra::getter::element(r1, i, j),
                        algebra::getter::element(r2, i, j), this->m_epsilon);
          }
        }
      }

      // Test the transposable_product method.
      {
        // Both untransposed
        {
          typename A::template matrix<ROWS, COLS> r1 = m1 * m2;
          typename A::template matrix<ROWS, COLS> r2 =
              algebra::matrix::transposed_product<false, false>(m1, m2);

          for (std::size_t i = 0; i < ROWS; ++i) {
            for (std::size_t j = 0; j < COLS; ++j) {
              ASSERT_NEAR(algebra::getter::element(r1, i, j),
                          algebra::getter::element(r2, i, j), this->m_epsilon);
            }
          }
        }
      }
    }
  }

  template <typename A, std::size_t N>
  void test_matrix_ops_square_matrix() {
    {
      typename A::template matrix<N, N> m1;
      typename A::template matrix<N, N> m2;

      for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
          algebra::getter::element(m1, i, j) =
              static_cast<algebra::traits::scalar_t<decltype(m1)>>(i * N + j);
          algebra::getter::element(m2, i, j) =
              static_cast<algebra::traits::scalar_t<decltype(m1)>>(
                  -1 * (i * N + j) + 42);
        }
      }

      // Test the set_product method.
      {
        typename A::template matrix<N, N> r1 = m1 * m2;
        typename A::template matrix<N, N> r2;
        algebra::matrix::set_product(r2, m1, m2);

        for (std::size_t i = 0; i < N; ++i) {
          for (std::size_t j = 0; j < N; ++j) {
            ASSERT_NEAR(algebra::getter::element(r1, i, j),
                        algebra::getter::element(r2, i, j), this->m_epsilon);
          }
        }
      }

      // Test the set_product_right_transpose method.
      {
        typename A::template matrix<N, N> r1 =
            m1 * algebra::matrix::transpose(m2);
        typename A::template matrix<N, N> r2;
        algebra::matrix::set_product_right_transpose(r2, m1, m2);

        for (std::size_t i = 0; i < N; ++i) {
          for (std::size_t j = 0; j < N; ++j) {
            ASSERT_NEAR(algebra::getter::element(r1, i, j),
                        algebra::getter::element(r2, i, j), this->m_epsilon);
          }
        }
      }

      // Test the set_product_left_transpose method.
      {
        typename A::template matrix<N, N> r1 =
            algebra::matrix::transpose(m1) * m2;
        typename A::template matrix<N, N> r2;
        algebra::matrix::set_product_left_transpose(r2, m1, m2);

        for (std::size_t i = 0; i < N; ++i) {
          for (std::size_t j = 0; j < N; ++j) {
            ASSERT_NEAR(algebra::getter::element(r1, i, j),
                        algebra::getter::element(r2, i, j), this->m_epsilon);
          }
        }
      }

      // Test the set_inplace_product_right method.
      {
        typename A::template matrix<N, N> r1 = m1 * m2;
        typename A::template matrix<N, N> r2 = m1;
        algebra::matrix::set_inplace_product_right(r2, m2);

        for (std::size_t i = 0; i < N; ++i) {
          for (std::size_t j = 0; j < N; ++j) {
            ASSERT_NEAR(algebra::getter::element(r1, i, j),
                        algebra::getter::element(r2, i, j), this->m_epsilon);
          }
        }
      }

      // Test the set_inplace_product_left method.
      {
        typename A::template matrix<N, N> r1 = m1 * m2;
        typename A::template matrix<N, N> r2 = m2;
        algebra::matrix::set_inplace_product_left(r2, m1);

        for (std::size_t i = 0; i < N; ++i) {
          for (std::size_t j = 0; j < N; ++j) {
            ASSERT_NEAR(algebra::getter::element(r1, i, j),
                        algebra::getter::element(r2, i, j), this->m_epsilon);
          }
        }
      }

      // Test the set_inplace_product_right_transpose method.
      {
        typename A::template matrix<N, N> r1 =
            m1 * algebra::matrix::transpose(m2);
        typename A::template matrix<N, N> r2 = m1;
        algebra::matrix::set_inplace_product_right_transpose(r2, m2);

        for (std::size_t i = 0; i < N; ++i) {
          for (std::size_t j = 0; j < N; ++j) {
            ASSERT_NEAR(algebra::getter::element(r1, i, j),
                        algebra::getter::element(r2, i, j), this->m_epsilon);
          }
        }
      }

      // Test the set_inplace_product_left_transpose method.
      {
        typename A::template matrix<N, N> r1 =
            algebra::matrix::transpose(m1) * m2;
        typename A::template matrix<N, N> r2 = m2;
        algebra::matrix::set_inplace_product_left_transpose(r2, m1);

        for (std::size_t i = 0; i < N; ++i) {
          for (std::size_t j = 0; j < N; ++j) {
            ASSERT_NEAR(algebra::getter::element(r1, i, j),
                        algebra::getter::element(r2, i, j), this->m_epsilon);
          }
        }
      }

      // Test the transposable_product method.
      {
        // Only left transposed
        {
          typename A::template matrix<N, N> r1 =
              algebra::matrix::transpose(m1) * m2;
          typename A::template matrix<N, N> r2 =
              algebra::matrix::transposed_product<true, false>(m1, m2);

          for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
              ASSERT_NEAR(algebra::getter::element(r1, i, j),
                          algebra::getter::element(r2, i, j), this->m_epsilon);
            }
          }
        }

        // Only right transposed
        {
          typename A::template matrix<N, N> r1 =
              m1 * algebra::matrix::transpose(m2);
          typename A::template matrix<N, N> r2 =
              algebra::matrix::transposed_product<false, true>(m1, m2);

          for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
              ASSERT_NEAR(algebra::getter::element(r1, i, j),
                          algebra::getter::element(r2, i, j), this->m_epsilon);
            }
          }
        }

        // Both transposed
        {
          typename A::template matrix<N, N> r1 =
              algebra::matrix::transpose(m1) * algebra::matrix::transpose(m2);
          typename A::template matrix<N, N> r2 =
              algebra::matrix::transposed_product<true, true>(m1, m2);

          for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
              ASSERT_NEAR(algebra::getter::element(r1, i, j),
                          algebra::getter::element(r2, i, j), this->m_epsilon);
            }
          }
        }
      }
    }

    this->template test_matrix_ops_any_matrix<A, N, N>();
  }

  template <typename A, std::size_t M, std::size_t N, std::size_t O>
  void test_matrix_ops_inhomogeneous_multipliable_matrices() {
    // Test NxM and MxO matrix multiplication
    {
      typename A::template matrix<N, M> m1;
      typename A::template matrix<M, O> m2;

      for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < M; ++j) {
          algebra::getter::element(m1, i, j) =
              static_cast<algebra::traits::scalar_t<decltype(m1)>>(i * N + j);
        }
      }

      for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t j = 0; j < O; ++j) {
          algebra::getter::element(m2, i, j) =
              static_cast<algebra::traits::scalar_t<decltype(m1)>>(i * M + j);
        }
      }

      {
        typename A::template matrix<N, O> r1 = m1 * m2;
        typename A::template matrix<N, O> r2 =
            algebra::matrix::transposed_product<false, false>(m1, m2);

        for (std::size_t i = 0; i < N; ++i) {
          for (std::size_t j = 0; j < O; ++j) {
            ASSERT_NEAR(algebra::getter::element(r1, i, j),
                        algebra::getter::element(r2, i, j), this->m_epsilon);
          }
        }
      }
    }

    // Test NxM and (OxM)^T matrix multiplication
    {
      typename A::template matrix<N, M> m1;
      typename A::template matrix<O, M> m2;

      for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < M; ++j) {
          algebra::getter::element(m1, i, j) =
              static_cast<algebra::traits::scalar_t<decltype(m1)>>(i * N + j);
        }
      }

      for (std::size_t i = 0; i < O; ++i) {
        for (std::size_t j = 0; j < M; ++j) {
          algebra::getter::element(m2, i, j) =
              static_cast<algebra::traits::scalar_t<decltype(m1)>>(i * O + j);
        }
      }

      {
        typename A::template matrix<N, O> r1 =
            m1 * algebra::matrix::transpose(m2);
        typename A::template matrix<N, O> r2 =
            algebra::matrix::transposed_product<false, true>(m1, m2);

        for (std::size_t i = 0; i < N; ++i) {
          for (std::size_t j = 0; j < O; ++j) {
            ASSERT_NEAR(algebra::getter::element(r1, i, j),
                        algebra::getter::element(r2, i, j), this->m_epsilon);
          }
        }
      }
    }

    // Test (MxN)^T and MxO matrix multiplication
    {
      typename A::template matrix<M, N> m1;
      typename A::template matrix<M, O> m2;

      for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
          algebra::getter::element(m1, i, j) =
              static_cast<algebra::traits::scalar_t<decltype(m1)>>(i * M + j);
        }
      }

      for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t j = 0; j < O; ++j) {
          algebra::getter::element(m2, i, j) =
              static_cast<algebra::traits::scalar_t<decltype(m1)>>(i * M + j);
        }
      }

      {
        typename A::template matrix<N, O> r1 =
            algebra::matrix::transpose(m1) * m2;
        typename A::template matrix<N, O> r2 =
            algebra::matrix::transposed_product<true, false>(m1, m2);

        for (std::size_t i = 0; i < N; ++i) {
          for (std::size_t j = 0; j < O; ++j) {
            ASSERT_NEAR(algebra::getter::element(r1, i, j),
                        algebra::getter::element(r2, i, j), this->m_epsilon);
          }
        }
      }
    }

    // Test (MxN)^T and (OxM)^T matrix multiplication
    {
      typename A::template matrix<M, N> m1;
      typename A::template matrix<O, M> m2;

      for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
          algebra::getter::element(m1, i, j) =
              static_cast<algebra::traits::scalar_t<decltype(m1)>>(i * M + j);
        }
      }

      for (std::size_t i = 0; i < O; ++i) {
        for (std::size_t j = 0; j < M; ++j) {
          algebra::getter::element(m2, i, j) =
              static_cast<algebra::traits::scalar_t<decltype(m1)>>(i * O + j);
        }
      }

      {
        typename A::template matrix<N, O> r1 =
            algebra::matrix::transpose(m1) * algebra::matrix::transpose(m2);
        typename A::template matrix<N, O> r2 =
            algebra::matrix::transposed_product<true, true>(m1, m2);

        for (std::size_t i = 0; i < N; ++i) {
          for (std::size_t j = 0; j < O; ++j) {
            ASSERT_NEAR(algebra::getter::element(r1, i, j),
                        algebra::getter::element(r2, i, j), this->m_epsilon);
          }
        }
      }
    }
  }
};

TYPED_TEST_SUITE_P(test_host_basics_matrix);

/// Test case class, to be specialised for the different plugins - transforms
template <typename T>
class test_host_basics_transform : public testing::Test, public test_base<T> {};
TYPED_TEST_SUITE_P(test_host_basics_transform);

// This defines the local frame test suite
TYPED_TEST_P(test_host_basics_vector, local_vectors) {
  // Print the linear algebra types of this backend
  using algebra::operator<<;

  // Construction
  typename TypeParam::point2 vA{0.f, 1.f};
  ASSERT_EQ(vA[0], 0.f);
  ASSERT_EQ(vA[1], 1.f);

  // Test printing
  std::cout << vA << std::endl;

  // Test comparison
  constexpr auto epsilon{
      std::numeric_limits<typename TypeParam::scalar>::epsilon()};
  ASSERT_TRUE(algebra::approx_equal(vA, vA));
  ASSERT_TRUE(algebra::approx_equal(vA, vA, epsilon));

  typename TypeParam::scalar rel_err{1.f + 10.f * epsilon};
  typename TypeParam::point2 vA_err = rel_err * vA;
  ASSERT_TRUE(algebra::approx_equal(vA, vA_err, 11.f * epsilon));
  ASSERT_FALSE(algebra::approx_equal(vA, vA_err, 9.f * epsilon));

  rel_err = 1.f + 17.f * epsilon;
  vA_err = rel_err * vA;
  ASSERT_TRUE(algebra::approx_equal(vA, vA_err, 18.f * epsilon));
  ASSERT_FALSE(algebra::approx_equal(vA, vA_err, 16.f * epsilon));
  // Cast to (different) precision
  const auto vA_cast_f = algebra::cast_to<float>(vA);

  for (typename TypeParam::size_type i = 0; i < 2; ++i) {
    auto elem_i = vA_cast_f[i];

    static_assert(std::same_as<decltype(elem_i), float>);
    ASSERT_FLOAT_EQ(elem_i, static_cast<float>(vA[i]));
  }

  const auto vA_cast_d = algebra::cast_to<double>(vA);

  for (typename TypeParam::size_type i = 0; i < 2; ++i) {
    auto elem_i = vA_cast_d[i];

    static_assert(std::same_as<decltype(elem_i), double>);
    ASSERT_DOUBLE_EQ(elem_i, static_cast<double>(vA[i]));
  }

  const auto vA_cast_i = algebra::cast_to<int>(vA);

  for (typename TypeParam::size_type i = 0; i < 2; ++i) {
    auto elem_i = vA_cast_i[i];

    static_assert(std::same_as<decltype(elem_i), int>);
    ASSERT_EQ(elem_i, static_cast<int>(vA[i]));
  }

  // Assignment
  typename TypeParam::point2 vB = vA;
  ASSERT_EQ(vB[0], 0.f);
  ASSERT_EQ(vB[1], 1.f);

  // Addition
  typename TypeParam::point2 vC = vA + vB;
  ASSERT_EQ(vC[0], 0.f);
  ASSERT_EQ(vC[1], 2.f);

  // Multiplication by scalar
  typename TypeParam::point2 vC2 = vC * 2.f;
  ASSERT_EQ(vC2[0], 0.f);
  ASSERT_EQ(vC2[1], 4.f);

  // Cast operations to phi, theta, eta, perp
  typename TypeParam::vector2 vD{1.f, 1.f};
  typename TypeParam::scalar phi = algebra::vector::phi(vD);
  ASSERT_NEAR(phi, M_PI_4, this->m_epsilon);

  typename TypeParam::scalar perp = algebra::vector::perp(vD);
  ASSERT_NEAR(perp, std::sqrt(2.), this->m_epsilon);

  typename TypeParam::scalar norm = algebra::vector::norm(vD);
  ASSERT_NEAR(norm, std::sqrt(2.), this->m_epsilon);

  typename TypeParam::vector2 vDnorm = algebra::vector::normalize(vD);
  ASSERT_NEAR(vDnorm[0], 1. / std::sqrt(2.), this->m_epsilon);
  ASSERT_NEAR(vDnorm[1], 1. / std::sqrt(2.), this->m_epsilon);

  // Test comparison
  typename TypeParam::point2 vE{100.f, 462809.f};

  ASSERT_TRUE(algebra::approx_equal(vE, vE));
  ASSERT_TRUE(algebra::approx_equal(vE, vE, epsilon));

  rel_err = 1.f + 10.f * epsilon;
  typename TypeParam::point2 vE_err = rel_err * vE;
  ASSERT_TRUE(algebra::approx_equal(vE, vE_err, 11.f * epsilon));
  ASSERT_FALSE(algebra::approx_equal(vE, vE_err, 9.f * epsilon));

  typename TypeParam::point2 vE_abs_err{100.00001f, 462809.05f};
  ASSERT_TRUE(algebra::approx_equal(vE, vE_abs_err, 0.00001f));
  ASSERT_FALSE(algebra::approx_equal(vE, vE_abs_err, epsilon));

  rel_err = 1.f + 17.f * epsilon;
  vE_err = rel_err * vE;
  ASSERT_TRUE(algebra::approx_equal(vE, vE_err, 18.f * epsilon));
  ASSERT_FALSE(algebra::approx_equal(vE, vE_err, 16.f * epsilon));
}

// This defines the vector3 test suite
TYPED_TEST_P(test_host_basics_vector, vector3) {
  // Print the linear algebra types of this backend
  using algebra::operator<<;

  // Test concepts
  static_assert(algebra::concepts::scalar<typename TypeParam::scalar>);
  static_assert(!algebra::concepts::vector<typename TypeParam::scalar>);

  static_assert(!algebra::concepts::scalar<typename TypeParam::vector3>);
  static_assert(algebra::concepts::vector<typename TypeParam::vector3>);
  static_assert(algebra::concepts::vector3D<typename TypeParam::vector3>);
  static_assert(!algebra::concepts::vector2D<typename TypeParam::vector3>);

  static_assert(!algebra::concepts::scalar<typename TypeParam::vector2>);
  static_assert(algebra::concepts::vector<typename TypeParam::vector2>);
  static_assert(algebra::concepts::vector2D<typename TypeParam::vector2>);
  static_assert(!algebra::concepts::vector3D<typename TypeParam::vector2>);

  // Construction
  typename TypeParam::vector3 vA{0.f, 1.f, 2.f};
  ASSERT_EQ(vA[0], 0.f);
  ASSERT_EQ(vA[1], 1.f);
  ASSERT_EQ(vA[2], 2.f);

  // Test printing
  std::cout << vA << std::endl;

  // Cast to (different) precision
  const auto vA_cast_f = algebra::cast_to<float>(vA);

  for (typename TypeParam::size_type i = 0; i < 3; ++i) {
    auto elem_i = vA_cast_f[i];

    static_assert(std::same_as<decltype(elem_i), float>);
    ASSERT_FLOAT_EQ(elem_i, static_cast<float>(vA[i]));
  }

  const auto vA_cast_d = algebra::cast_to<double>(vA);

  for (typename TypeParam::size_type i = 0; i < 3; ++i) {
    auto elem_i = vA_cast_d[i];

    static_assert(std::same_as<decltype(elem_i), double>);
    ASSERT_DOUBLE_EQ(elem_i, static_cast<double>(vA[i]));
  }

  const auto vA_cast_i = algebra::cast_to<int>(vA);

  for (typename TypeParam::size_type i = 0; i < 3; ++i) {
    auto elem_i = vA_cast_i[i];

    static_assert(std::same_as<decltype(elem_i), int>);
    ASSERT_EQ(elem_i, static_cast<int>(vA[i]));
  }

  // Assignment
  typename TypeParam::vector3 vB = vA;
  ASSERT_EQ(vB[0], 0.f);
  ASSERT_EQ(vB[1], 1.f);
  ASSERT_EQ(vB[2], 2.f);

  // Addition
  typename TypeParam::vector3 vC = vA + vB;
  ASSERT_EQ(vC[0], 0.f);
  ASSERT_EQ(vC[1], 2.f);
  ASSERT_EQ(vC[2], 4.f);

  // Multiplication by scalar
  typename TypeParam::vector3 vC2 = vC * 2.0f;
  ASSERT_EQ(vC2[0], 0.f);
  ASSERT_EQ(vC2[1], 4.f);
  ASSERT_EQ(vC2[2], 8.f);

  // Cast operations to phi, theta, eta, perp
  typename TypeParam::vector3 vD{1.f, 1.f, 1.f};
  typename TypeParam::scalar phi = algebra::vector::phi(vD);
  ASSERT_NEAR(phi, M_PI_4, this->m_epsilon);

  typename TypeParam::scalar theta = algebra::vector::theta(vD);
  ASSERT_NEAR(theta, std::atan2(std::sqrt(2.), 1.), this->m_epsilon);

  typename TypeParam::scalar eta = algebra::vector::eta(vD);
  ASSERT_NEAR(eta, 0.65847891569137573, this->m_isclose);

  typename TypeParam::scalar perp = algebra::vector::perp(vD);
  ASSERT_NEAR(perp, std::sqrt(2.), this->m_epsilon);

  typename TypeParam::scalar norm = algebra::vector::norm(vD);
  ASSERT_NEAR(norm, std::sqrt(3.), this->m_epsilon);
}

// This defines the vector operation test suite
TYPED_TEST_P(test_host_basics_vector, getter) {

  typename TypeParam::vector3 v3{1.f, 1.f, 1.f};

  // Normalization
  typename TypeParam::vector3 v3n = algebra::vector::normalize(v3);
  ASSERT_NEAR(v3n[0], 1. / std::sqrt(3.), this->m_epsilon);
  ASSERT_NEAR(v3n[1], 1. / std::sqrt(3.), this->m_epsilon);
  ASSERT_NEAR(v3n[2], 1. / std::sqrt(3.), this->m_epsilon);

  // Cross product
  typename TypeParam::vector3 z =
      algebra::vector::normalize(typename TypeParam::vector3{3.f, 2.f, 1.f});
  typename TypeParam::vector3 x =
      algebra::vector::normalize(typename TypeParam::vector3{2.f, -3.f, 0.f});
  typename TypeParam::vector3 y =
      algebra::vector::cross(z, algebra::vector::normalize(x));

  // Check with dot product
  ASSERT_NEAR(algebra::vector::dot(x, y), 0.f, this->m_epsilon);
  ASSERT_NEAR(algebra::vector::dot(y, z), 0.f, this->m_epsilon);
  ASSERT_NEAR(algebra::vector::dot(z, x), 0.f, this->m_epsilon);
}

TYPED_TEST_P(test_host_basics_matrix, matrix_2x3) {
  static constexpr typename TypeParam::size_type ROWS = 2;
  static constexpr typename TypeParam::size_type COLS = 3;

  using matrix_2x3_t = typename TypeParam::template matrix<2, 3>;

  // Test type traits
  static_assert(std::is_same_v<algebra::traits::index_t<matrix_2x3_t>,
                               typename TypeParam::size_type>);
  static_assert(std::is_same_v<algebra::traits::value_t<matrix_2x3_t>,
                               typename TypeParam::scalar>);
  static_assert(std::is_same_v<algebra::traits::scalar_t<matrix_2x3_t>,
                               typename TypeParam::scalar>);
  static_assert(std::is_same_v<algebra::traits::vector_t<matrix_2x3_t>,
                               typename TypeParam::vector2>);

  static_assert(algebra::traits::rows<matrix_2x3_t> == 2);
  static_assert(algebra::traits::columns<matrix_2x3_t> == 3);
  static_assert(algebra::traits::rank<matrix_2x3_t> == 2);
  static_assert(algebra::traits::size<matrix_2x3_t> == 6);
  static_assert(!algebra::traits::is_square<matrix_2x3_t>);
  static_assert(
      algebra::traits::is_square<typename TypeParam::template matrix<2, 2>>);
  static_assert(
      algebra::traits::is_square<typename TypeParam::template matrix<3, 3>>);

  // Test concepts
  static_assert(algebra::concepts::matrix<matrix_2x3_t>);
  static_assert(!algebra::concepts::scalar<matrix_2x3_t>);
  static_assert(!algebra::concepts::vector<matrix_2x3_t>);
  static_assert(!algebra::concepts::square_matrix<matrix_2x3_t>);

  static_assert(
      algebra::concepts::index<algebra::traits::index_t<matrix_2x3_t>>);
  static_assert(
      algebra::concepts::value<algebra::traits::value_t<matrix_2x3_t>>);
  static_assert(
      algebra::concepts::scalar<algebra::traits::scalar_t<matrix_2x3_t>>);
  static_assert(
      algebra::concepts::vector<algebra::traits::vector_t<matrix_2x3_t>>);

  // Test on matrix - vector operations
  typename TypeParam::vector3 vE{1.f, 2.f, 3.f};

  matrix_2x3_t m23;

  algebra::getter::element(m23, 0, 0) = 1.f;
  algebra::getter::element(m23, 0, 1) = 2.f;
  algebra::getter::element(m23, 0, 2) = 3.f;
  algebra::getter::element(m23, 1, 0) = 4.f;
  algebra::getter::element(m23, 1, 1) = 5.f;
  algebra::getter::element(m23, 1, 2) = 6.f;

  // Cast to (different) precision
  const auto m23_cast_f = algebra::cast_to<float>(m23);

  for (typename TypeParam::size_type j = 0; j < 3; ++j) {
    for (typename TypeParam::size_type i = 0; i < 2; ++i) {
      auto elem_i = algebra::getter::element(m23_cast_f, i, j);

      static_assert(std::same_as<decltype(elem_i), float>);
      ASSERT_FLOAT_EQ(elem_i,
                      static_cast<float>(algebra::getter::element(m23, i, j)));
    }
  }

  const auto m23_cast_d = algebra::cast_to<double>(m23);

  for (typename TypeParam::size_type j = 0; j < 3; ++j) {
    for (typename TypeParam::size_type i = 0; i < 2; ++i) {
      auto elem_i = algebra::getter::element(m23_cast_d, i, j);

      static_assert(std::same_as<decltype(elem_i), double>);
      ASSERT_DOUBLE_EQ(
          elem_i, static_cast<double>(algebra::getter::element(m23, i, j)));
    }
  }

  const auto m23_cast_i = algebra::cast_to<int>(m23);

  for (typename TypeParam::size_type j = 0; j < 3; ++j) {
    for (typename TypeParam::size_type i = 0; i < 2; ++i) {
      auto elem_i = algebra::getter::element(m23_cast_i, i, j);

      static_assert(std::same_as<decltype(elem_i), int>);
      ASSERT_EQ(elem_i, static_cast<int>(algebra::getter::element(m23, i, j)));
    }
  }

  typename TypeParam::vector2 v2 = m23 * vE;

  ASSERT_NEAR(v2[0], 14, this->m_epsilon);
  ASSERT_NEAR(v2[1], 32, this->m_epsilon);

  this->template test_matrix_ops_any_matrix<TypeParam, ROWS, COLS>();
}

TYPED_TEST_P(test_host_basics_matrix, matrix_3x1) {
  // Print the linear algebra types of this backend
  using algebra::operator<<;

  static constexpr typename TypeParam::size_type ROWS = 3;
  static constexpr typename TypeParam::size_type COLS = 1;

  // Cross product on vector3 and matrix<3,1>
  typename TypeParam::template matrix<3, 1> vF;
  algebra::getter::element(vF, 0, 0) = 5.f;
  algebra::getter::element(vF, 1, 0) = 6.f;
  algebra::getter::element(vF, 2, 0) = 13.f;

  // Test printing
  std::cout << vF << std::endl;

  // Cast to (different) precision
  const auto vF_cast_f = algebra::cast_to<float>(vF);

  for (typename TypeParam::size_type i = 0; i < 3; ++i) {
    auto elem_i = algebra::getter::element(vF_cast_f, i, 0);

    static_assert(std::same_as<decltype(elem_i), float>);
    ASSERT_FLOAT_EQ(elem_i,
                    static_cast<float>(algebra::getter::element(vF, i, 0)));
  }

  const auto vF_cast_d = algebra::cast_to<double>(vF);

  for (typename TypeParam::size_type i = 0; i < 3; ++i) {
    auto elem_i = algebra::getter::element(vF_cast_d, i, 0);

    static_assert(std::same_as<decltype(elem_i), double>);
    ASSERT_DOUBLE_EQ(elem_i,
                     static_cast<double>(algebra::getter::element(vF, i, 0)));
  }

  const auto vF_cast_i = algebra::cast_to<int>(vF);

  for (typename TypeParam::size_type i = 0; i < 3; ++i) {
    auto elem_i = algebra::getter::element(vF_cast_i, i, 0);

    static_assert(std::same_as<decltype(elem_i), int>);
    ASSERT_EQ(elem_i, static_cast<int>(algebra::getter::element(vF, i, 0)));
  }

  typename TypeParam::vector3 vD{1.f, 1.f, 1.f};
  typename TypeParam::vector3 vG = algebra::vector::cross(vD, vF);
  ASSERT_NEAR(vG[0], 7.f, this->m_epsilon);
  ASSERT_NEAR(vG[1], -8.f, this->m_epsilon);
  ASSERT_NEAR(vG[2], 1.f, this->m_epsilon);

  // Dot product on vector3 and matrix<3,1>
  auto dot = algebra::vector::dot(vG, vF);
  ASSERT_NEAR(dot, 0.f, this->m_epsilon);

  this->template test_matrix_ops_any_matrix<TypeParam, ROWS, COLS>();
}

TYPED_TEST_P(test_host_basics_matrix, matrix_6x4) {
  // Print the linear algebra types of this backend
  using algebra::operator<<;

  // Create the matrix.
  static constexpr typename TypeParam::size_type ROWS = 6;
  static constexpr typename TypeParam::size_type COLS = 4;
  using matrix_6x4_t = typename TypeParam::template matrix<ROWS, COLS>;
  matrix_6x4_t m;

  // Test type traits
  static_assert(std::is_same_v<algebra::traits::index_t<matrix_6x4_t>,
                               typename TypeParam::size_type>);
  static_assert(std::is_same_v<algebra::traits::value_t<matrix_6x4_t>,
                               typename TypeParam::scalar>);
  static_assert(std::is_same_v<algebra::traits::scalar_t<matrix_6x4_t>,
                               typename TypeParam::scalar>);

  static_assert(algebra::traits::rows<matrix_6x4_t> == 6);
  static_assert(algebra::traits::columns<matrix_6x4_t> == 4);
  static_assert(algebra::traits::rank<matrix_6x4_t> == 4);
  static_assert(algebra::traits::size<matrix_6x4_t> == 24);
  static_assert(!algebra::traits::is_square<matrix_6x4_t>);
  static_assert(
      algebra::traits::is_square<typename TypeParam::template matrix<4, 4>>);
  static_assert(
      algebra::traits::is_square<typename TypeParam::template matrix<6, 6>>);

  // Test concepts
  static_assert(algebra::concepts::matrix<matrix_6x4_t>);
  static_assert(!algebra::concepts::scalar<matrix_6x4_t>);
  static_assert(!algebra::concepts::vector<matrix_6x4_t>);
  static_assert(!algebra::concepts::square_matrix<matrix_6x4_t>);

  // Fill it.
  for (typename TypeParam::size_type i = 0; i < ROWS; ++i) {
    for (typename TypeParam::size_type j = 0; j < COLS; ++j) {
      algebra::getter::element(m, i, j) =
          0.5f * static_cast<typename TypeParam::scalar>(i + j);
    }
  }

  // Check its content.
  const typename TypeParam::template matrix<ROWS, COLS>& m_const_ref = m;
  for (typename TypeParam::size_type i = 0; i < ROWS; ++i) {
    for (typename TypeParam::size_type j = 0; j < COLS; ++j) {
      const typename TypeParam::scalar ref =
          0.5f * static_cast<typename TypeParam::scalar>(i + j);
      ASSERT_NEAR(algebra::getter::element(m, i, j), ref, this->m_epsilon);
      ASSERT_NEAR(algebra::getter::element(m_const_ref, i, j), ref,
                  this->m_epsilon);
    }
  }

  // Test set_zero
  algebra::matrix::set_zero(m);
  for (typename TypeParam::size_type i = 0; i < ROWS; ++i) {
    for (typename TypeParam::size_type j = 0; j < COLS; ++j) {
      ASSERT_NEAR(algebra::getter::element(m, i, j), 0., this->m_epsilon);
    }
  }

  // Test set_identity
  algebra::matrix::set_identity(m);
  for (typename TypeParam::size_type i = 0; i < ROWS; ++i) {
    for (typename TypeParam::size_type j = 0; j < COLS; ++j) {
      if (i == j) {
        ASSERT_NEAR(algebra::getter::element(m, i, j), 1.f, this->m_epsilon);
      } else {
        ASSERT_NEAR(algebra::getter::element(m, i, j), 0.f, this->m_epsilon);
      }
    }
  }

  // Test block operations
  auto b13 = algebra::getter::block<1, 3>(m, 0, 0);
  ASSERT_NEAR(algebra::getter::element(b13, 0, 0), 1.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(b13, 0, 1), 0.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(b13, 0, 2), 0.f, this->m_epsilon);

  auto b13_tp = algebra::matrix::transpose(b13);
  ASSERT_NEAR(algebra::getter::element(b13_tp, 0, 0), 1.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(b13_tp, 1, 0), 0.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(b13_tp, 2, 0), 0.f, this->m_epsilon);

  auto b32 = algebra::getter::block<3, 2>(m, 2, 2);
  ASSERT_NEAR(algebra::getter::element(b32, 0, 0), 1.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(b32, 0, 1), 0.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(b32, 1, 0), 0.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(b32, 1, 1), 1.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(b32, 2, 0), 0.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(b32, 2, 1), 0.f, this->m_epsilon);

  algebra::getter::element(b32, 0, 0) = 4.f;
  algebra::getter::element(b32, 0, 1) = 3.f;
  algebra::getter::element(b32, 1, 0) = 12.f;
  algebra::getter::element(b32, 1, 1) = 13.f;
  algebra::getter::element(b32, 2, 0) = 5.f;
  algebra::getter::element(b32, 2, 1) = 6.f;

  algebra::getter::set_block(m, b32, 2, 2);
  ASSERT_NEAR(algebra::getter::element(m, 2, 2), 4.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(m, 2, 3), 3.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(m, 3, 2), 12.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(m, 3, 3), 13.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(m, 4, 2), 5.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(m, 4, 3), 6.f, this->m_epsilon);

  typename TypeParam::vector3 v = {10.f, 20.f, 30.f};
  algebra::getter::set_block(m, v, 0, 2);
  ASSERT_NEAR(algebra::getter::element(m, 0, 2), 10., this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(m, 1, 2), 20., this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(m, 2, 2), 30., this->m_epsilon);

  // Test printing
  std::cout << m << std::endl;

  this->template test_matrix_ops_any_matrix<TypeParam, ROWS, COLS>();
}

TYPED_TEST_P(test_host_basics_matrix, matrix_3x3) {
  static constexpr typename TypeParam::size_type N = 3;

  {
    typename TypeParam::vector3 v = {10.f, 20.f, 30.f};
    typename TypeParam::template matrix<3, 3> m33;
    algebra::getter::element(m33, 0, 0) = 1;
    algebra::getter::element(m33, 1, 0) = 2;
    algebra::getter::element(m33, 2, 0) = 3;
    algebra::getter::element(m33, 0, 1) = 5;
    algebra::getter::element(m33, 1, 1) = 6;
    algebra::getter::element(m33, 2, 1) = 7;
    algebra::getter::element(m33, 0, 2) = 9;
    algebra::getter::element(m33, 1, 2) = 10;
    algebra::getter::element(m33, 2, 2) = 11;

    const typename TypeParam::vector3 v2 = m33 * v;
    ASSERT_NEAR(v2[0], 380., this->m_epsilon);
    ASSERT_NEAR(v2[1], 440., this->m_epsilon);
    ASSERT_NEAR(v2[2], 500., this->m_epsilon);
  }

  {
    typename TypeParam::template matrix<3, 3> m33;
    algebra::getter::element(m33, 0, 0) = 1.f;
    algebra::getter::element(m33, 0, 1) = 5.f;
    algebra::getter::element(m33, 0, 2) = 7.f;
    algebra::getter::element(m33, 1, 0) = 3.f;
    algebra::getter::element(m33, 1, 1) = 5.f;
    algebra::getter::element(m33, 1, 2) = 6.f;
    algebra::getter::element(m33, 2, 0) = 2.f;
    algebra::getter::element(m33, 2, 1) = 8.f;
    algebra::getter::element(m33, 2, 2) = 9.f;

    // Test 3 X 3 matrix determinant
    auto m33_det = algebra::matrix::determinant(m33);
    ASSERT_NEAR(m33_det, 20.f, this->m_isclose);

    // Test 3 X 3 matrix inverse
    auto m33_inv = algebra::matrix::inverse(m33);
    ASSERT_NEAR(algebra::getter::element(m33_inv, 0, 0), -3.f / 20.f,
                this->m_isclose);
    ASSERT_NEAR(algebra::getter::element(m33_inv, 0, 1), 11.f / 20.f,
                this->m_isclose);
    ASSERT_NEAR(algebra::getter::element(m33_inv, 0, 2), -5.f / 20.f,
                this->m_isclose);
    ASSERT_NEAR(algebra::getter::element(m33_inv, 1, 0), -15.f / 20.f,
                this->m_isclose);
    ASSERT_NEAR(algebra::getter::element(m33_inv, 1, 1), -5.f / 20.f,
                this->m_isclose);
    ASSERT_NEAR(algebra::getter::element(m33_inv, 1, 2), 15.f / 20.f,
                this->m_isclose);
    ASSERT_NEAR(algebra::getter::element(m33_inv, 2, 0), 14.f / 20.f,
                this->m_isclose);
    ASSERT_NEAR(algebra::getter::element(m33_inv, 2, 1), 2.f / 20.f,
                this->m_isclose);
    ASSERT_NEAR(algebra::getter::element(m33_inv, 2, 2), -10.f / 20.f,
                this->m_isclose);
  }

  this->template test_matrix_ops_square_matrix<TypeParam, N>();
}

TYPED_TEST_P(test_host_basics_matrix, matrix_2x2) {
  static constexpr typename TypeParam::size_type N = 2;

  typename TypeParam::template matrix<2, 2> m22;
  algebra::getter::element(m22, 0, 0) = 4.f;
  algebra::getter::element(m22, 0, 1) = 3.f;
  algebra::getter::element(m22, 1, 0) = 12.f;
  algebra::getter::element(m22, 1, 1) = 13.f;

  // Test 2 X 2 matrix determinant
  auto m22_det = algebra::matrix::determinant(m22);
  ASSERT_NEAR(m22_det, 16.f, this->m_isclose);

  // Test 2 X 2 matrix inverse
  auto m22_inv = algebra::matrix::inverse(m22);
  ASSERT_NEAR(algebra::getter::element(m22_inv, 0, 0), 13.f / 16.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m22_inv, 0, 1), -3.f / 16.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m22_inv, 1, 0), -12.f / 16.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m22_inv, 1, 1), 4.f / 16.f,
              this->m_isclose);

  this->template test_matrix_ops_square_matrix<TypeParam, N>();
}

TYPED_TEST_P(test_host_basics_matrix, matrix_5x5) {

  // Test 5 X 5 matrix
  typename TypeParam::template matrix<5, 5> m55;
  algebra::getter::element(m55, 0, 0) = 1.f;
  algebra::getter::element(m55, 0, 1) = 3.f;
  algebra::getter::element(m55, 0, 2) = -9.f;
  algebra::getter::element(m55, 0, 3) = -5.f;
  algebra::getter::element(m55, 0, 4) = -2.f;

  algebra::getter::element(m55, 1, 0) = -6.f;
  algebra::getter::element(m55, 1, 1) = -3.f;
  algebra::getter::element(m55, 1, 2) = 1.f;
  algebra::getter::element(m55, 1, 3) = 0.f;
  algebra::getter::element(m55, 1, 4) = 2.f;

  algebra::getter::element(m55, 2, 0) = 12.f;
  algebra::getter::element(m55, 2, 1) = 7.f;
  algebra::getter::element(m55, 2, 2) = -9.f;
  algebra::getter::element(m55, 2, 3) = 11.f;
  algebra::getter::element(m55, 2, 4) = 2.f;

  algebra::getter::element(m55, 3, 0) = -3.f;
  algebra::getter::element(m55, 3, 1) = 4.f;
  algebra::getter::element(m55, 3, 2) = 5.f;
  algebra::getter::element(m55, 3, 3) = -6.f;
  algebra::getter::element(m55, 3, 4) = 7.f;

  algebra::getter::element(m55, 4, 0) = 9.f;
  algebra::getter::element(m55, 4, 1) = 6.f;
  algebra::getter::element(m55, 4, 2) = 3.f;
  algebra::getter::element(m55, 4, 3) = 0.f;
  algebra::getter::element(m55, 4, 4) = -3.f;

  auto m55_det = algebra::matrix::determinant(m55);
  ASSERT_NEAR((m55_det - 17334.f) / 17334.f, 0.f, this->m_isclose);

  auto m55_inv = algebra::matrix::inverse(m55);

  ASSERT_NEAR(algebra::getter::element(m55_inv, 0, 0), -2106.f / 17334.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m55_inv, 0, 1), -12312.f / 17334.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m55_inv, 0, 2), -486.f / 17334.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m55_inv, 0, 3), 864.f / 17334.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m55_inv, 0, 4), -5112.f / 17334.f,
              this->m_isclose);

  ASSERT_NEAR(algebra::getter::element(m55_inv, 1, 0), 2754.f / 17334.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m55_inv, 1, 1), 13878.f / 17334.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m55_inv, 1, 2), 1080.f / 17334.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m55_inv, 1, 3), -315.f / 17334.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m55_inv, 1, 4), 7401.f / 17334.f,
              this->m_isclose);

  ASSERT_NEAR(algebra::getter::element(m55_inv, 2, 0), -918.f / 17334.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m55_inv, 2, 1), 1152.f / 17334.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m55_inv, 2, 2), -360.f / 17334.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m55_inv, 2, 3), 105.f / 17334.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m55_inv, 2, 4), 1385.f / 17334.f,
              this->m_isclose);

  ASSERT_NEAR(algebra::getter::element(m55_inv, 3, 0), 108.f / 17334.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m55_inv, 3, 1), 7002.f / 17334.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m55_inv, 3, 2), 1062.f / 17334.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m55_inv, 3, 3), -1032.f / 17334.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m55_inv, 3, 4), 2896.f / 17334.f,
              this->m_isclose);

  ASSERT_NEAR(algebra::getter::element(m55_inv, 4, 0), -1728.f / 17334.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m55_inv, 4, 1), -8028.f / 17334.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m55_inv, 4, 2), 342.f / 17334.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m55_inv, 4, 3), 2067.f / 17334.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m55_inv, 4, 4), -4927.f / 17334.f,
              this->m_isclose);
}

TYPED_TEST_P(test_host_basics_matrix, matrix_6x6) {
  static constexpr typename TypeParam::size_type N = 6;

  // Test 6 X 6 big matrix determinant
  typename TypeParam::template matrix<6, 6> m66_big;
  algebra::getter::element(m66_big, 0, 0) = 1.f;
  algebra::getter::element(m66_big, 0, 1) = 0.f;
  algebra::getter::element(m66_big, 0, 2) = 3.f;
  algebra::getter::element(m66_big, 0, 3) = 0.f;
  algebra::getter::element(m66_big, 0, 4) = 0.f;
  algebra::getter::element(m66_big, 0, 5) = 0.f;

  algebra::getter::element(m66_big, 1, 0) = 0.f;
  algebra::getter::element(m66_big, 1, 1) = -2.f;
  algebra::getter::element(m66_big, 1, 2) = 4.f;
  algebra::getter::element(m66_big, 1, 3) = 0.f;
  algebra::getter::element(m66_big, 1, 4) = 5.f;
  algebra::getter::element(m66_big, 1, 5) = 0.f;

  algebra::getter::element(m66_big, 2, 0) = 0.f;
  algebra::getter::element(m66_big, 2, 1) = 0.f;
  algebra::getter::element(m66_big, 2, 2) = 3.f;
  algebra::getter::element(m66_big, 2, 3) = 0.f;
  algebra::getter::element(m66_big, 2, 4) = 0.f;
  algebra::getter::element(m66_big, 2, 5) = 0.f;

  algebra::getter::element(m66_big, 3, 0) = 0.f;
  algebra::getter::element(m66_big, 3, 1) = 0.f;
  algebra::getter::element(m66_big, 3, 2) = 0.f;
  algebra::getter::element(m66_big, 3, 3) = 4.f;
  algebra::getter::element(m66_big, 3, 4) = 0.f;
  algebra::getter::element(m66_big, 3, 5) = 0.f;

  algebra::getter::element(m66_big, 4, 0) = 0.f;
  algebra::getter::element(m66_big, 4, 1) = 0.f;
  algebra::getter::element(m66_big, 4, 2) = 0.f;
  algebra::getter::element(m66_big, 4, 3) = 0.f;
  algebra::getter::element(m66_big, 4, 4) = 9.f;
  algebra::getter::element(m66_big, 4, 5) = 0.f;

  algebra::getter::element(m66_big, 5, 0) = -1.f;
  algebra::getter::element(m66_big, 5, 1) = -1.f;
  algebra::getter::element(m66_big, 5, 2) = -1.f;
  algebra::getter::element(m66_big, 5, 3) = -1.f;
  algebra::getter::element(m66_big, 5, 4) = -1.f;
  algebra::getter::element(m66_big, 5, 5) = -1.f;

  auto m66_big_det = algebra::matrix::determinant(m66_big);
  ASSERT_NEAR(m66_big_det, 216.f, 2.f * this->m_isclose);

  auto m66_big_inv = algebra::matrix::inverse(m66_big);

  // Test comparison
  constexpr auto epsilon{
      std::numeric_limits<typename TypeParam::scalar>::epsilon()};

  ASSERT_TRUE(algebra::approx_equal(m66_big, m66_big));
  ASSERT_TRUE(algebra::approx_equal(m66_big, m66_big, epsilon));
  ASSERT_TRUE(algebra::approx_equal(m66_big_inv, m66_big_inv));
  ASSERT_TRUE(algebra::approx_equal(m66_big_inv, m66_big_inv, epsilon));

  typename TypeParam::scalar rel_err{1.f + 10.f * epsilon};
  typename TypeParam::template matrix<6, 6> m66_big_err = rel_err * m66_big;
  ASSERT_TRUE(algebra::approx_equal(m66_big, m66_big_err, 11.f * epsilon));
  ASSERT_FALSE(algebra::approx_equal(m66_big, m66_big_err, 9.f * epsilon));

  rel_err = 1.f + 17.f * epsilon;
  m66_big_err = rel_err * m66_big;
  ASSERT_TRUE(algebra::approx_equal(m66_big, m66_big_err, 18.f * epsilon));
  ASSERT_FALSE(algebra::approx_equal(m66_big, m66_big_err, 16.f * epsilon));

  typename TypeParam::template matrix<6, 6> prod66 = m66_big_inv * m66_big;
  auto I66 = algebra::matrix::identity<decltype(m66_big)>();
  // Set the max error for Fastor and SMatrix plugins
  ASSERT_TRUE(algebra::approx_equal(prod66, I66, epsilon, 2.f * epsilon));

  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 0, 0), 36.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 0, 1), 0.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 0, 2), -36.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 0, 3), 0.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 0, 4), 0.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 0, 5), 0.f / 36.f,
              this->m_isclose);

  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 1, 0), 0.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 1, 1), -18.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 1, 2), 24.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 1, 3), 0.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 1, 4), 10.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 1, 5), 0.f / 36.f,
              this->m_isclose);

  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 2, 0), 0.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 2, 1), 0.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 2, 2), 12.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 2, 3), 0.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 2, 4), 0.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 2, 5), 0.f / 36.f,
              this->m_isclose);

  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 3, 0), 0.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 3, 1), 0.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 3, 2), 0.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 3, 3), 9.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 3, 4), 0.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 3, 5), 0.f / 36.f,
              this->m_isclose);

  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 4, 0), 0.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 4, 1), 0.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 4, 2), 0.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 4, 3), 0.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 4, 4), 4.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 4, 5), 0.f / 36.f,
              this->m_isclose);

  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 5, 0), -36.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 5, 1), 18.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 5, 2), 0.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 5, 3), -9.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 5, 4), -14.f / 36.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m66_big_inv, 5, 5), -36.f / 36.f,
              this->m_isclose);

  // Cast to (different) precision
  const auto m66_cast_f = algebra::cast_to<float>(m66_big_inv);

  for (typename TypeParam::size_type j = 0; j < 6; ++j) {
    for (typename TypeParam::size_type i = 0; i < 6; ++i) {
      auto elem_i = algebra::getter::element(m66_cast_f, i, j);

      static_assert(std::same_as<decltype(elem_i), float>);
      ASSERT_FLOAT_EQ(elem_i, static_cast<float>(
                                  algebra::getter::element(m66_big_inv, i, j)));
    }
  }

  const auto m66_cast_d = algebra::cast_to<double>(m66_big_inv);

  for (typename TypeParam::size_type j = 0; j < 6; ++j) {
    for (typename TypeParam::size_type i = 0; i < 6; ++i) {
      auto elem_i = algebra::getter::element(m66_cast_d, i, j);

      static_assert(std::same_as<decltype(elem_i), double>);
      ASSERT_DOUBLE_EQ(elem_i, static_cast<double>(algebra::getter::element(
                                   m66_big_inv, i, j)));
    }
  }

  const auto m66_cast_i = algebra::cast_to<int>(m66_big_inv);

  for (typename TypeParam::size_type j = 0; j < 6; ++j) {
    for (typename TypeParam::size_type i = 0; i < 6; ++i) {
      auto elem_i = algebra::getter::element(m66_cast_i, i, j);

      static_assert(std::same_as<decltype(elem_i), int>);
      ASSERT_EQ(elem_i,
                static_cast<int>(algebra::getter::element(m66_big_inv, i, j)));
    }
  }

  // Test 6 X 6 small matrix determinant
  typename TypeParam::template matrix<6, 6> m66_small;

  algebra::getter::element(m66_small, 0, 0) =
      static_cast<typename TypeParam::scalar>(10.792386);
  algebra::getter::element(m66_small, 0, 1) =
      static_cast<typename TypeParam::scalar>(0.216181);
  algebra::getter::element(m66_small, 0, 2) =
      static_cast<typename TypeParam::scalar>(0.057650);
  algebra::getter::element(m66_small, 0, 3) =
      static_cast<typename TypeParam::scalar>(-0.002764);
  algebra::getter::element(m66_small, 0, 4) =
      static_cast<typename TypeParam::scalar>(0.000001);
  algebra::getter::element(m66_small, 0, 5) =
      static_cast<typename TypeParam::scalar>(0);

  algebra::getter::element(m66_small, 1, 0) =
      static_cast<typename TypeParam::scalar>(43.909368);
  algebra::getter::element(m66_small, 1, 1) =
      static_cast<typename TypeParam::scalar>(10.372997);
  algebra::getter::element(m66_small, 1, 2) =
      static_cast<typename TypeParam::scalar>(0.231496);
  algebra::getter::element(m66_small, 1, 3) =
      static_cast<typename TypeParam::scalar>(-0.065972);
  algebra::getter::element(m66_small, 1, 4) =
      static_cast<typename TypeParam::scalar>(-0.000002);
  algebra::getter::element(m66_small, 1, 5) =
      static_cast<typename TypeParam::scalar>(0);

  algebra::getter::element(m66_small, 2, 0) =
      static_cast<typename TypeParam::scalar>(0.045474);
  algebra::getter::element(m66_small, 2, 1) =
      static_cast<typename TypeParam::scalar>(-0.001730);
  algebra::getter::element(m66_small, 2, 2) =
      static_cast<typename TypeParam::scalar>(0.000246);
  algebra::getter::element(m66_small, 2, 3) =
      static_cast<typename TypeParam::scalar>(0.000004);
  algebra::getter::element(m66_small, 2, 4) =
      static_cast<typename TypeParam::scalar>(0);
  algebra::getter::element(m66_small, 2, 5) =
      static_cast<typename TypeParam::scalar>(0);

  algebra::getter::element(m66_small, 3, 0) =
      static_cast<typename TypeParam::scalar>(-0.255134);
  algebra::getter::element(m66_small, 3, 1) =
      static_cast<typename TypeParam::scalar>(-0.059846);
  algebra::getter::element(m66_small, 3, 2) =
      static_cast<typename TypeParam::scalar>(-0.001345);
  algebra::getter::element(m66_small, 3, 3) =
      static_cast<typename TypeParam::scalar>(0.000383);
  algebra::getter::element(m66_small, 3, 4) =
      static_cast<typename TypeParam::scalar>(0);
  algebra::getter::element(m66_small, 3, 5) =
      static_cast<typename TypeParam::scalar>(0);

  algebra::getter::element(m66_small, 4, 0) =
      static_cast<typename TypeParam::scalar>(-0.001490);
  algebra::getter::element(m66_small, 4, 1) =
      static_cast<typename TypeParam::scalar>(-0.000057);
  algebra::getter::element(m66_small, 4, 2) =
      static_cast<typename TypeParam::scalar>(-0.000008);
  algebra::getter::element(m66_small, 4, 3) =
      static_cast<typename TypeParam::scalar>(0.000001);
  algebra::getter::element(m66_small, 4, 4) =
      static_cast<typename TypeParam::scalar>(0.000001);
  algebra::getter::element(m66_small, 4, 5) =
      static_cast<typename TypeParam::scalar>(0);

  algebra::getter::element(m66_small, 5, 0) =
      static_cast<typename TypeParam::scalar>(0);
  algebra::getter::element(m66_small, 5, 1) =
      static_cast<typename TypeParam::scalar>(0);
  algebra::getter::element(m66_small, 5, 2) =
      static_cast<typename TypeParam::scalar>(0);
  algebra::getter::element(m66_small, 5, 3) =
      static_cast<typename TypeParam::scalar>(0);
  algebra::getter::element(m66_small, 5, 4) =
      static_cast<typename TypeParam::scalar>(0);
  algebra::getter::element(m66_small, 5, 5) =
      static_cast<typename TypeParam::scalar>(89875.517874);

  auto m66_small_det = algebra::matrix::determinant(m66_small);
  ASSERT_NEAR((m66_small_det - 4.30636e-11f) / 4.30636e-11f, 0.f,
              2.f * this->m_isclose);

  this->template test_matrix_ops_square_matrix<TypeParam, N>();
}

TYPED_TEST_P(test_host_basics_matrix, matrix_7x4_4x12) {
  this->template test_matrix_ops_inhomogeneous_multipliable_matrices<
      TypeParam, 7, 4, 12>();
}

TYPED_TEST_P(test_host_basics_matrix, matrix_17x9_9x4) {
  this->template test_matrix_ops_inhomogeneous_multipliable_matrices<
      TypeParam, 17, 9, 4>();
}

TYPED_TEST_P(test_host_basics_matrix, matrix_5x2_2x3) {
  this->template test_matrix_ops_inhomogeneous_multipliable_matrices<TypeParam,
                                                                     5, 2, 3>();
}

TYPED_TEST_P(test_host_basics_matrix, matrix_small_mixed) {

  typename TypeParam::template matrix<2, 2> m22;
  algebra::getter::element(m22, 0, 0) = 4.f;
  algebra::getter::element(m22, 0, 1) = 3.f;
  algebra::getter::element(m22, 1, 0) = 12.f;
  algebra::getter::element(m22, 1, 1) = 13.f;

  // Test 2 X 2 matrix inverse
  auto m22_inv = algebra::matrix::inverse(m22);
  ASSERT_NEAR(algebra::getter::element(m22_inv, 0, 0), 13.f / 16.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m22_inv, 0, 1), -3.f / 16.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m22_inv, 1, 0), -12.f / 16.f,
              this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m22_inv, 1, 1), 4.f / 16.f,
              this->m_isclose);

  typename TypeParam::template matrix<3, 3> m33;
  algebra::getter::element(m33, 0, 0) = 1.f;
  algebra::getter::element(m33, 0, 1) = 5.f;
  algebra::getter::element(m33, 0, 2) = 7.f;
  algebra::getter::element(m33, 1, 0) = 3.f;
  algebra::getter::element(m33, 1, 1) = 5.f;
  algebra::getter::element(m33, 1, 2) = 6.f;
  algebra::getter::element(m33, 2, 0) = 2.f;
  algebra::getter::element(m33, 2, 1) = 8.f;
  algebra::getter::element(m33, 2, 2) = 9.f;

  auto m33_inv = algebra::matrix::inverse(m33);

  // Test Zero
  auto m23 = algebra::matrix::zero<typename TypeParam::template matrix<2, 3>>();
  algebra::getter::element(m23, 0, 0) += 2.f;
  algebra::getter::element(m23, 0, 1) += 3.f;
  algebra::getter::element(m23, 0, 2) += 4.f;
  algebra::getter::element(m23, 1, 0) += 5.f;
  algebra::getter::element(m23, 1, 1) += 6.f;
  algebra::getter::element(m23, 1, 2) += 7.f;

  // Test scalar X Matrix
  m23 = 2. * m23;
  ASSERT_NEAR(algebra::getter::element(m23, 0, 0), 4.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(m23, 0, 1), 6.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(m23, 0, 2), 8.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(m23, 1, 0), 10.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(m23, 1, 1), 12.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(m23, 1, 2), 14.f, this->m_epsilon);

  // Test Transpose
  auto m32 = algebra::matrix::transpose(m23);

  // Test Identity and (Matrix + Matrix)
  m32 = m32 +
        algebra::matrix::identity<typename TypeParam::template matrix<3, 2>>();

  // Test Matrix X scalar
  m32 = m32 * 2.f;

  ASSERT_NEAR(algebra::getter::element(m32, 0, 0), 10.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(m32, 0, 1), 20.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(m32, 1, 0), 12.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(m32, 1, 1), 26.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(m32, 2, 0), 16.f, this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(m32, 2, 1), 28.f, this->m_epsilon);

  // Test Matrix multiplication
  m22 = m22_inv * m23 * m33_inv * m32;

  ASSERT_NEAR(algebra::getter::element(m22, 0, 0), 6.225f, this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m22, 0, 1), 14.675f, this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m22, 1, 0), -3.3f, this->m_isclose);
  ASSERT_NEAR(algebra::getter::element(m22, 1, 1), -7.9f, this->m_isclose);
}

// clang-format off
#define TEST_HOST_BASICS_MATRIX_TESTS(...) \
  REGISTER_TYPED_TEST_SUITE_P(test_host_basics_matrix \
    , matrix_2x2 \
    , matrix_2x3 \
    , matrix_3x1 \
    , matrix_3x3 \
    , matrix_6x4 \
    , matrix_5x5 \
    , matrix_6x6 \
    , matrix_7x4_4x12 \
    , matrix_17x9_9x4 \
    , matrix_5x2_2x3 \
    , matrix_small_mixed \
    )
// clang-format on

// This defines the transform3 test suite
TYPED_TEST_P(test_host_basics_transform, transform3) {
  // Print the linear algebra types of this backend
  using algebra::operator<<;

  static_assert(algebra::concepts::transform3D<typename TypeParam::transform3>);

  // Preparation work
  typename TypeParam::vector3 z =
      algebra::vector::normalize(typename TypeParam::vector3{3.f, 2.f, 1.f});
  typename TypeParam::vector3 x =
      algebra::vector::normalize(typename TypeParam::vector3{2.f, -3.f, 0.f});
  typename TypeParam::vector3 y = algebra::vector::cross(z, x);
  typename TypeParam::point3 t = {2.f, 3.f, 4.f};

  // Test constructor from t, z, x
  typename TypeParam::transform3 trf1(t, z, x);
  ASSERT_TRUE(trf1 == trf1);
  typename TypeParam::transform3 trf2;
  trf2 = trf1;

  // Test printing
  std::cout << trf1 << std::endl;

  // Test comparison
  constexpr auto epsilon{
      std::numeric_limits<typename TypeParam::scalar>::epsilon()};

  ASSERT_TRUE(algebra::approx_equal(trf1, trf1));
  ASSERT_TRUE(algebra::approx_equal(trf1, trf1, epsilon));

  typename TypeParam::scalar rel_err{1.f + 10.f * epsilon};
  typename TypeParam::transform3 trf1_err(rel_err * t, rel_err * z,
                                          rel_err * x);
  ASSERT_TRUE(algebra::approx_equal(trf1, trf1_err, 200.f * epsilon));
  ASSERT_FALSE(algebra::approx_equal(trf1, trf1_err, 10.f * epsilon));
  // Cast to (different) precision
  const auto trf1_cast_f = algebra::cast_to<float>(trf1);
  const auto& mat_f = trf1_cast_f.matrix();
  const auto& mat_inv_f = trf1_cast_f.matrix_inverse();

  for (typename TypeParam::size_type j = 0; j < 3; ++j) {
    for (typename TypeParam::size_type i = 0; i < 2; ++i) {

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

  for (typename TypeParam::size_type j = 0; j < 3; ++j) {
    for (typename TypeParam::size_type i = 0; i < 2; ++i) {

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

  for (typename TypeParam::size_type j = 0; j < 3; ++j) {
    for (typename TypeParam::size_type i = 0; i < 2; ++i) {

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
  typename TypeParam::transform3 trfm(m44);

  // Make sure that algebra::getter:vector can be called.
  [[maybe_unused]] typename TypeParam::vector3
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
  std::array<typename TypeParam::scalar, 16> matray_helper = {
      1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f,
      0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 0.f};
  typename TypeParam::transform3::template array_type<16> matray;
  for (unsigned int i = 0u; i < 16u; ++i) {
    matray[i] = matray_helper[i];
  }
  typename TypeParam::transform3 trfma(matray);

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
TYPED_TEST_P(test_host_basics_transform, global_transformations) {

  // Preparation work
  typename TypeParam::vector3 z =
      algebra::vector::normalize(typename TypeParam::vector3{3.f, 2.f, 1.f});
  typename TypeParam::vector3 x =
      algebra::vector::normalize(typename TypeParam::vector3{2.f, -3.f, 0.f});
  [[maybe_unused]] typename TypeParam::vector3 y = algebra::vector::cross(z, x);
  typename TypeParam::point3 t = {2.f, 3.f, 4.f};
  typename TypeParam::transform3 trf(t, z, x);

  // Check that local origin translates into global translation
  typename TypeParam::point3 lzero = {0.f, 0.f, 0.f};
  typename TypeParam::point3 gzero = trf.point_to_global(lzero);
  ASSERT_NEAR(gzero[0], t[0], this->m_epsilon);
  ASSERT_NEAR(gzero[1], t[1], this->m_epsilon);
  ASSERT_NEAR(gzero[2], t[2], this->m_epsilon);

  // Check a round trip for point
  typename TypeParam::point3 lpoint = {3.f, 4.f, 5.f};
  typename TypeParam::point3 gpoint = trf.point_to_global(lpoint);
  typename TypeParam::point3 lpoint_r = trf.point_to_local(gpoint);
  ASSERT_NEAR(lpoint[0], lpoint_r[0], this->m_isclose);
  ASSERT_NEAR(lpoint[1], lpoint_r[1], this->m_isclose);
  ASSERT_NEAR(lpoint[2], lpoint_r[2], this->m_isclose);

  // Check a point versus vector transform
  // vector should not change if transformed by a pure translation
  typename TypeParam::transform3 ttrf(t);

  typename TypeParam::vector3 gvector = {1.f, 1.f, 1.f};
  typename TypeParam::vector3 lvector = ttrf.vector_to_local(gvector);
  ASSERT_NEAR(gvector[0], lvector[0], this->m_isclose);
  ASSERT_NEAR(gvector[1], lvector[1], this->m_isclose);
  ASSERT_NEAR(gvector[2], lvector[2], this->m_isclose);

  // Check a round trip for vector
  typename TypeParam::vector3 lvectorB = {7.f, 8.f, 9.f};
  typename TypeParam::vector3 gvectorB = trf.vector_to_local(lvectorB);
  typename TypeParam::vector3 lvectorC = trf.vector_to_global(gvectorB);
  ASSERT_NEAR(lvectorB[0], lvectorC[0], this->m_isclose);
  ASSERT_NEAR(lvectorB[1], lvectorC[1], this->m_isclose);
  ASSERT_NEAR(lvectorB[2], lvectorC[2], this->m_isclose);
}
