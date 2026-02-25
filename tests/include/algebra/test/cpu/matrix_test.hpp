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

/// Test case class, to be specialised for the different plugins - matrices
template <algebra::concepts::algebra A>
class matrix_test : public testing::Test, public test_base<A> {
 protected:
  template <std::size_t ROWS, std::size_t COLS>
  void matrix_test_ops_any_matrix() {
    // Test the set_product method.
    {
      algebra::get_matrix_t<A, ROWS, ROWS> m1;
      algebra::get_matrix_t<A, ROWS, COLS> m2;

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
        algebra::get_matrix_t<A, ROWS, COLS> r1 = m1 * m2;
        algebra::get_matrix_t<A, ROWS, COLS> r2;
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
      algebra::get_matrix_t<A, ROWS, ROWS> m1;
      algebra::get_matrix_t<A, COLS, ROWS> m2;

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
        algebra::get_matrix_t<A, ROWS, COLS> r1 =
            m1 * algebra::matrix::transpose(m2);
        algebra::get_matrix_t<A, ROWS, COLS> r2;
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
      algebra::get_matrix_t<A, ROWS, ROWS> m1;
      algebra::get_matrix_t<A, ROWS, COLS> m2;

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
        algebra::get_matrix_t<A, ROWS, COLS> r1 =
            algebra::matrix::transpose(m1) * m2;
        algebra::get_matrix_t<A, ROWS, COLS> r2;
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
          algebra::get_matrix_t<A, ROWS, COLS> r1 = m1 * m2;
          algebra::get_matrix_t<A, ROWS, COLS> r2 =
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

  template <std::size_t N>
  void matrix_test_ops_square_matrix() {
    {
      algebra::get_matrix_t<A, N, N> m1;
      algebra::get_matrix_t<A, N, N> m2;

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
        algebra::get_matrix_t<A, N, N> r1 = m1 * m2;
        algebra::get_matrix_t<A, N, N> r2;
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
        algebra::get_matrix_t<A, N, N> r1 = m1 * algebra::matrix::transpose(m2);
        algebra::get_matrix_t<A, N, N> r2;
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
        algebra::get_matrix_t<A, N, N> r1 = algebra::matrix::transpose(m1) * m2;
        algebra::get_matrix_t<A, N, N> r2;
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
        algebra::get_matrix_t<A, N, N> r1 = m1 * m2;
        algebra::get_matrix_t<A, N, N> r2 = m1;
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
        algebra::get_matrix_t<A, N, N> r1 = m1 * m2;
        algebra::get_matrix_t<A, N, N> r2 = m2;
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
        algebra::get_matrix_t<A, N, N> r1 = m1 * algebra::matrix::transpose(m2);
        algebra::get_matrix_t<A, N, N> r2 = m1;
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
        algebra::get_matrix_t<A, N, N> r1 = algebra::matrix::transpose(m1) * m2;
        algebra::get_matrix_t<A, N, N> r2 = m2;
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
          algebra::get_matrix_t<A, N, N> r1 =
              algebra::matrix::transpose(m1) * m2;
          algebra::get_matrix_t<A, N, N> r2 =
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
          algebra::get_matrix_t<A, N, N> r1 =
              m1 * algebra::matrix::transpose(m2);
          algebra::get_matrix_t<A, N, N> r2 =
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
          algebra::get_matrix_t<A, N, N> r1 =
              algebra::matrix::transpose(m1) * algebra::matrix::transpose(m2);
          algebra::get_matrix_t<A, N, N> r2 =
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

    this->template matrix_test_ops_any_matrix<N, N>();
  }

  template <std::size_t M, std::size_t N, std::size_t O>
  void matrix_test_ops_inhomogeneous_multipliable_matrices() {
    // Test NxM and MxO matrix multiplication
    {
      algebra::get_matrix_t<A, N, M> m1;
      algebra::get_matrix_t<A, M, O> m2;

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
        algebra::get_matrix_t<A, N, O> r1 = m1 * m2;
        algebra::get_matrix_t<A, N, O> r2 =
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
      algebra::get_matrix_t<A, N, M> m1;
      algebra::get_matrix_t<A, O, M> m2;

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
        algebra::get_matrix_t<A, N, O> r1 = m1 * algebra::matrix::transpose(m2);
        algebra::get_matrix_t<A, N, O> r2 =
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
      algebra::get_matrix_t<A, M, N> m1;
      algebra::get_matrix_t<A, M, O> m2;

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
        algebra::get_matrix_t<A, N, O> r1 = algebra::matrix::transpose(m1) * m2;
        algebra::get_matrix_t<A, N, O> r2 =
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
      algebra::get_matrix_t<A, M, N> m1;
      algebra::get_matrix_t<A, O, M> m2;

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
        algebra::get_matrix_t<A, N, O> r1 =
            algebra::matrix::transpose(m1) * algebra::matrix::transpose(m2);
        algebra::get_matrix_t<A, N, O> r2 =
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

TYPED_TEST_SUITE_P(matrix_test);

TYPED_TEST_P(matrix_test, matrix_2x3) {
  static constexpr algebra::get_size_t<TypeParam> ROWS = 2;
  static constexpr algebra::get_size_t<TypeParam> COLS = 3;

  using matrix_2x3_t = algebra::get_matrix_t<TypeParam, 2, 3>;

  // Test type traits
  static_assert(std::is_same_v<algebra::traits::index_t<matrix_2x3_t>,
                               algebra::get_size_t<TypeParam>>);
  static_assert(std::is_same_v<algebra::traits::value_t<matrix_2x3_t>,
                               algebra::get_scalar_t<TypeParam>>);
  static_assert(std::is_same_v<algebra::traits::scalar_t<matrix_2x3_t>,
                               algebra::get_scalar_t<TypeParam>>);
  static_assert(std::is_same_v<algebra::traits::vector_t<matrix_2x3_t>,
                               algebra::get_vector2D_t<TypeParam>>);

  static_assert(algebra::traits::rows<matrix_2x3_t> == 2);
  static_assert(algebra::traits::columns<matrix_2x3_t> == 3);
  static_assert(algebra::traits::rank<matrix_2x3_t> == 2);
  static_assert(algebra::traits::size<matrix_2x3_t> == 6);
  static_assert(!algebra::traits::is_square<matrix_2x3_t>);
  static_assert(
      algebra::traits::is_square<algebra::get_matrix_t<TypeParam, 2, 2>>);
  static_assert(
      algebra::traits::is_square<algebra::get_matrix_t<TypeParam, 3, 3>>);

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
  algebra::get_vector3D_t<TypeParam> vE{1.f, 2.f, 3.f};

  matrix_2x3_t m23;

  algebra::getter::element(m23, 0, 0) = 1.f;
  algebra::getter::element(m23, 0, 1) = 2.f;
  algebra::getter::element(m23, 0, 2) = 3.f;
  algebra::getter::element(m23, 1, 0) = 4.f;
  algebra::getter::element(m23, 1, 1) = 5.f;
  algebra::getter::element(m23, 1, 2) = 6.f;

  // Cast to (different) precision
  const auto m23_cast_f = algebra::cast_to<float>(m23);

  for (algebra::get_size_t<TypeParam> j = 0; j < 3; ++j) {
    for (algebra::get_size_t<TypeParam> i = 0; i < 2; ++i) {
      auto elem_i = algebra::getter::element(m23_cast_f, i, j);

      static_assert(std::same_as<decltype(elem_i), float>);
      ASSERT_FLOAT_EQ(elem_i,
                      static_cast<float>(algebra::getter::element(m23, i, j)));
    }
  }

  const auto m23_cast_d = algebra::cast_to<double>(m23);

  for (algebra::get_size_t<TypeParam> j = 0; j < 3; ++j) {
    for (algebra::get_size_t<TypeParam> i = 0; i < 2; ++i) {
      auto elem_i = algebra::getter::element(m23_cast_d, i, j);

      static_assert(std::same_as<decltype(elem_i), double>);
      ASSERT_DOUBLE_EQ(
          elem_i, static_cast<double>(algebra::getter::element(m23, i, j)));
    }
  }

  const auto m23_cast_i = algebra::cast_to<int>(m23);

  for (algebra::get_size_t<TypeParam> j = 0; j < 3; ++j) {
    for (algebra::get_size_t<TypeParam> i = 0; i < 2; ++i) {
      auto elem_i = algebra::getter::element(m23_cast_i, i, j);

      static_assert(std::same_as<decltype(elem_i), int>);
      ASSERT_EQ(elem_i, static_cast<int>(algebra::getter::element(m23, i, j)));
    }
  }

  algebra::get_vector2D_t<TypeParam> v2 = m23 * vE;

  ASSERT_NEAR(v2[0], 14, this->m_epsilon);
  ASSERT_NEAR(v2[1], 32, this->m_epsilon);

  this->template matrix_test_ops_any_matrix<ROWS, COLS>();
}

TYPED_TEST_P(matrix_test, matrix_3x1) {
  // Print the linear algebra types of this backend
  using algebra::operator<<;

  static constexpr algebra::get_size_t<TypeParam> ROWS = 3;
  static constexpr algebra::get_size_t<TypeParam> COLS = 1;

  // Cross product on vector3 and matrix<3,1>
  algebra::get_matrix_t<TypeParam, 3, 1> vF;
  algebra::getter::element(vF, 0, 0) = 5.f;
  algebra::getter::element(vF, 1, 0) = 6.f;
  algebra::getter::element(vF, 2, 0) = 13.f;

  // Test printing
  std::cout << vF << std::endl;

  // Cast to (different) precision
  const auto vF_cast_f = algebra::cast_to<float>(vF);

  for (algebra::get_size_t<TypeParam> i = 0; i < 3; ++i) {
    auto elem_i = algebra::getter::element(vF_cast_f, i, 0);

    static_assert(std::same_as<decltype(elem_i), float>);
    ASSERT_FLOAT_EQ(elem_i,
                    static_cast<float>(algebra::getter::element(vF, i, 0)));
  }

  const auto vF_cast_d = algebra::cast_to<double>(vF);

  for (algebra::get_size_t<TypeParam> i = 0; i < 3; ++i) {
    auto elem_i = algebra::getter::element(vF_cast_d, i, 0);

    static_assert(std::same_as<decltype(elem_i), double>);
    ASSERT_DOUBLE_EQ(elem_i,
                     static_cast<double>(algebra::getter::element(vF, i, 0)));
  }

  const auto vF_cast_i = algebra::cast_to<int>(vF);

  for (algebra::get_size_t<TypeParam> i = 0; i < 3; ++i) {
    auto elem_i = algebra::getter::element(vF_cast_i, i, 0);

    static_assert(std::same_as<decltype(elem_i), int>);
    ASSERT_EQ(elem_i, static_cast<int>(algebra::getter::element(vF, i, 0)));
  }

  algebra::get_vector3D_t<TypeParam> vD{1.f, 1.f, 1.f};
  algebra::get_vector3D_t<TypeParam> vG = algebra::vector::cross(vD, vF);
  ASSERT_NEAR(vG[0], 7.f, this->m_epsilon);
  ASSERT_NEAR(vG[1], -8.f, this->m_epsilon);
  ASSERT_NEAR(vG[2], 1.f, this->m_epsilon);

  // Dot product on vector3 and matrix<3,1>
  auto dot = algebra::vector::dot(vG, vF);
  ASSERT_NEAR(dot, 0.f, this->m_epsilon);

  this->template matrix_test_ops_any_matrix<ROWS, COLS>();
}

TYPED_TEST_P(matrix_test, matrix_6x4) {
  // Print the linear algebra types of this backend
  using algebra::operator<<;

  // Create the matrix.
  static constexpr algebra::get_size_t<TypeParam> ROWS = 6;
  static constexpr algebra::get_size_t<TypeParam> COLS = 4;
  using matrix_6x4_t = algebra::get_matrix_t<TypeParam, ROWS, COLS>;
  matrix_6x4_t m;

  // Test type traits
  static_assert(std::is_same_v<algebra::traits::index_t<matrix_6x4_t>,
                               algebra::get_size_t<TypeParam>>);
  static_assert(std::is_same_v<algebra::traits::value_t<matrix_6x4_t>,
                               algebra::get_scalar_t<TypeParam>>);
  static_assert(std::is_same_v<algebra::traits::scalar_t<matrix_6x4_t>,
                               algebra::get_scalar_t<TypeParam>>);

  static_assert(algebra::traits::rows<matrix_6x4_t> == 6);
  static_assert(algebra::traits::columns<matrix_6x4_t> == 4);
  static_assert(algebra::traits::rank<matrix_6x4_t> == 4);
  static_assert(algebra::traits::size<matrix_6x4_t> == 24);
  static_assert(!algebra::traits::is_square<matrix_6x4_t>);
  static_assert(
      algebra::traits::is_square<algebra::get_matrix_t<TypeParam, 4, 4>>);
  static_assert(
      algebra::traits::is_square<algebra::get_matrix_t<TypeParam, 6, 6>>);

  // Test concepts
  static_assert(algebra::concepts::matrix<matrix_6x4_t>);
  static_assert(!algebra::concepts::scalar<matrix_6x4_t>);
  static_assert(!algebra::concepts::vector<matrix_6x4_t>);
  static_assert(!algebra::concepts::square_matrix<matrix_6x4_t>);

  // Fill it.
  for (algebra::get_size_t<TypeParam> i = 0; i < ROWS; ++i) {
    for (algebra::get_size_t<TypeParam> j = 0; j < COLS; ++j) {
      algebra::getter::element(m, i, j) =
          0.5f * static_cast<algebra::get_scalar_t<TypeParam>>(i + j);
    }
  }

  // Check its content.
  const algebra::get_matrix_t<TypeParam, ROWS, COLS>& m_const_ref = m;
  for (algebra::get_size_t<TypeParam> i = 0; i < ROWS; ++i) {
    for (algebra::get_size_t<TypeParam> j = 0; j < COLS; ++j) {
      const algebra::get_scalar_t<TypeParam> ref =
          0.5f * static_cast<algebra::get_scalar_t<TypeParam>>(i + j);
      ASSERT_NEAR(algebra::getter::element(m, i, j), ref, this->m_epsilon);
      ASSERT_NEAR(algebra::getter::element(m_const_ref, i, j), ref,
                  this->m_epsilon);
    }
  }

  // Test set_zero
  algebra::matrix::set_zero(m);
  for (algebra::get_size_t<TypeParam> i = 0; i < ROWS; ++i) {
    for (algebra::get_size_t<TypeParam> j = 0; j < COLS; ++j) {
      ASSERT_NEAR(algebra::getter::element(m, i, j), 0., this->m_epsilon);
    }
  }

  // Test set_identity
  algebra::matrix::set_identity(m);
  for (algebra::get_size_t<TypeParam> i = 0; i < ROWS; ++i) {
    for (algebra::get_size_t<TypeParam> j = 0; j < COLS; ++j) {
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

  algebra::get_vector3D_t<TypeParam> v = {10.f, 20.f, 30.f};
  algebra::getter::set_block(m, v, 0, 2);
  ASSERT_NEAR(algebra::getter::element(m, 0, 2), 10., this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(m, 1, 2), 20., this->m_epsilon);
  ASSERT_NEAR(algebra::getter::element(m, 2, 2), 30., this->m_epsilon);

  // Test printing
  std::cout << m << std::endl;

  this->template matrix_test_ops_any_matrix<ROWS, COLS>();
}

TYPED_TEST_P(matrix_test, matrix_3x3) {
  static constexpr algebra::get_size_t<TypeParam> N = 3;

  {
    algebra::get_vector3D_t<TypeParam> v = {10.f, 20.f, 30.f};
    algebra::get_matrix_t<TypeParam, 3, 3> m33;
    algebra::getter::element(m33, 0, 0) = 1;
    algebra::getter::element(m33, 1, 0) = 2;
    algebra::getter::element(m33, 2, 0) = 3;
    algebra::getter::element(m33, 0, 1) = 5;
    algebra::getter::element(m33, 1, 1) = 6;
    algebra::getter::element(m33, 2, 1) = 7;
    algebra::getter::element(m33, 0, 2) = 9;
    algebra::getter::element(m33, 1, 2) = 10;
    algebra::getter::element(m33, 2, 2) = 11;

    const algebra::get_vector3D_t<TypeParam> v2 = m33 * v;
    ASSERT_NEAR(v2[0], 380., this->m_epsilon);
    ASSERT_NEAR(v2[1], 440., this->m_epsilon);
    ASSERT_NEAR(v2[2], 500., this->m_epsilon);
  }

  {
    algebra::get_matrix_t<TypeParam, 3, 3> m33;
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

  this->template matrix_test_ops_square_matrix<N>();
}

TYPED_TEST_P(matrix_test, matrix_2x2) {
  static constexpr algebra::get_size_t<TypeParam> N = 2;

  algebra::get_matrix_t<TypeParam, 2, 2> m22;
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

  this->template matrix_test_ops_square_matrix<N>();
}

TYPED_TEST_P(matrix_test, matrix_5x5) {

  // Test 5 X 5 matrix
  algebra::get_matrix_t<TypeParam, 5, 5> m55;
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

TYPED_TEST_P(matrix_test, matrix_6x6) {
  static constexpr algebra::get_size_t<TypeParam> N = 6;

  // Test 6 X 6 big matrix determinant
  algebra::get_matrix_t<TypeParam, 6, 6> m66_big;
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
      std::numeric_limits<algebra::get_scalar_t<TypeParam>>::epsilon()};

  ASSERT_TRUE(algebra::approx_equal(m66_big, m66_big));
  ASSERT_TRUE(algebra::approx_equal(m66_big, m66_big, epsilon));
  ASSERT_TRUE(algebra::approx_equal(m66_big_inv, m66_big_inv));
  ASSERT_TRUE(algebra::approx_equal(m66_big_inv, m66_big_inv, epsilon));

  algebra::get_scalar_t<TypeParam> rel_err{1.f + 10.f * epsilon};
  algebra::get_matrix_t<TypeParam, 6, 6> m66_big_err = rel_err * m66_big;
  ASSERT_TRUE(algebra::approx_equal(m66_big, m66_big_err, 11.f * epsilon));
  ASSERT_FALSE(algebra::approx_equal(m66_big, m66_big_err, 9.f * epsilon));

  rel_err = 1.f + 17.f * epsilon;
  m66_big_err = rel_err * m66_big;
  ASSERT_TRUE(algebra::approx_equal(m66_big, m66_big_err, 18.f * epsilon));
  ASSERT_FALSE(algebra::approx_equal(m66_big, m66_big_err, 16.f * epsilon));

  algebra::get_matrix_t<TypeParam, 6, 6> prod66 = m66_big_inv * m66_big;
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

  for (algebra::get_size_t<TypeParam> j = 0; j < 6; ++j) {
    for (algebra::get_size_t<TypeParam> i = 0; i < 6; ++i) {
      auto elem_i = algebra::getter::element(m66_cast_f, i, j);

      static_assert(std::same_as<decltype(elem_i), float>);
      ASSERT_FLOAT_EQ(elem_i, static_cast<float>(
                                  algebra::getter::element(m66_big_inv, i, j)));
    }
  }

  const auto m66_cast_d = algebra::cast_to<double>(m66_big_inv);

  for (algebra::get_size_t<TypeParam> j = 0; j < 6; ++j) {
    for (algebra::get_size_t<TypeParam> i = 0; i < 6; ++i) {
      auto elem_i = algebra::getter::element(m66_cast_d, i, j);

      static_assert(std::same_as<decltype(elem_i), double>);
      ASSERT_DOUBLE_EQ(elem_i, static_cast<double>(algebra::getter::element(
                                   m66_big_inv, i, j)));
    }
  }

  const auto m66_cast_i = algebra::cast_to<int>(m66_big_inv);

  for (algebra::get_size_t<TypeParam> j = 0; j < 6; ++j) {
    for (algebra::get_size_t<TypeParam> i = 0; i < 6; ++i) {
      auto elem_i = algebra::getter::element(m66_cast_i, i, j);

      static_assert(std::same_as<decltype(elem_i), int>);
      ASSERT_EQ(elem_i,
                static_cast<int>(algebra::getter::element(m66_big_inv, i, j)));
    }
  }

  // Test 6 X 6 small matrix determinant
  algebra::get_matrix_t<TypeParam, 6, 6> m66_small;

  algebra::getter::element(m66_small, 0, 0) =
      static_cast<algebra::get_scalar_t<TypeParam>>(10.792386);
  algebra::getter::element(m66_small, 0, 1) =
      static_cast<algebra::get_scalar_t<TypeParam>>(0.216181);
  algebra::getter::element(m66_small, 0, 2) =
      static_cast<algebra::get_scalar_t<TypeParam>>(0.057650);
  algebra::getter::element(m66_small, 0, 3) =
      static_cast<algebra::get_scalar_t<TypeParam>>(-0.002764);
  algebra::getter::element(m66_small, 0, 4) =
      static_cast<algebra::get_scalar_t<TypeParam>>(0.000001);
  algebra::getter::element(m66_small, 0, 5) =
      static_cast<algebra::get_scalar_t<TypeParam>>(0);

  algebra::getter::element(m66_small, 1, 0) =
      static_cast<algebra::get_scalar_t<TypeParam>>(43.909368);
  algebra::getter::element(m66_small, 1, 1) =
      static_cast<algebra::get_scalar_t<TypeParam>>(10.372997);
  algebra::getter::element(m66_small, 1, 2) =
      static_cast<algebra::get_scalar_t<TypeParam>>(0.231496);
  algebra::getter::element(m66_small, 1, 3) =
      static_cast<algebra::get_scalar_t<TypeParam>>(-0.065972);
  algebra::getter::element(m66_small, 1, 4) =
      static_cast<algebra::get_scalar_t<TypeParam>>(-0.000002);
  algebra::getter::element(m66_small, 1, 5) =
      static_cast<algebra::get_scalar_t<TypeParam>>(0);

  algebra::getter::element(m66_small, 2, 0) =
      static_cast<algebra::get_scalar_t<TypeParam>>(0.045474);
  algebra::getter::element(m66_small, 2, 1) =
      static_cast<algebra::get_scalar_t<TypeParam>>(-0.001730);
  algebra::getter::element(m66_small, 2, 2) =
      static_cast<algebra::get_scalar_t<TypeParam>>(0.000246);
  algebra::getter::element(m66_small, 2, 3) =
      static_cast<algebra::get_scalar_t<TypeParam>>(0.000004);
  algebra::getter::element(m66_small, 2, 4) =
      static_cast<algebra::get_scalar_t<TypeParam>>(0);
  algebra::getter::element(m66_small, 2, 5) =
      static_cast<algebra::get_scalar_t<TypeParam>>(0);

  algebra::getter::element(m66_small, 3, 0) =
      static_cast<algebra::get_scalar_t<TypeParam>>(-0.255134);
  algebra::getter::element(m66_small, 3, 1) =
      static_cast<algebra::get_scalar_t<TypeParam>>(-0.059846);
  algebra::getter::element(m66_small, 3, 2) =
      static_cast<algebra::get_scalar_t<TypeParam>>(-0.001345);
  algebra::getter::element(m66_small, 3, 3) =
      static_cast<algebra::get_scalar_t<TypeParam>>(0.000383);
  algebra::getter::element(m66_small, 3, 4) =
      static_cast<algebra::get_scalar_t<TypeParam>>(0);
  algebra::getter::element(m66_small, 3, 5) =
      static_cast<algebra::get_scalar_t<TypeParam>>(0);

  algebra::getter::element(m66_small, 4, 0) =
      static_cast<algebra::get_scalar_t<TypeParam>>(-0.001490);
  algebra::getter::element(m66_small, 4, 1) =
      static_cast<algebra::get_scalar_t<TypeParam>>(-0.000057);
  algebra::getter::element(m66_small, 4, 2) =
      static_cast<algebra::get_scalar_t<TypeParam>>(-0.000008);
  algebra::getter::element(m66_small, 4, 3) =
      static_cast<algebra::get_scalar_t<TypeParam>>(0.000001);
  algebra::getter::element(m66_small, 4, 4) =
      static_cast<algebra::get_scalar_t<TypeParam>>(0.000001);
  algebra::getter::element(m66_small, 4, 5) =
      static_cast<algebra::get_scalar_t<TypeParam>>(0);

  algebra::getter::element(m66_small, 5, 0) =
      static_cast<algebra::get_scalar_t<TypeParam>>(0);
  algebra::getter::element(m66_small, 5, 1) =
      static_cast<algebra::get_scalar_t<TypeParam>>(0);
  algebra::getter::element(m66_small, 5, 2) =
      static_cast<algebra::get_scalar_t<TypeParam>>(0);
  algebra::getter::element(m66_small, 5, 3) =
      static_cast<algebra::get_scalar_t<TypeParam>>(0);
  algebra::getter::element(m66_small, 5, 4) =
      static_cast<algebra::get_scalar_t<TypeParam>>(0);
  algebra::getter::element(m66_small, 5, 5) =
      static_cast<algebra::get_scalar_t<TypeParam>>(89875.517874);

  auto m66_small_det = algebra::matrix::determinant(m66_small);
  ASSERT_NEAR((m66_small_det - 4.30636e-11f) / 4.30636e-11f, 0.f,
              2.f * this->m_isclose);

  this->template matrix_test_ops_square_matrix<N>();
}

TYPED_TEST_P(matrix_test, matrix_7x4_4x12) {
  this->template matrix_test_ops_inhomogeneous_multipliable_matrices<7, 4,
                                                                     12>();
}

TYPED_TEST_P(matrix_test, matrix_17x9_9x4) {
  this->template matrix_test_ops_inhomogeneous_multipliable_matrices<17, 9,
                                                                     4>();
}

TYPED_TEST_P(matrix_test, matrix_5x2_2x3) {
  this->template matrix_test_ops_inhomogeneous_multipliable_matrices<5, 2, 3>();
}

TYPED_TEST_P(matrix_test, matrix_small_mixed) {

  algebra::get_matrix_t<TypeParam, 2, 2> m22;
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

  algebra::get_matrix_t<TypeParam, 3, 3> m33;
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
  auto m23 = algebra::matrix::zero<algebra::get_matrix_t<TypeParam, 2, 3>>();
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
  m32 =
      m32 + algebra::matrix::identity<algebra::get_matrix_t<TypeParam, 3, 2>>();

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
#define ALGEBRA_REGISTER_MATRIX_TESTS(...) \
  REGISTER_TYPED_TEST_SUITE_P(matrix_test \
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

}  // namespace algebra::test
