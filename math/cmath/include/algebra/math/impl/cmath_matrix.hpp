/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/algorithms/utils/algorithm_finder.hpp"
#include "algebra/qualifiers.hpp"

namespace algebra::cmath::matrix {

/// "Matrix actor", assuming a simple 2D matrix
template <typename size_type, template <typename, size_type> class array_t,
          template <typename, size_type, size_type> class matrix_t,
          typename scalar_t, class determinant_actor_t, class inverse_actor_t,
          class element_getter_t, class block_getter_t>
struct actor {

  /// Size type
  using size_ty = size_type;

  /// Scalar type
  using scalar_type = scalar_t;

  /// Function (object) used for accessing a matrix element
  using element_getter = element_getter_t;

  /// Function (object) used for accessing a matrix block
  using block_getter = block_getter_t;

  /// 2D matrix type
  template <size_type ROWS, size_type COLS>
  using matrix_type = matrix_t<scalar_t, ROWS, COLS>;

  /// vector type
  template <size_type N>
  using array_type = array_t<scalar_t, N>;

  /// Operator getting a reference to one element of a non-const matrix
  template <size_type ROWS, size_type COLS>
  ALGEBRA_HOST_DEVICE inline scalar_t &element(matrix_type<ROWS, COLS> &m,
                                               size_type row,
                                               size_type col) const {
    return element_getter()(m, row, col);
  }

  /// Operator getting one value of a const matrix
  template <size_type ROWS, size_type COLS>
  ALGEBRA_HOST_DEVICE inline scalar_t element(const matrix_type<ROWS, COLS> &m,
                                              size_type row,
                                              size_type col) const {
    return element_getter()(m, row, col);
  }

  /// Operator getting a block of a const matrix
  template <size_type ROWS, size_type COLS, class input_matrix_type>
  ALGEBRA_HOST_DEVICE matrix_type<ROWS, COLS> block(const input_matrix_type &m,
                                                    size_type row,
                                                    size_type col) const {
    return block_getter().template operator()<ROWS, COLS>(m, row, col);
  }

  /// Operator setting a block with a vector matrix
  template <size_type ROWS, size_type COLS, class input_matrix_type>
  ALGEBRA_HOST_DEVICE void set_block(input_matrix_type &m,
                                     const matrix_type<ROWS, COLS> &b,
                                     size_type row, size_type col) const {
    for (size_type j = 0; j < COLS; ++j) {
      for (size_type i = 0; i < ROWS; ++i) {
        element_getter()(m, i + row, j + col) = element_getter()(b, i, j);
      }
    }
  }

  /// Operator setting a block with a vector
  template <size_type ROWS, template <typename, size_type> class vector_t,
            class input_matrix_type>
  ALGEBRA_HOST_DEVICE void set_block(input_matrix_type &m,
                                     const vector_t<scalar_t, ROWS> &b,
                                     size_type row, size_type col) const {
    for (size_type i = 0; i < ROWS; ++i) {
      element_getter()(m, i + row, col) = b[i];
    }
  }

  // Create zero matrix
  template <size_type ROWS, size_type COLS>
  ALGEBRA_HOST_DEVICE inline matrix_type<ROWS, COLS> zero() const {
    matrix_type<ROWS, COLS> ret;

    for (size_type j = 0; j < COLS; ++j) {
      for (size_type i = 0; i < ROWS; ++i) {
        element_getter()(ret, i, j) = 0;
      }
    }

    return ret;
  }

  // Create identity matrix
  template <size_type ROWS, size_type COLS>
  ALGEBRA_HOST_DEVICE inline matrix_type<ROWS, COLS> identity() const {
    matrix_type<ROWS, COLS> ret;

    for (size_type j = 0; j < COLS; ++j) {
      for (size_type i = 0; i < ROWS; ++i) {
        if (i == j) {
          element_getter()(ret, i, j) = 1;
        } else {
          element_getter()(ret, i, j) = 0;
        }
      }
    }

    return ret;
  }

  // Set input matrix as zero matrix
  template <size_type ROWS, size_type COLS>
  ALGEBRA_HOST_DEVICE inline void set_zero(matrix_type<ROWS, COLS> &m) const {

    for (size_type j = 0; j < COLS; ++j) {
      for (size_type i = 0; i < ROWS; ++i) {
        element_getter()(m, i, j) = 0;
      }
    }
  }

  // Set input matrix as identity matrix
  template <size_type ROWS, size_type COLS>
  ALGEBRA_HOST_DEVICE inline void set_identity(
      matrix_type<ROWS, COLS> &m) const {

    for (size_type j = 0; j < COLS; ++j) {
      for (size_type i = 0; i < ROWS; ++i) {
        if (i == j) {
          element_getter()(m, i, j) = 1;
        } else {
          element_getter()(m, i, j) = 0;
        }
      }
    }
  }

  // Create transpose matrix
  template <size_type ROWS, size_type COLS>
  ALGEBRA_HOST_DEVICE inline matrix_type<COLS, ROWS> transpose(
      const matrix_type<ROWS, COLS> &m) const {

    matrix_type<COLS, ROWS> ret;

    for (size_type i = 0; i < ROWS; ++i) {
      for (size_type j = 0; j < COLS; ++j) {
        element_getter()(ret, j, i) = element_getter()(m, i, j);
      }
    }

    return ret;
  }

  // Get determinant
  template <size_type N>
  ALGEBRA_HOST_DEVICE inline scalar_t determinant(
      const matrix_type<N, N> &m) const {

    return determinant_actor_t()(m);
  }

  // Create inverse matrix
  template <size_type N>
  ALGEBRA_HOST_DEVICE inline matrix_type<N, N> inverse(
      const matrix_type<N, N> &m) const {

    return inverse_actor_t()(m);
  }

  // Set matrix C to the product AB
  template <size_type M, size_type N, size_type O>
  ALGEBRA_HOST_DEVICE inline void set_product(
      matrix_type<M, O> &C, const matrix_type<M, N> &A,
      const matrix_type<N, O> &B) const {

    for (size_type i = 0; i < M; ++i) {
      for (size_type j = 0; j < O; ++j) {
        scalar_t T = 0.f;

        for (size_type k = 0; k < N; ++k) {
          T += element_getter()(A, i, k) * element_getter()(B, k, j);
        }

        element_getter()(C, i, j) = T;
      }
    }
  }

  // Set matrix C to the product A^TB
  template <size_type M, size_type N, size_type O>
  ALGEBRA_HOST_DEVICE inline void set_product_left_transpose(
      matrix_type<M, O> &C, const matrix_type<N, M> &A,
      const matrix_type<N, O> &B) const {

    for (size_type i = 0; i < M; ++i) {
      for (size_type j = 0; j < O; ++j) {
        scalar_t T = 0.f;

        for (size_type k = 0; k < N; ++k) {
          T += element_getter()(A, k, i) * element_getter()(B, k, j);
        }

        element_getter()(C, i, j) = T;
      }
    }
  }

  // Set matrix C to the product AB^T
  template <size_type M, size_type N, size_type O>
  ALGEBRA_HOST_DEVICE inline void set_product_right_transpose(
      matrix_type<M, O> &C, const matrix_type<M, N> &A,
      const matrix_type<O, N> &B) const {

    for (size_type i = 0; i < M; ++i) {
      for (size_type j = 0; j < O; ++j) {
        scalar_t T = 0.f;

        for (size_type k = 0; k < N; ++k) {
          T += element_getter()(A, i, k) * element_getter()(B, j, k);
        }

        element_getter()(C, i, j) = T;
      }
    }
  }

  // Set matrix A to the product AB in place
  template <size_type M>
  ALGEBRA_HOST_DEVICE inline void set_inplace_product_right(
      matrix_type<M, M> &A, const matrix_type<M, M> &B) const {

    for (size_type i = 0; i < M; ++i) {
      matrix_type<1, M> Q;

      for (size_type j = 0; j < M; ++j) {
        element_getter()(Q, 0, j) = element_getter()(A, i, j);
      }

      for (size_type j = 0; j < M; ++j) {
        scalar_t T = 0.f;

        for (size_type k = 0; k < M; ++k) {
          T += element_getter()(Q, 0, k) * element_getter()(B, k, j);
        }

        element_getter()(A, i, j) = T;
      }
    }
  }

  // Set matrix A to the product BA in place
  template <size_type M>
  ALGEBRA_HOST_DEVICE inline void set_inplace_product_left(
      matrix_type<M, M> &A, const matrix_type<M, M> &B) const {

    for (size_type j = 0; j < M; ++j) {
      matrix_type<1, M> Q;

      for (size_type i = 0; i < M; ++i) {
        element_getter()(Q, 0, i) = element_getter()(A, i, j);
      }

      for (size_type i = 0; i < M; ++i) {
        scalar_t T = 0.f;

        for (size_type k = 0; k < M; ++k) {
          T += element_getter()(B, i, k) * element_getter()(Q, 0, k);
        }

        element_getter()(A, i, j) = T;
      }
    }
  }

  // Set matrix A to the product AB^T in place
  template <size_type M>
  ALGEBRA_HOST_DEVICE inline void set_inplace_product_right_transpose(
      matrix_type<M, M> &A, const matrix_type<M, M> &B) const {

    for (size_type i = 0; i < M; ++i) {
      matrix_type<1, M> Q;

      for (size_type j = 0; j < M; ++j) {
        element_getter()(Q, 0, j) = element_getter()(A, i, j);
      }

      for (size_type j = 0; j < M; ++j) {
        scalar_t T = 0.f;

        for (size_type k = 0; k < M; ++k) {
          T += element_getter()(Q, 0, k) * element_getter()(B, j, k);
        }

        element_getter()(A, i, j) = T;
      }
    }
  }

  // Set matrix A to the product B^TA in place
  template <size_type M>
  ALGEBRA_HOST_DEVICE inline void set_inplace_product_left_transpose(
      matrix_type<M, M> &A, const matrix_type<M, M> &B) const {

    for (size_type j = 0; j < M; ++j) {
      matrix_type<1, M> Q;

      for (size_type i = 0; i < M; ++i) {
        element_getter()(Q, 0, i) = element_getter()(A, i, j);
      }

      for (size_type i = 0; i < M; ++i) {
        scalar_t T = 0.f;

        for (size_type k = 0; k < M; ++k) {
          T += element_getter()(B, k, i) * element_getter()(Q, 0, k);
        }

        element_getter()(A, i, j) = T;
      }
    }
  }
};

namespace determinant {

template <typename size_type,
          template <typename, size_type, size_type> class matrix_t,
          typename scalar_t, class... As>
struct actor {

  /// 2D matrix type
  template <size_type ROWS, size_type COLS>
  using matrix_type = matrix_t<scalar_t, ROWS, COLS>;

  template <size_type N>
  ALGEBRA_HOST_DEVICE inline scalar_t operator()(
      const matrix_type<N, N> &m) const {

    return typename find_algorithm<size_type, N, As...>::algorithm_type()(m);
  }
};

}  // namespace determinant

namespace inverse {

template <typename size_type,
          template <typename, size_type, size_type> class matrix_t,
          typename scalar_t, class... As>
struct actor {

  /// 2D matrix type
  template <size_type ROWS, size_type COLS>
  using matrix_type = matrix_t<scalar_t, ROWS, COLS>;

  template <size_type N>
  ALGEBRA_HOST_DEVICE inline matrix_type<N, N> operator()(
      const matrix_type<N, N> &m) const {

    return typename find_algorithm<size_type, N, As...>::algorithm_type()(m);
  }
};

}  // namespace inverse

}  // namespace algebra::cmath::matrix
