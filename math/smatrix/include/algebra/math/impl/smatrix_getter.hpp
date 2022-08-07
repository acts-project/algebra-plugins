/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/qualifiers.hpp"

// ROOT/Smatrix include(s).
#include <Math/Expression.h>
#include <Math/Functions.h>
#include <Math/SMatrix.h>
#include <Math/SVector.h>
#include <TMath.h>

// System include(s).
#include <cassert>

namespace algebra::smatrix::math {

/// Functor used to access elements of Vc matrices
template <typename scalar_t>
struct element_getter {

  template <unsigned int ROWS, unsigned int COLS>
  using matrix_type = ROOT::Math::SMatrix<scalar_t, ROWS, COLS>;

  template <unsigned int ROWS, unsigned int COLS>
  ALGEBRA_HOST_DEVICE inline scalar_t &operator()(matrix_type<ROWS, COLS> &m,
                                                  unsigned int row,
                                                  unsigned int col) const {

    assert(row < ROWS);
    assert(col < COLS);
    return m(row, col);
  }

  template <unsigned int ROWS, unsigned int COLS>
  ALGEBRA_HOST_DEVICE inline scalar_t operator()(
      const matrix_type<ROWS, COLS> &m, unsigned int row,
      unsigned int col) const {

    assert(row < ROWS);
    assert(col < COLS);
    return m(row, col);
  }
};  // element_getter

/// Function extracting an element from a matrix (const)
template <typename scalar_t, unsigned int ROWS, unsigned int COLS>
ALGEBRA_HOST_DEVICE inline scalar_t element(
    const ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &m, unsigned int row,
    unsigned int col) {

  return element_getter<scalar_t>()(m, row, col);
}

/// Function extracting an element from a matrix (non-const)
template <typename scalar_t, unsigned int ROWS, unsigned int COLS>
ALGEBRA_HOST_DEVICE inline scalar_t &element(
    ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &m, unsigned int row,
    unsigned int col) {

  return element_getter<scalar_t>()(m, row, col);
}

/// Functor used to extract a block from SMatrix matrices
template <typename scalar_t>
struct block_getter {

  template <unsigned int ROWS, unsigned int COLS>
  using matrix_type = ROOT::Math::SMatrix<scalar_t, ROWS, COLS>;

  template <unsigned int ROWS, unsigned int COLS, class input_matrix_type>
  ALGEBRA_HOST_DEVICE matrix_type<ROWS, COLS> operator()(
      const input_matrix_type &m, unsigned int row, unsigned int col) const {

    return m.template Sub<matrix_type<ROWS, COLS> >(row, col);
  }
};  // struct block_getter

}  // namespace algebra::smatrix::math
