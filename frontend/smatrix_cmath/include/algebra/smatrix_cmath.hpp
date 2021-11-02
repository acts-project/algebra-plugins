/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/scalar.hpp"
#include "algebra/math/cmath.hpp"
#include "algebra/math/smatrix.hpp"
#include "algebra/storage/smatrix.hpp"

// ROOT/Smatrix include(s).
#include <Math/SMatrix.h>

namespace algebra {
namespace smatrix {

/// Functor used to access elements of Vc matrices
struct element_getter {

  template <unsigned int ROWS, unsigned int COLS>
  using matrix_type = ROOT::Math::SMatrix<scalar, ROWS, COLS>;

  template <unsigned int ROWS, unsigned int COLS>
  ALGEBRA_HOST_DEVICE inline scalar& operator()(matrix_type<ROWS, COLS>& m,
                                                unsigned int row,
                                                unsigned int col) const {

    return m(col, row);
  }

  template <unsigned int ROWS, unsigned int COLS>
  ALGEBRA_HOST_DEVICE inline scalar operator()(const matrix_type<ROWS, COLS>& m,
                                               unsigned int row,
                                               unsigned int col) const {

    return m(col, row);
  }
};  // element_getter

/// Functor used to extract a block from Vc matrices
struct block_getter {

  template <unsigned int ROWS, unsigned int COLS>
  using matrix_type = ROOT::Math::SMatrix<scalar, ROWS, COLS>;

  template <unsigned int ROWS, unsigned int COLS, class input_matrix_type>
  ALGEBRA_HOST_DEVICE matrix_type<ROWS, COLS> operator()(
      const input_matrix_type& m, unsigned int row, unsigned int col) const {

    return m.template Sub<matrix_type<ROWS, COLS> >(row, col);
  }
};  // struct block_getter

/// @name cmath based transforms on @c algebra::smatrix::storage_type
/// @{

using transform3 = cmath::transform3<unsigned int, smatrix::storage_type,
                                     scalar, ROOT::Math::SMatrix<scalar, 4, 4>,
                                     element_getter, block_getter>;
using cartesian2 = cmath::cartesian2<transform3>;
using polar2 = cmath::polar2<transform3>;
using cylindrical2 = cmath::cylindrical2<transform3>;

/// @}

}  // namespace smatrix

namespace getter {

/// @name Getter functions on @c algebra::smatrix::storage_type
/// @{

using smatrix::math::eta;
using smatrix::math::norm;
using smatrix::math::perp;
using smatrix::math::phi;
using smatrix::math::theta;

/// @}

/// Function extracting a slice from the matrix used by
/// @c algebra::smatrix::transform3
template <unsigned int SIZE, unsigned int ROWS, unsigned int COLS>
ALGEBRA_HOST_DEVICE inline auto vector(
    const ROOT::Math::SMatrix<scalar, ROWS, COLS>& m, unsigned int row,
    unsigned int col) {

  return m.template SubCol<smatrix::storage_type<scalar, SIZE> >(col, row);
}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::smatrix::storage_type
/// @{

using smatrix::math::cross;
using smatrix::math::dot;
using smatrix::math::normalize;

/// @}

}  // namespace vector
}  // namespace algebra
