/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/cmath.hpp"
#include "algebra/math/smatrix.hpp"
#include "algebra/storage/smatrix.hpp"

// ROOT/Smatrix include(s).
#include <Math/SMatrix.h>

namespace algebra {
namespace smatrix {

/// Functor used to access elements of Vc matrices
template <typename scalar_t>
struct element_getter {

  template <unsigned int ROWS, unsigned int COLS>
  using matrix_type = ROOT::Math::SMatrix<scalar_t, ROWS, COLS>;

  template <unsigned int ROWS, unsigned int COLS>
  ALGEBRA_HOST_DEVICE inline scalar_t& operator()(matrix_type<ROWS, COLS>& m,
                                                  unsigned int row,
                                                  unsigned int col) const {

    return m(col, row);
  }

  template <unsigned int ROWS, unsigned int COLS>
  ALGEBRA_HOST_DEVICE inline scalar_t operator()(
      const matrix_type<ROWS, COLS>& m, unsigned int row,
      unsigned int col) const {

    return m(col, row);
  }
};  // element_getter

/// Functor used to extract a block from Vc matrices
template <typename scalar_t>
struct block_getter {

  template <unsigned int ROWS, unsigned int COLS>
  using matrix_type = ROOT::Math::SMatrix<scalar_t, ROWS, COLS>;

  template <unsigned int ROWS, unsigned int COLS, class input_matrix_type>
  ALGEBRA_HOST_DEVICE matrix_type<ROWS, COLS> operator()(
      const input_matrix_type& m, unsigned int row, unsigned int col) const {

    return m.template Sub<matrix_type<ROWS, COLS> >(row, col);
  }
};  // struct block_getter

/// @name cmath based transforms on @c algebra::smatrix::storage_type
/// @{

template <typename T>
using transform3 = cmath::transform3<unsigned int, smatrix::storage_type, T,
                                     ROOT::Math::SMatrix<T, 4, 4>,
                                     element_getter<T>, block_getter<T> >;
template <typename T>
using cartesian2 = cmath::cartesian2<transform3<T> >;
template <typename T>
using polar2 = cmath::polar2<transform3<T> >;
template <typename T>
using cylindrical2 = cmath::cylindrical2<transform3<T> >;
template <typename T, std::size_t ROWS, std::size_t COLS>
using matrix = math::matrix<matrix_type, T, ROWS, COLS>;

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
template <unsigned int SIZE, unsigned int ROWS, unsigned int COLS,
          typename scalar_t>
ALGEBRA_HOST_DEVICE inline auto vector(
    const ROOT::Math::SMatrix<scalar_t, ROWS, COLS>& m, unsigned int row,
    unsigned int col) {

  return m.template SubCol<smatrix::storage_type<scalar_t, SIZE> >(col, row);
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
