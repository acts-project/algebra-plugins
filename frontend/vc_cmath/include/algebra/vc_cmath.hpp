/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/cmath.hpp"
#include "algebra/math/vc.hpp"
#include "algebra/storage/vc.hpp"

/// @name Operators on @c algebra::vc types
/// @{

using algebra::cmath::operator*;
using algebra::cmath::operator-;
using algebra::cmath::operator+;

/// @}

namespace algebra {
namespace vc {

/// Functor used to access elements of Vc matrices
template <typename scalar_t>
struct element_getter {

  template <std::size_t ROWS, std::size_t COLS>
  using matrix_type = Vc::array<Vc::array<scalar_t, ROWS>, COLS>;

  template <std::size_t ROWS, std::size_t COLS>
  ALGEBRA_HOST_DEVICE inline scalar_t& operator()(matrix_type<ROWS, COLS>& m,
                                                  std::size_t row,
                                                  std::size_t col) const {

    return m[col][row];
  }

  template <std::size_t ROWS, std::size_t COLS>
  ALGEBRA_HOST_DEVICE inline scalar_t operator()(
      const matrix_type<ROWS, COLS>& m, std::size_t row,
      std::size_t col) const {

    return m[col][row];
  }
};  // element_getter

/// Functor used to extract a block from Vc matrices
template <typename scalar_t>
struct block_getter {

  template <std::size_t ROWS, std::size_t COLS>
  using matrix_type = Vc::array<Vc::array<scalar_t, ROWS>, COLS>;

  template <std::size_t ROWS, std::size_t COLS, class input_matrix_type>
  ALGEBRA_HOST_DEVICE matrix_type<ROWS, COLS> operator()(
      const input_matrix_type& m, std::size_t row, std::size_t col) const {

    matrix_type<ROWS, COLS> submatrix{};
    for (std::size_t icol = col; icol < col + COLS; ++icol) {
      for (std::size_t irow = row; irow < row + ROWS; ++irow) {
        submatrix[icol - col][irow - row] = m[icol][irow];
      }
    }
    return submatrix;
  }
};  // struct block_getter

/// @name cmath based transforms on @c algebra::vc types
/// @{

// Pull in the definitions needed by the cmath transforms, into this namespace.
using math::cross;
using math::perp;
using math::phi;

template <typename T>
using transform3 =
    cmath::transform3<std::size_t, vc::storage_type, T,
                      Vc::array<Vc::array<T, 4>, 4>, element_getter<T>,
                      block_getter<T>, vc::vector3<T>, vc::point2<T> >;
template <typename T>
using cartesian2 = cmath::cartesian2<transform3<T> >;
template <typename T>
using polar2 = cmath::polar2<transform3<T> >;
template <typename T>
using cylindrical2 = cmath::cylindrical2<transform3<T> >;
template <typename T, std::size_t ROWS, std::size_t COLS>
using matrix = math::matrix<matrix_type, T, ROWS, COLS>;

/// @}

}  // namespace vc

namespace getter {

/// @name Getter functions on @c algebra::vc types
/// @{

using cmath::eta;
using cmath::norm;
using cmath::perp;
using cmath::phi;
using cmath::theta;

using vc::math::eta;
using vc::math::norm;
using vc::math::perp;
using vc::math::phi;
using vc::math::theta;

/// @|

/// Function extracting a slice from the matrix used by
/// @c algebra::vc::transform3
template <std::size_t SIZE, std::size_t ROWS, std::size_t COLS,
          typename scalar_t>
ALGEBRA_HOST_DEVICE inline vc::storage_type<scalar_t, SIZE> vector(
    const Vc::array<Vc::array<scalar_t, ROWS>, COLS>& m, std::size_t row,
    std::size_t col) {

  return cmath::vector_getter<std::size_t, Vc::array, scalar_t, SIZE,
                              vc::storage_type<scalar_t, SIZE> >()(m, row, col);
}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::vc::storage_type
/// @{

using cmath::cross;
using cmath::dot;
using cmath::normalize;

using vc::math::cross;
using vc::math::dot;
using vc::math::normalize;

/// @}

}  // namespace vector
}  // namespace algebra
