/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/scalar.hpp"
#include "algebra/math/cmath.hpp"
#include "algebra/storage/vc.hpp"

namespace algebra {

using cmath::operator*;
using cmath::operator-;
using cmath::operator+;

namespace vc {

/// Functor used to access elements of Vc matrices
struct element_getter {

  template <std::size_t ROWS, std::size_t COLS>
  using matrix_type = Vc::array<Vc::array<scalar, ROWS>, COLS>;

  template <std::size_t ROWS, std::size_t COLS>
  ALGEBRA_HOST_DEVICE inline scalar& operator()(matrix_type<ROWS, COLS>& m,
                                                std::size_t row,
                                                std::size_t col) const {

    return m[col][row];
  }

  template <std::size_t ROWS, std::size_t COLS>
  ALGEBRA_HOST_DEVICE inline scalar operator()(const matrix_type<ROWS, COLS>& m,
                                               std::size_t row,
                                               std::size_t col) const {

    return m[col][row];
  }
};  // element_getter

/// Functor used to extract a block from Vc matrices
struct block_getter {

  template <std::size_t ROWS, std::size_t COLS>
  using matrix_type = Vc::array<Vc::array<scalar, ROWS>, COLS>;

  template <std::size_t ROWS, std::size_t COLS, class input_matrix_type>
  ALGEBRA_HOST_DEVICE matrix_type<ROWS, COLS> operator()(
      const input_matrix_type& m, std::size_t row, std::size_t col) const {

    matrix_type<ROWS, COLS> submatrix;
    for (std::size_t icol = col; icol < col + COLS; ++icol) {
      for (std::size_t irow = row; irow < row + ROWS; ++irow) {
        submatrix[icol - col][irow - row] = m[icol][irow];
      }
    }
    return submatrix;
  }
};  // struct block_getter

using transform3 =
    cmath::transform3<vc::storage_type, scalar,
                      Vc::array<Vc::array<scalar, 4>, 4>, element_getter,
                      block_getter, vc::vector3, vc::point2>;
using cartesian2 = cmath::cartesian2<transform3>;
using polar2 = cmath::polar2<transform3>;
using cylindrical2 = cmath::cylindrical2<transform3>;

}  // namespace vc

namespace getter {

using cmath::eta;
using cmath::norm;
using cmath::perp;
using cmath::phi;
using cmath::theta;

template <auto SIZE, auto ROWS, auto COLS>
ALGEBRA_HOST_DEVICE inline vc::storage_type<scalar, SIZE> vector(
    const Vc::array<Vc::array<scalar, ROWS>, COLS>& m, std::size_t row,
    std::size_t col) {

  return cmath::vector_getter<vc::storage_type, scalar>()
      .template operator()<SIZE>(m, row, col);
}

}  // namespace getter

namespace vector {

using cmath::cross;
using cmath::dot;
using cmath::normalize;

}  // namespace vector
}  // namespace algebra
