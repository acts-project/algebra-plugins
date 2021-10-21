/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/storage/vc.hpp"
#include "algebra/math/cmath.hpp"

namespace algebra {
namespace vc {

/// Functor used to access elements of Vc matrices
struct element_getter {

  template<std::size_t ROWS, std::size_t COLS>
  using matrix_type = Vc::array<Vc::array<scalar, ROWS>, COLS>;

  template<std::size_t ROWS, std::size_t COLS>
  ALGEBRA_HOST_DEVICE inline scalar&
  operator()(matrix_type<ROWS, COLS> &m, std::size_t row, std::size_t col) const {

    return m[col][row];
  }

  template<std::size_t ROWS, std::size_t COLS>
  ALGEBRA_HOST_DEVICE inline scalar
  operator()(const matrix_type<ROWS, COLS> &m, std::size_t row, std::size_t col) const {

    return m[col][row];
  }
};  // element_getter

/// Functor used to extract a block from Vc matrices
struct block_getter {

  template<std::size_t ROWS, std::size_t COLS>
  using matrix_type = Vc::array<Vc::array<scalar, ROWS>, COLS>;

  template <std::size_t ROWS, std::size_t COLS, class input_matrix_type>
  ALGEBRA_HOST_DEVICE
  matrix_type<ROWS, COLS> operator()(const input_matrix_type &m,
                                     std::size_t row, std::size_t col) const {

    matrix_type< ROWS, COLS> submatrix;
    for (std::size_t icol = col; icol < col + COLS; ++icol) {
      for (std::size_t irow = row; irow < row + ROWS; ++irow) {
        submatrix[icol - col][irow - row] = m[icol][irow];
      }
    }
    return submatrix;
  }
};  // struct block_getter

using transform3 = cmath::transform3<storage_type, scalar, Vc::array<Vc::array<scalar, 4>, 4 >, element_getter, block_getter >;
using cartesian2 = cmath::cartesian2<storage_type, scalar, transform3>;
using polar2 = cmath::polar2<storage_type, scalar, transform3>;
using cylindrical2 = cmath::cylindrical2<storage_type, scalar, transform3>;

}  // namespace vc

namespace getter {

auto phi = [](const auto& a) { return cmath::phi<vc::storage_type>(a); };
auto theta = [](const auto& a) { return cmath::theta<vc::storage_type>(a); };
auto perp = [](const auto& a) { return cmath::perp<vc::storage_type>(a); };
auto norm = [](const auto& a) { return cmath::norm<vc::storage_type>(a); };
auto eta = [](const auto& a) { return cmath::eta<vc::storage_type>(a); };

}  // namespace getter

namespace vector {

auto cross = [](const auto& a, const auto& b) { return cmath::cross<vc::storage_type>(a, b); };
auto dot = [](const auto& a, const auto& b) { return cmath::dot<vc::storage_type>(a, b); };
auto normalize = [](const auto& a) { return cmath::normalize<vc::storage_type>(a); };

}  // namespace vector
}  // namespace algebra
