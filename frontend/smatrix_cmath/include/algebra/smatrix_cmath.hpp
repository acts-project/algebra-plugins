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

using transform3 =
    cmath::transform3<storage_type, scalar, ROOT::Math::SMatrix<scalar, 4, 4>,
                      element_getter, block_getter>;
using cartesian2 = cmath::cartesian2<storage_type, scalar, transform3>;
using polar2 = cmath::polar2<storage_type, scalar, transform3>;
using cylindrical2 = cmath::cylindrical2<storage_type, scalar, transform3>;

}  // namespace smatrix

namespace getter {

auto phi = [](const auto& a) { return cmath::phi<smatrix::storage_type>(a); };
auto theta = [](const auto& a) {
  return cmath::theta<smatrix::storage_type>(a);
};
auto perp = [](const auto& a) { return cmath::perp<smatrix::storage_type>(a); };
auto norm = [](const auto& a) { return cmath::norm<smatrix::storage_type>(a); };
auto eta = [](const auto& a) { return cmath::eta<smatrix::storage_type>(a); };

}  // namespace getter

namespace vector {

auto cross = [](const auto& a, const auto& b) {
  return cmath::cross<smatrix::storage_type>(a, b);
};
auto dot = [](const auto& a, const auto& b) {
  return cmath::dot<smatrix::storage_type>(a, b);
};
auto normalize = [](const auto& a) {
  return cmath::normalize<smatrix::storage_type>(a);
};

}  // namespace vector
}  // namespace algebra
