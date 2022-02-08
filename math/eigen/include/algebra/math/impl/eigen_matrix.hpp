/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/qualifiers.hpp"

namespace algebra::eigen::math {

template <template <typename T, std::size_t, std::size_t> class matrix_type,
          typename T, std::size_t ROWS, std::size_t COLS>
struct matrix {

  matrix_type<T, ROWS, COLS> _data;

  ALGEBRA_HOST_DEVICE
  T& operator()(unsigned int row, unsigned int col) { return _data(row, col); }

  ALGEBRA_HOST_DEVICE
  const T& operator()(unsigned int row, unsigned int col) const {
    return _data(row, col);
  }
};

}  // namespace algebra::eigen::math