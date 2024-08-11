/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/algorithms/utils/algorithm_finder.hpp"
#include "algebra/math/common.hpp"
#include "algebra/qualifiers.hpp"

namespace algebra::generic::math {

/// Create zero matrix - generic transform3
template <class matrix_t>
ALGEBRA_HOST_DEVICE inline matrix_t zero() {

  using index_t = algebra::trait::index_t<matrix_t>;
  using element_getter_t = algebra::trait::element_getter_t<matrix_t>;

  matrix_t ret;

  for (index_t j = 0; j < algebra::trait::columns<matrix_t>; ++j) {
    for (index_t i = 0; i < algebra::trait::rows<matrix_t>; ++i) {
      element_getter_t{}(ret, i, j) = 0;
    }
  }

  return ret;
}

/// Create identity matrix - generic transform3
template <class matrix_t>
ALGEBRA_HOST_DEVICE inline matrix_t identity() {

  using index_t = algebra::trait::index_t<matrix_t>;
  using element_getter_t = algebra::trait::element_getter_t<matrix_t>;

  auto ret{zero<matrix_t>()};

  for (index_t i = 0; i < algebra::trait::rank<matrix_t>; ++i) {
    element_getter_t{}(ret, i, i) = 1;
  }

  return ret;
}

/// Set @param m as zero matrix
template <class matrix_t>
ALGEBRA_HOST_DEVICE inline void set_zero(matrix_t &m) {
  m = zero<matrix_t>();
}

/// Set @param m as identity matrix
template <class matrix_t>
ALGEBRA_HOST_DEVICE inline void set_identity(matrix_t &m) {
  m = identity<matrix_t>();
}

/// @returns the transpose matrix of @param m
template <class matrix_t>
ALGEBRA_HOST_DEVICE inline auto transpose(const matrix_t &m) {

  using index_t = algebra::trait::index_t<matrix_t>;
  using value_t = algebra::trait::value_t<matrix_t>;
  using element_getter_t = algebra::trait::element_getter_t<matrix_t>;

  constexpr index_t rows{algebra::trait::rows<matrix_t>};
  constexpr index_t columns{algebra::trait::columns<matrix_t>};

  algebra::trait::get_matrix_t<matrix_t, columns, rows, value_t> ret;

  for (index_t i = 0; i < rows; ++i) {
    for (index_t j = 0; j < columns; ++j) {
      element_getter_t{}(ret, j, i) = element_getter_t{}(m, i, j);
    }
  }

  return ret;
}

/// @returns the determinant of @param m
template <class matrix_t>
ALGEBRA_HOST_DEVICE inline algebra::trait::scalar_t<matrix_t> determinant(
    const matrix_t &m) {

  return determinant_t<matrix_t>{}(m);
}

/// @returns the determinant of @param m
template <class matrix_t>
ALGEBRA_HOST_DEVICE inline matrix_t inverse(const matrix_t &m) {

  return inversion_t<matrix_t>{}(m);
}

}  // namespace algebra::generic::math
