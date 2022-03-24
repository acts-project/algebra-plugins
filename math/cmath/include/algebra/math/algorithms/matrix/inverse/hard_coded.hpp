/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/algorithms/matrix/determinant/hard_coded.hpp"
#include "algebra/qualifiers.hpp"

// System include(s)
#include <type_traits>

namespace algebra::cmath::matrix::inverse {

/// "inverse getter", assuming a N X N matrix
template <typename size_type,
          template <typename, size_type, size_type> class matrix_t,
          typename scalar_t, class element_getter_t, size_type... Ds>
struct hard_coded {

  using _dims = std::integer_sequence<size_type, Ds...>;

  /// Function (object) used for accessing a matrix element
  using element_getter = element_getter_t;

  /// 2D matrix type
  template <size_type ROWS, size_type COLS>
  using matrix_type = matrix_t<scalar_t, ROWS, COLS>;

  using determinant_getter_type =
      determinant::hard_coded<size_type, matrix_t, scalar_t, element_getter_t>;

  // 2 X 2 matrix inverse
  template <size_type N, std::enable_if_t<N == 2, bool> = true>
  ALGEBRA_HOST_DEVICE inline matrix_type<N, N> operator()(
      const matrix_type<N, N> &m) const {

    matrix_type<N, N> ret;

    scalar_t det = determinant_getter_type()(m);

    element_getter()(ret, 0, 0) = element_getter()(m, 1, 1) / det;
    element_getter()(ret, 0, 1) = -1 * element_getter()(m, 0, 1) / det;
    element_getter()(ret, 1, 0) = -1 * element_getter()(m, 1, 0) / det;
    element_getter()(ret, 1, 1) = element_getter()(m, 0, 0) / det;

    return ret;
  }
};

}  // namespace algebra::cmath::matrix::inverse