/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/qualifiers.hpp"
#include "algebra/type_traits.hpp"

// System include(s)
#include <type_traits>

namespace algebra::cmath::matrix::determinant {

/// "Determinant getter", assuming a N X N matrix
template <class matrix_t, class element_getter_t>
struct hard_coded {

  using scalar_type = algebra::trait::value_t<matrix_t>;
  using size_type = algebra::trait::index_t<matrix_t>;

  /// Function (object) used for accessing a matrix element
  using element_getter = element_getter_t;

  // 2 X 2 matrix determinant
  template <typename M = matrix_t,
            std::enable_if_t<algebra::trait::rank<M> == 2, bool> = true>
  ALGEBRA_HOST_DEVICE inline scalar_type operator()(const matrix_t &m) const {

    return element_getter()(m, 0, 0) * element_getter()(m, 1, 1) -
           element_getter()(m, 0, 1) * element_getter()(m, 1, 0);
  }

  // 4 X 4 matrix determinant
  template <typename M = matrix_t,
            std::enable_if_t<algebra::trait::rank<M> == 4, bool> = true>
  ALGEBRA_HOST_DEVICE inline scalar_type operator()(const matrix_t &m) const {

    return element_getter()(m, 0, 3) * element_getter()(m, 1, 2) *
               element_getter()(m, 2, 1) * element_getter()(m, 3, 0) -
           element_getter()(m, 0, 2) * element_getter()(m, 1, 3) *
               element_getter()(m, 2, 1) * element_getter()(m, 3, 0) -
           element_getter()(m, 0, 3) * element_getter()(m, 1, 1) *
               element_getter()(m, 2, 2) * element_getter()(m, 3, 0) +
           element_getter()(m, 0, 1) * element_getter()(m, 1, 3) *
               element_getter()(m, 2, 2) * element_getter()(m, 3, 0) +
           element_getter()(m, 0, 2) * element_getter()(m, 1, 1) *
               element_getter()(m, 2, 3) * element_getter()(m, 3, 0) -
           element_getter()(m, 0, 1) * element_getter()(m, 1, 2) *
               element_getter()(m, 2, 3) * element_getter()(m, 3, 0) -
           element_getter()(m, 0, 3) * element_getter()(m, 1, 2) *
               element_getter()(m, 2, 0) * element_getter()(m, 3, 1) +
           element_getter()(m, 0, 2) * element_getter()(m, 1, 3) *
               element_getter()(m, 2, 0) * element_getter()(m, 3, 1) +
           element_getter()(m, 0, 3) * element_getter()(m, 1, 0) *
               element_getter()(m, 2, 2) * element_getter()(m, 3, 1) -
           element_getter()(m, 0, 0) * element_getter()(m, 1, 3) *
               element_getter()(m, 2, 2) * element_getter()(m, 3, 1) -
           element_getter()(m, 0, 2) * element_getter()(m, 1, 0) *
               element_getter()(m, 2, 3) * element_getter()(m, 3, 1) +
           element_getter()(m, 0, 0) * element_getter()(m, 1, 2) *
               element_getter()(m, 2, 3) * element_getter()(m, 3, 1) +
           element_getter()(m, 0, 3) * element_getter()(m, 1, 1) *
               element_getter()(m, 2, 0) * element_getter()(m, 3, 2) -
           element_getter()(m, 0, 1) * element_getter()(m, 1, 3) *
               element_getter()(m, 2, 0) * element_getter()(m, 3, 2) -
           element_getter()(m, 0, 3) * element_getter()(m, 1, 0) *
               element_getter()(m, 2, 1) * element_getter()(m, 3, 2) +
           element_getter()(m, 0, 0) * element_getter()(m, 1, 3) *
               element_getter()(m, 2, 1) * element_getter()(m, 3, 2) +
           element_getter()(m, 0, 1) * element_getter()(m, 1, 0) *
               element_getter()(m, 2, 3) * element_getter()(m, 3, 2) -
           element_getter()(m, 0, 0) * element_getter()(m, 1, 1) *
               element_getter()(m, 2, 3) * element_getter()(m, 3, 2) -
           element_getter()(m, 0, 2) * element_getter()(m, 1, 1) *
               element_getter()(m, 2, 0) * element_getter()(m, 3, 3) +
           element_getter()(m, 0, 1) * element_getter()(m, 1, 2) *
               element_getter()(m, 2, 0) * element_getter()(m, 3, 3) +
           element_getter()(m, 0, 2) * element_getter()(m, 1, 0) *
               element_getter()(m, 2, 1) * element_getter()(m, 3, 3) -
           element_getter()(m, 0, 0) * element_getter()(m, 1, 2) *
               element_getter()(m, 2, 1) * element_getter()(m, 3, 3) -
           element_getter()(m, 0, 1) * element_getter()(m, 1, 0) *
               element_getter()(m, 2, 2) * element_getter()(m, 3, 3) +
           element_getter()(m, 0, 0) * element_getter()(m, 1, 1) *
               element_getter()(m, 2, 2) * element_getter()(m, 3, 3);
  }
};

}  // namespace algebra::cmath::matrix::determinant