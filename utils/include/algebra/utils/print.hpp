/** Algebra plugins, part of the ACTS project
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "algebra/concepts.hpp"
#include "algebra/qualifiers.hpp"

// System include(s).
#include <iostream>

namespace algebra {

/// Print a generic vector or point @param v
template <typename vector_t>
requires(concepts::vector<vector_t> ||
         concepts::point<vector_t>) ALGEBRA_HOST std::ostream&
operator<<(std::ostream& out, const vector_t& v) {

  using index_t = algebra::traits::index_t<vector_t>;

  constexpr index_t size{algebra::traits::size<vector_t>};

  out << "[";
  for (index_t i = 0; i < size; ++i) {
    out << v[i];
    if (i != size - 1) {
      out << ", ";
    }
  }
  out << "]";

  return out;
}

/// Print a column matrix @param v
template <concepts::column_matrix vector_t>
ALGEBRA_HOST std::ostream& operator<<(std::ostream& out, const vector_t& v) {

  using index_t = algebra::traits::index_t<vector_t>;
  using element_getter_t = algebra::traits::element_getter_t<vector_t>;

  constexpr index_t rows{algebra::traits::rows<vector_t>};

  out << "[";
  for (index_t i = 0; i < rows; ++i) {
    out << element_getter_t{}(v, i, 0);
    if (i != rows - 1) {
      out << ", ";
    }
  }
  out << "]";

  return out;
}

/// Print a generic matrix @param m
template <concepts::matrix matrix_t>
ALGEBRA_HOST std::ostream& operator<<(std::ostream& out, const matrix_t& m) {

  using index_t = algebra::traits::index_t<matrix_t>;
  using element_getter_t = algebra::traits::element_getter_t<matrix_t>;

  constexpr index_t rows{algebra::traits::rows<matrix_t>};
  constexpr index_t columns{algebra::traits::columns<matrix_t>};

  out << "[";
  for (index_t i = 0; i < rows; ++i) {
    out << "[";
    for (index_t j = 0; j < columns; ++j) {
      out << element_getter_t{}(m, i, j);
      if (j != columns - 1) {
        out << ", ";
      }
    }
    out << "]";
    if (i != rows - 1) {
      out << ", ";
    }
  }
  out << "]";

  return out;
}

/// Print a 3D transform @param trf
template <concepts::transform3D transform_t>
ALGEBRA_HOST std::ostream& operator<<(std::ostream& out,
                                      const transform_t& trf) {
  out << trf.matrix();

  return out;
}

}  // namespace algebra
