/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/concepts.hpp"
#include "algebra/math/common.hpp"
#include "algebra/qualifiers.hpp"
#include "algebra/type_traits.hpp"

namespace algebra::generic::math {

/// This method retrieves phi from a vector with rows >= 2
///
/// @param v the input vector
template <concepts::vector vector_t>
ALGEBRA_HOST_DEVICE inline algebra::traits::scalar_t<vector_t> phi(
    const vector_t &v) noexcept {

  using element_getter_t = algebra::traits::element_getter_t<vector_t>;

  return algebra::math::atan2(element_getter_t{}(v, 1),
                              element_getter_t{}(v, 0));
}

/// This method retrieves the perpendicular magnitude of a vector with rows >= 2
///
/// @param v the input vector
template <concepts::vector vector_t>
ALGEBRA_HOST_DEVICE inline algebra::traits::scalar_t<vector_t> perp(
    const vector_t &v) noexcept {

  using element_getter_t = algebra::traits::element_getter_t<vector_t>;

  return algebra::math::sqrt(
      algebra::math::fma(element_getter_t{}(v, 0), element_getter_t{}(v, 0),
                         element_getter_t{}(v, 1) * element_getter_t{}(v, 1)));
}

/// This method retrieves theta from a vector with rows >= 3
///
/// @param v the input vector
template <concepts::vector vector_t>
ALGEBRA_HOST_DEVICE inline algebra::traits::scalar_t<vector_t> theta(
    const vector_t &v) noexcept {

  using element_getter_t = algebra::traits::element_getter_t<vector_t>;

  return algebra::math::atan2(perp(v), element_getter_t{}(v, 2));
}

/// Cross product between two input vectors - 3 Dim
///
/// @tparam vector1_t first vector or column matrix type
/// @tparam vector2_t second vector or column matrix type
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return a vector representing the cross product
template <typename vector1_t, typename vector2_t>
requires(
    (concepts::vector3D<vector1_t> || concepts::column_matrix3D<vector1_t>)&&(
        concepts::vector3D<vector2_t> ||
        concepts::column_matrix3D<vector2_t>)) ALGEBRA_HOST_DEVICE
    inline algebra::traits::vector_t<vector1_t> cross(const vector1_t &a,
                                                      const vector2_t &b) {

  using element_getter_t = algebra::traits::element_getter_t<vector1_t>;

  return {
      algebra::math::fma(element_getter_t{}(a, 1), element_getter_t{}(b, 2),
                         -element_getter_t{}(b, 1) * element_getter_t{}(a, 2)),
      algebra::math::fma(element_getter_t{}(a, 2), element_getter_t{}(b, 0),
                         -element_getter_t{}(b, 2) * element_getter_t{}(a, 0)),
      algebra::math::fma(element_getter_t{}(a, 0), element_getter_t{}(b, 1),
                         -element_getter_t{}(b, 0) * element_getter_t{}(a, 1))};
}

/// Dot product between two input vectors
///
/// @tparam vector1_t first vector or column matrix type
/// @tparam vector2_t second vector or column matrix type
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return the scalar dot product value
template <typename vector1_t, typename vector2_t>
requires((concepts::vector<vector1_t> || concepts::column_matrix<vector1_t>)&&(
    concepts::vector<vector2_t> ||
    concepts::column_matrix<vector2_t>)) ALGEBRA_HOST_DEVICE
    inline algebra::traits::scalar_t<vector1_t> dot(const vector1_t &a,
                                                    const vector2_t &b) {

  using scalar_t = algebra::traits::scalar_t<vector1_t>;
  using index_t = algebra::traits::index_t<vector1_t>;
  using element_getter_t = algebra::traits::element_getter_t<vector1_t>;

  scalar_t ret{element_getter_t{}(a, 0) * element_getter_t{}(b, 0)};

  for (index_t i = 1; i < algebra::traits::size<vector1_t>; i++) {
    ret = std::fma(element_getter_t{}(a, i), element_getter_t{}(b, i), ret);
  }

  return ret;
}

/// This method retrieves the norm of a vector with rows >= 3
///
/// @param v the input vector
template <concepts::vector vector_t>
ALGEBRA_HOST_DEVICE inline algebra::traits::scalar_t<vector_t> norm(
    const vector_t &v) {

  return algebra::math::sqrt(dot(v, v));
}

/// This method retrieves the pseudo-rapidity from a vector or vector base with
/// rows >= 3
///
/// @param v the input vector
template <concepts::vector vector_t>
ALGEBRA_HOST_DEVICE inline algebra::traits::scalar_t<vector_t> eta(
    const vector_t &v) noexcept {

  using element_getter_t = algebra::traits::element_getter_t<vector_t>;

  return algebra::math::atanh(element_getter_t{}(v, 2) / norm(v));
}

/// Get a normalized version of the input vector
///
/// @tparam vector_t vector or column matrix type
///
/// @param v the input vector
///
/// @returns the normalized vector
template <concepts::vector vector_t>
ALGEBRA_HOST_DEVICE inline vector_t normalize(const vector_t &v) {

  using scalar_t = algebra::traits::scalar_t<vector_t>;

  return v * static_cast<scalar_t>(1.) / norm(v);
}

}  // namespace algebra::generic::math
