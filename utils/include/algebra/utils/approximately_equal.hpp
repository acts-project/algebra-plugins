/** Algebra plugins, part of the ACTS project
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "algebra/concepts.hpp"
#include "algebra/math/common.hpp"
#include "algebra/qualifiers.hpp"

// System include(s)
#include <concepts>
#include <limits>

namespace algebra {

/// Compare two scalars according to a max relative error tolerance
/// @see
/// https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
///
/// @note This is by no means safe for all comparisons. Use with caution!
///
/// @param a first value
/// @param b second value
/// @param rel_error maximal relative error
///
/// @returns true if the two values are approximaterly equal
template <concepts::scalar scalar_t, concepts::value value_t>
requires std::convertible_to<scalar_t, value_t>
    ALGEBRA_HOST_DEVICE constexpr auto approx_equal(
        const scalar_t a, const scalar_t b,
        const value_t rel_error = 16.f *
                                  std::numeric_limits<value_t>::epsilon(),
        const value_t max_error = std::numeric_limits<value_t>::epsilon()) {

  if constexpr (std::integral<scalar_t>) {
    return a == b;
  } else {
    // Calculate the difference.
    const scalar_t diff{math::fabs(a - b)};

    // If the numbers are is close to zero
    if (diff <= max_error) {
      return true;
    }

    // Find the largest value to scale the epsilon
    const scalar_t largest = math::max(math::fabs(a), math::fabs(b));

    return (diff <= largest * rel_error);
  }
}

/// Elementwise compare two vectors according to a max relative error tolerance
///
/// @note This is by no means safe for all comparisons. Use with caution!
///
/// @param v1 first vector
/// @param v2 second vector
/// @param rel_error maximal relative error
///
/// @returns true if the two vectors are approximaterly equal
template <typename vector1_t, typename vector2_t>
requires((concepts::vector<vector1_t> || concepts::point<vector1_t>)&&(
             concepts::vector<vector2_t> ||
             concepts::point<vector2_t>)&&!concepts::simd_scalar<vector1_t> &&
         !concepts::simd_scalar<vector2_t>) ALGEBRA_HOST_DEVICE
    constexpr auto approx_equal(
        const vector1_t& v1, const vector2_t& v2,
        const algebra::traits::value_t<vector1_t> rel_error =
            16.f *
            std::numeric_limits<algebra::traits::value_t<vector1_t>>::epsilon(),
        const algebra::traits::value_t<vector1_t> max_error = std::
            numeric_limits<algebra::traits::value_t<vector1_t>>::epsilon()) {

  static_assert(std::same_as<algebra::traits::value_t<vector1_t>,
                             algebra::traits::value_t<vector2_t>>);

  using index_t = algebra::traits::index_t<vector1_t>;

  constexpr index_t size{algebra::traits::size<vector1_t>};

  bool ret{true};
  for (index_t i = 0; i < size; ++i) {
    ret &= ::algebra::approx_equal(v1[i], v2[i], rel_error, max_error);
  }

  return ret;
}

/// Elementwise compare two column matrices according to a max relative error
/// tolerance
///
/// @note This is by no means safe for all comparisons. Use with caution!
///
/// @param v1 first column matrix
/// @param v2 second column matrix
/// @param rel_error maximal relative error
///
/// @returns true if the two column matrices are approximaterly equal
template <concepts::column_matrix vector1_t, concepts::column_matrix vector2_t>
ALGEBRA_HOST_DEVICE constexpr auto approx_equal(
    const vector1_t& v1, const vector2_t& v2,
    const algebra::traits::value_t<vector1_t> rel_error =
        16.f *
        std::numeric_limits<algebra::traits::value_t<vector1_t>>::epsilon(),
    const algebra::traits::value_t<vector1_t> max_error =
        std::numeric_limits<algebra::traits::value_t<vector1_t>>::epsilon()) {

  static_assert(std::same_as<algebra::traits::value_t<vector1_t>,
                             algebra::traits::value_t<vector2_t>>);

  using index_t = algebra::traits::index_t<vector1_t>;
  using element_getter_t = algebra::traits::element_getter_t<vector1_t>;

  constexpr index_t rows{algebra::traits::rows<vector1_t>};

  bool ret{true};
  for (index_t i = 0; i < rows; ++i) {
    ret &= ::algebra::approx_equal(element_getter_t{}(v1, i, 0),
                                   element_getter_t{}(v2, i, 0), rel_error,
                                   max_error);
  }

  return ret;
}

/// Elementwise compare two matrices according to a max relative error tolerance
///
/// @note This is by no means safe for all comparisons. Use with caution!
///
/// @param m1 first matrix
/// @param m2 second matrix
/// @param rel_error maximal relative error
///
/// @returns true if the two matrices are approximaterly equal
template <concepts::matrix matrix1_t, concepts::matrix matrix2_t>
ALGEBRA_HOST_DEVICE constexpr auto approx_equal(
    const matrix1_t& m1, const matrix2_t& m2,
    const algebra::traits::value_t<matrix1_t> rel_error =
        16.f *
        std::numeric_limits<algebra::traits::value_t<matrix1_t>>::epsilon(),
    const algebra::traits::value_t<matrix1_t> max_error =
        std::numeric_limits<algebra::traits::value_t<matrix1_t>>::epsilon()) {

  static_assert(std::same_as<algebra::traits::value_t<matrix1_t>,
                             algebra::traits::value_t<matrix2_t>>);

  using index_t = algebra::traits::index_t<matrix1_t>;
  using element_getter_t = algebra::traits::element_getter_t<matrix1_t>;

  constexpr index_t rows{algebra::traits::rows<matrix1_t>};
  constexpr index_t columns{algebra::traits::columns<matrix1_t>};

  bool ret{true};
  for (index_t j = 0; j < columns; ++j) {
    for (index_t i = 0; i < rows; ++i) {
      ret &= ::algebra::approx_equal(element_getter_t{}(m1, i, j),
                                     element_getter_t{}(m2, i, j), rel_error,
                                     max_error);
    }
  }

  return ret;
}

/// Elementwise compare two transforms according to a max relative error
/// tolerance
///
/// @note This is by no means safe for all comparisons. Use with caution!
///
/// @param trf1 first transform
/// @param trf2 second transform
/// @param rel_error maximal relative error
///
/// @returns true if the two transforms are approximaterly equal
template <concepts::transform3D transform_t>
ALGEBRA_HOST_DEVICE constexpr auto approx_equal(
    const transform_t& trf1, const transform_t& trf2,
    const algebra::traits::value_t<typename transform_t::scalar_type>
        rel_error = 16.f * std::numeric_limits<algebra::traits::value_t<
                               typename transform_t::scalar_type>>::epsilon(),
    const algebra::traits::value_t<typename transform_t::scalar_type>
        max_error = std::numeric_limits<algebra::traits::value_t<
            typename transform_t::scalar_type>>::epsilon()) {

  return (::algebra::approx_equal(trf1.matrix(), trf2.matrix(), rel_error,
                                  max_error) &&
          ::algebra::approx_equal(trf1.matrix_inverse(), trf2.matrix_inverse(),
                                  rel_error, max_error));
}

}  // namespace algebra
