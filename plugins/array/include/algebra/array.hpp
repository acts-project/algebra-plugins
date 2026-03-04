/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/generic.hpp"
#include "algebra/impl/array_getter.hpp"
#include "algebra/impl/array_matrix.hpp"
#include "algebra/impl/array_operators.hpp"
#include "algebra/impl/array_types.hpp"
#include "algebra/impl/array_vector.hpp"
#include "algebra/impl/generic_matrix.hpp"

namespace algebra {

/// @name Operators on @c algebra::array::storage_type
/// @{

using algebra::array::operator*;
using algebra::array::operator-;
using algebra::array::operator+;

/// @}

namespace getter {

/// @name Getter functions on @c algebra::array::storage_type
/// @{

using array::storage::block;
using array::storage::element;
using array::storage::set_block;
using array::storage::vector;

/// @}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::array::storage_type
/// @{

// array specific implementations
using array::dot;
using array::normalize;

// generic implementations
using array::cross;
using array::eta;
using array::norm;
using array::perp;
using array::phi;
using array::theta;

/// @}

}  // namespace vector

// Use special algorithms for 4 dimensional matrices
namespace generic {

// Determinant algorithms
template <concepts::scalar T, auto ROWS, auto COLS>
struct determinant_selector<4, array::matrix_type<T, ROWS, COLS>,
                            array::element_getter> {
  using type =
      matrix::determinant::hard_coded<array::matrix_type<T, ROWS, COLS>>;
};

// Inversion algorithms
template <concepts::scalar T, auto ROWS, auto COLS>
struct inversion_selector<4, array::matrix_type<T, ROWS, COLS>,
                          array::element_getter> {
  using type = matrix::inverse::hard_coded<array::matrix_type<T, ROWS, COLS>>;
};

}  // namespace generic

namespace matrix {

/// @name Matrix functions on @c algebra::array::storage_type
/// @{

using array::identity;
using array::set_identity;
using array::set_zero;
using array::zero;

// Uses generic implementation in the background
using array::determinant;
using array::inverse;
using array::transpose;

using generic::math::set_inplace_product_left;
using generic::math::set_inplace_product_left_transpose;
using generic::math::set_inplace_product_right;
using generic::math::set_inplace_product_right_transpose;
using generic::math::set_product;
using generic::math::set_product_left_transpose;
using generic::math::set_product_right_transpose;
using generic::math::transposed_product;

/// @}

}  // namespace matrix

namespace array {

/// @name array based transforms on @c algebra::array
/// @{

template <concepts::scalar T>
using transform3 =
    generic::math::transform3<array::index_type, T, array::matrix_type,
                              array::storage_type>;

/// @}

}  // namespace array

namespace plugin {

/// Define the plugin types
/// @{
template <concepts::value V>
struct array {
  /// Define scalar type
  using value_type = V;

  template <concepts::element T>
  using simd = T;

  using boolean = bool;
  using scalar = value_type;
  using index_type = algebra::array::index_type;
  using transform3D = algebra::array::transform3<value_type>;
  using point2D = algebra::array::point2<value_type>;
  using point3D = algebra::array::point3<value_type>;
  using vector2D = algebra::array::vector2<value_type>;
  using vector3D = algebra::array::vector3<value_type>;

  template <std::size_t ROWS, std::size_t COLS>
  using matrix = algebra::array::matrix_type<value_type, ROWS, COLS>;
};
/// @}

}  // namespace plugin

}  // namespace algebra
