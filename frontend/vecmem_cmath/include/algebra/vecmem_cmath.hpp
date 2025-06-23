/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/cmath.hpp"
#include "algebra/math/generic.hpp"
#include "algebra/storage/vecmem.hpp"

namespace algebra {

/// @name Operators on @c algebra::vecmem::storage_type
/// @{

using algebra::cmath::operator*;
using algebra::cmath::operator-;
using algebra::cmath::operator+;

/// @}

namespace getter {

/// @name Getter functions on @c algebra::vecmem::matrix_type
/// @{

using cmath::storage::block;
using cmath::storage::element;
using cmath::storage::set_block;
using cmath::storage::vector;

/// @}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::vecmem::storage_type
/// @{

// vecmem::static_array specific implementations
using cmath::dot;
using cmath::normalize;

// generic implementations
using generic::math::cross;
using generic::math::eta;
using generic::math::norm;
using generic::math::perp;
using generic::math::phi;
using generic::math::theta;

/// @}

}  // namespace vector

// Use special algorithms for 4 dimensional matrices
namespace generic {

// Determinant algorithms
template <concepts::scalar T, auto ROWS, auto COLS>
struct determinant_selector<4, vecmem::matrix_type<T, ROWS, COLS>> {
  using type =
      matrix::determinant::hard_coded<vecmem::matrix_type<T, ROWS, COLS>,
                                      vecmem::element_getter>;
};

// Inversion algorithms
template <concepts::scalar T, auto ROWS, auto COLS>
struct inversion_selector<4, vecmem::matrix_type<T, ROWS, COLS>> {
  using type = matrix::inverse::hard_coded<vecmem::matrix_type<T, ROWS, COLS>,
                                           vecmem::element_getter>;
};

}  // namespace generic

namespace matrix {

/// @name Matrix functions on @c algebra::vecmem::storage_type
/// @{

using cmath::identity;
using cmath::set_identity;
using cmath::set_zero;
using cmath::zero;

using generic::math::determinant;
using generic::math::inverse;
using generic::math::transpose;

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

namespace vecmem {

/// @name cmath based transforms on @c algebra::vecmem
/// @{

template <concepts::scalar T>
using transform3 =
    generic::math::transform3<std::size_t, T, vecmem::matrix_type,
                              vecmem::storage_type>;

/// @}

}  // namespace vecmem

namespace plugin {

/// Define the plugin types
/// @{
template <concepts::value V>
struct vecmem {
  /// Define scalar type
  using value_type = V;

  template <concepts::value T>
  using simd = T;

  using boolean = bool;
  using scalar = value_type;
  using size_type = algebra::vecmem::size_type;
  using transform3D = algebra::vecmem::transform3<value_type>;
  using point2D = algebra::vecmem::point2<value_type>;
  using point3D = algebra::vecmem::point3<value_type>;
  using vector3D = algebra::vecmem::vector3<value_type>;

  template <std::size_t ROWS, std::size_t COLS>
  using matrix = algebra::vecmem::matrix_type<value_type, ROWS, COLS>;
};
/// @}

}  // namespace plugin

}  // namespace algebra
