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
#include "algebra/math/vc_aos.hpp"
#include "algebra/storage/array.hpp"
#include "algebra/storage/vc_aos.hpp"

/// @name Operators on @c algebra::vc_aos types
/// @{

using algebra::cmath::operator*;
using algebra::cmath::operator-;
using algebra::cmath::operator+;

/// @}

namespace algebra {

namespace getter {

/// @name Getter functions on @c algebra::vc_aos::matrix_type
/// @{

using vc_aos::storage::block;
using vc_aos::storage::element;
using vc_aos::storage::set_block;
using vc_aos::storage::vector;

/// @}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::vc_aos::storage_type
/// @{

using generic::math::cross;
using generic::math::dot;
using generic::math::eta;
using generic::math::norm;
using generic::math::normalize;
using generic::math::perp;
using generic::math::phi;
using generic::math::theta;

/// @}

}  // namespace vector

namespace matrix {

/// @name Matrix functions on @c algebra::vc_aos::storage_type
/// @{

using generic::math::determinant;
using generic::math::identity;
using generic::math::inverse;
using generic::math::set_identity;
using generic::math::set_zero;
using generic::math::transpose;
using generic::math::zero;

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

namespace vc_aos {

/// @name generic based transforms on @c algebra::vc_aos
/// @{

template <concepts::value T>
using transform3 =
    generic::math::transform3<vc_aos::size_type, T, vc_aos::matrix_type,
                              vc_aos::storage_type>;

/// @}

}  // namespace vc_aos

namespace plugin {

/// Define the plugin types
/// @{
template <concepts::value V>
struct vc_aos_generic {
  /// Define scalar type
  using value_type = V;

  template <concepts::value T>
  using simd = T;

  using boolean = bool;
  using scalar = value_type;
  using size_type = algebra::vc_aos::size_type;
  using transform3D = algebra::vc_aos::transform3<value_type>;
  using point2D = algebra::vc_aos::point2<value_type>;
  using point3D = algebra::vc_aos::point3<value_type>;
  using vector3D = algebra::vc_aos::vector3<value_type>;

  template <std::size_t ROWS, std::size_t COLS>
  using matrix = algebra::vc_aos::matrix_type<value_type, ROWS, COLS>;
};
/// @}

}  // namespace plugin

}  // namespace algebra
