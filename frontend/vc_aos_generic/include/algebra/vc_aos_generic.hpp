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

/// @}

}  // namespace matrix

namespace vc_aos {

/// @name generic based transforms on @c algebra::vc_aos
/// @{

template <typename T>
using transform3 =
    generic::math::transform3<vc_aos::size_type, T, vc_aos::matrix_type,
                              vc_aos::storage_type>;

/// @}

}  // namespace vc_aos

}  // namespace algebra
