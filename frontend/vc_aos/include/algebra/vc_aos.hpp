/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/vc_aos.hpp"
#include "algebra/print.hpp"
#include "algebra/storage/vc_aos.hpp"

// System include(s).
#include <cassert>
#include <type_traits>

/// Print the linear algebra types of this backend
using algebra::operator<<;

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

/// @name Vector functions on @c algebra::vc_aos types
/// @{

// Vc array specific
using vc_aos::math::cross;
using vc_aos::math::dot;
using vc_aos::math::eta;
using vc_aos::math::norm;
using vc_aos::math::normalize;
using vc_aos::math::perp;
using vc_aos::math::phi;
using vc_aos::math::theta;

/// @}

}  // namespace vector

namespace matrix {

/// @name Matrix functions on @c algebra::vc_aos types
/// @{

using vc_aos::math::determinant;
using vc_aos::math::identity;
using vc_aos::math::inverse;
using vc_aos::math::set_identity;
using vc_aos::math::set_zero;
using vc_aos::math::transpose;
using vc_aos::math::zero;

/// @}

}  // namespace matrix

namespace vc_aos {

/// @name Vc based transforms on @c algebra::vc_aos::storage_type
/// @{

template <concepts::value T>
using transform3 = math::transform3<vc_aos::storage_type, T>;

/// @}

}  // namespace vc_aos

}  // namespace algebra
