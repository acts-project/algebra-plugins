/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/impl/vc_aos_transform3.hpp"
#include "algebra/math/vc_soa.hpp"
#include "algebra/print.hpp"
#include "algebra/storage/vc_soa.hpp"

// System include(s).
#include <cassert>
#include <type_traits>

/// @name Operators on @c algebra::storage::vector types
/// @{

using algebra::storage::operator*;
using algebra::storage::operator/;
using algebra::storage::operator-;
using algebra::storage::operator+;

/// @}

namespace algebra {

namespace getter {

/// @name Getter functions on @c algebra::vc_soa types
/// @{

using vc_soa::storage::block;
using vc_soa::storage::element;
using vc_soa::storage::set_block;
using vc_soa::storage::vector;

/// Print the linear algebra types of this backend
using algebra::operator<<;

/// @}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::vc_soa types
/// @{

using vc_soa::math::cross;
using vc_soa::math::dot;
using vc_soa::math::eta;
using vc_soa::math::norm;
using vc_soa::math::normalize;
using vc_soa::math::perp;
using vc_soa::math::phi;
using vc_soa::math::theta;

/// @}

}  // namespace vector

// Produces clash with matrix typedefs in other plugins
namespace matrix {

using vc_soa::math::determinant;
using vc_soa::math::identity;
using vc_soa::math::inverse;
using vc_soa::math::set_identity;
using vc_soa::math::set_zero;
using vc_soa::math::transpose;
using vc_soa::math::zero;

}  // namespace matrix

namespace vc_soa {

/// @name Vc based transforms on @c algebra::vc_soa types
/// @{

template <concepts::value T>
using transform3 =
    algebra::vc_aos::math::transform3<algebra::vc_soa::storage_type, T>;

/// @}

}  // namespace vc_soa

}  // namespace algebra
