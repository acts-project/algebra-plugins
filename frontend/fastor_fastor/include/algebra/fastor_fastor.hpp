/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/fastor.hpp"
#include "algebra/print.hpp"
#include "algebra/storage/fastor.hpp"

// Fastor include(s).
#ifdef _MSC_VER
#pragma warning(disable : 4244 4701 4702)
#endif  // MSVC
#include <Fastor/Fastor.h>
#ifdef _MSC_VER
#pragma warning(default : 4244 4701 4702)
#endif  // MSVC

/// Print the linear algebra types of this backend
using algebra::operator<<;

namespace algebra {

namespace getter {

/// @name Getter functions on @c algebra::fastor::matrix_type
/// @{

using fastor::storage::block;
using fastor::storage::element;
using fastor::storage::set_block;
using fastor::storage::vector;

/// @}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::fastor::storage_type
/// @{

using fastor::math::cross;
using fastor::math::dot;
using fastor::math::eta;
using fastor::math::norm;
using fastor::math::normalize;
using fastor::math::perp;
using fastor::math::phi;
using fastor::math::theta;

/// @}

}  // namespace vector

namespace matrix {

/// @name Matrix functions on @c algebra::fastor::storage_type
/// @{

using fastor::math::determinant;
using fastor::math::identity;
using fastor::math::inverse;
using fastor::math::set_identity;
using fastor::math::set_zero;
using fastor::math::transpose;
using fastor::math::zero;

/// @}

}  // namespace matrix

namespace fastor {

/// @name Transform on @c algebra::fastor::storage_type
/// @{

template <concepts::scalar T>
using transform3 = math::transform3<T>;

/// @}

}  // namespace fastor

}  // namespace algebra
