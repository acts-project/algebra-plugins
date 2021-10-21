/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/scalar.hpp"

// System include(s).
#include <array>

/// Main algebra namespace
namespace algebra {
/// @c std::array algebra definitions
namespace array {

/// Array type used in this storage model
template <typename T, std::size_t N>
using storage_type = std::array<T, N>;

/// 3-element "vector" type, using @c std::array
using vector3 = storage_type<scalar, 3>;
/// Point in 3D space, using @c std::array
using point3 = vector3;
/// Point in 2D space, using @c std::array
using point2 = storage_type<scalar, 2>;

}  // namespace array
}  // namespace algebra
