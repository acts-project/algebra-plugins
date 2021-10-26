/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/scalar.hpp"
#include "algebra/storage/impl/vc_array.hpp"

// System include(s).
#include <array>
#include <cstddef>

namespace algebra::vc {

/// Array type used in the Vc storage model
template <typename T, std::size_t N>
using storage_type = array<T, N>;

/// 3-element "vector" type, using @c algebra::vc::array
using vector3 = storage_type<scalar, 3>;
/// Point in 3D space, using @c algebra::vc::array
using point3 = vector3;
/// 2-element "vector" type, using @c std::array
using vector2 = std::array<scalar, 2>;
/// Point in 2D space, using @c std::array
using point2 = vector2;

}  // namespace algebra::vc
