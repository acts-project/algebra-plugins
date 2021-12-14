/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/storage/impl/vc_array4.hpp"

// System include(s).
#include <array>
#include <cstddef>

// Vc include(s).
#include <Vc/Vc>

namespace algebra::vc {

/// Array type used in the Vc storage model
template <typename T, std::size_t N>
using storage_type = Vc::SimdArray<T, N>;

/// 3-element "vector" type, using @c algebra::vc::array4
template <typename T>
using vector3 = array4<T>;
/// Point in 3D space, using @c algebra::vc::array4
template <typename T>
using point3 = vector3<T>;
/// 2-element "vector" type, using @c std::array
template <typename T>
using vector2 = std::array<T, 2>;
/// Point in 2D space, using @c std::array
template <typename T>
using point2 = vector2<T>;

}  // namespace algebra::vc
