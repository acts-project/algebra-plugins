/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s).
#include <array>
#include <cstddef>

namespace algebra::array {

/// Array type used in the Array storage model
template <typename T, std::size_t N>
using storage_type = std::array<T, N>;

/// 3-element "vector" type, using @c std::array
template <typename T>
using vector3 = storage_type<T, 3>;
/// Point in 3D space, using @c std::array
template <typename T>
using point3 = vector3<T>;
/// 2-element "vector" type, using @c std::array
template <typename T>
using vector2 = storage_type<T, 2>;
/// Point in 2D space, using @c std::array
template <typename T>
using point2 = vector2<T>;
/// Matrix, using @c std::array
template <typename T, std::size_t ROWS, std::size_t COLS>
using matrix_type = std::array<std::array<T, ROWS>, COLS>;

}  // namespace algebra::array
