/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "algebra/type_traits.hpp"

// System include(s).
#include <array>
#include <cstddef>

namespace algebra {

namespace array {

/// size type for Array storage model
using size_type = std::size_t;
/// Array type used in the Array storage model
template <typename T, size_type N>
using storage_type = std::array<T, N>;
/// Vector type used in the Array storage model
template <typename T, std::size_t N>
using vector_type = storage_type<T, N>;
/// Matrix type used in the Array storage model
template <typename T, size_type ROWS, size_type COLS>
using matrix_type = storage_type<storage_type<T, ROWS>, COLS>;

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

}  // namespace array

namespace trait {

/// Type trait specializations
/// @{
template <typename T, std::size_t ROWS, std::size_t COLS>
struct index<std::array<std::array<T, ROWS>, COLS>> {
  using type = algebra::array::size_type;
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct dimensions<std::array<std::array<T, ROWS>, COLS>> {

  using size_type = index_t<std::array<std::array<T, ROWS>, COLS>>;

  static constexpr size_type rows{ROWS};
  static constexpr size_type columns{COLS};
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct value<std::array<std::array<T, ROWS>, COLS>> {
  using type = T;
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct vector<std::array<std::array<T, ROWS>, COLS>> {
  using type = std::array<T, ROWS>;
};
/// @}

}  // namespace trait

}  // namespace algebra
