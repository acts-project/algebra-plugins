/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "algebra/type_traits.hpp"

// VecMem include(s).
#include <vecmem/containers/static_array.hpp>

// System include(s).
#include <cstddef>

namespace algebra {

namespace vecmem {

/// Array type used in the VecMem storage model
template <typename T, std::size_t N>
using storage_type = ::vecmem::static_array<T, N>;
/// Matrix type used in the VecMem storage model
template <typename T, std::size_t ROWS, std::size_t COLS>
using matrix_type = storage_type<storage_type<T, ROWS>, COLS>;

/// 3-element "vector" type, using @c vecmem::static_array
template <typename T>
using vector3 = storage_type<T, 3>;
/// Point in 3D space, using @c vecmem::static_array
template <typename T>
using point3 = vector3<T>;
/// 2-element "vector" type, using @c vecmem::static_array
template <typename T>
using vector2 = storage_type<T, 2>;
/// Point in 2D space, using @c vecmem::static_array
template <typename T>
using point2 = vector2<T>;

}  // namespace vecmem

namespace trait {

/// Type trait specializations
/// @{
template <typename T, std::size_t ROWS, std::size_t COLS>
struct index<::vecmem::static_array<::vecmem::static_array<T, ROWS>, COLS>> {
  using type = std::size_t;
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct dimensions<
    ::vecmem::static_array<::vecmem::static_array<T, ROWS>, COLS>> {

  using size_type =
      index_t<::vecmem::static_array<::vecmem::static_array<T, ROWS>, COLS>>;

  static constexpr size_type rows{ROWS};
  static constexpr size_type columns{COLS};
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct value<::vecmem::static_array<::vecmem::static_array<T, ROWS>, COLS>> {
  using type = T;
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct vector<::vecmem::static_array<::vecmem::static_array<T, ROWS>, COLS>> {
  using type = ::vecmem::static_array<T, ROWS>;
};
/// @}

}  // namespace trait

}  // namespace algebra
