/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "algebra/storage/impl/cmath_getter.hpp"
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

/// Index
/// @{
template <typename T, std::size_t N>
struct index<std::array<T, N>> {
  using type = algebra::array::size_type;
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct index<std::array<std::array<T, ROWS>, COLS>> {
  using type = algebra::array::size_type;
};
/// @}

/// Dimension
/// @{
template <typename T, std::size_t N>
struct dimensions<std::array<T, N>> {

  using size_type = index_t<std::array<T, N>>;

  static constexpr size_type dim{1};
  static constexpr size_type rows{N};
  static constexpr size_type columns{1};
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct dimensions<std::array<std::array<T, ROWS>, COLS>> {

  using size_type = index_t<std::array<std::array<T, ROWS>, COLS>>;

  static constexpr size_type dim{2};
  static constexpr size_type rows{ROWS};
  static constexpr size_type columns{COLS};
};
/// @}

/// Value
/// @{
template <typename T, std::size_t N>
struct value<std::array<T, N>> {
  using type = T;
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct value<std::array<std::array<T, ROWS>, COLS>> {
  using type = T;
};
/// @}

/// Vector
/// @{
template <typename T, std::size_t N>
struct vector<std::array<T, N>> {

  template <typename other_T, std::size_t other_N>
  using other_type = std::array<other_T, other_N>;

  using type = other_type<T, N>;
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct vector<std::array<std::array<T, ROWS>, COLS>> {

  template <typename other_T, std::size_t other_N>
  using other_type = std::array<other_T, other_N>;

  using type = other_type<T, ROWS>;
};
/// @}

/// Matrix
/// @{
template <typename T, std::size_t ROWS, std::size_t COLS>
struct matrix<std::array<std::array<T, ROWS>, COLS>> {
  template <typename other_T, std::size_t other_ROWS, std::size_t other_COLS>
  using other_type = std::array<std::array<other_T, other_ROWS>, other_COLS>;

  using type = std::array<std::array<T, ROWS>, COLS>;
};

template <typename T, int N>
struct matrix<std::array<T, N>> {
  template <typename other_T, int other_ROWS, int other_COLS>
  using other_type = std::array<std::array<other_T, other_ROWS>, other_COLS>;

  using type = other_type<T, N, 1>;
};
/// @}

/// Elemet/Block Getter
/// @{
template <typename T, std::size_t N>
struct element_getter<std::array<T, N>> {
  using type = cmath::storage::element_getter;
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct element_getter<std::array<std::array<T, ROWS>, COLS>> {
  using type = cmath::storage::element_getter;
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct block_getter<std::array<std::array<T, ROWS>, COLS>> {
  using type = cmath::storage::block_getter;
};
/// @}

}  // namespace trait

}  // namespace algebra
