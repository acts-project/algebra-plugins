/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "algebra/storage/impl/smatrix_getter.hpp"
#include "algebra/type_traits.hpp"

// ROOT/Smatrix include(s).
#include <Math/SMatrix.h>
#include <Math/SVector.h>

// System include(s).
#include <cstddef>

namespace algebra {

namespace smatrix {

/// size type for SMatrix storage model
using size_type = unsigned int;
/// Array type used in the SMatrix storage model
template <typename T, size_type N>
using storage_type = ROOT::Math::SVector<T, N>;
/// Matrix type used in the SMatrix storage model
template <typename T, size_type ROWS, size_type COLS>
using matrix_type = ROOT::Math::SMatrix<T, ROWS, COLS>;

/// 3-element "vector" type, using @c ROOT::Math::SVector
template <typename T>
using vector3 = storage_type<T, 3>;
/// Point in 3D space, using @c ROOT::Math::SVector
template <typename T>
using point3 = vector3<T>;
/// 2-element "vector" type, using @c ROOT::Math::SVector
template <typename T>
using vector2 = storage_type<T, 2>;
/// Point in 2D space, using @c ROOT::Math::SVector
template <typename T>
using point2 = vector2<T>;

}  // namespace smatrix

namespace trait {

/// Type trait specializations
/// @{

/// Index
/// @{
template <typename T, std::size_t ROWS, std::size_t COLS>
struct index<ROOT::Math::SMatrix<T, ROWS, COLS>> {
  using type = unsigned int;
};
/// @}

/// Dimension
/// @{
template <typename T, std::size_t ROWS, std::size_t COLS>
struct dimensions<ROOT::Math::SMatrix<T, ROWS, COLS>> {

  using size_type = std::size_t;

  static constexpr size_type dim{2};
  static constexpr size_type rows{ROWS};
  static constexpr size_type columns{COLS};
};
/// @}

/// Value
/// @{
template <typename T, std::size_t ROWS, std::size_t COLS>
struct value<ROOT::Math::SMatrix<T, ROWS, COLS>> {
  using type = T;
};
/// @}

/// Vector
/// @{
template <typename T, std::size_t ROWS, std::size_t COLS>
struct vector<ROOT::Math::SMatrix<T, ROWS, COLS>> {

  template <typename other_T, std::size_t other_N>
  using other_type = ROOT::Math::SVector<other_T, other_N>;

  using type = other_type<T, ROWS>;
};

template <typename T, std::size_t N>
struct vector<ROOT::Math::SVector<T, N>> {
  using type = ROOT::Math::SVector<T, N>;

  template <typename other_T, std::size_t other_N>
  using other_type = ROOT::Math::SVector<other_T, other_N>;

  using type = other_type<T, N>;
};
/// @}

/// Matrix
/// @{
template <typename T, std::size_t ROWS, std::size_t COLS>
struct matrix<ROOT::Math::SMatrix<T, ROWS, COLS>> {
  template <typename other_T, std::size_t other_ROWS, std::size_t other_COLS>
  using other_type = ROOT::Math::SMatrix<other_T, other_ROWS, other_COLS>;

  using type = ROOT::Math::SMatrix<T, ROWS, COLS>;
};

template <typename T, int N>
struct matrix<ROOT::Math::SVector<T, N>> {
  template <typename other_T, int other_ROWS, int other_COLS>
  using other_type = ROOT::Math::SMatrix<other_T, other_ROWS, other_COLS>;

  using type = other_type<T, N, 1>;
};
/// @}

/// Elemet/Block Getter
/// @{
template <typename T, std::size_t ROWS, std::size_t COLS>
struct element_getter<ROOT::Math::SMatrix<T, ROWS, COLS>> {
  using type = smatrix::storage::element_getter;
};

template <typename T, std::size_t N>
struct element_getter<ROOT::Math::SVector<T, N>> {
  using type = smatrix::storage::element_getter;
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct block_getter<ROOT::Math::SMatrix<T, ROWS, COLS>> {
  using type = smatrix::storage::block_getter;
};
/// @}

}  // namespace trait

}  // namespace algebra
