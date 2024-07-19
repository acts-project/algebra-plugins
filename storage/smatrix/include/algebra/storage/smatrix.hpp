/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
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
template <typename T, std::size_t ROWS, std::size_t COLS>
struct dimensions<ROOT::Math::SMatrix<T, ROWS, COLS>> {

  using size_type = std::size_t;

  static constexpr size_type rows{ROWS};
  static constexpr size_type columns{COLS};
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct value<ROOT::Math::SMatrix<T, ROWS, COLS>> {
  using type = T;
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct index<ROOT::Math::SMatrix<T, ROWS, COLS>> {
  using type = unsigned int;
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct vector<ROOT::Math::SMatrix<T, ROWS, COLS>> {
  using type = ROOT::Math::SVector<T, N>;
};
/// @}

}  // namespace trait

}  // namespace algebra
