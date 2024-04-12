/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/storage/impl/fastor_matrix.hpp"
#include "algebra/type_traits.hpp"

// System include(s).
#include <cstddef>

namespace algebra {

namespace fastor {

/// size type for Fastor storage model
using size_type = std::size_t;
/// Array type used in the Fastor storage model
template <typename T, size_type N>
using storage_type = Fastor::Tensor<T, N>;
/// Matrix type used in the Fastor storage model
template <typename T, size_type ROWS, size_type COLS>
using matrix_type = algebra::fastor::Matrix<T, ROWS, COLS>;

/// 3-element "vector" type, using @c Fastor::Tensor
template <typename T>
using vector3 = storage_type<T, 3>;
/// Point in 3D space, using @c Fastor::Tensor
template <typename T>
using point3 = vector3<T>;
/// 2-element "vector" type, using @c Fastor::Tensor
template <typename T>
using vector2 = storage_type<T, 2>;
/// Point in 2D space, using @c Fastor::Tensor
template <typename T>
using point2 = vector2<T>;

}  // namespace fastor

namespace trait {

/// Type trait specializations
/// @{
template <typename T, std::size_t ROWS, std::size_t COLS>
struct dimensions<Fastor::Tensor<T, ROWS, COLS>> {

  using size_type = std::size_t;

  static constexpr size_type rows{ROWS};
  static constexpr size_type columns{COLS};
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct dimensions<algebra::fastor::Matrix<T, ROWS, COLS>> {

  using size_type = std::size_t;

  static constexpr size_type rows{ROWS};
  static constexpr size_type columns{COLS};
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct value<Fastor::Tensor<T, ROWS, COLS>> {
  using type = T;
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct value<algebra::fastor::Matrix<T, ROWS, COLS>> {
  using type = T;
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct index<Fastor::Tensor<T, ROWS, COLS>> {
  using type = std::size_t;
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct index<algebra::fastor::Matrix<T, ROWS, COLS>> {
  using type = std::size_t;
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct vector<Fastor::Tensor<T, ROWS, COLS>> {
  using type = Fastor::Tensor<T, ROWS>;
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct vector<algebra::fastor::Matrix<T, ROWS, COLS>> {
  using type = Fastor::Tensor<T, ROWS>;
};
/// @}

}  // namespace trait

}  // namespace algebra
