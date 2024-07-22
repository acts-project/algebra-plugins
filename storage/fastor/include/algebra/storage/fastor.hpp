/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/storage/impl/fastor_getter.hpp"
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

/// Index
/// @{
template <typename T, std::size_t ROWS, std::size_t COLS>
struct index<Fastor::Tensor<T, ROWS, COLS>> {
  using type = std::size_t;
};
/// @}

/// Dimension
/// @{
template <typename T, std::size_t ROWS, std::size_t COLS>
struct dimensions<Fastor::Tensor<T, ROWS, COLS>> {

  using size_type = std::size_t;

  static constexpr size_type dim{2};
  static constexpr size_type rows{ROWS};
  static constexpr size_type columns{COLS};
};
/// @}

/// Value
/// @{
template <typename T, std::size_t ROWS, std::size_t COLS>
struct value<Fastor::Tensor<T, ROWS, COLS>> {
  using type = T;
};
/// @}

/// Vector
/// @{
template <typename T, std::size_t N>
struct vector<Fastor::Tensor<T, N>> {

  template <typename other_T, std::size_t other_N>
  using other_type = Fastor::Tensor<other_T, other_ROWS>;

  using type = other_type<T, N>;
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct vector<Fastor::Tensor<T, ROWS, COLS>> {

  template <typename other_T, std::size_t other_N>
  using other_type = Fastor::Tensor<other_T, other_ROWS>;

  using type = other_type<T, ROWS>;
};
/// @}

/// Matrix
/// @{
template <typename T, std::size_t ROWS, std::size_t COLS>
struct matrix<Fastor::Tensor<T, ROWS, COLS>> {
  template <typename other_T, std::size_t other_ROWS, std::size_t other_COLS>
  using other_type = Fastor::Tensor<other_T, other_ROWS, other_COLS>;

  using type = Fastor::Tensor<T, ROWS, COLS>;
};

template <typename T, int N>
struct matrix<Fastor::Tensor<T, N>> {
  template <typename other_T, int other_ROWS, int other_COLS>
  using other_type = Fastor::Tensor<other_T, other_ROWS, other_COLS>;

  using type = other_type<T, N, 1>;
};
/// @}

/// Element/Block Getter
/// @{
template <typename T, std::size_t ROWS, std::size_t COLS>
struct element_getter<Fastor::Tensor<T, ROWS, COLS>> {
  using type = fastor::storage::element_getter;
};

template <typename T, std::size_t N>
struct element_getter<Fastor::Tensor<T, N>> {
  using type = fastor::storage::element_getter;
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct block_getter<Fastor::Tensor<T, ROWS, COLS>> {
  using type = fastor::storage::block_getter;
};
/// @}

}  // namespace trait

}  // namespace algebra
