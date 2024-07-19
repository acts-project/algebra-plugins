/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/storage/impl/eigen_array.hpp"
#include "algebra/type_traits.hpp"

// System include(s).
#include <cstddef>

namespace algebra {

namespace eigen {

/// size type for Eigen storage model
using size_type = int;
/// Array type used in the Eigen storage model
template <typename T, size_type N>
using storage_type = array<T, N>;
/// Vector type used in the Eigen storage model
template <typename T, std::size_t N>
using vector_type = storage_type<T, N>;
/// Matrix type used in the Eigen storage model
/// If the number of rows is 1, make it RowMajor
template <typename T, size_type ROWS, size_type COLS>
using matrix_type = Eigen::Matrix<T, ROWS, COLS, (ROWS == 1), ROWS, COLS>;

/// 3-element "vector" type, using @c algebra::eigen::array
template <typename T>
using vector3 = storage_type<T, 3>;
/// Point in 3D space, using @c algebra::eigen::array
template <typename T>
using point3 = vector3<T>;
/// 2-element "vector" type, using @c algebra::eigen::array
template <typename T>
using vector2 = storage_type<T, 2>;
/// Point in 2D space, using @c algebra::eigen::array
template <typename T>
using point2 = vector2<T>;

}  // namespace eigen

namespace trait {

/// Type trait specializations
/// @{
template <typename derived_t>
struct index<Eigen::MatrixBase<derived_t>> {
  using type = algebra::eigen::size_type;
};

template <typename derived_t>
struct dimensions<Eigen::MatrixBase<derived_t>> {

  using size_type = index_t<Eigen::MatrixBase<derived_t>>;

  static constexpr size_type rows{
      Eigen::MatrixBase<derived_t>::RowsAtCompileTime};
  static constexpr size_type columns{
      Eigen::MatrixBase<derived_t>::ColsAtCompileTime};
};

template <typename derived_t>
struct value<Eigen::MatrixBase<derived_t>> {
  using type = typename Eigen::MatrixBase<derived_t>::value_type;
};

template <typename derived_t>
struct vector<Eigen::MatrixBase<derived_t>> {
  using type = algebra::eigen::array<value_t<Eigen::MatrixBase<derived_t>>,
                                     rows<Eigen::MatrixBase<derived_t>>>;
};

template <typename T, int ROWS, int COLS>
struct index<Eigen::Matrix<T, ROWS, COLS, (ROWS == 1), ROWS, COLS>> {
  using type = algebra::eigen::size_type;
};

template <typename T, int ROWS, int COLS>
struct dimensions<Eigen::Matrix<T, ROWS, COLS, (ROWS == 1), ROWS, COLS>> {

  using size_type =
      index_t<Eigen::Matrix<T, ROWS, COLS, (ROWS == 1), ROWS, COLS>>;

  static constexpr size_type rows{ROWS};
  static constexpr size_type columns{COLS};
};

template <typename T, int ROWS, int COLS>
struct value<Eigen::Matrix<T, ROWS, COLS, (ROWS == 1), ROWS, COLS>> {
  using type = T;
};

template <typename T, int ROWS, int COLS>
struct vector<Eigen::Matrix<T, ROWS, COLS, (ROWS == 1), ROWS, COLS>> {
  using type = algebra::eigen::array<T, ROWS>;
};
/// @}

}  // namespace trait

}  // namespace algebra
