/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/storage/impl/eigen_array.hpp"
#include "algebra/storage/impl/eigen_getter.hpp"
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
template <typename T, size_type N>
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

/// Index
/// @{
template <typename T, int N>
struct index<algebra::eigen::array<T, N>> {
  using type = algebra::eigen::size_type;
};

template <typename derived_t>
struct index<Eigen::MatrixBase<derived_t>> {
  using type = algebra::eigen::size_type;
};

template <typename T, int ROWS, int COLS>
struct index<Eigen::Matrix<T, ROWS, COLS, (ROWS == 1), ROWS, COLS>> {
  using type = algebra::eigen::size_type;
};

template <typename T, int ROWS, int COLS, int bROWS, int bCOLS>
struct index<Eigen::Block<Eigen::Matrix<T, ROWS, COLS, (ROWS == 1), ROWS, COLS>,
                          bROWS, bCOLS, false>> {
  using type = algebra::eigen::size_type;
};
/// @}

/// Dimensions
/// @{
template <typename T, int N>
struct dimensions<algebra::eigen::array<T, N>> {

  using size_type = index_t<algebra::eigen::array<T, N>>;

  static constexpr size_type dim{1};
  static constexpr size_type rows{N};
  static constexpr size_type columns{1};
};

template <typename derived_t>
struct dimensions<Eigen::MatrixBase<derived_t>> {

  using size_type = index_t<Eigen::MatrixBase<derived_t>>;

  static constexpr size_type dim{2};
  static constexpr size_type rows{
      Eigen::MatrixBase<derived_t>::RowsAtCompileTime};
  static constexpr size_type columns{
      Eigen::MatrixBase<derived_t>::ColsAtCompileTime};
};

template <typename T, int ROWS, int COLS>
struct dimensions<Eigen::Matrix<T, ROWS, COLS, (ROWS == 1), ROWS, COLS>> {

  using size_type =
      index_t<Eigen::Matrix<T, ROWS, COLS, (ROWS == 1), ROWS, COLS>>;

  static constexpr size_type dim{2};
  static constexpr size_type rows{ROWS};
  static constexpr size_type columns{COLS};
};

template <typename T, int ROWS, int COLS, int bROWS, int bCOLS>
struct dimensions<
    Eigen::Block<Eigen::Matrix<T, ROWS, COLS, (ROWS == 1), ROWS, COLS>, bROWS,
                 bCOLS, false>> {

  using size_type = index_t<
      Eigen::Block<Eigen::Matrix<T, ROWS, COLS, (ROWS == 1), ROWS, COLS>, bROWS,
                   bCOLS, false>>;

  static constexpr size_type dim{2};
  static constexpr size_type rows{ROWS};
  static constexpr size_type columns{COLS};
};
/// @}

/// Value
/// @{
template <typename T, int N>
struct value<algebra::eigen::array<T, N>> {
  using type = T;
};

template <typename derived_t>
struct value<Eigen::MatrixBase<derived_t>> {
  using type = typename Eigen::MatrixBase<derived_t>::value_type;
};

template <typename T, int ROWS, int COLS>
struct value<Eigen::Matrix<T, ROWS, COLS, (ROWS == 1), ROWS, COLS>> {
  using type = T;
};

template <typename T, int ROWS, int COLS, int bROWS, int bCOLS>
struct value<Eigen::Block<Eigen::Matrix<T, ROWS, COLS, (ROWS == 1), ROWS, COLS>,
                          bROWS, bCOLS, false>> {
  using type = T;
};
/// @}

/// Vector
/// @{
template <typename T, int N>
struct vector<algebra::eigen::array<T, N>> {

  template <typename other_T, int other_N>
  using other_type = algebra::eigen::array<other_T, other_N>;

  using type = other_type<T, N>;
};

template <typename derived_t>
struct vector<Eigen::MatrixBase<derived_t>> {

  template <typename other_T, int other_N>
  using other_type = algebra::eigen::array<other_T, other_N>;

  using type = other_type<value_t<Eigen::MatrixBase<derived_t>>,
                          rows<Eigen::MatrixBase<derived_t>>>;
};

template <typename T, int ROWS, int COLS>
struct vector<Eigen::Matrix<T, ROWS, COLS, (ROWS == 1), ROWS, COLS>> {

  template <typename other_T, int other_N>
  using other_type = algebra::eigen::array<other_T, other_N>;

  using type = other_type<T, ROWS>;
};
/// @}

/// Matrix
/// @{
template <typename derived_t>
struct matrix<Eigen::MatrixBase<derived_t>> {
  template <typename other_T, int other_ROWS, int other_COLS>
  using other_type =
      eigen::matrix_type<typename Eigen::MatrixBase<derived_t>::value_type,
                         Eigen::MatrixBase<derived_t>::RowsAtCompileTime,
                         Eigen::MatrixBase<derived_t>::ColsAtCompileTime>;

  using type = Eigen::MatrixBase<derived_t>;
};

template <typename T, int ROWS, int COLS>
struct matrix<Eigen::Matrix<T, ROWS, COLS, (ROWS == 1), ROWS, COLS>> {
  template <typename other_T, int other_ROWS, int other_COLS>
  using other_type = eigen::matrix_type<other_T, other_ROWS, other_COLS>;

  using type = Eigen::Matrix<T, ROWS, COLS, (ROWS == 1), ROWS, COLS>;
};

template <typename T, int N>
struct matrix<algebra::eigen::array<T, N>> {
  template <typename other_T, int other_ROWS, int other_COLS>
  using other_type = eigen::matrix_type<other_T, other_ROWS, other_COLS>;

  using type = other_type<T, N, 1>;
};

template <typename T, int ROWS, int COLS, int bROWS, int bCOLS>
struct matrix<
    Eigen::Block<Eigen::Matrix<T, ROWS, COLS, (ROWS == 1), ROWS, COLS>, bROWS,
                 bCOLS, false>> {
  template <typename other_T, int other_ROWS, int other_COLS>
  using other_type = eigen::matrix_type<other_T, other_ROWS, other_COLS>;

  using type =
      Eigen::Block<Eigen::Matrix<T, ROWS, COLS, (ROWS == 1), ROWS, COLS>, bROWS,
                   bCOLS, false>;
};
/// @}

/// Element getter
/// @{
template <typename T, int N>
struct element_getter<algebra::eigen::array<T, N>> {
  using type = eigen::storage::element_getter;
};

template <typename derived_t>
struct element_getter<Eigen::MatrixBase<derived_t>> {
  using type = eigen::storage::element_getter;
};

template <typename T, int ROWS, int COLS>
struct element_getter<Eigen::Matrix<T, ROWS, COLS, (ROWS == 1), ROWS, COLS>> {
  using type = eigen::storage::element_getter;
};

template <typename T, int ROWS, int COLS, int bROWS, int bCOLS>
struct element_getter<
    Eigen::Block<Eigen::Matrix<T, ROWS, COLS, (ROWS == 1), ROWS, COLS>, bROWS,
                 bCOLS, false>> {
  using type = eigen::storage::element_getter;
};
/// @}

/// Block getter
/// @{
template <typename derived_t>
struct block_getter<Eigen::MatrixBase<derived_t>> {
  using type = eigen::storage::block_getter;
};

template <typename T, int ROWS, int COLS>
struct block_getter<Eigen::Matrix<T, ROWS, COLS, (ROWS == 1), ROWS, COLS>> {
  using type = eigen::storage::block_getter;
};
/// @}
/// @}

}  // namespace trait

}  // namespace algebra
