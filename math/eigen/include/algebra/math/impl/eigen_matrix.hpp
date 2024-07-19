/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/qualifiers.hpp"
#include "algebra/storage/eigen.hpp"

// Eigen include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diagnostic push
#pragma nv_diag_suppress 20012
#endif  // __NVCC_DIAG_PRAGMA_SUPPORT__
#include <Eigen/Core>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diagnostic pop
#endif  // __NVCC_DIAG_PRAGMA_SUPPORT__

namespace algebra::eigen::math {

/// Operator getting a block of a const matrix
template <
    int ROWS, int COLS, class derived_type, typename size_type_1,
    typename size_type_2,
    std::enable_if_t<std::is_convertible<size_type_1, Eigen::Index>::value &&
                         std::is_convertible<size_type_2, Eigen::Index>::value,
                     bool> = true>
ALGEBRA_HOST_DEVICE auto block(const Eigen::MatrixBase<derived_type> &m,
                               size_type_1 row, size_type_2 col) {
  return m.template block<ROWS, COLS>(static_cast<Eigen::Index>(row),
                                      static_cast<Eigen::Index>(col));
}

/// Operator getting a block of a const matrix
template <
    int ROWS, int COLS, class derived_type, typename size_type_1,
    typename size_type_2,
    std::enable_if_t<std::is_convertible<size_type_1, Eigen::Index>::value &&
                         std::is_convertible<size_type_2, Eigen::Index>::value,
                     bool> = true>
ALGEBRA_HOST_DEVICE auto block(Eigen::MatrixBase<derived_type> &m,
                               size_type_1 row, size_type_2 col) {
  return m.template block<ROWS, COLS>(static_cast<Eigen::Index>(row),
                                      static_cast<Eigen::Index>(col));
}

/// Operator setting a block
template <
    typename derived_type1, typename derived_type2, typename size_type_1,
    typename size_type_2,
    std::enable_if_t<std::is_convertible<size_type_1, Eigen::Index>::value &&
                         std::is_convertible<size_type_2, Eigen::Index>::value,
                     bool> = true>
ALGEBRA_HOST_DEVICE void set_block(Eigen::MatrixBase<derived_type1> &m,
                                   const Eigen::MatrixBase<derived_type2> &b,
                                   size_type_1 row, size_type_2 col) {
  using block_t = Eigen::MatrixBase<derived_type2>;
  constexpr auto R{block_t::RowsAtCompileTime};
  constexpr auto C{block_t::ColsAtCompileTime};
  m.template block<R, C>(static_cast<Eigen::Index>(row),
                         static_cast<Eigen::Index>(col)) = b;
}

// Create zero matrix
template <typename matrix_t>
ALGEBRA_HOST_DEVICE inline matrix_t zero() {
  return matrix_t::Zero();
}

// Create identity matrix
template <typename matrix_t>
ALGEBRA_HOST_DEVICE inline matrix_t identity() {
  return matrix_t::Identity();
}

// Set input matrix as zero matrix
template <typename derived_type>
ALGEBRA_HOST_DEVICE inline void set_zero(Eigen::MatrixBase<derived_type> &m) {
  m.setZero();
}

// Set input matrix as identity matrix
template <typename derived_type>
ALGEBRA_HOST_DEVICE inline void set_identity(
    Eigen::MatrixBase<derived_type> &m) {
  m.setIdentity();
}

// Create transpose matrix
template <typename derived_type>
ALGEBRA_HOST_DEVICE inline matrix_type<
    typename Eigen::MatrixBase<derived_type>::value_type,
    Eigen::MatrixBase<derived_type>::ColsAtCompileTime,
    Eigen::MatrixBase<derived_type>::RowsAtCompileTime>
transpose(const Eigen::MatrixBase<derived_type> &m) {
  return m.transpose();
}

/// @returns the determinant of @param m
template <typename derived_type>
ALGEBRA_HOST_DEVICE inline typename Eigen::MatrixBase<derived_type>::value_type
determinant(const Eigen::MatrixBase<derived_type> &m) {
  return m.determinant();
}

/// @returns the inverse of @param m
template <typename derived_type>
ALGEBRA_HOST_DEVICE inline auto inverse(
    const Eigen::MatrixBase<derived_type> &m) {
  return m.inverse();
}

}  // namespace algebra::eigen::math
