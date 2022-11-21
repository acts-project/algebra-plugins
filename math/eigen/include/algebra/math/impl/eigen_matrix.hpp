/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
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
#include <Eigen/Core>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC

namespace algebra::eigen::matrix {

/// "Matrix actor", assuming an Eigen matrix
template <typename scalar_t>
struct actor {

  /// Size type
  using size_ty = int;

  /// Scalar type
  using scalar_type = scalar_t;

  /// 2D matrix type
  template <int ROWS, int COLS>
  using matrix_type = algebra::eigen::matrix_type<scalar_t, ROWS, COLS>;

  /// Array type
  template <size_ty N>
  using array_type = storage_type<scalar_type, N>;

  /// 3-element "vector" type
  using vector3 = array_type<3>;

  /// Operator getting a reference to one element of a non-const matrix
  template <class Derived>
  ALGEBRA_HOST_DEVICE inline constexpr scalar_t &element(
      Eigen::MatrixBase<Derived> &m, int row, int col) const {
    return m(row, col);
  }

  /// Operator getting one value of a const matrix
  template <class Derived>
  ALGEBRA_HOST_DEVICE inline scalar_t constexpr element(
      const Eigen::MatrixBase<Derived> &m, int row, int col) const {
    return m(row, col);
  }

  /// Operator getting a reference to one element of a non-const matrix
  template <typename Derived>
  ALGEBRA_HOST_DEVICE inline scalar_t constexpr element(
      Eigen::MatrixBase<Derived> &&m, int row, int col) const {
    return m(row, col);
  }

  /// Operator getting a block of a const matrix
  template <
      int ROWS, int COLS, class Derived,
      std::enable_if_t<ROWS <= Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                       bool> = true,
      std::enable_if_t<COLS <= Eigen::MatrixBase<Derived>::ColsAtCompileTime,
                       bool> = true>
  ALGEBRA_HOST_DEVICE inline constexpr matrix_type<ROWS, COLS> block(
      const Eigen::MatrixBase<Derived> &m, int row, int col) const {
    return m.template block<ROWS, COLS>(row, col);
  }

  /// Operator setting a block
  template <class DerivedA, class DerivedB,
            std::enable_if_t<Eigen::MatrixBase<DerivedB>::RowsAtCompileTime <=
                                 Eigen::MatrixBase<DerivedA>::RowsAtCompileTime,
                             bool> = true,
            std::enable_if_t<Eigen::MatrixBase<DerivedB>::ColsAtCompileTime <=
                                 Eigen::MatrixBase<DerivedA>::ColsAtCompileTime,
                             bool> = true>
  ALGEBRA_HOST_DEVICE inline constexpr void set_block(
      Eigen::MatrixBase<DerivedA> &m, const Eigen::MatrixBase<DerivedB> &b,
      int row, int col) const {
    m.template block<Eigen::MatrixBase<DerivedB>::RowsAtCompileTime,
                     Eigen::MatrixBase<DerivedB>::ColsAtCompileTime>(row, col) =
        b;
  }

  /// Create zero matrix
  template <int ROWS, int COLS>
  ALGEBRA_HOST_DEVICE inline constexpr matrix_type<ROWS, COLS> zero() const {
    return matrix_type<ROWS, COLS>::Zero();
  }

  /// Create identity matrix
  template <int ROWS, int COLS>
  ALGEBRA_HOST_DEVICE inline constexpr matrix_type<ROWS, COLS> identity()
      const {
    return matrix_type<ROWS, COLS>::Identity();
  }

  /// Set input matrix as zero matrix
  template <int ROWS, int COLS>
  ALGEBRA_HOST_DEVICE inline constexpr void set_zero(
      matrix_type<ROWS, COLS> &m) const {
    m.setZero();
  }

  /// Set input matrix as identity matrix
  template <int ROWS, int COLS>
  ALGEBRA_HOST_DEVICE inline constexpr void set_identity(
      matrix_type<ROWS, COLS> &m) const {
    m.setIdentity();
  }

  /// Multiply two matrices/expressions
  template <class DerivedA, class DerivedB,
            std::enable_if_t<Eigen::MatrixBase<DerivedA>::ColsAtCompileTime ==
                                 Eigen::MatrixBase<DerivedB>::RowsAtCompileTime,
                             bool> = true>
  ALGEBRA_HOST_DEVICE inline constexpr decltype(auto) mat_mul(
      const Eigen::MatrixBase<DerivedA> &a,
      const Eigen::MatrixBase<DerivedB> &b) const {
    return a * b;
  }

  /// Create transpose matrix
  template <class Derived>
  ALGEBRA_HOST_DEVICE inline constexpr matrix_type<
      Eigen::MatrixBase<Derived>::ColsAtCompileTime,
      Eigen::MatrixBase<Derived>::RowsAtCompileTime>
  transpose(const Eigen::MatrixBase<Derived> &m) const {
    return m.transpose();
  }

  /// Get determinant
  template <class Derived,
            std::enable_if_t<Eigen::MatrixBase<Derived>::RowsAtCompileTime ==
                                 Eigen::MatrixBase<Derived>::ColsAtCompileTime,
                             bool> = true>
  ALGEBRA_HOST_DEVICE inline constexpr scalar_t determinant(
      const Eigen::MatrixBase<Derived> &m) {
    return m.determinant();
  }

  /// Create inverse matrix
  template <class Derived,
            std::enable_if_t<Eigen::MatrixBase<Derived>::RowsAtCompileTime ==
                                 Eigen::MatrixBase<Derived>::ColsAtCompileTime,
                             bool> = true>
  ALGEBRA_HOST_DEVICE inline constexpr decltype(auto) inverse(
      const Eigen::MatrixBase<Derived> &m) const {
    return m.inverse();
  }
};

}  // namespace algebra::eigen::matrix