/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/qualifiers.hpp"

// Eigen include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#include <Eigen/Core>
#include <Eigen/Geometry>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC

// System include(s).
#include <type_traits>

namespace algebra::eigen::math {

namespace internal {

/// Functor used to access elements of Eigen matrices
template <typename scalar_t>
struct element_getter {
  /// Get non-const access to a matrix element
  template <
      typename derived_type,
      std::enable_if_t<std::is_base_of<Eigen::DenseCoeffsBase<
                                           derived_type, Eigen::WriteAccessors>,
                                       Eigen::MatrixBase<derived_type> >::value,
                       bool> = true>
  ALGEBRA_HOST_DEVICE inline scalar_t &operator()(
      Eigen::MatrixBase<derived_type> &m, unsigned int row,
      unsigned int col) const {

    return m(row, col);
  }
  /// Get const access to a matrix element
  template <typename derived_type>
  ALGEBRA_HOST_DEVICE inline scalar_t operator()(
      const Eigen::MatrixBase<derived_type> &m, unsigned int row,
      unsigned int col) const {

    return m(row, col);
  }
};  // struct element_getter

}  // namespace internal

/** Transform wrapper class to ensure standard API within differnt plugins */
template <typename scalar_t>
struct transform3 {

  /// @name Type definitions for the struct
  /// @{

  /// Array type used by the transform
  template <typename T, std::size_t N>
  using array_type = Eigen::Matrix<T, N, 1>;
  /// Scalar type used by the transform
  using scalar_type = scalar_t;

  /// 3-element "vector" type
  using vector3 = array_type<scalar_type, 3>;
  /// Point in 3D space
  using point3 = vector3;
  /// Point in 2D space
  using point2 = array_type<scalar_type, 2>;

  /// 4x4 matrix type
  using matrix44 =
      typename Eigen::Transform<scalar_type, 3, Eigen::Affine>::MatrixType;

  /// Function (object) used for accessing a matrix element
  using element_getter = internal::element_getter<scalar_type>;

  /// @}

  /// @name Data objects
  /// @{

  Eigen::Transform<scalar_type, 3, Eigen::Affine> _data;

  Eigen::Transform<scalar_type, 3, Eigen::Affine> _data_inv;

  /// @}

  /** Contructor with arguments: t, z, x
   *
   * @param t the translation (or origin of the new frame)
   * @param z the z axis of the new frame, normal vector for planes
   * @param x the x axis of the new frame
   *
   **/
  ALGEBRA_HOST_DEVICE
  transform3(const vector3 &t, const vector3 &z, const vector3 &x) {

    auto y = z.cross(x);

    auto &matrix = _data.matrix();
    matrix.template block<3, 1>(0, 0) = x;
    matrix.template block<3, 1>(0, 1) = y;
    matrix.template block<3, 1>(0, 2) = z;
    matrix.template block<3, 1>(0, 3) = t;

    _data_inv = _data.inverse();
  }

  /** Constructor with arguments: translation
   *
   * @param t is the transform
   **/
  ALGEBRA_HOST_DEVICE
  transform3(const vector3 &t) {
    _data.setIdentity();

    auto &matrix = _data.matrix();
    matrix.template block<3, 1>(0, 3) = t;

    _data_inv = _data.inverse();
  }

  /** Constructor with arguments: matrix
   *
   * @param m is the full 4x4 matrix
   **/
  ALGEBRA_HOST_DEVICE
  transform3(const matrix44 &m) {

    _data.matrix() = m;

    _data_inv = _data.inverse();
  }

  /** Constructor with arguments: matrix as std::aray of scalar
   *
   * @param ma is the full 4x4 matrix as a 16 array
   **/
  ALGEBRA_HOST_DEVICE
  transform3(const array_type<scalar_type, 16> &ma) {

    _data.matrix() << ma[0], ma[1], ma[2], ma[3], ma[4], ma[5], ma[6], ma[7],
        ma[8], ma[9], ma[10], ma[11], ma[12], ma[13], ma[14], ma[15];

    _data_inv = _data.inverse();
  }

  /** Default constructor: set contents to identity matrices */
  ALGEBRA_HOST_DEVICE
  transform3() {
    _data.setIdentity();
    _data_inv.setIdentity();
  }

  /** Default contructors */
  transform3(const transform3 &rhs) = default;
  ~transform3() = default;

  /** Equality operator */
  ALGEBRA_HOST_DEVICE
  inline bool operator==(const transform3 &rhs) const {

    return (_data.isApprox(rhs._data));
  }

  /** Rotate a vector into / from a frame
   *
   * @param m is the rotation matrix
   * @param v is the vector to be rotated
   */
  template <
      typename derived_type,
      std::enable_if_t<Eigen::MatrixBase<derived_type>::RowsAtCompileTime == 3,
                       bool> = true,
      std::enable_if_t<Eigen::MatrixBase<derived_type>::ColsAtCompileTime == 1,
                       bool> = true>
  ALGEBRA_HOST_DEVICE static inline auto rotate(
      const Eigen::Transform<scalar_type, 3, Eigen::Affine> &m,
      const Eigen::MatrixBase<derived_type> &v) {

    return m.matrix().template block<3, 3>(0, 0) * v;
  }

  /** This method retrieves the rotation of a transform  **/
  ALGEBRA_HOST_DEVICE
  inline auto rotation() const {

    return _data.matrix().template block<3, 3>(0, 0);
  }

  /** This method retrieves the translation of a transform **/
  ALGEBRA_HOST_DEVICE
  inline auto translation() const {

    return _data.matrix().template block<3, 1>(0, 3);
  }

  /** This method retrieves the 4x4 matrix of a transform */
  ALGEBRA_HOST_DEVICE
  inline const matrix44 &matrix() const { return _data.matrix(); }

  /** This method transform from a point from the local 3D cartesian frame to
   * the global 3D cartesian frame */
  template <
      typename derived_type,
      std::enable_if_t<Eigen::MatrixBase<derived_type>::RowsAtCompileTime == 3,
                       bool> = true,
      std::enable_if_t<Eigen::MatrixBase<derived_type>::ColsAtCompileTime == 1,
                       bool> = true>
  ALGEBRA_HOST_DEVICE inline auto point_to_global(
      const Eigen::MatrixBase<derived_type> &v) const {

    return (_data * v);
  }

  /** This method transform from a vector from the global 3D cartesian frame
   * into the local 3D cartesian frame */
  template <
      typename derived_type,
      std::enable_if_t<Eigen::MatrixBase<derived_type>::RowsAtCompileTime == 3,
                       bool> = true,
      std::enable_if_t<Eigen::MatrixBase<derived_type>::ColsAtCompileTime == 1,
                       bool> = true>
  ALGEBRA_HOST_DEVICE inline auto point_to_local(
      const Eigen::MatrixBase<derived_type> &v) const {

    return (_data_inv * v);
  }

  /** This method transform from a vector from the local 3D cartesian frame to
   * the global 3D cartesian frame */
  template <
      typename derived_type,
      std::enable_if_t<Eigen::MatrixBase<derived_type>::RowsAtCompileTime == 3,
                       bool> = true,
      std::enable_if_t<Eigen::MatrixBase<derived_type>::ColsAtCompileTime == 1,
                       bool> = true>
  ALGEBRA_HOST_DEVICE inline auto vector_to_global(
      const Eigen::MatrixBase<derived_type> &v) const {

    return (_data.linear() * v);
  }

  /** This method transform from a vector from the global 3D cartesian frame
   * into the local 3D cartesian frame */
  template <
      typename derived_type,
      std::enable_if_t<Eigen::MatrixBase<derived_type>::RowsAtCompileTime == 3,
                       bool> = true,
      std::enable_if_t<Eigen::MatrixBase<derived_type>::ColsAtCompileTime == 1,
                       bool> = true>
  ALGEBRA_HOST_DEVICE inline auto vector_to_local(
      const Eigen::MatrixBase<derived_type> &v) const {

    return (_data_inv.linear() * v);
  }
};  // struct transform3

}  // namespace algebra::eigen::math
