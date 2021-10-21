/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/algebra_qualifiers.hpp"

// System include(s).
#include <cstddef>

namespace algebra {
namespace cmath {

/** Transform wrapper class to ensure standard API within differnt plugins
 **/
template<template <typename,std::size_t> class array_t,
         typename scalar_t>
struct transform3 {

  /// @name Type definitions for the struct
  /// @{

  /// Array type used by the transform
  template <typename T, std::size_t N>
  using array_type = array_t<T, N>;
  /// Scalar type used by the transform
  using scalar_type = scalar_t;

  /// 3-element "vector" type
  using vector3 = array_type<scalar_type, 3>;
  /// Point in 3D space
  using point3 = vector3;
  /// Point in 2D space
  using point2 = array_type<scalar_type, 2>;

  /// 4x4 matrix type
  using matrix44 = array_type<array_type<scalar_type, 4>, 4>;

  /// @}

  matrix44 _data;
  matrix44 _data_inv;

  /** Contructor with arguments: t, z, x
   *
   * @param t the translation (or origin of the new frame)
   * @param z the z axis of the new frame, normal vector for planes
   * @param x the x axis of the new frame
   *
   * @note y will be constructed by cross product
   *
   **/
  ALGEBRA_HOST_DEVICE
  transform3(const vector3 &t, const vector3 &z, const vector3 &x) {

    auto y = vector::cross(z, x);
    _data[0][0] = x[0];
    _data[0][1] = x[1];
    _data[0][2] = x[2];
    _data[0][3] = 0.;
    _data[1][0] = y[0];
    _data[1][1] = y[1];
    _data[1][2] = y[2];
    _data[1][3] = 0.;
    _data[2][0] = z[0];
    _data[2][1] = z[1];
    _data[2][2] = z[2];
    _data[2][3] = 0.;
    _data[3][0] = t[0];
    _data[3][1] = t[1];
    _data[3][2] = t[2];
    _data[3][3] = 1.;

    _data_inv = invert(_data);
  }

  /** Constructor with arguments: translation
   *
   * @param t is the transform
   **/
  ALGEBRA_HOST_DEVICE
  transform3(const vector3 &t) {

    _data[0][0] = 1.;
    _data[0][1] = 0.;
    _data[0][2] = 0.;
    _data[0][3] = 0.;
    _data[1][0] = 0.;
    _data[1][1] = 1.;
    _data[1][2] = 0.;
    _data[1][3] = 0.;
    _data[2][0] = 0.;
    _data[2][1] = 0.;
    _data[2][2] = 1.;
    _data[2][3] = 0.;
    _data[3][0] = t[0];
    _data[3][1] = t[1];
    _data[3][2] = t[2];
    _data[3][3] = 1.;

    _data_inv = invert(_data);
  }

  /** Constructor with arguments: matrix
   *
   * @param m is the full 4x4 matrix
   **/
  ALGEBRA_HOST_DEVICE
  transform3(const matrix44 &m) { _data = m; }

  /** Constructor with arguments: matrix as std::aray of scalar
   *
   * @param ma is the full 4x4 matrix 16 array
   **/
  ALGEBRA_HOST_DEVICE
  transform3(const array_type<scalar, 16> &ma) {

    _data[0][0] = ma[0];
    _data[0][1] = ma[4];
    _data[0][2] = ma[8];
    _data[0][3] = ma[12];
    _data[1][0] = ma[1];
    _data[1][1] = ma[5];
    _data[1][2] = ma[9];
    _data[1][3] = ma[13];
    _data[2][0] = ma[2];
    _data[2][1] = ma[6];
    _data[2][2] = ma[10];
    _data[2][3] = ma[14];
    _data[3][0] = ma[3];
    _data[3][1] = ma[7];
    _data[3][2] = ma[11];
    _data[3][3] = ma[15];

    _data_inv = invert(_data);
  }

  /** Constructor with arguments: identity
   *
   **/
  ALGEBRA_HOST_DEVICE
  transform3() {

    _data[0][0] = 1.;
    _data[0][1] = 0.;
    _data[0][2] = 0.;
    _data[0][3] = 0.;
    _data[1][0] = 0.;
    _data[1][1] = 1.;
    _data[1][2] = 0.;
    _data[1][3] = 0.;
    _data[2][0] = 0.;
    _data[2][1] = 0.;
    _data[2][2] = 1.;
    _data[2][3] = 0.;
    _data[3][0] = 0.;
    _data[3][1] = 0.;
    _data[3][2] = 0.;
    _data[3][3] = 1.;

    _data_inv = _data;
  }

  /** Default contructors */
  transform3(const transform3 &rhs) = default;

  /** Equality operator */
  ALGEBRA_HOST_DEVICE
  inline bool operator==(const transform3 &rhs) const {

    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        if (_data[i][j] != rhs._data[i][j]) {
          return false;
        }
      }
    }

    return true;
  }

  /** The determinant of a 4x4 matrix
   *
   * @param m is the matrix
   *
   * @return a sacalar determinant - no checking done
   */
  ALGEBRA_HOST_DEVICE
  static inline scalar determinant(const matrix44 &m) {

    return m[3][0] * m[2][1] * m[1][2] * m[0][3] -
           m[2][0] * m[3][1] * m[1][2] * m[0][3] -
           m[3][0] * m[1][1] * m[2][2] * m[0][3] +
           m[1][0] * m[3][1] * m[2][2] * m[0][3] +
           m[2][0] * m[1][1] * m[3][2] * m[0][3] -
           m[1][0] * m[2][1] * m[3][2] * m[0][3] -
           m[3][0] * m[2][1] * m[0][2] * m[1][3] +
           m[2][0] * m[3][1] * m[0][2] * m[1][3] +
           m[3][0] * m[0][1] * m[2][2] * m[1][3] -
           m[0][0] * m[3][1] * m[2][2] * m[1][3] -
           m[2][0] * m[0][1] * m[3][2] * m[1][3] +
           m[0][0] * m[2][1] * m[3][2] * m[1][3] +
           m[3][0] * m[1][1] * m[0][2] * m[2][3] -
           m[1][0] * m[3][1] * m[0][2] * m[2][3] -
           m[3][0] * m[0][1] * m[1][2] * m[2][3] +
           m[0][0] * m[3][1] * m[1][2] * m[2][3] +
           m[1][0] * m[0][1] * m[3][2] * m[2][3] -
           m[0][0] * m[1][1] * m[3][2] * m[2][3] -
           m[2][0] * m[1][1] * m[0][2] * m[3][3] +
           m[1][0] * m[2][1] * m[0][2] * m[3][3] +
           m[2][0] * m[0][1] * m[1][2] * m[3][3] -
           m[0][0] * m[2][1] * m[1][2] * m[3][3] -
           m[1][0] * m[0][1] * m[2][2] * m[3][3] +
           m[0][0] * m[1][1] * m[2][2] * m[3][3];
  }

  /** The inverse of a 4x4 matrix
   *
   * @param m is the matrix
   *
   * @return an inverse matrix
   */
  ALGEBRA_HOST_DEVICE
  static inline matrix44 invert(const matrix44 &m) {

    matrix44 i;
    i[0][0] = m[2][1] * m[3][2] * m[1][3] - m[3][1] * m[2][2] * m[1][3] +
              m[3][1] * m[1][2] * m[2][3] - m[1][1] * m[3][2] * m[2][3] -
              m[2][1] * m[1][2] * m[3][3] + m[1][1] * m[2][2] * m[3][3];
    i[1][0] = m[3][0] * m[2][2] * m[1][3] - m[2][0] * m[3][2] * m[1][3] -
              m[3][0] * m[1][2] * m[2][3] + m[1][0] * m[3][2] * m[2][3] +
              m[2][0] * m[1][2] * m[3][3] - m[1][0] * m[2][2] * m[3][3];
    i[2][0] = m[2][0] * m[3][1] * m[1][3] - m[3][0] * m[2][1] * m[1][3] +
              m[3][0] * m[1][1] * m[2][3] - m[1][0] * m[3][1] * m[2][3] -
              m[2][0] * m[1][1] * m[3][3] + m[1][0] * m[2][1] * m[3][3];
    i[3][0] = m[3][0] * m[2][1] * m[1][2] - m[2][0] * m[3][1] * m[1][2] -
              m[3][0] * m[1][1] * m[2][2] + m[1][0] * m[3][1] * m[2][2] +
              m[2][0] * m[1][1] * m[3][2] - m[1][0] * m[2][1] * m[3][2];
    i[0][1] = m[3][1] * m[2][2] * m[0][3] - m[2][1] * m[3][2] * m[0][3] -
              m[3][1] * m[0][2] * m[2][3] + m[0][1] * m[3][2] * m[2][3] +
              m[2][1] * m[0][2] * m[3][3] - m[0][1] * m[2][2] * m[3][3];
    i[1][1] = m[2][0] * m[3][2] * m[0][3] - m[3][0] * m[2][2] * m[0][3] +
              m[3][0] * m[0][2] * m[2][3] - m[0][0] * m[3][2] * m[2][3] -
              m[2][0] * m[0][2] * m[3][3] + m[0][0] * m[2][2] * m[3][3];
    i[2][1] = m[3][0] * m[2][1] * m[0][3] - m[2][0] * m[3][1] * m[0][3] -
              m[3][0] * m[0][1] * m[2][3] + m[0][0] * m[3][1] * m[2][3] +
              m[2][0] * m[0][1] * m[3][3] - m[0][0] * m[2][1] * m[3][3];
    i[3][1] = m[2][0] * m[3][1] * m[0][2] - m[3][0] * m[2][1] * m[0][2] +
              m[3][0] * m[0][1] * m[2][2] - m[0][0] * m[3][1] * m[2][2] -
              m[2][0] * m[0][1] * m[3][2] + m[0][0] * m[2][1] * m[3][2];
    i[0][2] = m[1][1] * m[3][2] * m[0][3] - m[3][1] * m[1][2] * m[0][3] +
              m[3][1] * m[0][2] * m[1][3] - m[0][1] * m[3][2] * m[1][3] -
              m[1][1] * m[0][2] * m[3][3] + m[0][1] * m[1][2] * m[3][3];
    i[1][2] = m[3][0] * m[1][2] * m[0][3] - m[1][0] * m[3][2] * m[0][3] -
              m[3][0] * m[0][2] * m[1][3] + m[0][0] * m[3][2] * m[1][3] +
              m[1][0] * m[0][2] * m[3][3] - m[0][0] * m[1][2] * m[3][3];
    i[2][2] = m[1][0] * m[3][1] * m[0][3] - m[3][0] * m[1][1] * m[0][3] +
              m[3][0] * m[0][1] * m[1][3] - m[0][0] * m[3][1] * m[1][3] -
              m[1][0] * m[0][1] * m[3][3] + m[0][0] * m[1][1] * m[3][3];
    i[3][2] = m[3][0] * m[1][1] * m[0][2] - m[1][0] * m[3][1] * m[0][2] -
              m[3][0] * m[0][1] * m[1][2] + m[0][0] * m[3][1] * m[1][2] +
              m[1][0] * m[0][1] * m[3][2] - m[0][0] * m[1][1] * m[3][2];
    i[0][3] = m[2][1] * m[1][2] * m[0][3] - m[1][1] * m[2][2] * m[0][3] -
              m[2][1] * m[0][2] * m[1][3] + m[0][1] * m[2][2] * m[1][3] +
              m[1][1] * m[0][2] * m[2][3] - m[0][1] * m[1][2] * m[2][3];
    i[1][3] = m[1][0] * m[2][2] * m[0][3] - m[2][0] * m[1][2] * m[0][3] +
              m[2][0] * m[0][2] * m[1][3] - m[0][0] * m[2][2] * m[1][3] -
              m[1][0] * m[0][2] * m[2][3] + m[0][0] * m[1][2] * m[2][3];
    i[2][3] = m[2][0] * m[1][1] * m[0][3] - m[1][0] * m[2][1] * m[0][3] -
              m[2][0] * m[0][1] * m[1][3] + m[0][0] * m[2][1] * m[1][3] +
              m[1][0] * m[0][1] * m[2][3] - m[0][0] * m[1][1] * m[2][3];
    i[3][3] = m[1][0] * m[2][1] * m[0][2] - m[2][0] * m[1][1] * m[0][2] +
              m[2][0] * m[0][1] * m[1][2] - m[0][0] * m[2][1] * m[1][2] -
              m[1][0] * m[0][1] * m[2][2] + m[0][0] * m[1][1] * m[2][2];
    scalar idet = 1. / determinant(i);
    for (unsigned int c = 0; c < 4; ++c) {
      for (unsigned int r = 0; r < 4; ++r) {
        i[c][r] *= idet;
      }
    }
    return i;
  }

  /** Rotate a vector into / from a frame
   *
   * @param m is the rotation matrix
   * @param v is the vector to be rotated
   */
  ALGEBRA_HOST_DEVICE
  static inline vector3 rotate(const matrix44 &m, const vector3 &v) {

    return vector3{m[0][0] * v[0] + m[1][0] * v[1] + m[2][0] * v[2],
                   m[0][1] * v[0] + m[1][1] * v[1] + m[2][1] * v[2],
                   m[0][2] * v[0] + m[1][2] * v[1] + m[2][2] * v[2]};
  }

  /** This method retrieves the rotation of a transform */
  ALGEBRA_HOST_DEVICE
  auto inline rotation() const { return getter::block<3, 3>(_data, 0, 0); }

  /** This method retrieves the translation of a transform */
  ALGEBRA_HOST_DEVICE
  inline point3 translation() const {
    return point3{_data[3][0], _data[3][1], _data[3][2]};
  }

  /** This method retrieves the 4x4 matrix of a transform */
  ALGEBRA_HOST_DEVICE
  inline const matrix44 &matrix() const { return _data; }

  /** This method transform from a point from the local 3D cartesian frame to
   * the global 3D cartesian frame */
  template <typename point_type>
  ALGEBRA_HOST_DEVICE inline const point_type point_to_global(
      const point_type &v) const {
    vector3 rg = rotate(_data, v);
    return point3{rg[0] + _data[3][0], rg[1] + _data[3][1],
                  rg[2] + _data[3][2]};
  }

  /** This method transform from a vector from the global 3D cartesian frame
   * into the local 3D cartesian frame */
  template <typename point_type>
  ALGEBRA_HOST_DEVICE inline const point_type point_to_local(
      const point_type &v) const {
    vector3 rg = rotate(_data_inv, v);
    return point3{rg[0] + _data_inv[3][0], rg[1] + _data_inv[3][1],
                  rg[2] + _data_inv[3][2]};
  }

  /** This method transform from a vector from the local 3D cartesian frame to
   * the global 3D cartesian frame */
  template <typename vector_type>
  ALGEBRA_HOST_DEVICE inline const vector_type vector_to_global(
      const vector_type &v) const {
    return rotate(_data, v);
  }

  /** This method transform from a vector from the global 3D cartesian frame
   * into the local 3D cartesian frame */
  template <typename vector_type>
  ALGEBRA_HOST_DEVICE inline const auto vector_to_local(
      const vector_type &v) const {
    return rotate(_data_inv, v);
  }

};  // struct transform3

}  // namespace cmath
}  // namespace algebra
