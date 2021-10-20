/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/algebra_qualifiers.hpp"
#include "algebra/common/types.hpp"

// System include(s).
#include <array>
#include <cmath>

#ifdef ALGEBRA_PLUGIN_CUSTOM_SCALARTYPE
using algebra_scalar = ALGEBRA_PLUGIN_CUSTOM_SCALARTYPE;
#else
using algebra_scalar = double;
#endif

#define __plugin_without_matrix_element_accessor 1

namespace algebra {

using scalar = algebra_scalar;

ALGEBRA_HOST_DEVICE
inline __plugin_array<scalar, 2> operator*(const __plugin_array<scalar, 2> &a,
                                           scalar s) {
  return {a[0] * s, a[1] * s};
}

ALGEBRA_HOST_DEVICE
inline __plugin_array<scalar, 2> operator*(scalar s,
                                           const __plugin_array<scalar, 2> &a) {
  return {s * a[0], s * a[1]};
}

ALGEBRA_HOST_DEVICE
inline __plugin_array<scalar, 2> operator-(const __plugin_array<scalar, 2> &a,
                                           const __plugin_array<scalar, 2> &b) {
  return {a[0] - b[0], a[1] - b[1]};
}

ALGEBRA_HOST_DEVICE
inline __plugin_array<scalar, 2> operator+(const __plugin_array<scalar, 2> &a,
                                           const __plugin_array<scalar, 2> &b) {
  return {a[0] + b[0], a[1] + b[1]};
}

ALGEBRA_HOST_DEVICE
inline __plugin_array<scalar, 3> operator*(const __plugin_array<scalar, 3> &a,
                                           scalar s) {
  return {a[0] * s, a[1] * s, a[2] * s};
}

ALGEBRA_HOST_DEVICE
inline __plugin_array<scalar, 3> operator*(scalar s,
                                           const __plugin_array<scalar, 3> &a) {
  return {s * a[0], s * a[1], s * a[2]};
}

ALGEBRA_HOST_DEVICE
inline __plugin_array<scalar, 3> operator-(const __plugin_array<scalar, 3> &a,
                                           const __plugin_array<scalar, 3> &b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

ALGEBRA_HOST_DEVICE
inline __plugin_array<scalar, 3> operator+(const __plugin_array<scalar, 3> &a,
                                           const __plugin_array<scalar, 3> &b) {
  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

namespace vector {

/** Cross product between two input vectors - 3 Dim
 *
 * @tparam derived_type_lhs is the first matrix (epresseion) template
 * @tparam derived_type_rhs is the second matrix (epresseion) template
 *
 * @param a the first input vector
 * @param b the second input vector
 *
 * @return a vector (expression) representing the cross product
 **/
ALGEBRA_HOST_DEVICE
inline __plugin_array<scalar, 3> cross(const __plugin_array<scalar, 3> &a,
                                       const __plugin_array<scalar, 3> &b) {
  return {a[1] * b[2] - b[1] * a[2], a[2] * b[0] - b[2] * a[0],
          a[0] * b[1] - b[0] * a[1]};
}
}  // namespace vector

// array getter methdos
namespace getter {
/** This method retrieves phi from a vector, vector base with rows > 2
 *
 * @param v the input vector
 **/
template <typename vector_type>
ALGEBRA_HOST_DEVICE inline auto phi(const vector_type &v) noexcept {
  return std::atan2(v[1], v[0]);
}

/** This method retrieves theta from a vector, vector base with rows >= 3
 *
 * @param v the input vector
 **/
template <typename vector_type>
ALGEBRA_HOST_DEVICE inline auto theta(const vector_type &v) noexcept {
  return std::atan2(std::sqrt(v[0] * v[0] + v[1] * v[1]), v[2]);
}

/** This method retrieves the perpenticular magnitude of a vector with rows >= 2
 *
 * @param v the input vector
 **/
template <typename vector_type>
ALGEBRA_HOST_DEVICE inline auto perp(const vector_type &v) noexcept {
  return std::sqrt(v[0] * v[0] + v[1] * v[1]);
}

/** This method retrieves the norm of a vector, no dimension restriction
 *
 * @param v the input vector
 **/
ALGEBRA_HOST_DEVICE
inline auto norm(const __plugin_array<scalar, 2> &v) {
  return perp(v);
}

/** This method retrieves the norm of a vector, no dimension restriction
 *
 * @param v the input vector
 **/
ALGEBRA_HOST_DEVICE
inline auto norm(const __plugin_array<scalar, 3> &v) {
  return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

/** This method retrieves the pseudo-rapidity from a vector or vector base with
 *rows >= 3
 *
 * @param v the input vector
 **/
template <typename vector_type>
ALGEBRA_HOST_DEVICE inline auto eta(const vector_type &v) noexcept {
  return std::atanh(v[2] / norm(v));
}

/** This method retrieves a column from a matrix
 *
 * @param m the input matrix
 **/
template <unsigned int kROWS, typename matrix_type>
ALGEBRA_HOST_DEVICE inline auto vector(const matrix_type &m, unsigned int row,
                                       unsigned int col) noexcept {
  __plugin_array<scalar, kROWS> subvector;
  for (unsigned int irow = row; irow < row + kROWS; ++irow) {
    subvector[irow - row] = m[col][irow];
  }
  return subvector;
}

/** This method retrieves a column from a matrix
 *
 * @param m the input matrix
 **/
template <unsigned int kROWS, unsigned int kCOLS, typename matrix_type>
ALGEBRA_HOST_DEVICE inline auto block(const matrix_type &m, unsigned int row,
                                      unsigned int col) noexcept {
  __plugin_array<__plugin_array<scalar, kROWS>, kCOLS> submatrix;
  for (unsigned int icol = col; icol < col + kCOLS; ++icol) {
    for (unsigned int irow = row; irow < row + kROWS; ++irow) {
      submatrix[icol - col][irow - row] = m[icol][irow];
    }
  }
  return submatrix;
}

}  // namespace getter

// array definitions
namespace array {

using vector3 = __plugin_array<scalar, 3>;
using point3 = vector3;
using point2 = __plugin_array<scalar, 2>;

/** Transform wrapper class to ensure standard API within differnt plugins
 **/
struct transform3 {
  using matrix44 = __plugin_array<__plugin_array<scalar, 4>, 4>;

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
  transform3(const array_s<scalar, 16> &ma) {
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
  ~transform3() = default;

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
};

/** Frame projection into a cartesian coordinate frame
 */
struct cartesian2 {
  /** This method transform from a point from the global 3D cartesian frame to
   *the local 2D cartesian frame
   *
   * @param trf the transform from global to local thredimensional frame
   * @param p the point in global frame
   *
   * @return a local point2
   **/
  ALGEBRA_HOST_DEVICE
  inline point2 operator()(const transform3 &trf, const point3 &p) const {
    return operator()(trf.point_to_local(p));
  }

  /** This method transform from a point from the global 3D cartesian frame to
   * the local 2D cartesian frame
   *
   * @param v the point in local frame
   *
   * @return a local point2
   */
  ALGEBRA_HOST_DEVICE
  inline point2 operator()(const point3 &v) const { return {v[0], v[1]}; }
};

/** Local frame projection into a polar coordinate frame */
struct polar2 {
  /** This method transform from a point from the global 3D cartesian frame to
   *the local 2D cartesian frame
   *
   * @param trf the transform from global to local thredimensional frame
   * @param p the point in global frame
   *
   * @return a local point2
   **/
  ALGEBRA_HOST_DEVICE
  inline point2 operator()(const transform3 &trf, const point3 &p) const {
    return operator()(trf.point_to_local(p));
  }

  /** This method transform from a point from 2D or 3D cartesian frame to a 2D
   * polar point */
  template <typename point_type>
  ALGEBRA_HOST_DEVICE inline point2 operator()(const point_type &v) const {
    return point2{getter::perp(v), getter::phi(v)};
  }
};

/** Local frame projection into a polar coordinate frame */
struct cylindrical2 {
  /** This method transform from a point from the global 3D cartesian frame to
   *the local 2D cartesian frame
   *
   * @param trf the transform from global to local thredimensional frame
   * @param p the point in global frame
   *
   * @return a local point2
   **/
  ALGEBRA_HOST_DEVICE
  inline point2 operator()(const transform3 &trf, const point3 &p) const {
    return operator()(trf.point_to_local(p));
  }

  /** This method transform from a point from 2 3D cartesian frame to a 2D
   * cylindrical point */
  ALGEBRA_HOST_DEVICE
  inline point2 operator()(const point3 &v) const {
    return point2{getter::perp(v) * getter::phi(v), v[2]};
  }
};

}  // namespace array

// Vector transfroms
namespace vector {

/** Dot product between two input vectors - 2 Dim
 *
 * @param a the first input vector
 * @param b the second input vector
 *
 * @return the scalar dot product value
 **/
ALGEBRA_HOST_DEVICE
inline scalar dot(const __plugin_array<scalar, 2> &a,
                  const __plugin_array<scalar, 2> &b) {
  return a[0] * b[0] + a[1] * b[1];
}

/** Get a normalized version of the input vector
 *
 * @param v the input vector
 **/
ALGEBRA_HOST_DEVICE
inline __plugin_array<scalar, 2> normalize(const __plugin_array<scalar, 2> &v) {
  scalar oon = 1. / std::sqrt(dot(v, v));
  return {v[0] * oon, v[1] * oon};
}

/** Dot product between two input vectors - 3 Dim
 *
 * @param a the first input vector
 * @param b the second input vector
 *
 * @return the scalar dot product value
 **/
ALGEBRA_HOST_DEVICE
inline scalar dot(const __plugin_array<scalar, 3> &a,
                  const __plugin_array<scalar, 3> &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/** Get a normalized version of the input vector
 *
 * @param v the input vector
 **/
ALGEBRA_HOST_DEVICE
inline __plugin_array<scalar, 3> normalize(const __plugin_array<scalar, 3> &v) {
  scalar oon = 1. / std::sqrt(dot(v, v));
  return {v[0] * oon, v[1] * oon, v[2] * oon};
}

}  // namespace vector

}  // namespace algebra
