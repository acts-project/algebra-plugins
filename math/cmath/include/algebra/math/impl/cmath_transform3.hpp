/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/algebra_qualifiers.hpp"
#include "algebra/math/impl/cmath_getter.hpp"
#include "algebra/math/impl/cmath_vector.hpp"

namespace algebra::cmath {

/** Transform wrapper class to ensure standard API within differnt plugins
 **/
template <template <typename, auto> class array_t, typename scalar_t,
          typename matrix44_t = array_t<array_t<scalar_t, 4>, 4>,
          class element_getter_t = element_getter<array_t, scalar_t>,
          class block_getter_t = block_getter<array_t, scalar_t>,
          typename vector3_t = array_t<scalar_t, 3>,
          typename point2_t = array_t<scalar_t, 2> >
struct transform3 {

  /// @name Type definitions for the struct
  /// @{

  /// Array type used by the transform
  template <typename T, auto N>
  using array_type = array_t<T, N>;
  /// Scalar type used by the transform
  using scalar_type = scalar_t;

  /// 3-element "vector" type
  using vector3 = vector3_t;
  /// Point in 3D space
  using point3 = vector3;
  /// Point in 2D space
  using point2 = point2_t;

  /// 4x4 matrix type
  using matrix44 = matrix44_t;

  /// Function (object) used for accessing a matrix element
  using element_getter = element_getter_t;
  /// Function (object) used for accessing a sub-matrix of a matrix
  using block_getter = block_getter_t;

  /// @}

  /// @name Data objects
  /// @{

  matrix44 _data;
  matrix44 _data_inv;

  /// @}

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

    auto y = cross(z, x);
    element_getter()(_data, 0, 0) = x[0];
    element_getter()(_data, 1, 0) = x[1];
    element_getter()(_data, 2, 0) = x[2];
    element_getter()(_data, 3, 0) = 0.;
    element_getter()(_data, 0, 1) = y[0];
    element_getter()(_data, 1, 1) = y[1];
    element_getter()(_data, 2, 1) = y[2];
    element_getter()(_data, 3, 1) = 0.;
    element_getter()(_data, 0, 2) = z[0];
    element_getter()(_data, 1, 2) = z[1];
    element_getter()(_data, 2, 2) = z[2];
    element_getter()(_data, 3, 2) = 0.;
    element_getter()(_data, 0, 3) = t[0];
    element_getter()(_data, 1, 3) = t[1];
    element_getter()(_data, 2, 3) = t[2];
    element_getter()(_data, 3, 3) = 1.;

    _data_inv = invert(_data);
  }

  /** Constructor with arguments: translation
   *
   * @param t is the transform
   **/
  ALGEBRA_HOST_DEVICE
  transform3(const vector3 &t) {

    element_getter()(_data, 0, 0) = 1.;
    element_getter()(_data, 1, 0) = 0.;
    element_getter()(_data, 2, 0) = 0.;
    element_getter()(_data, 3, 0) = 0.;
    element_getter()(_data, 0, 1) = 0.;
    element_getter()(_data, 1, 1) = 1.;
    element_getter()(_data, 2, 1) = 0.;
    element_getter()(_data, 3, 1) = 0.;
    element_getter()(_data, 0, 2) = 0.;
    element_getter()(_data, 1, 2) = 0.;
    element_getter()(_data, 2, 2) = 1.;
    element_getter()(_data, 3, 2) = 0.;
    element_getter()(_data, 0, 3) = t[0];
    element_getter()(_data, 1, 3) = t[1];
    element_getter()(_data, 2, 3) = t[2];
    element_getter()(_data, 3, 3) = 1.;

    _data_inv = invert(_data);
  }

  /** Constructor with arguments: matrix
   *
   * @param m is the full 4x4 matrix
   **/
  ALGEBRA_HOST_DEVICE
  transform3(const matrix44 &m) {

    _data = m;
    _data_inv = invert(_data);
  }

  /** Constructor with arguments: matrix as array of scalar
   *
   * @param ma is the full 4x4 matrix 16 array
   **/
  ALGEBRA_HOST_DEVICE
  transform3(const array_type<scalar_type, 16> &ma) {

    element_getter()(_data, 0, 0) = ma[0];
    element_getter()(_data, 1, 0) = ma[4];
    element_getter()(_data, 2, 0) = ma[8];
    element_getter()(_data, 3, 0) = ma[12];
    element_getter()(_data, 0, 1) = ma[1];
    element_getter()(_data, 1, 1) = ma[5];
    element_getter()(_data, 2, 1) = ma[9];
    element_getter()(_data, 3, 1) = ma[13];
    element_getter()(_data, 0, 2) = ma[2];
    element_getter()(_data, 1, 2) = ma[6];
    element_getter()(_data, 2, 2) = ma[10];
    element_getter()(_data, 3, 2) = ma[14];
    element_getter()(_data, 0, 3) = ma[3];
    element_getter()(_data, 1, 3) = ma[7];
    element_getter()(_data, 2, 3) = ma[11];
    element_getter()(_data, 3, 3) = ma[15];

    _data_inv = invert(_data);
  }

  /** Constructor with arguments: identity
   *
   **/
  ALGEBRA_HOST_DEVICE
  transform3() {

    element_getter()(_data, 0, 0) = 1.;
    element_getter()(_data, 0, 1) = 0.;
    element_getter()(_data, 0, 2) = 0.;
    element_getter()(_data, 0, 3) = 0.;
    element_getter()(_data, 1, 0) = 0.;
    element_getter()(_data, 1, 1) = 1.;
    element_getter()(_data, 1, 2) = 0.;
    element_getter()(_data, 1, 3) = 0.;
    element_getter()(_data, 2, 0) = 0.;
    element_getter()(_data, 2, 1) = 0.;
    element_getter()(_data, 2, 2) = 1.;
    element_getter()(_data, 2, 3) = 0.;
    element_getter()(_data, 3, 0) = 0.;
    element_getter()(_data, 3, 1) = 0.;
    element_getter()(_data, 3, 2) = 0.;
    element_getter()(_data, 3, 3) = 1.;

    _data_inv = _data;
  }

  /** Default contructors */
  transform3(const transform3 &rhs) = default;

  /** Equality operator */
  ALGEBRA_HOST_DEVICE
  inline bool operator==(const transform3 &rhs) const {

    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        if (element_getter()(_data, i, j) !=
            element_getter()(rhs._data, i, j)) {
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
  static inline scalar_type determinant(const matrix44 &m) {

    return element_getter()(m, 0, 3) * element_getter()(m, 1, 2) *
               element_getter()(m, 2, 1) * element_getter()(m, 3, 0) -
           element_getter()(m, 0, 2) * element_getter()(m, 1, 3) *
               element_getter()(m, 2, 1) * element_getter()(m, 3, 0) -
           element_getter()(m, 0, 3) * element_getter()(m, 1, 1) *
               element_getter()(m, 2, 2) * element_getter()(m, 3, 0) +
           element_getter()(m, 0, 1) * element_getter()(m, 1, 3) *
               element_getter()(m, 2, 2) * element_getter()(m, 3, 0) +
           element_getter()(m, 0, 2) * element_getter()(m, 1, 1) *
               element_getter()(m, 2, 3) * element_getter()(m, 3, 0) -
           element_getter()(m, 0, 1) * element_getter()(m, 1, 2) *
               element_getter()(m, 2, 3) * element_getter()(m, 3, 0) -
           element_getter()(m, 0, 3) * element_getter()(m, 1, 2) *
               element_getter()(m, 2, 0) * element_getter()(m, 3, 1) +
           element_getter()(m, 0, 2) * element_getter()(m, 1, 3) *
               element_getter()(m, 2, 0) * element_getter()(m, 3, 1) +
           element_getter()(m, 0, 3) * element_getter()(m, 1, 0) *
               element_getter()(m, 2, 2) * element_getter()(m, 3, 1) -
           element_getter()(m, 0, 0) * element_getter()(m, 1, 3) *
               element_getter()(m, 2, 2) * element_getter()(m, 3, 1) -
           element_getter()(m, 0, 2) * element_getter()(m, 1, 0) *
               element_getter()(m, 2, 3) * element_getter()(m, 3, 1) +
           element_getter()(m, 0, 0) * element_getter()(m, 1, 2) *
               element_getter()(m, 2, 3) * element_getter()(m, 3, 1) +
           element_getter()(m, 0, 3) * element_getter()(m, 1, 1) *
               element_getter()(m, 2, 0) * element_getter()(m, 3, 2) -
           element_getter()(m, 0, 1) * element_getter()(m, 1, 3) *
               element_getter()(m, 2, 0) * element_getter()(m, 3, 2) -
           element_getter()(m, 0, 3) * element_getter()(m, 1, 0) *
               element_getter()(m, 2, 1) * element_getter()(m, 3, 2) +
           element_getter()(m, 0, 0) * element_getter()(m, 1, 3) *
               element_getter()(m, 2, 1) * element_getter()(m, 3, 2) +
           element_getter()(m, 0, 1) * element_getter()(m, 1, 0) *
               element_getter()(m, 2, 3) * element_getter()(m, 3, 2) -
           element_getter()(m, 0, 0) * element_getter()(m, 1, 1) *
               element_getter()(m, 2, 3) * element_getter()(m, 3, 2) -
           element_getter()(m, 0, 2) * element_getter()(m, 1, 1) *
               element_getter()(m, 2, 0) * element_getter()(m, 3, 3) +
           element_getter()(m, 0, 1) * element_getter()(m, 1, 2) *
               element_getter()(m, 2, 0) * element_getter()(m, 3, 3) +
           element_getter()(m, 0, 2) * element_getter()(m, 1, 0) *
               element_getter()(m, 2, 1) * element_getter()(m, 3, 3) -
           element_getter()(m, 0, 0) * element_getter()(m, 1, 2) *
               element_getter()(m, 2, 1) * element_getter()(m, 3, 3) -
           element_getter()(m, 0, 1) * element_getter()(m, 1, 0) *
               element_getter()(m, 2, 2) * element_getter()(m, 3, 3) +
           element_getter()(m, 0, 0) * element_getter()(m, 1, 1) *
               element_getter()(m, 2, 2) * element_getter()(m, 3, 3);
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
    element_getter()(i, 0, 0) =
        element_getter()(m, 1, 2) * element_getter()(m, 2, 3) *
            element_getter()(m, 3, 1) -
        element_getter()(m, 1, 3) * element_getter()(m, 2, 2) *
            element_getter()(m, 3, 1) +
        element_getter()(m, 1, 3) * element_getter()(m, 2, 1) *
            element_getter()(m, 3, 2) -
        element_getter()(m, 1, 1) * element_getter()(m, 2, 3) *
            element_getter()(m, 3, 2) -
        element_getter()(m, 1, 2) * element_getter()(m, 2, 1) *
            element_getter()(m, 3, 3) +
        element_getter()(m, 1, 1) * element_getter()(m, 2, 2) *
            element_getter()(m, 3, 3);
    element_getter()(i, 0, 1) =
        element_getter()(m, 0, 3) * element_getter()(m, 2, 2) *
            element_getter()(m, 3, 1) -
        element_getter()(m, 0, 2) * element_getter()(m, 2, 3) *
            element_getter()(m, 3, 1) -
        element_getter()(m, 0, 3) * element_getter()(m, 2, 1) *
            element_getter()(m, 3, 2) +
        element_getter()(m, 0, 1) * element_getter()(m, 2, 3) *
            element_getter()(m, 3, 2) +
        element_getter()(m, 0, 2) * element_getter()(m, 2, 1) *
            element_getter()(m, 3, 3) -
        element_getter()(m, 0, 1) * element_getter()(m, 2, 2) *
            element_getter()(m, 3, 3);
    element_getter()(i, 0, 2) =
        element_getter()(m, 0, 2) * element_getter()(m, 1, 3) *
            element_getter()(m, 3, 1) -
        element_getter()(m, 0, 3) * element_getter()(m, 1, 2) *
            element_getter()(m, 3, 1) +
        element_getter()(m, 0, 3) * element_getter()(m, 1, 1) *
            element_getter()(m, 3, 2) -
        element_getter()(m, 0, 1) * element_getter()(m, 1, 3) *
            element_getter()(m, 3, 2) -
        element_getter()(m, 0, 2) * element_getter()(m, 1, 1) *
            element_getter()(m, 3, 3) +
        element_getter()(m, 0, 1) * element_getter()(m, 1, 2) *
            element_getter()(m, 3, 3);
    element_getter()(i, 0, 3) =
        element_getter()(m, 0, 3) * element_getter()(m, 1, 2) *
            element_getter()(m, 2, 1) -
        element_getter()(m, 0, 2) * element_getter()(m, 1, 3) *
            element_getter()(m, 2, 1) -
        element_getter()(m, 0, 3) * element_getter()(m, 1, 1) *
            element_getter()(m, 2, 2) +
        element_getter()(m, 0, 1) * element_getter()(m, 1, 3) *
            element_getter()(m, 2, 2) +
        element_getter()(m, 0, 2) * element_getter()(m, 1, 1) *
            element_getter()(m, 2, 3) -
        element_getter()(m, 0, 1) * element_getter()(m, 1, 2) *
            element_getter()(m, 2, 3);
    element_getter()(i, 1, 0) =
        element_getter()(m, 1, 3) * element_getter()(m, 2, 2) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 1, 2) * element_getter()(m, 2, 3) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 1, 3) * element_getter()(m, 2, 0) *
            element_getter()(m, 3, 2) +
        element_getter()(m, 1, 0) * element_getter()(m, 2, 3) *
            element_getter()(m, 3, 2) +
        element_getter()(m, 1, 2) * element_getter()(m, 2, 0) *
            element_getter()(m, 3, 3) -
        element_getter()(m, 1, 0) * element_getter()(m, 2, 2) *
            element_getter()(m, 3, 3);
    element_getter()(i, 1, 1) =
        element_getter()(m, 0, 2) * element_getter()(m, 2, 3) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 0, 3) * element_getter()(m, 2, 2) *
            element_getter()(m, 3, 0) +
        element_getter()(m, 0, 3) * element_getter()(m, 2, 0) *
            element_getter()(m, 3, 2) -
        element_getter()(m, 0, 0) * element_getter()(m, 2, 3) *
            element_getter()(m, 3, 2) -
        element_getter()(m, 0, 2) * element_getter()(m, 2, 0) *
            element_getter()(m, 3, 3) +
        element_getter()(m, 0, 0) * element_getter()(m, 2, 2) *
            element_getter()(m, 3, 3);
    element_getter()(i, 1, 2) =
        element_getter()(m, 0, 3) * element_getter()(m, 1, 2) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 0, 2) * element_getter()(m, 1, 3) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 0, 3) * element_getter()(m, 1, 0) *
            element_getter()(m, 3, 2) +
        element_getter()(m, 0, 0) * element_getter()(m, 1, 3) *
            element_getter()(m, 3, 2) +
        element_getter()(m, 0, 2) * element_getter()(m, 1, 0) *
            element_getter()(m, 3, 3) -
        element_getter()(m, 0, 0) * element_getter()(m, 1, 2) *
            element_getter()(m, 3, 3);
    element_getter()(i, 1, 3) =
        element_getter()(m, 0, 2) * element_getter()(m, 1, 3) *
            element_getter()(m, 2, 0) -
        element_getter()(m, 0, 3) * element_getter()(m, 1, 2) *
            element_getter()(m, 2, 0) +
        element_getter()(m, 0, 3) * element_getter()(m, 1, 0) *
            element_getter()(m, 2, 2) -
        element_getter()(m, 0, 0) * element_getter()(m, 1, 3) *
            element_getter()(m, 2, 2) -
        element_getter()(m, 0, 2) * element_getter()(m, 1, 0) *
            element_getter()(m, 2, 3) +
        element_getter()(m, 0, 0) * element_getter()(m, 1, 2) *
            element_getter()(m, 2, 3);
    element_getter()(i, 2, 0) =
        element_getter()(m, 1, 1) * element_getter()(m, 2, 3) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 1, 3) * element_getter()(m, 2, 1) *
            element_getter()(m, 3, 0) +
        element_getter()(m, 1, 3) * element_getter()(m, 2, 0) *
            element_getter()(m, 3, 1) -
        element_getter()(m, 1, 0) * element_getter()(m, 2, 3) *
            element_getter()(m, 3, 1) -
        element_getter()(m, 1, 1) * element_getter()(m, 2, 0) *
            element_getter()(m, 3, 3) +
        element_getter()(m, 1, 0) * element_getter()(m, 2, 1) *
            element_getter()(m, 3, 3);
    element_getter()(i, 2, 1) =
        element_getter()(m, 0, 3) * element_getter()(m, 2, 1) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 0, 1) * element_getter()(m, 2, 3) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 0, 3) * element_getter()(m, 2, 0) *
            element_getter()(m, 3, 1) +
        element_getter()(m, 0, 0) * element_getter()(m, 2, 3) *
            element_getter()(m, 3, 1) +
        element_getter()(m, 0, 1) * element_getter()(m, 2, 0) *
            element_getter()(m, 3, 3) -
        element_getter()(m, 0, 0) * element_getter()(m, 2, 1) *
            element_getter()(m, 3, 3);
    element_getter()(i, 2, 2) =
        element_getter()(m, 0, 1) * element_getter()(m, 1, 3) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 0, 3) * element_getter()(m, 1, 1) *
            element_getter()(m, 3, 0) +
        element_getter()(m, 0, 3) * element_getter()(m, 1, 0) *
            element_getter()(m, 3, 1) -
        element_getter()(m, 0, 0) * element_getter()(m, 1, 3) *
            element_getter()(m, 3, 1) -
        element_getter()(m, 0, 1) * element_getter()(m, 1, 0) *
            element_getter()(m, 3, 3) +
        element_getter()(m, 0, 0) * element_getter()(m, 1, 1) *
            element_getter()(m, 3, 3);
    element_getter()(i, 2, 3) =
        element_getter()(m, 0, 3) * element_getter()(m, 1, 1) *
            element_getter()(m, 2, 0) -
        element_getter()(m, 0, 1) * element_getter()(m, 1, 3) *
            element_getter()(m, 2, 0) -
        element_getter()(m, 0, 3) * element_getter()(m, 1, 0) *
            element_getter()(m, 2, 1) +
        element_getter()(m, 0, 0) * element_getter()(m, 1, 3) *
            element_getter()(m, 2, 1) +
        element_getter()(m, 0, 1) * element_getter()(m, 1, 0) *
            element_getter()(m, 2, 3) -
        element_getter()(m, 0, 0) * element_getter()(m, 1, 1) *
            element_getter()(m, 2, 3);
    element_getter()(i, 3, 0) =
        element_getter()(m, 1, 2) * element_getter()(m, 2, 1) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 1, 1) * element_getter()(m, 2, 2) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 1, 2) * element_getter()(m, 2, 0) *
            element_getter()(m, 3, 1) +
        element_getter()(m, 1, 0) * element_getter()(m, 2, 2) *
            element_getter()(m, 3, 1) +
        element_getter()(m, 1, 1) * element_getter()(m, 2, 0) *
            element_getter()(m, 3, 2) -
        element_getter()(m, 1, 0) * element_getter()(m, 2, 1) *
            element_getter()(m, 3, 2);
    element_getter()(i, 3, 1) =
        element_getter()(m, 0, 1) * element_getter()(m, 2, 2) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 0, 2) * element_getter()(m, 2, 1) *
            element_getter()(m, 3, 0) +
        element_getter()(m, 0, 2) * element_getter()(m, 2, 0) *
            element_getter()(m, 3, 1) -
        element_getter()(m, 0, 0) * element_getter()(m, 2, 2) *
            element_getter()(m, 3, 1) -
        element_getter()(m, 0, 1) * element_getter()(m, 2, 0) *
            element_getter()(m, 3, 2) +
        element_getter()(m, 0, 0) * element_getter()(m, 2, 1) *
            element_getter()(m, 3, 2);
    element_getter()(i, 3, 2) =
        element_getter()(m, 0, 2) * element_getter()(m, 1, 1) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 0, 1) * element_getter()(m, 1, 2) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 0, 2) * element_getter()(m, 1, 0) *
            element_getter()(m, 3, 1) +
        element_getter()(m, 0, 0) * element_getter()(m, 1, 2) *
            element_getter()(m, 3, 1) +
        element_getter()(m, 0, 1) * element_getter()(m, 1, 0) *
            element_getter()(m, 3, 2) -
        element_getter()(m, 0, 0) * element_getter()(m, 1, 1) *
            element_getter()(m, 3, 2);
    element_getter()(i, 3, 3) =
        element_getter()(m, 0, 1) * element_getter()(m, 1, 2) *
            element_getter()(m, 2, 0) -
        element_getter()(m, 0, 2) * element_getter()(m, 1, 1) *
            element_getter()(m, 2, 0) +
        element_getter()(m, 0, 2) * element_getter()(m, 1, 0) *
            element_getter()(m, 2, 1) -
        element_getter()(m, 0, 0) * element_getter()(m, 1, 2) *
            element_getter()(m, 2, 1) -
        element_getter()(m, 0, 1) * element_getter()(m, 1, 0) *
            element_getter()(m, 2, 2) +
        element_getter()(m, 0, 0) * element_getter()(m, 1, 1) *
            element_getter()(m, 2, 2);
    scalar_type idet = 1. / determinant(i);
    for (unsigned int c = 0; c < 4; ++c) {
      for (unsigned int r = 0; r < 4; ++r) {
        element_getter()(i, c, r) *= idet;
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

    return {
        element_getter()(m, 0, 0) * v[0] + element_getter()(m, 0, 1) * v[1] +
            element_getter()(m, 0, 2) * v[2],
        element_getter()(m, 1, 0) * v[0] + element_getter()(m, 1, 1) * v[1] +
            element_getter()(m, 1, 2) * v[2],
        element_getter()(m, 2, 0) * v[0] + element_getter()(m, 2, 1) * v[1] +
            element_getter()(m, 2, 2) * v[2]};
  }

  /** This method retrieves the rotation of a transform */
  ALGEBRA_HOST_DEVICE
  auto inline rotation() const {

    return block_getter().template operator()<3, 3>(_data, 0, 0);
  }

  /** This method retrieves the translation of a transform */
  ALGEBRA_HOST_DEVICE
  inline point3 translation() const {

    return {element_getter()(_data, 0, 3), element_getter()(_data, 1, 3),
            element_getter()(_data, 2, 3)};
  }

  /** This method retrieves the 4x4 matrix of a transform */
  ALGEBRA_HOST_DEVICE
  inline const matrix44 &matrix() const { return _data; }

  /** This method transform from a point from the local 3D cartesian frame to
   * the global 3D cartesian frame */
  ALGEBRA_HOST_DEVICE inline point3 point_to_global(const point3 &v) const {

    const vector3 rg = rotate(_data, v);
    return {rg[0] + element_getter()(_data, 0, 3),
            rg[1] + element_getter()(_data, 1, 3),
            rg[2] + element_getter()(_data, 2, 3)};
  }

  /** This method transform from a vector from the global 3D cartesian frame
   * into the local 3D cartesian frame */

  ALGEBRA_HOST_DEVICE inline point3 point_to_local(const point3 &v) const {

    const vector3 rg = rotate(_data_inv, v);
    return {rg[0] + element_getter()(_data_inv, 0, 3),
            rg[1] + element_getter()(_data_inv, 1, 3),
            rg[2] + element_getter()(_data_inv, 2, 3)};
  }

  /** This method transform from a vector from the local 3D cartesian frame to
   * the global 3D cartesian frame */
  ALGEBRA_HOST_DEVICE inline vector3 vector_to_global(const vector3 &v) const {

    return rotate(_data, v);
  }

  /** This method transform from a vector from the global 3D cartesian frame
   * into the local 3D cartesian frame */
  ALGEBRA_HOST_DEVICE inline vector3 vector_to_local(const vector3 &v) const {

    return rotate(_data_inv, v);
  }

};  // struct transform3

}  // namespace algebra::cmath
