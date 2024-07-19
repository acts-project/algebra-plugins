/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/algorithms/matrix/inverse/hard_coded.hpp"
#include "algebra/math/impl/cmath_getter.hpp"
#include "algebra/math/impl/cmath_matrix.hpp"
#include "algebra/math/impl/cmath_vector.hpp"
#include "algebra/qualifiers.hpp"

namespace algebra::cmath {

/** Transform wrapper class to ensure standard API within differnt plugins
 **/
template <typename size_ty, typename scalar_t,
          template <typename, size_ty, size_ty> class matrix_t,
          template <typename, size_ty> class array_t, typename element_getter_t,
          typename block_getter_t>
struct transform3 {

  /// @name Type definitions for the struct
  /// @{

  /// Size type
  using size_type = size_ty;

  /// Scalar type
  using scalar_type = scalar_t;

  /// 2D Matrix type
  template <size_type ROWS, size_type COLS>
  using matrix_type = matrix_t<scalar_t, ROWS, COLS>;

  /// Array type
  template <size_type N>
  using array_type = array_t<scalar_t, N>;

  // 4 x 4 Matrix
  using matrix44 = matrix_type<4, 4>;

  /// 3-element "vector" type
  using vector3 = array_type<3>;
  /// Point in 3D space
  using point3 = vector3;
  /// Point in 2D space
  using point2 = array_type<2>;

  /// Function (object) used for accessing a matrix element
  using element_getter = element_getter_t;

  /// Function (object) used for accessing a matrix element
  using block_getter = block_getter_t;

  /// Matrix inversion algorithm
  using matrix_inversion =
      cmath::matrix::inverse::hard_coded<matrix44, element_getter_t>;

  /// @}

  /// @name Data objects
  /// @{
  matrix44 _data{cmath::identity<matrix44, element_getter>()};
  matrix44 _data_inv{cmath::identity<matrix44, element_getter>()};

  /// @}

  /** Default constructor: identity
   **/
  constexpr transform3() = default;

  /** Contructor with arguments: t, x, y, z
   *
   * @param t the translation (or origin of the new frame)
   * @param x the x axis of the new frame
   * @param y the y axis of the new frame
   * @param z the z axis of the new frame, normal vector for planes
   *
   **/
  ALGEBRA_HOST_DEVICE
  transform3(const vector3 &t, const vector3 &x, const vector3 &y,
             const vector3 &z, bool get_inverse = true) {

    element_getter{}(_data, 0, 0) = x[0];
    element_getter{}(_data, 1, 0) = x[1];
    element_getter{}(_data, 2, 0) = x[2];
    element_getter{}(_data, 3, 0) = 0.;
    element_getter{}(_data, 0, 1) = y[0];
    element_getter{}(_data, 1, 1) = y[1];
    element_getter{}(_data, 2, 1) = y[2];
    element_getter{}(_data, 3, 1) = 0.;
    element_getter{}(_data, 0, 2) = z[0];
    element_getter{}(_data, 1, 2) = z[1];
    element_getter{}(_data, 2, 2) = z[2];
    element_getter{}(_data, 3, 2) = 0.;
    element_getter{}(_data, 0, 3) = t[0];
    element_getter{}(_data, 1, 3) = t[1];
    element_getter{}(_data, 2, 3) = t[2];
    element_getter{}(_data, 3, 3) = 1.;

    if (get_inverse) {
      _data_inv = matrix_inversion{}(_data);
    }
  }

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
  transform3(const vector3 &t, const vector3 &z, const vector3 &x,
             bool get_inverse = true)
      : transform3(t, x, cross(z, x), z, get_inverse) {}

  /** Constructor with arguments: translation
   *
   * @param t is the transform
   **/
  ALGEBRA_HOST_DEVICE
  explicit transform3(const vector3 &t) {

    element_getter{}(_data, 0, 0) = 1.;
    element_getter{}(_data, 1, 0) = 0.;
    element_getter{}(_data, 2, 0) = 0.;
    element_getter{}(_data, 3, 0) = 0.;
    element_getter{}(_data, 0, 1) = 0.;
    element_getter{}(_data, 1, 1) = 1.;
    element_getter{}(_data, 2, 1) = 0.;
    element_getter{}(_data, 3, 1) = 0.;
    element_getter{}(_data, 0, 2) = 0.;
    element_getter{}(_data, 1, 2) = 0.;
    element_getter{}(_data, 2, 2) = 1.;
    element_getter{}(_data, 3, 2) = 0.;
    element_getter{}(_data, 0, 3) = t[0];
    element_getter{}(_data, 1, 3) = t[1];
    element_getter{}(_data, 2, 3) = t[2];
    element_getter{}(_data, 3, 3) = 1.;

    _data_inv = matrix_inversion{}(_data);
  }

  /** Constructor with arguments: matrix
   *
   * @param m is the full 4x4 matrix
   **/
  ALGEBRA_HOST_DEVICE
  explicit transform3(const matrix44 &m) {

    _data = m;
    _data_inv = matrix_inversion{}(_data);
  }

  /** Constructor with arguments: matrix as array of scalar
   *
   * @param ma is the full 4x4 matrix 16 array
   **/
  ALGEBRA_HOST_DEVICE
  explicit transform3(const array_type<16> &ma) {

    element_getter{}(_data, 0, 0) = ma[0];
    element_getter{}(_data, 1, 0) = ma[4];
    element_getter{}(_data, 2, 0) = ma[8];
    element_getter{}(_data, 3, 0) = ma[12];
    element_getter{}(_data, 0, 1) = ma[1];
    element_getter{}(_data, 1, 1) = ma[5];
    element_getter{}(_data, 2, 1) = ma[9];
    element_getter{}(_data, 3, 1) = ma[13];
    element_getter{}(_data, 0, 2) = ma[2];
    element_getter{}(_data, 1, 2) = ma[6];
    element_getter{}(_data, 2, 2) = ma[10];
    element_getter{}(_data, 3, 2) = ma[14];
    element_getter{}(_data, 0, 3) = ma[3];
    element_getter{}(_data, 1, 3) = ma[7];
    element_getter{}(_data, 2, 3) = ma[11];
    element_getter{}(_data, 3, 3) = ma[15];

    _data_inv = matrix_inversion{}(_data);
  }

  /** Equality operator */
  ALGEBRA_HOST_DEVICE
  inline bool operator==(const transform3 &rhs) const {

    for (size_type j = 0; j < 4; j++) {
      for (size_type i = 0; i < 4; i++) {
        if (element_getter{}(_data, i, j) !=
            element_getter{}(rhs._data, i, j)) {
          return false;
        }
      }
    }

    return true;
  }

  /** Rotate a vector into / from a frame
   *
   * @param m is the rotation matrix
   * @param v is the vector to be rotated
   */
  ALGEBRA_HOST_DEVICE
  static inline vector3 rotate(const matrix44 &m, const vector3 &v) {

    vector3 ret{0.f, 0.f, 0.f};

    ret[0] += element_getter{}(m, 0, 0) * v[0];
    ret[1] += element_getter{}(m, 1, 0) * v[0];
    ret[2] += element_getter{}(m, 2, 0) * v[0];

    ret[0] += element_getter{}(m, 0, 1) * v[1];
    ret[1] += element_getter{}(m, 1, 1) * v[1];
    ret[2] += element_getter{}(m, 2, 1) * v[1];

    ret[0] += element_getter{}(m, 0, 2) * v[2];
    ret[1] += element_getter{}(m, 1, 2) * v[2];
    ret[2] += element_getter{}(m, 2, 2) * v[2];

    return ret;
  }

  /** This method retrieves the rotation of a transform */
  ALGEBRA_HOST_DEVICE
  auto inline rotation() const {

    return block_getter{}.template operator()<3, 3>(_data, 0, 0);
  }

  /** This method retrieves x axis */
  ALGEBRA_HOST_DEVICE
  inline point3 x() const {
    return {element_getter{}(_data, 0, 0), element_getter{}(_data, 1, 0),
            element_getter{}(_data, 2, 0)};
  }

  /** This method retrieves y axis */
  ALGEBRA_HOST_DEVICE
  inline point3 y() const {
    return {element_getter{}(_data, 0, 1), element_getter{}(_data, 1, 1),
            element_getter{}(_data, 2, 1)};
  }

  /** This method retrieves z axis */
  ALGEBRA_HOST_DEVICE
  inline point3 z() const {
    return {element_getter{}(_data, 0, 2), element_getter{}(_data, 1, 2),
            element_getter{}(_data, 2, 2)};
  }

  /** This method retrieves the translation of a transform */
  ALGEBRA_HOST_DEVICE
  inline point3 translation() const {

    return {element_getter{}(_data, 0, 3), element_getter{}(_data, 1, 3),
            element_getter{}(_data, 2, 3)};
  }

  /** This method retrieves the 4x4 matrix of a transform */
  ALGEBRA_HOST_DEVICE
  inline const matrix44 &matrix() const { return _data; }

  /** This method retrieves the 4x4 matrix of an inverse transform */
  ALGEBRA_HOST_DEVICE
  inline const matrix44 &matrix_inverse() const { return _data_inv; }

  /** This method transform from a point from the local 3D cartesian frame to
   * the global 3D cartesian frame */
  ALGEBRA_HOST_DEVICE inline point3 point_to_global(const point3 &v) const {

    const vector3 rg = rotate(_data, v);
    return {rg[0] + element_getter{}(_data, 0, 3),
            rg[1] + element_getter{}(_data, 1, 3),
            rg[2] + element_getter{}(_data, 2, 3)};
  }

  /** This method transform from a vector from the global 3D cartesian frame
   * into the local 3D cartesian frame */

  ALGEBRA_HOST_DEVICE inline point3 point_to_local(const point3 &v) const {

    const vector3 rg = rotate(_data_inv, v);
    return {rg[0] + element_getter{}(_data_inv, 0, 3),
            rg[1] + element_getter{}(_data_inv, 1, 3),
            rg[2] + element_getter{}(_data_inv, 2, 3)};
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
