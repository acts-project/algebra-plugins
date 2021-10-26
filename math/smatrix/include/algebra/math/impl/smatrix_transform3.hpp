/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/algebra_qualifiers.hpp"
#include "algebra/math/impl/smatrix_errorcheck.hpp"

// ROOT/Smatrix include(s).
#include "Math/SMatrix.h"
#include "Math/SVector.h"

namespace algebra::smatrix::math {

namespace internal {

/// Functor used to access elements of SMatrix matrices
template <typename scalar_t>
struct element_getter {
  /// The matrix type used by this struct
  template <unsigned int ROWS, unsigned int COLS>
  using matrix_type = ROOT::Math::SMatrix<scalar_t, ROWS, COLS>;
  /// Get const access to a matrix element
  template <unsigned int ROWS, unsigned int COLS>
  ALGEBRA_HOST inline scalar_t operator()(const matrix_type<ROWS, COLS> &m,
                                          unsigned int row,
                                          unsigned int col) const {
    return m(row, col);
  }
};  // element_getter

}  // namespace internal

/** Transform wrapper class to ensure standard API within differnt plugins
 *
 **/
template <typename scalar_t>
struct transform3 {

  /// @name Type definitions for the struct
  /// @{

  /// Array type used by the transform
  template <typename T, std::size_t N>
  using array_type = ROOT::Math::SVector<T, N>;
  /// Scalar type used by the transform
  using scalar_type = scalar_t;

  /// 3-element "vector" type
  using vector3 = array_type<scalar_type, 3>;
  /// Point in 3D space
  using point3 = vector3;
  /// Point in 2D space
  using point2 = array_type<scalar_type, 2>;

  /// 4x4 matrix type
  using matrix44 = ROOT::Math::SMatrix<scalar_type, 4, 4>;

  /// Function (object) used for accessing a matrix element
  using element_getter = internal::element_getter<scalar_type>;

  /// @}

  /// @name Data objects
  /// @{

  matrix44 _data = ROOT::Math::SMatrixIdentity();
  matrix44 _data_inv = ROOT::Math::SMatrixIdentity();

  /// @}

  /** Contructor with arguments: t, z, x
   *
   * @param t the translation (or origin of the new frame)
   * @param z the z axis of the new frame, normal vector for planes
   * @param x the x axis of the new frame
   *
   **/
  ALGEBRA_HOST
  transform3(const vector3 &t, const vector3 &z, const vector3 &x) {

    auto y = ROOT::Math::Cross(z, x);

    _data(0, 0) = x[0];
    _data(1, 0) = x[1];
    _data(2, 0) = x[2];
    _data(0, 1) = y[0];
    _data(1, 1) = y[1];
    _data(2, 1) = y[2];
    _data(0, 2) = z[0];
    _data(1, 2) = z[1];
    _data(2, 2) = z[2];
    _data(0, 3) = t[0];
    _data(1, 3) = t[1];
    _data(2, 3) = t[2];

    int ifail = 0;
    _data_inv = _data.Inverse(ifail);
    SMATRIX_CHECK(ifail);
  }

  /** Constructor with arguments: translation
   *
   * @param t is the translation
   **/
  ALGEBRA_HOST
  transform3(const vector3 &t) {

    _data(0, 3) = t[0];
    _data(1, 3) = t[1];
    _data(2, 3) = t[2];

    int ifail = 0;
    _data_inv = _data.Inverse(ifail);
    SMATRIX_CHECK(ifail);
  }

  /** Constructor with arguments: matrix
   *
   * @param m is the full 4x4 matrix
   **/
  ALGEBRA_HOST
  transform3(const matrix44 &m) {
    _data = m;

    int ifail = 0;
    _data_inv = _data.Inverse(ifail);
    SMATRIX_CHECK(ifail);
  }

  /** Constructor with arguments: matrix as std::aray of scalar
   *
   * @param ma is the full 4x4 matrix asa 16 array
   **/
  ALGEBRA_HOST
  transform3(const array_type<scalar_type, 16> &ma) {

    _data(0, 0) = ma[0];
    _data(1, 0) = ma[4];
    _data(2, 0) = ma[8];
    _data(3, 0) = ma[12];
    _data(0, 1) = ma[1];
    _data(1, 1) = ma[5];
    _data(2, 1) = ma[9];
    _data(3, 1) = ma[13];
    _data(0, 2) = ma[2];
    _data(1, 2) = ma[6];
    _data(2, 2) = ma[10];
    _data(3, 2) = ma[14];
    _data(0, 3) = ma[3];
    _data(1, 3) = ma[7];
    _data(2, 3) = ma[11];
    _data(3, 3) = ma[15];

    int ifail = 0;
    _data_inv = _data.Inverse(ifail);
    // Ignore failures here, since the unit test does manage to trigger an error
    // from ROOT in this place...
  }

  /** Default contructors */
  transform3() = default;
  transform3(const transform3 &rhs) = default;
  ~transform3() = default;

  /** Equality operator */
  ALGEBRA_HOST
  inline bool operator==(const transform3 &rhs) const {

    return _data == rhs._data;
  }

  /** This method retrieves the rotation of a transform */
  ALGEBRA_HOST
  inline auto rotation() const {

    return (_data.template Sub<ROOT::Math::SMatrix<scalar_type, 3, 3> >(0, 0));
  }

  /** This method retrieves the translation of a transform */
  ALGEBRA_HOST
  inline vector3 translation() const {

    return (_data.template SubCol<vector3>(3, 0));
  }

  /** This method retrieves the 4x4 matrix of a transform */
  ALGEBRA_HOST
  inline matrix44 matrix() const { return _data; }

  /** This method transform from a point from the local 3D cartesian frame to
   * the global 3D cartesian frame */
  ALGEBRA_HOST
  inline const point3 point_to_global(const point3 &v) const {

    ROOT::Math::SVector<scalar, 4> vector_4 =
        ROOT::Math::SVector<scalar_type, 4>();
    vector_4.Place_at(v, 0);
    vector_4[3] = static_cast<scalar_type>(1);
    return ROOT::Math::SVector<scalar_type, 4>(_data * vector_4)
        .template Sub<point3>(0);
  }

  /** This method transform from a vector from the global 3D cartesian frame
   * into the local 3D cartesian frame */
  ALGEBRA_HOST
  inline const point3 point_to_local(const point3 &v) const {

    ROOT::Math::SVector<scalar, 4> vector_4 =
        ROOT::Math::SVector<scalar_type, 4>();
    vector_4.Place_at(v, 0);
    vector_4[3] = static_cast<scalar_type>(1);
    return ROOT::Math::SVector<scalar_type, 4>(_data_inv * vector_4)
        .template Sub<point3>(0);
  }

  /** This method transform from a vector from the local 3D cartesian frame to
   * the global 3D cartesian frame */
  ALGEBRA_HOST
  inline const point3 vector_to_global(const vector3 &v) const {

    ROOT::Math::SVector<scalar, 4> vector_4 =
        ROOT::Math::SVector<scalar_type, 4>();
    vector_4.Place_at(v, 0);
    return ROOT::Math::SVector<scalar_type, 4>(_data * vector_4)
        .template Sub<point3>(0);
  }

  /** This method transform from a vector from the global 3D cartesian frame
   * into the local 3D cartesian frame */
  ALGEBRA_HOST
  inline const point3 vector_to_local(const vector3 &v) const {

    ROOT::Math::SVector<scalar, 4> vector_4 =
        ROOT::Math::SVector<scalar_type, 4>();
    vector_4.Place_at(v, 0);
    return ROOT::Math::SVector<scalar_type, 4>(_data_inv * vector_4)
        .template Sub<point3>(0);
  }
};  // struct transform3

}  // namespace algebra::smatrix::math
