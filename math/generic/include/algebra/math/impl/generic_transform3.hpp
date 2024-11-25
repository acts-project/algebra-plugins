/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/algorithms/matrix/inverse/hard_coded.hpp"
#include "algebra/math/impl/generic_matrix.hpp"
#include "algebra/math/impl/generic_vector.hpp"
#include "algebra/qualifiers.hpp"
#include "algebra/type_traits.hpp"

namespace algebra::generic::math {

/// Transform wrapper class to ensure standard API within differnt plugins
template <concepts::index index_t, concepts::scalar scalar_t,
          template <typename, index_t, index_t> class matrix_t,
          template <typename, index_t> class array_t>
struct transform3 {

  /// @name Type definitions for the struct
  /// @{

  /// Scalar type
  using scalar_type = scalar_t;

  /// Array type used by the transform
  template <index_t N>
  using array_type = array_t<scalar_t, N>;

  // 4 x 4 Matrix (keep four rows for better alignment)
  using matrix44 = matrix_t<scalar_t, 4, 4>;
  static_assert(concepts::square_matrix<matrix44>);

  /// 3-element "vector" type
  using vector3 = algebra::traits::get_vector_t<matrix44, 3, scalar_t>;
  /// Point in 3D space
  using point3 = vector3;
  /// Point in 2D space
  using point2 = algebra::traits::get_vector_t<matrix44, 2, scalar_t>;

  /// Function (object) used for accessing a matrix element
  using element_getter = algebra::traits::element_getter_t<matrix44>;

  /// Function (object) used for accessing a matrix element
  using block_getter = algebra::traits::block_getter_t<matrix44>;

  /// Matrix inversion algorithm
  using matrix_inversion =
      generic::matrix::inverse::hard_coded<matrix44, element_getter>;

  /// @}

  /// @name Data objects
  /// @{
  matrix44 _data{generic::math::identity<matrix44>()};
  matrix44 _data_inv{generic::math::identity<matrix44>()};

  /// @}

  /// Default constructor: identity
  constexpr transform3() = default;

  /// Contructor with arguments: t, x, y, z
  ///
  /// @param t the translation (or origin of the new frame)
  /// @param x the x axis of the new frame
  /// @param y the y axis of the new frame
  /// @param z the z axis of the new frame, normal vector for planes
  ALGEBRA_HOST_DEVICE
  transform3(const vector3 &t, const vector3 &x, const vector3 &y,
             const vector3 &z, bool get_inverse = true) {

    element_getter{}(_data, 0, 0) = element_getter{}(x, 0);
    element_getter{}(_data, 1, 0) = element_getter{}(x, 1);
    element_getter{}(_data, 2, 0) = element_getter{}(x, 2);
    element_getter{}(_data, 3, 0) = 0.f;
    element_getter{}(_data, 0, 1) = element_getter{}(y, 0);
    element_getter{}(_data, 1, 1) = element_getter{}(y, 1);
    element_getter{}(_data, 2, 1) = element_getter{}(y, 2);
    element_getter{}(_data, 3, 1) = 0.f;
    element_getter{}(_data, 0, 2) = element_getter{}(z, 0);
    element_getter{}(_data, 1, 2) = element_getter{}(z, 1);
    element_getter{}(_data, 2, 2) = element_getter{}(z, 2);
    element_getter{}(_data, 3, 2) = 0.f;
    element_getter{}(_data, 0, 3) = element_getter{}(t, 0);
    element_getter{}(_data, 1, 3) = element_getter{}(t, 1);
    element_getter{}(_data, 2, 3) = element_getter{}(t, 2);
    element_getter{}(_data, 3, 3) = 1.f;

    if (get_inverse) {
      _data_inv = matrix_inversion{}(_data);
    }
  }

  /// Contructor with arguments: t, z, x
  ///
  /// @param t the translation (or origin of the new frame)
  /// @param z the z axis of the new frame, normal vector for planes
  /// @param x the x axis of the new frame
  ///
  /// @note y will be constructed by cross product
  ALGEBRA_HOST_DEVICE
  transform3(const vector3 &t, const vector3 &z, const vector3 &x,
             bool get_inverse = true)
      : transform3(t, x, cross(z, x), z, get_inverse) {}

  /// Constructor with arguments: translation
  ///
  /// @param t is the transform
  ALGEBRA_HOST_DEVICE
  explicit transform3(const vector3 &t) {

    element_getter{}(_data, 0, 0) = 1.f;
    element_getter{}(_data, 1, 0) = 0.f;
    element_getter{}(_data, 2, 0) = 0.f;
    element_getter{}(_data, 3, 0) = 0.f;
    element_getter{}(_data, 0, 1) = 0.f;
    element_getter{}(_data, 1, 1) = 1.f;
    element_getter{}(_data, 2, 1) = 0.f;
    element_getter{}(_data, 3, 1) = 0.f;
    element_getter{}(_data, 0, 2) = 0.f;
    element_getter{}(_data, 1, 2) = 0.f;
    element_getter{}(_data, 2, 2) = 1.f;
    element_getter{}(_data, 3, 2) = 0.f;
    element_getter{}(_data, 0, 3) = element_getter{}(t, 0);
    element_getter{}(_data, 1, 3) = element_getter{}(t, 1);
    element_getter{}(_data, 2, 3) = element_getter{}(t, 2);
    element_getter{}(_data, 3, 3) = 1.f;

    _data_inv = matrix_inversion{}(_data);
  }

  /// Constructor with arguments: matrix
  ///
  /// @param m is the full 4x4 matrix
  ALGEBRA_HOST_DEVICE
  explicit transform3(const matrix44 &m) : _data{m} {
    _data_inv = matrix_inversion{}(_data);
  }

  /// Constructor with arguments: matrix as array of scalar
  ///
  /// @param ma is the full 4x4 matrix 16 array
  ALGEBRA_HOST_DEVICE
  explicit transform3(const array_type<16> &ma) {

    element_getter{}(_data, 0, 0) = ma[0];
    element_getter{}(_data, 1, 0) = ma[4];
    element_getter{}(_data, 2, 0) = ma[8];
    element_getter{}(_data, 3, 0) = 0.f;
    element_getter{}(_data, 0, 1) = ma[1];
    element_getter{}(_data, 1, 1) = ma[5];
    element_getter{}(_data, 2, 1) = ma[9];
    element_getter{}(_data, 3, 1) = 0.f;
    element_getter{}(_data, 0, 2) = ma[2];
    element_getter{}(_data, 1, 2) = ma[6];
    element_getter{}(_data, 2, 2) = ma[10];
    element_getter{}(_data, 3, 2) = 0.f;
    element_getter{}(_data, 0, 3) = ma[3];
    element_getter{}(_data, 1, 3) = ma[7];
    element_getter{}(_data, 2, 3) = ma[11];
    element_getter{}(_data, 3, 3) = 1.f;

    _data_inv = matrix_inversion{}(_data);
  }

  /// Equality operator
  ALGEBRA_HOST_DEVICE
  inline bool operator==(const transform3 &rhs) const {

    for (index_t j = 0; j < 4; j++) {
      // Check only the rows that can differ
      for (index_t i = 0; i < 3; i++) {
        if (element_getter{}(_data, i, j) !=
            element_getter{}(rhs._data, i, j)) {
          return false;
        }
      }
    }

    return true;
  }

  /// Rotate a vector into / from a frame
  ///
  /// @param m is the rotation matrix
  /// @param v is the vector to be rotated
  ALGEBRA_HOST_DEVICE
  static inline vector3 rotate(const matrix44 &m, const vector3 &v) {

    vector3 ret{0.f, 0.f, 0.f};

    element_getter{}(ret, 0) +=
        element_getter{}(m, 0, 0) * element_getter{}(v, 0);
    element_getter{}(ret, 1) +=
        element_getter{}(m, 1, 0) * element_getter{}(v, 0);
    element_getter{}(ret, 2) +=
        element_getter{}(m, 2, 0) * element_getter{}(v, 0);

    element_getter{}(ret, 0) +=
        element_getter{}(m, 0, 1) * element_getter{}(v, 1);
    element_getter{}(ret, 1) +=
        element_getter{}(m, 1, 1) * element_getter{}(v, 1);
    element_getter{}(ret, 2) +=
        element_getter{}(m, 2, 1) * element_getter{}(v, 1);

    element_getter{}(ret, 0) +=
        element_getter{}(m, 0, 2) * element_getter{}(v, 2);
    element_getter{}(ret, 1) +=
        element_getter{}(m, 1, 2) * element_getter{}(v, 2);
    element_getter{}(ret, 2) +=
        element_getter{}(m, 2, 2) * element_getter{}(v, 2);

    return ret;
  }

  /// This method retrieves the rotation of a transform
  ALGEBRA_HOST_DEVICE
  auto inline rotation() const {
    return block_getter{}.template operator()<3, 3>(_data, 0, 0);
  }

  /// This method retrieves x axis
  ALGEBRA_HOST_DEVICE
  inline point3 x() const {
    return {element_getter{}(_data, 0, 0), element_getter{}(_data, 1, 0),
            element_getter{}(_data, 2, 0)};
  }

  /// This method retrieves y axis
  ALGEBRA_HOST_DEVICE
  inline point3 y() const {
    return {element_getter{}(_data, 0, 1), element_getter{}(_data, 1, 1),
            element_getter{}(_data, 2, 1)};
  }

  /// This method retrieves z axis
  ALGEBRA_HOST_DEVICE
  inline point3 z() const {
    return {element_getter{}(_data, 0, 2), element_getter{}(_data, 1, 2),
            element_getter{}(_data, 2, 2)};
  }

  /// This method retrieves the translation of a transform
  ALGEBRA_HOST_DEVICE
  inline point3 translation() const {
    return {element_getter{}(_data, 0, 3), element_getter{}(_data, 1, 3),
            element_getter{}(_data, 2, 3)};
  }

  /// This method retrieves the 4x4 matrix of a transform
  ALGEBRA_HOST_DEVICE
  inline const matrix44 &matrix() const { return _data; }

  /// This method retrieves the 4x4 matrix of an inverse transform
  ALGEBRA_HOST_DEVICE
  inline const matrix44 &matrix_inverse() const { return _data_inv; }

  /// This method transform from a point from the local 3D cartesian frame to
  /// the global 3D cartesian frame
  ALGEBRA_HOST_DEVICE inline point3 point_to_global(const point3 &v) const {

    const vector3 rg = rotate(_data, v);

    return {element_getter{}(rg, 0) + element_getter{}(_data, 0, 3),
            element_getter{}(rg, 1) + element_getter{}(_data, 1, 3),
            element_getter{}(rg, 2) + element_getter{}(_data, 2, 3)};
  }

  /// This method transform from a vector from the global 3D cartesian frame
  /// into the local 3D cartesian frame
  ALGEBRA_HOST_DEVICE inline point3 point_to_local(const point3 &v) const {

    const vector3 rg = rotate(_data_inv, v);

    return {element_getter{}(rg, 0) + element_getter{}(_data_inv, 0, 3),
            element_getter{}(rg, 1) + element_getter{}(_data_inv, 1, 3),
            element_getter{}(rg, 2) + element_getter{}(_data_inv, 2, 3)};
  }

  /// This method transform from a vector from the local 3D cartesian frame to
  /// the global 3D cartesian frame
  ALGEBRA_HOST_DEVICE inline vector3 vector_to_global(const vector3 &v) const {
    return rotate(_data, v);
  }

  /// This method transform from a vector from the global 3D cartesian frame
  /// into the local 3D cartesian frame
  ALGEBRA_HOST_DEVICE inline vector3 vector_to_local(const vector3 &v) const {
    return rotate(_data_inv, v);
  }

};  // struct transform3

}  // namespace algebra::generic::math
