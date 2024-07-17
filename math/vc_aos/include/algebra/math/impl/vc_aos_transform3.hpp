/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/common.hpp"
#include "algebra/qualifiers.hpp"
#include "algebra/storage/matrix.hpp"
#include "algebra/storage/vector.hpp"

// Vc include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#include <Vc/Vc>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC

// System include(s).
#include <cassert>

namespace algebra::vc_aos::math {

using algebra::storage::operator*;
using algebra::storage::operator/;
using algebra::storage::operator-;
using algebra::storage::operator+;

/// Transform wrapper class to ensure standard API within differnt plugins
template <template <typename, std::size_t> class array_t, typename scalar_t>
struct transform3 {

  // How to address matrix columns
  enum elem : std::size_t {
    e_x = 0u,
    e_y = 1u,
    e_z = 2u,
    e_t = 3u,
  };

  /// @name Type definitions for the struct
  /// @{

  /// Scalar type used by the transform
  using scalar_type = scalar_t;
  /// The type of the matrix elements (scalar for AoS, Vc::Vector for SoA)
  using value_type =
      std::conditional_t<Vc::is_simd_vector<array_t<scalar_t, 4>>::value,
                         scalar_t, Vc::Vector<scalar_t>>;

  template <std::size_t N>
  using array_type = array_t<value_type, N>;

  /// 3-element "vector" type (does not observe translations)
  using vector3 = storage::vector<3u, value_type, array_t>;
  /// Point in 3D space (does observe translations)
  using point3 = vector3;
  /// Point in 2D space
  using point2 = storage::vector<2u, value_type, array_t>;

  /// 4x4 matrix type
  using matrix44 = storage::matrix<array_t, value_type, 3u, 4u>;
  using column_t = typename matrix44::vector_type;

  /// Function (object) used for accessing a matrix element
  using element_getter = storage::element_getter;

  /// @}

  /// @name Data objects
  /// @{

  matrix44 _data;
  matrix44 _data_inv;

  /// @}
  /// Default constructor: identity
  ALGEBRA_HOST_DEVICE
  constexpr transform3()
      : _data{storage::identity<matrix44>()},
        _data_inv{storage::identity<matrix44>()} {}

  /// Contructor with arguments: t, x, y, z
  ///
  /// @param t the translation (or origin of the new frame)
  /// @param x the x axis of the new frame
  /// @param y the y axis of the new frame
  /// @param z the z axis of the new frame, normal vector for planes
  ALGEBRA_HOST_DEVICE
  transform3(const vector3 &t, const vector3 &x, const vector3 &y,
             const vector3 &z)
      : _data{x, y, z, t}, _data_inv{invert(_data)} {}

  /// Contructor with arguments: t, z, x
  ///
  /// @param t the translation (or origin of the new frame)
  /// @param z the z axis of the new frame, normal vector for planes
  /// @param x the x axis of the new frame
  ///
  /// @note y will be constructed by cross product
  ALGEBRA_HOST_DEVICE
  transform3(const vector3 &t, const vector3 &z, const vector3 &x)
      : transform3(
            t, x,
            column_t(z[1] * x[2] - x[1] * z[2], z[2] * x[0] - x[2] * z[0],
                     z[0] * x[1] - x[0] * z[1]),
            z) {}

  /// Constructor with arguments: translation
  ///
  /// @param t is the transform
  ALGEBRA_HOST_DEVICE
  explicit transform3(const vector3 &t)
      : _data{column_t{1.f, 0.f, 0.f}, column_t{0.f, 1.f, 0.f},
              column_t{0.f, 0.f, 1.f}, t},
        _data_inv{invert(_data)} {}

  /// Constructor with arguments: matrix
  ///
  /// @param m is the full 4x4 matrix with simd-vector elements
  ALGEBRA_HOST_DEVICE
  explicit transform3(const matrix44 &m) : _data{m}, _data_inv{invert(_data)} {}

  /// Constructor with arguments: matrix as std::aray of scalar
  ///
  /// @param ma is the full 4x4 matrix 16 array
  ALGEBRA_HOST_DEVICE
  explicit transform3(const array_type<16> &ma) {

    // The values that are not set here, are known to be zero or one
    // and never used explicitly
    _data[e_x] = column_t{ma[0], ma[4], ma[8]};
    _data[e_y] = column_t{ma[1], ma[5], ma[9]};
    _data[e_z] = column_t{ma[2], ma[6], ma[10]};
    _data[e_t] = column_t{ma[3], ma[7], ma[11]};

    _data_inv = invert(_data);
  }

  /// Defaults
  transform3(const transform3 &rhs) = default;
  ~transform3() = default;

  /// Equality operator
  ALGEBRA_HOST_DEVICE
  inline constexpr bool operator==(const transform3 &rhs) const {
    return (_data == rhs._data);
  }

  /// Matrix access operator
  ALGEBRA_HOST_DEVICE
  inline const value_type &operator()(std::size_t row, std::size_t col) const {
    return _data[col][row];
  }
  ALGEBRA_HOST_DEVICE
  inline value_type &operator()(std::size_t row, std::size_t col) {
    return _data[col][row];
  }

  /// The determinant of a 4x4 matrix
  ///
  /// @param m is the matrix
  ///
  /// @return a sacalar determinant - no checking done
  ALGEBRA_HOST_DEVICE
  inline constexpr value_type determinant(const matrix44 &m) const {
    return -m[e_z][0] * m[e_y][1] * m[e_x][2] +
           m[e_y][0] * m[e_z][1] * m[e_x][2] +
           m[e_z][0] * m[e_x][1] * m[e_y][2] -
           m[e_x][0] * m[e_z][1] * m[e_y][2] -
           m[e_y][0] * m[e_x][1] * m[e_z][2] +
           m[e_x][0] * m[e_y][1] * m[e_z][2];
  }

  /// The inverse of a 4x4 matrix
  ///
  /// @param m is the matrix
  ///
  /// @return an inverse matrix
  ALGEBRA_HOST_DEVICE
  inline constexpr matrix44 invert(const matrix44 &m) const {
    matrix44 i;
    i[e_x][0] = -m[e_z][1] * m[e_y][2] + m[e_y][1] * m[e_z][2];
    i[e_x][1] = m[e_z][1] * m[e_x][2] - m[e_x][1] * m[e_z][2];
    i[e_x][2] = -m[e_y][1] * m[e_x][2] + m[e_x][1] * m[e_y][2];
    // i[e_x][3] = 0;
    i[e_y][0] = m[e_z][0] * m[e_y][2] - m[e_y][0] * m[e_z][2];
    i[e_y][1] = -m[e_z][0] * m[e_x][2] + m[e_x][0] * m[e_z][2];
    i[e_y][2] = m[e_y][0] * m[e_x][2] - m[e_x][0] * m[e_y][2];
    // i[e_y][3] = 0;
    i[e_z][0] = -m[e_z][0] * m[e_y][1] + m[e_y][0] * m[e_z][1];
    i[e_z][1] = m[e_z][0] * m[e_x][1] - m[e_x][0] * m[e_z][1];
    i[e_z][2] = -m[e_y][0] * m[e_x][1] + m[e_x][0] * m[e_y][1];
    // i[e_z][3] = 0;
    i[e_t][0] =
        m[e_t][0] * m[e_z][1] * m[e_y][2] - m[e_z][0] * m[e_t][1] * m[e_y][2] -
        m[e_t][0] * m[e_y][1] * m[e_z][2] + m[e_y][0] * m[e_t][1] * m[e_z][2] +
        m[e_z][0] * m[e_y][1] * m[e_t][2] - m[e_y][0] * m[e_z][1] * m[e_t][2];
    i[e_t][1] =
        m[e_z][0] * m[e_t][1] * m[e_x][2] - m[e_t][0] * m[e_z][1] * m[e_x][2] +
        m[e_t][0] * m[e_x][1] * m[e_z][2] - m[e_x][0] * m[e_t][1] * m[e_z][2] -
        m[e_z][0] * m[e_x][1] * m[e_t][2] + m[e_x][0] * m[e_z][1] * m[e_t][2];
    i[e_t][2] =
        m[e_t][0] * m[e_y][1] * m[e_x][2] - m[e_y][0] * m[e_t][1] * m[e_x][2] -
        m[e_t][0] * m[e_x][1] * m[e_y][2] + m[e_x][0] * m[e_t][1] * m[e_y][2] +
        m[e_y][0] * m[e_x][1] * m[e_t][2] - m[e_x][0] * m[e_y][1] * m[e_t][2];
    // i[e_t][3] = 1;
    const value_type idet{value_type(1.f) / determinant(i)};

    i[e_x] = i[e_x] * idet;
    i[e_y] = i[e_y] * idet;
    i[e_z] = i[e_z] * idet;
    i[e_t] = i[e_t] * idet;

    return i;
  }

  /// Rotate a vector into / from a frame
  ///
  /// @param m is the rotation matrix
  /// @param v is the vector to be rotated
  ALGEBRA_HOST_DEVICE
  inline constexpr auto rotate(const matrix44 &m, const vector3 &v) const {

    return m[e_x] * v[0] + m[e_y] * v[1] + m[e_z] * v[2];
  }

  /// This method retrieves the rotation of a transform
  ALGEBRA_HOST_DEVICE
  inline array_type<16> rotation() const {

    array_type<16> submatrix;
    for (unsigned int irow = 0; irow < 3; ++irow) {
      for (unsigned int icol = 0; icol < 3; ++icol) {
        submatrix[icol + irow * 4] = _data[icol][irow];
      }
    }
    return submatrix;
  }

  /// This method retrieves the new x-axis
  ALGEBRA_HOST_DEVICE
  inline const auto &x() const { return _data[e_x]; }

  /// This method retrieves the new y-axis
  ALGEBRA_HOST_DEVICE
  inline const auto &y() const { return _data[e_y]; }

  /// This method retrieves the new z-axis
  ALGEBRA_HOST_DEVICE
  inline const auto &z() const { return _data[e_z]; }

  /// This method retrieves the translation
  ALGEBRA_HOST_DEVICE
  inline const auto &translation() const { return _data[e_t]; }

  /// This method retrieves the 4x4 matrix of a transform
  ALGEBRA_HOST_DEVICE
  inline const matrix44 &matrix() const { return _data; }

  /// This method retrieves the 4x4 matrix of an inverse transform
  ALGEBRA_HOST_DEVICE
  inline const matrix44 &matrix_inverse() const { return _data_inv; }

  /// This method transform from a point from the local 3D cartesian frame
  ///  to the global 3D cartesian frame
  ///
  /// @tparam point_type 3D point
  ///
  /// @param v is the point to be transformed
  ///
  /// @return a global point
  template <typename point3_type>
  ALGEBRA_HOST_DEVICE inline auto point_to_global(const point3_type &p) const {

    return _data[e_x] * p[0] + _data[e_y] * p[1] + _data[e_z] * p[2] +
           _data[e_t];
  }

  /// This method transform from a vector from the global 3D cartesian frame
  ///  into the local 3D cartesian frame
  ///
  /// @tparam point_type 3D point
  ///
  /// @param v is the point to be transformed
  ///
  /// @return a local point
  template <typename point3_type>
  ALGEBRA_HOST_DEVICE inline auto point_to_local(const point3_type &p) const {

    return _data_inv[e_x] * p[0] + _data_inv[e_y] * p[1] +
           _data_inv[e_z] * p[2] + _data_inv[e_t];
  }

  /// This method transform from a vector from the local 3D cartesian frame
  ///  to the global 3D cartesian frame
  ///
  /// @tparam vector3_type 3D vector
  ///
  /// @param v is the vector to be transformed
  ///
  /// @return a vector in global coordinates
  template <typename vector3_type>
  ALGEBRA_HOST_DEVICE inline auto vector_to_global(
      const vector3_type &v) const {

    return rotate(_data, v);
  }

  /// This method transform from a vector from the global 3D cartesian frame
  ///  into the local 3D cartesian frame
  ///
  /// @tparam vector3_type 3D vector
  ///
  /// @param v is the vector to be transformed
  ///
  /// @return a vector in global coordinates
  template <typename vector3_type>
  ALGEBRA_HOST_DEVICE inline auto vector_to_local(const vector3_type &v) const {

    return rotate(_data_inv, v);
  }
};  // struct transform3

}  // namespace algebra::vc_aos::math
