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
#include "algebra/storage/matrix44.hpp"
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
/// @note this is templated on the storage array type to be reused for Vc SoA
template <template <typename, std::size_t> class array_t, typename scalar_t>
struct transform3 {

  /// @name Type definitions for the struct
  /// @{

  /// Scalar type used by the transform
  using scalar_type = scalar_t;
  /// The type of the matrix elements (scalar for AoS, Vc::Vector for SoA)
  using value_type =
      std::conditional_t<Vc::is_simd_vector<array_t<scalar_t, 4>>::value,
                         scalar_t, Vc::Vector<scalar_t>>;

  // AoS stores four elements per vector for alignment
  using storage_dim =
      std::conditional_t<Vc::is_simd_vector<value_type>::value,
                         std::integral_constant<std::size_t, 3u>,
                         std::integral_constant<std::size_t, 4u>>;

  template <std::size_t N>
  using array_type = array_t<value_type, N>;

  /// 3-element "vector" type (does not observe translations)
  using vector3 = storage::vector<storage_dim::value, value_type, array_t>;
  /// Point in 3D space (does observe translations)
  using point3 = vector3;
  /// Point in 2D space
  using point2 = storage::vector<2, value_type, array_t>;

  /// 4x4 matrix type
  using matrix44 = storage::matrix44<array_t, value_type, storage_dim::value>;

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
  transform3() = default;

  /// Contructor with arguments: t, x, y, z
  ///
  /// @param t the translation (or origin of the new frame)
  /// @param x the x axis of the new frame
  /// @param y the y axis of the new frame
  /// @param z the z axis of the new frame, normal vector for planes
  ALGEBRA_HOST_DEVICE
  transform3(const vector3 &t, const vector3 &x, const vector3 &y,
             const vector3 &z, [[maybe_unused]] bool get_inverse = true)
      : _data{x, y, z, t}, _data_inv{invert(_data)} {}

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
      : transform3(t, x,
                   vector3{z[1] * x[2] - x[1] * z[2], z[2] * x[0] - x[2] * z[0],
                           z[0] * x[1] - x[0] * z[1]},
                   z, get_inverse) {}

  /// Constructor with arguments: translation
  ///
  /// @param t is the transform
  ALGEBRA_HOST_DEVICE
  transform3(const vector3 &t) : _data{t}, _data_inv{invert(_data)} {}

  /// Constructor with arguments: matrix
  ///
  /// @param m is the full 4x4 matrix with simd-vector elements
  ALGEBRA_HOST_DEVICE
  transform3(const matrix44 &m) : _data{m}, _data_inv{invert(_data)} {}

  /// Constructor with arguments: matrix as std::aray of scalar
  ///
  /// @param ma is the full 4x4 matrix 16 array
  ALGEBRA_HOST_DEVICE
  transform3(const array_type<16> &ma) {

    _data.x = typename matrix44::vector_type{ma[0], ma[4], ma[8], ma[12]};
    _data.y = typename matrix44::vector_type{ma[1], ma[5], ma[9], ma[13]};
    _data.z = typename matrix44::vector_type{ma[2], ma[6], ma[10], ma[14]};
    _data.t = typename matrix44::vector_type{ma[3], ma[7], ma[11], ma[15]};
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
  inline const value_type &operator()(std::size_t i, std::size_t j) const {
    return element_getter{}(_data, i, j);
  }
  ALGEBRA_HOST_DEVICE
  inline value_type &operator()(std::size_t i, std::size_t j) {
    return element_getter{}(_data, i, j);
  }

  /// The determinant of a 4x4 matrix
  ///
  /// @param m is the matrix
  ///
  /// @return a sacalar determinant - no checking done
  ALGEBRA_HOST_DEVICE
  inline constexpr value_type determinant(const matrix44 &m) const {
    return -m.z[0] * m.y[1] * m.x[2] + m.y[0] * m.z[1] * m.x[2] +
           m.z[0] * m.x[1] * m.y[2] - m.x[0] * m.z[1] * m.y[2] -
           m.y[0] * m.x[1] * m.z[2] + m.x[0] * m.y[1] * m.z[2];
  }

  /// The inverse of a 4x4 matrix
  ///
  /// @param m is the matrix
  ///
  /// @return an inverse matrix
  ALGEBRA_HOST_DEVICE
  inline constexpr matrix44 invert(const matrix44 &m) const {
    matrix44 i;
    i.x[0] = -m.z[1] * m.y[2] + m.y[1] * m.z[2];
    i.x[1] = m.z[1] * m.x[2] - m.x[1] * m.z[2];
    i.x[2] = -m.y[1] * m.x[2] + m.x[1] * m.y[2];
    // i.x[3] = 0;
    i.y[0] = m.z[0] * m.y[2] - m.y[0] * m.z[2];
    i.y[1] = -m.z[0] * m.x[2] + m.x[0] * m.z[2];
    i.y[2] = m.y[0] * m.x[2] - m.x[0] * m.y[2];
    // i.y[3] = 0;
    i.z[0] = -m.z[0] * m.y[1] + m.y[0] * m.z[1];
    i.z[1] = m.z[0] * m.x[1] - m.x[0] * m.z[1];
    i.z[2] = -m.y[0] * m.x[1] + m.x[0] * m.y[1];
    // i.z[3] = 0;
    i.t[0] = m.t[0] * m.z[1] * m.y[2] - m.z[0] * m.t[1] * m.y[2] -
             m.t[0] * m.y[1] * m.z[2] + m.y[0] * m.t[1] * m.z[2] +
             m.z[0] * m.y[1] * m.t[2] - m.y[0] * m.z[1] * m.t[2];
    i.t[1] = m.z[0] * m.t[1] * m.x[2] - m.t[0] * m.z[1] * m.x[2] +
             m.t[0] * m.x[1] * m.z[2] - m.x[0] * m.t[1] * m.z[2] -
             m.z[0] * m.x[1] * m.t[2] + m.x[0] * m.z[1] * m.t[2];
    i.t[2] = m.t[0] * m.y[1] * m.x[2] - m.y[0] * m.t[1] * m.x[2] -
             m.t[0] * m.x[1] * m.y[2] + m.x[0] * m.t[1] * m.y[2] +
             m.y[0] * m.x[1] * m.t[2] - m.x[0] * m.y[1] * m.t[2];
    // i.t[3] = 1;
    const value_type idet{value_type(1.f) / determinant(i)};

    i.x = i.x * idet;
    i.y = i.y * idet;
    i.z = i.z * idet;
    i.t = i.t * idet;

    return i;
  }

  /// Rotate a vector into / from a frame
  ///
  /// @param m is the rotation matrix
  /// @param v is the vector to be rotated
  ALGEBRA_HOST_DEVICE
  inline constexpr auto rotate(const matrix44 &m, const vector3 &v) const {

    return m.x * v[0] + m.y * v[1] + m.z * v[2];
  }

  /// This method retrieves the rotation of a transform
  ALGEBRA_HOST_DEVICE
  inline array_type<16> rotation() const {

    array_type<16> submatrix;
    for (unsigned int irow = 0; irow < 3; ++irow) {
      for (unsigned int icol = 0; icol < 3; ++icol) {
        submatrix[icol + irow * 4] = element_getter()(_data, irow, icol);
      }
    }
    return submatrix;
  }

  /// This method retrieves x axis
  ALGEBRA_HOST_DEVICE
  inline const auto &x() const { return _data.x; }

  /// This method retrieves y axis
  ALGEBRA_HOST_DEVICE
  inline const auto &y() const { return _data.y; }

  /// This method retrieves z axis
  ALGEBRA_HOST_DEVICE
  inline const auto &z() const { return _data.z; }

  /// This method retrieves the translation
  ALGEBRA_HOST_DEVICE
  inline const auto &translation() const { return _data.t; }

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

    return _data.x * p[0] + _data.y * p[1] + _data.z * p[2] + _data.t;
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

    return _data_inv.x * p[0] + _data_inv.y * p[1] + _data_inv.z * p[2] +
           _data_inv.t;
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
