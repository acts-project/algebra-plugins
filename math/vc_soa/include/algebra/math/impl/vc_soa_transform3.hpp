/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/cmath.hpp"
#include "algebra/qualifiers.hpp"
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
#include <type_traits>

namespace algebra::vc_soa::math {

using algebra::storage::operator*;
using algebra::storage::operator/;
using algebra::storage::operator-;
using algebra::storage::operator+;

namespace internal {

/// 4x4 matrix type used by @c algebra::vc_soa::math::transform3 that has simd
/// vectors as matrix elements
template <typename scalar_t, template <typename, std::size_t> class array_t>
struct matrix44 {

  using vector_type = storage::vector<3, Vc::Vector<scalar_t>, array_t>;
  using value_type = Vc::Vector<scalar_t>;
  using scalar_type = scalar_t;

  /// Default constructor: Identity with no translation
  matrix44()
      : x{value_type::One(), value_type::Zero(), value_type::Zero()},
        y{value_type::Zero(), value_type::One(), value_type::Zero()},
        z{value_type::Zero(), value_type::Zero(), value_type::One()},
        t{value_type::Zero(), value_type::Zero(), value_type::Zero()} {}

  /// Identity rotation with translation @param translation
  matrix44(const vector_type &v)
      : x{value_type::One(), value_type::Zero(), value_type::Zero()},
        y{value_type::Zero(), value_type::One(), value_type::Zero()},
        z{value_type::Zero(), value_type::Zero(), value_type::One()},
        t{v} {}

  /// Construct from given column vectors @param x, @param y, @param z, @param t
  matrix44(const vector_type &v_0, const vector_type &v_1,
           const vector_type &v_2, const vector_type &v_3)
      : x{v_0}, y{v_1}, z{v_2}, t{v_3} {}

  /// Identity rotation with translation from single elemenst @param t_0,
  /// @param t_1, @param t_2
  matrix44(const value_type &t_0, const value_type &t_1, const value_type &t_2)
      : matrix44{{t_0, t_1, t_2}} {}

  /// Construct from elements (simd vector) of matrix column vectors:
  /// - column 0: @param x_0, @param x_1, @param x_2
  /// - column 1: @param y_0, @param y_1, @param y_2
  /// - column 2: @param z_0, @param z_1, @param z_2
  /// - column 3: @param t_0, @param t_1, @param t_2
  matrix44(const value_type &x_0, const value_type &x_1, const value_type &x_2,
           const value_type &y_0, const value_type &y_1, const value_type &y_2,
           const value_type &z_0, const value_type &z_1, const value_type &z_2,
           const value_type &t_0, const value_type &t_1, const value_type &t_2)
      : x{{x_0, x_1, x_2}},
        y{{y_0, y_1, y_2}},
        z{{z_0, z_1, z_2}},
        t{{t_0, t_1, t_2}} {}

  /// Equality operator between two matrices
  bool operator==(const matrix44 &rhs) const {
    return ((x == rhs.x) && (y == rhs.y) && (z == rhs.z) && (t == rhs.t));
  }

  /// Data variables
  vector_type x, y, z, t;

};  // struct matrix44

/// Functor used to access elements of matrix44
template <typename scalar_t, template <typename, std::size_t> class array_t>
struct element_getter {

  /// Get const access to a matrix element
  ALGEBRA_HOST inline const auto &operator()(
      const matrix44<scalar_t, array_t> &m, std::size_t row,
      std::size_t col) const {

    // Make sure that the indices are valid.
    assert(row < 4u);
    assert(col < 4u);

    // Return the selected element.
    switch (col) {
      case 0u:
        return m.x[row];
      case 1u:
        return m.y[row];
      case 2u:
        return m.z[row];
      case 3u:
        return m.t[row];
      default:
        return m.x[0];
    }
  }

  /// Get const access to a matrix element
  ALGEBRA_HOST inline auto &operator()(matrix44<scalar_t, array_t> &m,
                                       std::size_t row, std::size_t col) const {

    // Make sure that the indices are valid.
    assert(row < 4u);
    assert(col < 4u);

    // Return the selected element.
    switch (col) {
      case 0u:
        return m.x[row];
      case 1u:
        return m.y[row];
      case 2u:
        return m.z[row];
      case 3u:
        return m.t[row];
      default:
        return m.x[0];
    }
  }

};  // struct element_getter

}  // namespace internal

/// Transform wrapper class to ensure standard API within differnt plugins
template <template <typename, std::size_t> class array_t, typename scalar_t>
struct transform3 {

  /// @name Type definitions for the struct
  /// @{

  /// Array type used by the transform
  template <std::size_t N>
  using array_type = array_t<scalar_t, N>;
  /// Scalar type used by the transform
  using scalar_type = scalar_t;
  /// The type of the matrix elements (in this case: Vc::Vector)
  using value_type = Vc::Vector<scalar_t>;

  /// 3-element "vector" type (does not observe translations)
  using vector3 = storage::vector<3, value_type, array_t>;
  /// Point in 3D space (does observe translations)
  using point3 = vector3;
  /// Point in 2D space
  using point2 = storage::vector<2, value_type, array_t>;

  /// 4x4 matrix type
  using matrix44 = internal::matrix44<scalar_type, array_t>;

  /// Function (object) used for accessing a matrix element
  using element_getter = internal::element_getter<scalar_type, array_t>;

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
      : transform3(t, x, cross(z, x), z, get_inverse) {}

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
    const value_type idet{value_type::One() / determinant(i)};

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

  /// @returns the translation of the transform
  ALGEBRA_HOST_DEVICE
  inline point3 translation() const { return _data.t; }

  /// @returns the 4x4 matrix of the transform
  ALGEBRA_HOST_DEVICE
  inline const matrix44 &matrix() const { return _data; }

  /// @returns the 4x4 matrix of the inverse transform
  ALGEBRA_HOST_DEVICE
  inline const matrix44 &matrix_inverse() const { return _data_inv; }

  /// This method transforms a point from the local 3D cartesian frame
  /// to the global 3D cartesian frame
  ///
  /// @tparam point3_t 3D point
  ///
  /// @param v is the point to be transformed
  ///
  /// @return a global point
  template <
      typename point3_t,
      std::enable_if_t<std::is_convertible_v<point3_t, point3>, bool> = true>
  ALGEBRA_HOST_DEVICE inline point3 point_to_global(const point3_t &v) const {

    return rotate(_data, v) + _data.t;
  }

  /// This method transforms a point from the global 3D cartesian frame
  /// into the local 3D cartesian frame
  ///
  /// @tparam point3_t 3D point
  ///
  /// @param v is the point to be transformed
  ///
  /// @return a local point
  template <
      typename point3_t,
      std::enable_if_t<std::is_convertible_v<point3_t, point3>, bool> = true>
  ALGEBRA_HOST_DEVICE inline point3 point_to_local(const point3_t &v) const {

    return rotate(_data_inv, v) + _data_inv.t;
  }

  /// This method transforms a vector from the local 3D cartesian frame
  /// to the global 3D cartesian frame
  ///
  /// @tparam vector3_t 3D vector
  ///
  /// @param v is the vector to be transformed
  ///
  /// @return a vector in global coordinates
  template <
      typename vector3_t,
      std::enable_if_t<std::is_convertible_v<vector3_t, vector3>, bool> = true>
  ALGEBRA_HOST_DEVICE inline vector3 vector_to_global(
      const vector3_t &v) const {

    return rotate(_data, v);
  }

  /// This method transforms a vector from the global 3D cartesian frame
  /// into the local 3D cartesian frame
  ///
  /// @tparam vector3_t 3D vector
  ///
  /// @param v is the vector to be transformed
  ///
  /// @return a vector in local coordinates
  template <
      typename vector3_t,
      std::enable_if_t<std::is_convertible_v<vector3_t, vector3>, bool> = true>
  ALGEBRA_HOST_DEVICE inline vector3 vector_to_local(const vector3_t &v) const {

    return rotate(_data_inv, v);
  }
};  // struct transform3

}  // namespace algebra::vc_soa::math
