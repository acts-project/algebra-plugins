/** Algebra plugins, part of the ACTS project
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

// Project include(s).
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
#include <array>
#include <cstddef>
#include <limits>

namespace algebra::storage {

/// 4x4 matrix type used by @c algebra::vc_soa::math::transform3 that has simd
/// vectors as matrix elements
template <template <typename, std::size_t> class array_t, typename value_t,
          std::size_t N>
struct matrix44 {

  using vector_type = storage::vector<N, value_t, array_t>;
  using value_type = value_t;

  /// Default constructor: Identity with no translation
  matrix44()
      : x{1.f, 0.f, 0.f},
        y{0.f, 1.f, 0.f},
        z{0.f, 0.f, 1.f},
        t{0.f, 0.f, 0.f} {}

  /// Identity rotation with translation @param translation
  matrix44(vector_type v)
      : x{1.f, 0.f, 0.f}, y{0.f, 1.f, 0.f}, z{0.f, 0.f, 1.f}, t{std::move(v)} {}

  /// Construct from given column vectors @param x, @param y, @param z, @param t
  matrix44(vector_type v_0, vector_type v_1, vector_type v_2, vector_type v_3)
      : x{std::move(v_0)},
        y{std::move(v_1)},
        z{std::move(v_2)},
        t{std::move(v_3)} {}

  /// Identity rotation with translation from single elemenst @param t_0,
  /// @param t_1, @param t_2
  matrix44(const value_type t_0, const value_type t_1, const value_type t_2)
      : matrix44{{t_0, t_1, t_2}} {}

  /// Construct from elements (simd vector) of matrix column vectors:
  /// - column 0: @param x_0, @param x_1, @param x_2
  /// - column 1: @param y_0, @param y_1, @param y_2
  /// - column 2: @param z_0, @param z_1, @param z_2
  /// - column 3: @param t_0, @param t_1, @param t_2
  matrix44(const value_type x_0, const value_type x_1, const value_type x_2,
           const value_type y_0, const value_type y_1, const value_type y_2,
           const value_type z_0, const value_type z_1, const value_type z_2,
           const value_type t_0, const value_type t_1, const value_type t_2)
      : x{{x_0, x_1, x_2}},
        y{{y_0, y_1, y_2}},
        z{{z_0, z_1, z_2}},
        t{{t_0, t_1, t_2}} {}

  /// Equality operator between two matrices
  /// @{
  /// SoA
  template <typename V = value_t,
            typename std::enable_if_t<!std::is_scalar_v<V>, bool> = true>
  bool operator==(const matrix44 &rhs) const {
    return ((x == rhs.x) && (y == rhs.y) && (z == rhs.z) && (t == rhs.t));
  }

  /// AoS
  template <typename V = value_t,
            typename std::enable_if_t<std::is_scalar_v<V>, bool> = true>
  bool operator==(const matrix44 &rhs) const {
    return ((x.get() == rhs.x.get()).isFull() &&
            (y.get() == rhs.y.get()).isFull() &&
            (z.get() == rhs.z.get()).isFull() &&
            (t.get() == rhs.t.get()).isFull());
  }
  /// @}

  /// Data variables
  vector_type x, y, z, t;

};  // struct matrix44

/// Functor used to access elements of matrix44
struct element_getter {

  /// Get const access to a matrix element
  template <template <typename, std::size_t> class array_t, typename value_t,
            std::size_t N>
  ALGEBRA_HOST inline decltype(auto) operator()(
      const matrix44<array_t, value_t, N> &m, std::size_t row,
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
#ifndef _MSC_VER
        __builtin_unreachable();
#else
        return m.x[0];
#endif
    }
  }

  /// Get non-const access to a matrix element
  template <template <typename, std::size_t> class array_t, typename value_t,
            std::size_t N>
  ALGEBRA_HOST inline decltype(auto) operator()(
      matrix44<array_t, value_t, N> &m, std::size_t row,
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
#ifndef _MSC_VER
        __builtin_unreachable();
#else
        return m.x[0];
#endif
    }
  }

  /// Get const access to a matrix element
  template <template <typename, std::size_t> class array_t, typename value_t>
  ALGEBRA_HOST inline decltype(auto) operator()(const array_t<value_t, 16> &m,
                                                unsigned int row,
                                                unsigned int col) const {
    // Make sure that the indices are valid.
    assert(row < 4);
    assert(col < 4);
    // Return the selected element.
    return m[row * 4 + col];
  }

};  // struct element_getter

}  // namespace algebra::storage
