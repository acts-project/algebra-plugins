/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "algebra/concepts.hpp"
#include "test_types.hpp"

// System include(s).
#include <limits>

/// Invalid implementation of the test base class
template <class T>
class test_base {};

/// Test base class, using a @c test_types type argument
template <
    algebra::concepts::scalar scalar_t, algebra::concepts::point2D point2_t,
    algebra::concepts::point3D point3_t, algebra::concepts::vector2D vector2_t,
    algebra::concepts::vector3D vector3_t,
    algebra::concepts::transform3D transform3_t,
    algebra::concepts::index size_ty,
    template <typename, size_ty, size_ty> class matrix_t>
class test_base<test_types<scalar_t, point2_t, point3_t, vector2_t, vector3_t,
                           transform3_t, size_ty, matrix_t> > {

 public:
  /// @name Type definitions
  /// @{

  using scalar = scalar_t;
  using point2 = point2_t;
  using point3 = point3_t;
  using vector2 = vector2_t;
  using vector3 = vector3_t;
  using transform3 = transform3_t;
  using size_type = size_ty;
  template <size_type ROWS, size_type COLS>
  using matrix = matrix_t<scalar, ROWS, COLS>;

  /// @}

 protected:
  /// Epsilon parameter for the floating point comparisons
  static constexpr scalar m_epsilon = std::numeric_limits<scalar>::epsilon();
  /// Variable defining when two floating point values are "close"
  static constexpr scalar m_isclose = static_cast<scalar>(1e-5);

};  // class test_base
