/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "test_types.hpp"

// System include(s).
#include <limits>

/// Invalid implementation of the test base class
template <class T>
class test_base {};

/// Test base class, using a @c test_types type argument
template <typename scalar_t, typename point2_t, typename point3_t,
          typename vector2_t, typename vector3_t, typename transform3_t,
          typename cartesian2_t, typename polar2_t, typename cylindrical2_t>
class test_base<
    test_types<scalar_t, point2_t, point3_t, vector2_t, vector3_t, transform3_t,
               cartesian2_t, polar2_t, cylindrical2_t> > {

 public:
  /// @name Type definitions
  /// @{

  using scalar = scalar_t;
  using point2 = point2_t;
  using point3 = point3_t;
  using vector2 = vector2_t;
  using vector3 = vector3_t;
  using transform3 = transform3_t;
  using cartesian2 = cartesian2_t;
  using polar2 = polar2_t;
  using cylindrical2 = cylindrical2_t;

  /// @}

 protected:
  /// Epsilon parameter for the floating point comparisons
  static constexpr scalar m_epsilon = std::numeric_limits<scalar>::epsilon();
  /// Variable defining when two floating point values are "close"
  static constexpr scalar m_isclose = static_cast<scalar>(1e-5);

};  // class test_base
