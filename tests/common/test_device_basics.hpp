/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/algebra_qualifiers.hpp"

// Local include(s).
#include "test_base.hpp"

/// Operations that can be executed on the host or a device
///
/// This class does not execute any tests itself. It performs "operations",
/// whose results could be tested in code that uses this type.
///
template <typename T>
class test_device_basics : public test_base<T> {

  /// @name Type definitions
  /// @{

  using scalar = typename test_base<T>::scalar;
  using point2 = typename test_base<T>::point2;
  using point3 = typename test_base<T>::point3;
  using vector2 = typename test_base<T>::vector2;
  using vector3 = typename test_base<T>::vector3;
  using transform3 = typename test_base<T>::transform3;
  using cartesian2 = typename test_base<T>::cartesian2;
  using polar2 = typename test_base<T>::polar2;
  using cylindrical2 = typename test_base<T>::cylindrical2;

  /// @}

 public:
  /// Perform various 2D vector operations, and produce a scalar output
  ALGEBRA_HOST_DEVICE
  scalar vector_2d_ops(point2 a, point2 b) const {

    point2 c = a + b;
    point2 c2 = c * 2.0;
    scalar phi = algebra::getter::phi(c2);
    scalar perp = algebra::getter::perp(c2);
    scalar norm1 = algebra::getter::norm(c2);

    scalar dot = algebra::vector::dot(a, b);
    point2 norm2 = algebra::vector::normalize(c);
    scalar norm3 = algebra::getter::norm(norm2);

    return (phi + perp + norm1 + dot + norm3);
  }

};  // class test_device_basics
