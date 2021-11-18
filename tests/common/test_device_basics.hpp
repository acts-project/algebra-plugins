/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/qualifiers.hpp"

// Local include(s).
#include "test_base.hpp"

/// Operations that can be executed on the host or a device
///
/// This class does not execute any tests itself. It performs "operations",
/// whose results could be tested in code that uses this type.
///
template <typename T>
class test_device_basics : public test_base<T> {

 public:
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

  /// Perform various 3D vector operations, and produce a scalar output
  ALGEBRA_HOST_DEVICE
  scalar vector_3d_ops(vector3 a, vector3 b) const {

    vector3 c = a + b;
    vector3 c2 = c * 2.0;

    scalar phi = algebra::getter::phi(c2);
    scalar perp = algebra::getter::perp(c2);
    scalar norm1 = algebra::getter::norm(c2);

    vector3 d = algebra::vector::cross(a, b);

    scalar dot = algebra::vector::dot(a, d);
    vector3 norm2 = algebra::vector::normalize(c);
    scalar norm3 = algebra::getter::norm(norm2);

    return (phi + perp + norm1 + dot + norm3);
  }

  /// Perform various operations using the @c transform3 type
  ALGEBRA_HOST_DEVICE
  scalar transform3_ops(vector3 t1, vector3 t2, vector3 t3, vector3 a,
                        vector3 b) const {

    transform3 tr(t1, t2, t3);

    point3 translation = tr.translation();

    point3 gpoint = tr.point_to_global(a);
    point3 lpoint = tr.point_to_local(b);

    vector3 gvec = tr.vector_to_global(a);
    vector3 lvec = tr.vector_to_local(b);

    return {algebra::getter::norm(translation) + algebra::getter::perp(gpoint) +
            algebra::getter::phi(lpoint) + algebra::vector::dot(gvec, lvec)};
  }

  /// Perform various operations using the @c cartesian2 type
  ALGEBRA_HOST_DEVICE
  scalar cartesian2_ops(vector3 t1, vector3 t2, vector3 t3, vector3 a,
                        vector3 b) const {

    transform3 tr(t1, t2, t3);
    cartesian2 ca;

    point2 p1 = ca(tr, a);
    point2 p2 = ca(b);

    return {algebra::getter::phi(p1) + algebra::getter::norm(p2)};
  }

  /// Perform various operations using the @c cylintridcal2 type
  ALGEBRA_HOST_DEVICE
  scalar cylindrical2_ops(vector3 t1, vector3 t2, vector3 t3, vector3 a,
                          vector3 b) const {

    transform3 tr(t1, t2, t3);
    cylindrical2 cy;

    point2 p1 = cy(tr, a);
    point2 p2 = cy(b);

    return {algebra::getter::phi(p1) + algebra::getter::norm(p2)};
  }

  /// Perform various operations using the @c polar2 type
  ALGEBRA_HOST_DEVICE
  scalar polar2_ops(vector3 t1, vector3 t2, vector3 t3, vector3 a,
                    vector3 b) const {

    transform3 tr(t1, t2, t3);
    polar2 po;

    point2 p1 = po(tr, a);
    point2 p2 = po(b);
    point2 p3 = po(p1);

    return {algebra::getter::phi(p2) + algebra::getter::norm(p3)};
  }

};  // class test_device_basics
