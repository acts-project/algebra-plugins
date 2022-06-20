/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/impl/eigen_getter.hpp"
#include "algebra/qualifiers.hpp"

namespace algebra::eigen::math {

/** Local frame projection into a polar coordinate frame */
template <typename transform3_t>
struct cylindrical2 {

  /// @name Type definitions for the struct
  /// @{

  /// Transformation matching this struct
  using transform3_type = transform3_t;
  /// Scalar type used by the transform
  using scalar_type = typename transform3_type::scalar_type;

  /// Point in 2D space
  using point2 = typename transform3_type::point2;
  /// Point in 3D space
  using point3 = typename transform3_type::point3;

  /// @}

  /** This method transform from a point from 2D or 3D cartesian frame to a 2D
   * polar point
   *
   * @param v the point in local frame
   *
   * @return a local point2
   */
  ALGEBRA_HOST_DEVICE inline point2 operator()(const point3 &v) const {

    return {perp(v) * phi(v), v[2]};
  }

  /** This method transform from a point from the global 3D cartesian frame to
   *the local 2D cartesian frame
   *
   * @param trf the transform from global to local thredimensional frame
   * @param p the point in global frame
   *
   * @return a local point2
   **/
  ALGEBRA_HOST_DEVICE
  inline point2 operator()(const transform3_type &trf, const point3 &p) const {

    return operator()(trf.point_to_local(p));
  }
};  // struct cylindrical2

}  // namespace algebra::eigen::math
