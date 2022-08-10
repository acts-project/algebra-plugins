/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/qualifiers.hpp"

namespace algebra::common {

/** Frame projection into a cartesian coordinate frame
 */
template <typename transform3_t>
struct cartesian2 {

  /// @name Type definitions for the struct
  /// @{

  /// Transformation matching this struct
  using transform3_type = transform3_t;

  /// Point in 2D space
  using point2 = typename transform3_type::point2;
  /// Point in 3D space
  using point3 = typename transform3_type::point3;
  /// Vector in 3D space
  using vector3 = typename transform3_type::vector3;

  /// @}

  /** This method transform from a point from the global 3D cartesian frame to
   *the local 2D cartesian frame
   *
   * @param trf the transform from global to local thredimensional frame
   * @param p the point in global frame
   * @param v unused vector
   *
   * @return a local point2
   **/
  ALGEBRA_HOST_DEVICE
  inline point2 operator()(const transform3_type &trf, const point3 &p,
                           const vector3 & /*v*/) const {

    return operator()(trf, p);
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

  /** This method transform from a point from the global 3D cartesian frame to
   * the local 2D cartesian frame
   *
   * @param v the point in local frame
   *
   * @return a local point2
   */
  ALGEBRA_HOST_DEVICE
  inline point2 operator()(const point3 &v) const { return {v[0], v[1]}; }

};  // struct cartesian2

}  // namespace algebra::common