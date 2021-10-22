/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/algebra_qualifiers.hpp"
#include "algebra/math/impl/cmath_getter.hpp"
#include "algebra/math/impl/cmath_transform3.hpp"

// System include(s).
#include <cstddef>

namespace algebra::cmath {

/** Local frame projection into a polar coordinate frame
 */
template <template <typename, std::size_t> class array_t, typename transform3_t>
struct cylindrical2 {

  /// @name Type definitions for the struct
  /// @{

  /// Transformation matching this struct
  using transform3_type = transform3_t;

  /// Point in 2D space
  using point2 = typename transform3_type::point2;
  /// Point in 3D space
  using point3 = typename transform3_type::point3;

  /// @}

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

  /** This method transform from a point from 3D cartesian frame to a 2D
   * cylindrical point */
  ALGEBRA_HOST_DEVICE
  inline point2 operator()(const point3 &v) const {
    return {perp<array_t>(v) * phi<array_t>(v), v[2]};
  }

};  // struct cylindrical2

}  // namespace algebra::cmath
