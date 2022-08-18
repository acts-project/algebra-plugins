/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/coordinates/coordinate_base.hpp"
#include "algebra/qualifiers.hpp"

// System include(s).
#include <cmath>

namespace algebra::common {

/** Local frame projection into a polar coordinate frame
 */
template <typename transform3_t, typename E>
struct cylindrical2 : public coordinate_base<cylindrical2, transform3_t, E> {

  /// @name Type definitions for the struct
  /// @{

  /// Base type
  using base_type = coordinate_base<cylindrical2, transform3_t, E>;
  /// Transformation matching this struct
  using transform3_type = typename base_type::transform3_type;
  /// Sclar type
  using scalar_type = typename transform3_type::scalar_type;
  /// Point in 2D space
  using point2 = typename transform3_type::point2;
  /// Point in 3D space
  using point3 = typename transform3_type::point3;
  /// Vector in 3D space
  using vector3 = typename transform3_type::vector3;
  /// Vector actor
  using vector_actor = typename transform3_type::vector_actor;

  /// @}

  /** This method transform from a point from 3D cartesian frame to a 2D
   * cylindrical point */
  ALGEBRA_HOST_DEVICE
  inline point2 operator()(const point3 &p) const {

    return {vector_actor().perp(p) * vector_actor().phi(p), p[2]};
  }

  /** This method transform from a point from global cartesian 3D frame to a
   * local 2D cylindrical point */
  ALGEBRA_HOST_DEVICE
  inline point2 global_to_local(const transform3_type &trf,
                                const point3 &p) const {
    const auto local3 = trf.point_to_local(p);
    return this->operator()(local3);
  }

  /** This method transform from a local 2D cylindrical point to a point global
   * cartesian 3D frame*/
  template <typename cylinder_mask_t>
  ALGEBRA_HOST_DEVICE inline point3 local_to_global(
      const transform3_type &trf, const point2 &p,
      const cylinder_mask_t &mask) const {
    const scalar_type r = mask.radius();
    const scalar_type phi = p[0] / r;
    const scalar_type x = r * std::cos(phi);
    const scalar_type y = r * std::sin(phi);
    const scalar_type z = p[1];

    return trf.point_to_global(point3{x, y, z});
  }

};  // struct cylindrical2

}  // namespace algebra::common