/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/coordinates/coordinate_base.hpp"
#include "algebra/qualifiers.hpp"

namespace algebra::common {

/** Frame projection into a cartesian coordinate frame
 */
template <typename transform3_t, typename E>
struct cylindrical3 final
    : public coordinate_base<cylindrical3, transform3_t, E> {

  /// @name Type definitions for the struct
  /// @{

  /// Base type
  using base_type = coordinate_base<cylindrical3, transform3_t, E>;
  /// Transformation matching this struct
  using transform3_type = typename base_type::transform3_type;
  /// Sclar type
  using scalar_type = typename base_type::scalar_type;
  /// Point in 2D space
  using point2 = typename base_type::point2;
  /// Point in 3D space
  using point3 = typename base_type::point3;
  /// Vector in 3D space
  using vector3 = typename base_type::vector3;
  /// Vector actor
  using vector_actor = typename base_type::vector_actor;

  /// @}

  /** This method transform from a point from 3D cartesian frame to a 3D
   * cartesian point */
  ALGEBRA_HOST_DEVICE
  inline point3 operator()(const point3 &p) const {

    return {vector_actor().perp(p), vector_actor().phi(p), p[2]};
  }

  /** This method transform from a point from global cartesian 3D frame to a
   * local 3D cylindrical point */
  ALGEBRA_HOST_DEVICE
  inline point3 global_to_local(const transform3_type &trf,
                                const point3 &p) const {
    const auto local3 = trf.point_to_local(p);
    return this->operator()(local3);
  }

  /** This method transform from a local 3D cylindrical point to a point global
   * cartesian 3D frame*/
  ALGEBRA_HOST_DEVICE
  inline point3 local_to_global(const transform3_type &trf,
                                const point3 &p) const {
    const scalar_type x = p[0] * std::cos(p[1]);
    const scalar_type y = p[0] * std::sin(p[1]);

    return trf.point_to_global(point3{x, y, p[2]});
  }

};  // struct cylindrical3

}  // namespace algebra::common