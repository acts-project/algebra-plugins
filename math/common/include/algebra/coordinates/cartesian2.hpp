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

namespace algebra {

/** Frame projection into a cartesian coordinate frame
 */
template <typename transform3_t>
struct cartesian2 final : public coordinate_base<transform3_t> {

  /// @name Type definitions for the struct
  /// @{

  /// Base type
  using base_type = coordinate_base<transform3_t>;
  /// Sclar type
  using scalar_type = typename base_type::scalar_type;
  /// Transformation matching this struct
  using transform3_type = typename base_type::transform3_type;
  /// Point in 2D space
  using point2 = typename base_type::point2;
  /// Point in 3D space
  using point3 = typename base_type::point3;
  /// Vector in 3D space
  using vector3 = typename base_type::vector3;
  /// Vector actor
  using vector_actor = typename base_type::vector_actor;

  /// @}

  /** This method transform from a point from 2D cartesian frame to a 2D
   * cartesian point */
  ALGEBRA_HOST_DEVICE
  inline point2 operator()(const point2 &local2) const {

    return {local2[0], local2[1]};
  }

  /** This method transform from a point from 3D cartesian frame to a 2D
   * cartesian point */
  ALGEBRA_HOST_DEVICE
  inline point2 operator()(const point3 &local3) const {

    return {local3[0], local3[1]};
  }

  ALGEBRA_HOST_DEVICE
  inline point2 global_to_local(const transform3_type &trf, const point3 &p) {
    const auto local3 = trf.point_to_local(p);
    return this->operator()(local3);
  }

  ALGEBRA_HOST_DEVICE
  inline point2 global_to_local(const transform3_type &trf, const point3 &p,
                                const vector3 & /*d*/) {
    return global_to_local(trf, p);
  }

};  // struct cartesian2

}  // namespace algebra