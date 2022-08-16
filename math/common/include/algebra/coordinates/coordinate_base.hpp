/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/qualifiers.hpp"

namespace algebra {

/** Coordinate base struct
 */
template <typename transform3_t>
struct coordinate_base {

  /// @name Type definitions for the struct
  /// @{

  /// Transformation matching this struct
  using transform3_type = transform3_t;
  /// Scalar type
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
};

}  // namespace algebra