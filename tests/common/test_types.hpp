/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

/// Simple struct holding the types that describe a given plugin
template <typename scalar_t, typename point2_t, typename point3_t,
          typename vector2_t, typename vector3_t, typename transform3_t,
          typename cartesian2_t, typename polar2_t, typename cylindrical2_t,
          typename matrix_t>
struct test_types {

  using scalar = scalar_t;
  using point2 = point2_t;
  using point3 = point3_t;
  using vector2 = vector2_t;
  using vector3 = vector3_t;
  using transform3 = transform3_t;
  using cartesian2 = cartesian2_t;
  using polar2 = polar2_t;
  using cylindrical2 = cylindrical2_t;
  using matrix = matrix_t;

};  // struct test_types
