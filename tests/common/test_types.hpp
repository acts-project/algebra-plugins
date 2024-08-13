/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "algebra/concepts.hpp"

/// Simple struct holding the types that describe a given plugin
template <
    algebra::concepts::scalar scalar_t, algebra::concepts::point2D point2_t,
    algebra::concepts::point3D point3_t, algebra::concepts::vector2D vector2_t,
    algebra::concepts::vector3D vector3_t,
    algebra::concepts::transform3D transform3_t,
    algebra::concepts::index size_ty,
    template <typename, size_ty, size_ty> class matrix_t>
struct test_types {

  using scalar = scalar_t;
  using point2 = point2_t;
  using point3 = point3_t;
  using vector2 = vector2_t;
  using vector3 = vector3_t;
  using transform3 = transform3_t;
  using size_type = size_ty;
  template <size_type ROWS, size_type COLS>
  using matrix = matrix_t<scalar, ROWS, COLS>;

};  // struct test_types
