/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/scalar.hpp"
#include "algebra/math/cmath.hpp"
#include "algebra/math/vc.hpp"
#include "algebra/storage/vc.hpp"

namespace algebra {

using cmath::operator*;
using cmath::operator-;
using cmath::operator+;

namespace vc {

using transform3 = math::transform3<storage_type, scalar, vector3, point2>;
using cartesian2 = math::cartesian2<transform3>;
using polar2 = math::polar2<storage_type, transform3>;
using cylindrical2 = math::cylindrical2<storage_type, transform3>;

}  // namespace vc

namespace getter {

using cmath::eta;
using cmath::norm;
using cmath::perp;
using cmath::phi;
using cmath::theta;

template <auto SIZE, typename input_matrix_type>
ALGEBRA_HOST_DEVICE inline vc::storage_type<scalar, SIZE> vector(
    const input_matrix_type& m, std::size_t row, std::size_t col) {

  switch (col) {
    case 0:
      return m.x;
    case 1:
      return m.y;
    case 2:
      return m.z;
    case 3:
      return m.t;
    default:
      return {};
  }
}

}  // namespace getter

namespace vector {

using cmath::dot;
using cmath::normalize;
using vc::math::cross;
using vc::math::dot;
using vc::math::normalize;

}  // namespace vector
}  // namespace algebra
