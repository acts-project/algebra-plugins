/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/cmath.hpp"
#include "algebra/math/vc.hpp"
#include "algebra/storage/vecmem.hpp"

namespace algebra {

using cmath::operator*;
using cmath::operator-;
using cmath::operator+;

namespace vecmem {

using transform3 = vc::math::transform3<vecmem::storage_type, scalar>;
using cartesian2 = vc::math::cartesian2<transform3>;
using polar2 = vc::math::polar2<vecmem::storage_type, transform3>;
using cylindrical2 = vc::math::cylindrical2<vecmem::storage_type, transform3>;

}  // namespace vecmem

namespace getter {

using cmath::eta;
using cmath::norm;
using cmath::perp;
using cmath::phi;
using cmath::theta;

template <auto SIZE, typename input_matrix_type>
ALGEBRA_HOST_DEVICE inline vecmem::storage_type<scalar, 4> vector(
    const input_matrix_type& m, std::size_t row, std::size_t col) {

  switch (row) {
    case 0:
      return m.x;
    case 1:
      return m.y;
    case 2:
      return m.z;
    case 3:
      return m.t;
    default:
      return vecmem::storage_type<scalar, 4>();
  }
}

}  // namespace getter

namespace vector {

using cmath::cross;
using cmath::dot;
using cmath::normalize;

}  // namespace vector
}  // namespace algebra
