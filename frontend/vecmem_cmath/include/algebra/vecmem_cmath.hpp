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
#include "algebra/storage/vecmem.hpp"

namespace algebra {

using cmath::operator*;
using cmath::operator-;
using cmath::operator+;

namespace vecmem {

using transform3 = cmath::transform3<vecmem::storage_type, scalar>;
using cartesian2 = cmath::cartesian2<transform3>;
using polar2 = cmath::polar2<transform3>;
using cylindrical2 = cmath::cylindrical2<transform3>;

}  // namespace vecmem

namespace getter {

using cmath::eta;
using cmath::norm;
using cmath::perp;
using cmath::phi;
using cmath::theta;

template <auto SIZE, auto ROWS, auto COLS>
ALGEBRA_HOST_DEVICE inline vecmem::storage_type<scalar, SIZE> vector(
    const vecmem::storage_type<vecmem::storage_type<scalar, ROWS>, COLS>& m,
    std::size_t row, std::size_t col) {

  return cmath::vector_getter<vecmem::storage_type, scalar>()
      .template operator()<SIZE>(m, row, col);
}

}  // namespace getter

namespace vector {

using cmath::cross;
using cmath::dot;
using cmath::normalize;

}  // namespace vector
}  // namespace algebra
