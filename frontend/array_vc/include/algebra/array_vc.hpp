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
#include "algebra/storage/array.hpp"

namespace algebra {

using cmath::operator*;
using cmath::operator-;
using cmath::operator+;

namespace array {

using transform3 = vc::math::transform3<array::storage_type, scalar>;
using cartesian2 = vc::math::cartesian2<transform3>;
using polar2 = vc::math::polar2<array::storage_type, transform3>;
using cylindrical2 = vc::math::cylindrical2<array::storage_type, transform3>;

}  // namespace array

namespace getter {

auto phi = [](const auto& a) { return cmath::phi<array::storage_type>(a); };
auto theta = [](const auto& a) { return cmath::theta<array::storage_type>(a); };
auto perp = [](const auto& a) { return cmath::perp<array::storage_type>(a); };
auto norm = [](const auto& a) { return cmath::norm<array::storage_type>(a); };
auto eta = [](const auto& a) { return cmath::eta<array::storage_type>(a); };

template <auto SIZE, typename input_matrix_type>
ALGEBRA_HOST_DEVICE inline array::storage_type<scalar, 4> vector(
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
      return array::storage_type<scalar, 4>();
  }
}

}  // namespace getter

namespace vector {

auto cross = [](const auto& a, const auto& b) {
  return cmath::cross<array::storage_type>(a, b);
};
auto dot = [](const auto& a, const auto& b) {
  return cmath::dot<array::storage_type>(a, b);
};
auto normalize = [](const auto& a) {
  return cmath::normalize<array::storage_type>(a);
};

}  // namespace vector
}  // namespace algebra
