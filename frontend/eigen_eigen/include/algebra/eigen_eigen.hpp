/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "algebra/math/eigen.hpp"
#include "algebra/storage/eigen.hpp"

namespace algebra {
namespace eigen {

using transform3 = math::transform3<scalar>;
using cartesian2 = math::cartesian2<transform3>;
using polar2 = math::polar2<transform3>;
using cylindrical2 = math::cylindrical2<transform3>;

}  // namespace eigen

namespace getter {

auto phi = [](const auto& a) { return cmath::phi<eigen::storage_type>(a); };
auto theta = [](const auto& a) { return cmath::theta<eigen::storage_type>(a); };
auto perp = [](const auto& a) { return cmath::perp<eigen::storage_type>(a); };
auto norm = [](const auto& a) { return cmath::norm<eigen::storage_type>(a); };
auto eta = [](const auto& a) { return cmath::eta<eigen::storage_type>(a); };

}  // namespace getter

namespace vector {

auto cross = [](const auto& a, const auto& b) {
  return eigen::math::cross(a, b);
};
auto dot = [](const auto& a, const auto& b) { return eigen::math::dot(a, b); };
auto normalize = [](const auto& a) { return eigen::math::normalize(a); };

}  // namespace vector
}  // namespace algebra
