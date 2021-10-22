/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/scalar.hpp"
#include "algebra/math/vc.hpp"
#include "algebra/storage/vc.hpp"

namespace algebra {
namespace vc {

using transform3 = math::transform3<scalar>;
using cartesian2 = math::cartesian2<transform3>;
using polar2 = math::polar2<transform3>;
using cylindrical2 = math::cylindrical2<transform3>;

}  // namespace vc

namespace getter {

auto phi = [](const auto& a) { return cmath::phi<vc::storage_type>(a); };
auto theta = [](const auto& a) { return cmath::theta<vc::storage_type>(a); };
auto perp = [](const auto& a) { return cmath::perp<vc::storage_type>(a); };
auto norm = [](const auto& a) { return cmath::norm<vc::storage_type>(a); };
auto eta = [](const auto& a) { return cmath::eta<vc::storage_type>(a); };

}  // namespace getter

namespace vector {

auto cross = [](const auto& a, const auto& b) { return vc::math::cross(a, b); };
auto dot = [](const auto& a, const auto& b) { return vc::math::dot(a, b); };
auto normalize = [](const auto& a) { return vc::math::normalize(a); };

}  // namespace vector
}  // namespace algebra
