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

auto phi = [](const auto& a) { return eigen::math::phi(a); };
auto theta = [](const auto& a) { return eigen::math::theta(a); };
auto perp = [](const auto& a) { return eigen::math::perp(a); };
auto norm = [](const auto& a) { return eigen::math::norm(a); };
auto eta = [](const auto& a) { return eigen::math::eta(a); };

template <unsigned int SIZE, typename derived_type>
ALGEBRA_HOST_DEVICE inline auto vector(const Eigen::MatrixBase<derived_type>& m,
                                       std::size_t row, std::size_t col) {

  return m.template block<SIZE, 1>(row, col);
}

}  // namespace getter

namespace vector {

auto cross = [](const auto& a, const auto& b) {
  return eigen::math::cross(a, b);
};
auto dot = [](const auto& a, const auto& b) { return eigen::math::dot(a, b); };
auto normalize = [](const auto& a) { return eigen::math::normalize(a); };

}  // namespace vector
}  // namespace algebra
