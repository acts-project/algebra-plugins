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

using eigen::math::eta;
using eigen::math::norm;
using eigen::math::perp;
using eigen::math::phi;
using eigen::math::theta;

template <unsigned int SIZE, typename derived_type>
ALGEBRA_HOST_DEVICE inline auto vector(const Eigen::MatrixBase<derived_type>& m,
                                       std::size_t row, std::size_t col) {

  return m.template block<SIZE, 1>(row, col);
}

}  // namespace getter

namespace vector {

using eigen::math::cross;
using eigen::math::dot;
using eigen::math::normalize;

}  // namespace vector
}  // namespace algebra
