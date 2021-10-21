/** Algebra plugins library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "algebra/storage/vc.hpp"
#include "algebra/math/cmath.hpp"

namespace algebra {
namespace vc {

using transform3 = cmath::transform3<storage_type, scalar, Vc::array<storage_type<scalar, 4>, 4 > >;
using cartesian2 = cmath::cartesian2<storage_type, scalar, transform3>;
using polar2 = cmath::polar2<storage_type, scalar, transform3>;
using cylindrical2 = cmath::cylindrical2<storage_type, scalar, transform3>;

}  // namespace vc
}  // namespace algebra

// Test include(s).
#define __plugin algebra::vc
#define __plugin_without_matrix_element_accessor 1
#include "test_plugin.inl"
