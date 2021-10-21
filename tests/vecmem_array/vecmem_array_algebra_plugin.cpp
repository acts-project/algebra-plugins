/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "algebra/storage/vecmem.hpp"
#include "algebra/math/cmath.hpp"

namespace algebra {
namespace vecmem {

using transform3 = cmath::transform3<storage_type, scalar>;
using cartesian2 = cmath::cartesian2<storage_type, scalar>;
using polar2 = cmath::polar2<storage_type, scalar>;
using cylindrical2 = cmath::cylindrical2<storage_type, scalar>;

}  // namespace vecmem
}  // namespace algebra

// Test include(s).
#define __plugin algebra::vecmem
#define __plugin_without_matrix_element_accessor 1
#include "test_plugin.inl"
