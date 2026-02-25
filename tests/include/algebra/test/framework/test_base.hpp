/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/concepts.hpp"
#include "algebra/type_traits.hpp"

// System include(s).
#include <limits>

namespace algebra::test {

/// Test base class, providing basic type definitions and test constants
template <algebra::concepts::algebra A>
class test_base {

  using scalar_t = algebra::get_scalar_t<A>;

 protected:
  /// Epsilon parameter for the floating point comparisons
  static constexpr scalar_t m_epsilon =
      std::numeric_limits<scalar_t>::epsilon();
  /// Variable defining when two floating point values are "close"
  static constexpr scalar_t m_isclose = static_cast<scalar_t>(1e-5);

};  // class test_base

}  // namespace algebra::test
