/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "algebra/vecmem_cmath.hpp"

// Test include(s).
#include "test_basics.hpp"

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <string>

/// Struct providing a readable name for the test
struct test_specialisation_name {
  template <typename T>
  static std::string GetName(int) {
    return "vecmem_cmath";
  }
};

// Instantiate the test(s).
typedef testing::Types<test_types<
    algebra::scalar, algebra::vecmem::point2, algebra::vecmem::point3,
    algebra::vecmem::vector2, algebra::vecmem::vector3,
    algebra::vecmem::transform3, algebra::vecmem::cartesian2,
    algebra::vecmem::polar2, algebra::vecmem::cylindrical2> >
    vecmem_cmath_types;
INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, test_basics, vecmem_cmath_types,
                               test_specialisation_name);
