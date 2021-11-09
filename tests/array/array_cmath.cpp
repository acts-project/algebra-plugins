/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "algebra/array_cmath.hpp"

// Test include(s).
#include "test_host_basics.hpp"

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <string>

/// Struct providing a readable name for the test
struct test_specialisation_name {
  template <typename T>
  static std::string GetName(int) {
    return "array_cmath";
  }
};

// Instantiate the test(s).
typedef testing::Types<
    test_types<algebra::scalar, algebra::array::point2, algebra::array::point3,
               algebra::array::vector2, algebra::array::vector3,
               algebra::array::transform3, algebra::array::cartesian2,
               algebra::array::polar2, algebra::array::cylindrical2> >
    array_cmath_types;
INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, test_host_basics,
                               array_cmath_types, test_specialisation_name);
