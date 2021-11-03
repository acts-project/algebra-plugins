/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "algebra/vc_cmath.hpp"

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
    return "vc_cmath";
  }
};

// This test needs the algebra operators.
using algebra::operator+;
using algebra::operator*;

// Instantiate the test(s).
typedef testing::Types<test_types<
    algebra::scalar, algebra::vc::point2, algebra::vc::point3,
    algebra::vc::vector2, algebra::vc::vector3, algebra::vc::transform3,
    algebra::vc::cartesian2, algebra::vc::polar2, algebra::vc::cylindrical2> >
    vc_cmath_types;
INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, test_basics, vc_cmath_types,
                               test_specialisation_name);
