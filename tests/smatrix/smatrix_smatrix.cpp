/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "algebra/smatrix_smatrix.hpp"

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
    return "smatrix_smatrix";
  }
};

// Instantiate the test(s).
typedef testing::Types<test_types<
    algebra::scalar, algebra::smatrix::point2, algebra::smatrix::point3,
    algebra::smatrix::vector2, algebra::smatrix::vector3,
    algebra::smatrix::transform3, algebra::smatrix::cartesian2,
    algebra::smatrix::polar2, algebra::smatrix::cylindrical2> >
    smatrix_smatrix_types;
INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, test_basics,
                               smatrix_smatrix_types, test_specialisation_name);
