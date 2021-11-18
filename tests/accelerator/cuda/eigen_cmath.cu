/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "algebra/eigen_cmath.hpp"

// Local include(s).
#include "test_cuda_basics.cuh"

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <string>

/// Struct providing a readable name for the test
struct test_specialisation_name {
  template <typename T>
  static std::string GetName(int) {
    return "cuda_eigen_cmath";
  }
};

// Instantiate the test(s).
typedef testing::Types<
    test_types<algebra::scalar, algebra::eigen::point2, algebra::eigen::point3,
               algebra::eigen::vector2, algebra::eigen::vector3,
               algebra::eigen::transform3, algebra::eigen::cartesian2,
               algebra::eigen::polar2, algebra::eigen::cylindrical2> >
    eigen_cmath_types;
INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, test_cuda_basics,
                               eigen_cmath_types, test_specialisation_name);
