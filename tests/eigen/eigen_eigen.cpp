/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "algebra/eigen_eigen.hpp"

// Test include(s).
#include "test_host_basics.hpp"

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <string>

/// Struct providing a readable name for the test
struct test_specialisation_name {
  template <typename T>
  static std::string GetName(int i) {
    switch (i) {
      case 0:
        return "eigen_eigen<float>";
      case 1:
        return "eigen_eigen<double>";
      default:
        return "unknown";
    }
  }
};

// Instantiate the test(s).
typedef testing::Types<
    test_types<
        float, algebra::eigen::point2<float>, algebra::eigen::point3<float>,
        algebra::eigen::vector2<float>, algebra::eigen::vector3<float>,
        algebra::eigen::transform3<float>, algebra::eigen::cartesian2<float>,
        algebra::eigen::polar2<float>, algebra::eigen::cylindrical2<float>,
        std::size_t, algebra::eigen::matrix_type>,
    test_types<
        double, algebra::eigen::point2<double>, algebra::eigen::point3<double>,
        algebra::eigen::vector2<double>, algebra::eigen::vector3<double>,
        algebra::eigen::transform3<double>, algebra::eigen::cartesian2<double>,
        algebra::eigen::polar2<double>, algebra::eigen::cylindrical2<double>,
        std::size_t, algebra::eigen::matrix_type> >
    eigen_eigen_types;
INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, test_host_basics,
                               eigen_eigen_types, test_specialisation_name);
