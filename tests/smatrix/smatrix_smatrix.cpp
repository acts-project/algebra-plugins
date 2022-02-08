/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "algebra/smatrix_smatrix.hpp"

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
        return "smatrix_smatrix<float>";
      case 1:
        return "smatrix_smatrix<double>";
      default:
        return "unknown";
    }
  }
};

// Instantiate the test(s).
typedef testing::Types<
    test_types<
        float, algebra::smatrix::point2<float>, algebra::smatrix::point3<float>,
        algebra::smatrix::vector2<float>, algebra::smatrix::vector3<float>,
        algebra::smatrix::transform3<float>,
        algebra::smatrix::cartesian2<float>, algebra::smatrix::polar2<float>,
        algebra::smatrix::cylindrical2<float>,
        algebra::smatrix::matrix<float, 4, 4> >,
    test_types<
        double, algebra::smatrix::point2<double>,
        algebra::smatrix::point3<double>, algebra::smatrix::vector2<double>,
        algebra::smatrix::vector3<double>, algebra::smatrix::transform3<double>,
        algebra::smatrix::cartesian2<double>, algebra::smatrix::polar2<double>,
        algebra::smatrix::cylindrical2<double>,
        algebra::smatrix::matrix<double, 4, 4> > >
    smatrix_smatrix_types;
INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, test_host_basics,
                               smatrix_smatrix_types, test_specialisation_name);
