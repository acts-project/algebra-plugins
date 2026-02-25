/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Algebra plugin selection
#if ALGEBRA_TEST_ARRAY
#include "algebra/array.hpp"
#endif
#if ALGEBRA_TEST_EIGEN
#include "algebra/eigen.hpp"
#endif
#if ALGEBRA_TEST_FASTOR
#include "algebra/fastor.hpp"
#endif
#if ALGEBRA_TEST_SMATRIX
#include "algebra/smatrix.hpp"
#endif
#if ALGEBRA_TEST_VC_AOS
#include "algebra/vc_aos.hpp"
#endif
#if ALGEBRA_TEST_VC_SOA
#include "algebra/vc_soa.hpp"
#endif

// Project include(s).
#include "algebra/concepts.hpp"
#include "algebra/type_traits.hpp"

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <limits>
#include <string>

namespace algebra::test {

#define ALGEBRA_DEFINE_TEST_TYPES(ALGEBRA)            \
  template <algebra::concepts::value T>               \
  using algebra_type = ::algebra::plugin::ALGEBRA<T>; \
                                                      \
  struct test_specialisation_name {                   \
    template <typename T>                             \
    static std::string GetName(int i) {               \
      switch (i) {                                    \
        case 0:                                       \
          return std::string(#ALGEBRA) + "<float>";   \
        case 1:                                       \
          return std::string(#ALGEBRA) + "<double>";  \
        default:                                      \
          return "unknown";                           \
      }                                               \
    }                                                 \
  };

// Select algebra-plugin to compile the test with
#if ALGEBRA_TEST_ARRAY
ALGEBRA_DEFINE_TEST_TYPES(array)
#elif ALGEBRA_TEST_EIGEN
ALGEBRA_DEFINE_TEST_TYPES(eigen)
#elif ALGEBRA_TEST_FASTOR
ALGEBRA_DEFINE_TEST_TYPES(fastor)
#elif ALGEBRA_TEST_SMATRIX
ALGEBRA_DEFINE_TEST_TYPES(smatrix)
#elif ALGEBRA_TEST_VC_AOS
ALGEBRA_DEFINE_TEST_TYPES(vc_aos)
#elif ALGEBRA_TEST_VC_SOA
ALGEBRA_DEFINE_TEST_TYPES(vc_soa)
#endif

// Instantiate the test(s).
using test_types = testing::Types<algebra_type<float>, algebra_type<double>>;

}  // namespace algebra::test
