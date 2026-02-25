/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Test include(s).
#include "algebra/test/cpu/matrix_test.hpp"
#include "algebra/test/cpu/transform_test.hpp"
#include "algebra/test/cpu/vector_test.hpp"

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <string>

namespace algebra::test {

// Register the tests
ALGEBRA_REGISTER_VECTOR_TESTS();
#if !ALGEBRA_TEST_VC_AOS
ALGEBRA_REGISTER_MATRIX_TESTS();
#endif
ALGEBRA_REGISTER_TRANSFORM_TESTS();

// The 'test_types' are defined in 'algebra/test/framework/types.hpp'
INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, vector_test, test_types,
                               test_specialisation_name);
// TODO: Implement final parts of Vc AoS plugin
#if !ALGEBRA_TEST_VC_AOS
INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, matrix_test, test_types,
                               test_specialisation_name);
#endif
INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, transform_test, test_types,
                               test_specialisation_name);

}  // namespace algebra::test
