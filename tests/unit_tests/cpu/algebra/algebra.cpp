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
REGISTER_TYPED_TEST_SUITE_P(vector_test, vector2D, vector3D, element_getter,
                            dot_product_with_ops, cross_product_add_sub);
#if !ALGEBRA_TEST_VC_AOS
REGISTER_TYPED_TEST_SUITE_P(matrix_test, matrix_2x2, matrix_2x3, matrix_3x1,
                            matrix_3x3, matrix_6x4, matrix_5x5, matrix_6x6,
                            matrix_7x4_4x12, matrix_17x9_9x4, matrix_5x2_2x3,
                            matrix_small_mixed);
#endif
REGISTER_TYPED_TEST_SUITE_P(transform_test, transform3D,
                            global_transformations);

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
