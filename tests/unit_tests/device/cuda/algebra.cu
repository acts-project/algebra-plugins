/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Test include(s).
#include "algebra/test/device/cuda/algebra_test_suite.cuh"

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <string>

namespace algebra::test::cuda {

ALGEBRA_REGISTER_CUDA_VECTOR_TESTS();
ALGEBRA_REGISTER_CUDA_MATRIX_TESTS();
ALGEBRA_REGISTER_CUDA_TRANSFORM_TESTS();

// The 'test_types' are defined in 'algebra/test/framework/types.hpp'
INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, cuda_vector_test, test_types,
                               test_specialisation_name);
INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, cuda_matrix_test, test_types,
                               test_specialisation_name);
INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, cuda_transform_test, test_types,
                               test_specialisation_name);

}  // namespace algebra::test::cuda
