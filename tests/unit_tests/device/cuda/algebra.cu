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

REGISTER_TYPED_TEST_SUITE_P(cuda_vector_test, vector_2d_ops, vector_3d_ops);
REGISTER_TYPED_TEST_SUITE_P(cuda_matrix_test, matrix64_ops, matrix22_ops);
REGISTER_TYPED_TEST_SUITE_P(cuda_transform_test, transform3D);

// The 'test_types' are defined in 'algebra/test/framework/types.hpp'
INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, cuda_vector_test, test_types,
                               test_specialisation_name);
INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, cuda_matrix_test, test_types,
                               test_specialisation_name);
INSTANTIATE_TYPED_TEST_SUITE_P(algebra_plugins, cuda_transform_test, test_types,
                               test_specialisation_name);

}  // namespace algebra::test::cuda
