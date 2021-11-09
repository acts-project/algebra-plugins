/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "algebra/array_cmath.hpp"

// Test include(s).
#include "array_cmath_kernels.hpp"

// VecMem include(s).
#include <vecmem/containers/vector.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>

// GoogleTest include(s).
#include <gtest/gtest.h>

/// Memory resource for all of the tests.
static vecmem::cuda::managed_memory_resource resource;

/// Test for some basic 2D "vector operations"
TEST(cuda_array_cmath, 2d_vector_ops) {

  // Array size for the test.
  static constexpr std::size_t ARRAY_SIZE = 100;

  // Set up the input vector(s).
  vecmem::vector<algebra::array::point2> input_a(ARRAY_SIZE, &resource),
      input_b(ARRAY_SIZE, &resource);
  for (std::size_t i = 0; i < ARRAY_SIZE; ++i) {
    input_a[i] = {i * 0.5, (i + 1) * 1.0};
    input_b[i] = {(i - 1) * 1.2, i * 0.6};
  }

  // Set up the output vector(s).
  vecmem::vector<algebra::scalar> output_host(ARRAY_SIZE, &resource),
      output_device(ARRAY_SIZE, &resource);

  // Run the test on the host, and on the/a device.
  array_cmath_2d_vector_ops(
      vecmem::get_data(input_a), vecmem::get_data(input_b),
      vecmem::get_data(output_host), vecmem::get_data(output_device));

  // Compare the outputs.
  for (std::size_t i = 0; i < ARRAY_SIZE; ++i) {
    ASSERT_NEAR(output_host[i], output_device[i], 0.001);
  }
}
