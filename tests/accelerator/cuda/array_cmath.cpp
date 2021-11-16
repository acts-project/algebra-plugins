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

/// CUDA test case for algebra::array_cmath
class cuda_array_cmath : public testing::Test {

 public:
  cuda_array_cmath() {

    // Initialise the vectors with some dummy values.
    for (std::size_t i = 0; i < s_arraySize; ++i) {

      m_t1[i] = {1., 2., 3.};
      m_t2[i] = {4., 5., 6.};
      m_t3[i] = {7., 8., 9.};

      m_p1[i] = {i * 0.5, (i + 1) * 1.0};
      m_p2[i] = {(i - 1) * 1.2, i * 0.6};

      m_v1[i] = {i * 0.6, (i + 1) * 1.2, (i + 2) * 1.3};
      m_v2[i] = {(i - 1) * 1.8, i * 2.3, (i - 2) * 3.4};
    }
  }

 protected:
  /// Memory resource for all of the tests.
  vecmem::cuda::managed_memory_resource m_resource;

  /// Size for the tested arrays.
  static constexpr std::size_t s_arraySize = 5000;

  /// @name Inputs for the tests
  /// @{

  vecmem::vector<algebra::array::vector3> m_t1{s_arraySize, &m_resource};
  vecmem::vector<algebra::array::vector3> m_t2{s_arraySize, &m_resource};
  vecmem::vector<algebra::array::vector3> m_t3{s_arraySize, &m_resource};

  vecmem::vector<algebra::array::point2> m_p1{s_arraySize, &m_resource};
  vecmem::vector<algebra::array::point2> m_p2{s_arraySize, &m_resource};

  vecmem::vector<algebra::array::vector3> m_v1{s_arraySize, &m_resource};
  vecmem::vector<algebra::array::vector3> m_v2{s_arraySize, &m_resource};

  /// @}

};  // class cuda_array_cmath

/// Test for some basic 2D "vector operations"
TEST_F(cuda_array_cmath, 2d_vector_ops) {

  // Set up the output vector(s).
  vecmem::vector<algebra::scalar> output_host(s_arraySize, &m_resource),
      output_device(s_arraySize, &m_resource);

  // Run the test on the host, and on the/a device.
  array_cmath_2d_vector_ops(vecmem::get_data(m_p1), vecmem::get_data(m_p2),
                            vecmem::get_data(output_host),
                            vecmem::get_data(output_device));

  // Compare the outputs.
  for (std::size_t i = 0; i < s_arraySize; ++i) {
    EXPECT_FLOAT_EQ(output_host[i], output_device[i]);
  }
}

/// Test for some basic 3D "vector operations"
TEST_F(cuda_array_cmath, 3d_vector_ops) {

  // Set up the output vector(s).
  vecmem::vector<algebra::scalar> output_host(s_arraySize, &m_resource),
      output_device(s_arraySize, &m_resource);

  // Run the test on the host, and on the/a device.
  array_cmath_3d_vector_ops(vecmem::get_data(m_v1), vecmem::get_data(m_v2),
                            vecmem::get_data(output_host),
                            vecmem::get_data(output_device));

  // Compare the outputs.
  for (std::size_t i = 0; i < s_arraySize; ++i) {
    EXPECT_FLOAT_EQ(output_host[i], output_device[i]);
  }
}

/// Test for some operations with @c algebra::cmath::transform3
TEST_F(cuda_array_cmath, transform3) {

  // Set up the output vector(s).
  vecmem::vector<algebra::scalar> output_host(s_arraySize, &m_resource),
      output_device(s_arraySize, &m_resource);

  // Run the test on the host, and on the/a device.
  array_cmath_transform3_ops(
      vecmem::get_data(m_t1), vecmem::get_data(m_t2), vecmem::get_data(m_t3),
      vecmem::get_data(m_v1), vecmem::get_data(m_v2),
      vecmem::get_data(output_host), vecmem::get_data(output_device));

  // Compare the outputs.
  for (std::size_t i = 0; i < s_arraySize; ++i) {
    EXPECT_FLOAT_EQ(output_host[i], output_device[i]);
  }
}

/// Test for some operations with @c algebra::cmath::cartesian2
TEST_F(cuda_array_cmath, cartesian2) {

  // Set up the output vector(s).
  vecmem::vector<algebra::scalar> output_host(s_arraySize, &m_resource),
      output_device(s_arraySize, &m_resource);

  // Run the test on the host, and on the/a device.
  array_cmath_cartesian2_ops(
      vecmem::get_data(m_t1), vecmem::get_data(m_t2), vecmem::get_data(m_t3),
      vecmem::get_data(m_v1), vecmem::get_data(m_v2),
      vecmem::get_data(output_host), vecmem::get_data(output_device));

  // Compare the outputs.
  for (std::size_t i = 0; i < s_arraySize; ++i) {
    EXPECT_FLOAT_EQ(output_host[i], output_device[i]);
  }
}

/// Test for some operations with @c algebra::cmath::cylindrical2
TEST_F(cuda_array_cmath, cylindrical2) {

  // Set up the output vector(s).
  vecmem::vector<algebra::scalar> output_host(s_arraySize, &m_resource),
      output_device(s_arraySize, &m_resource);

  // Run the test on the host, and on the/a device.
  array_cmath_cylindrical2_ops(
      vecmem::get_data(m_t1), vecmem::get_data(m_t2), vecmem::get_data(m_t3),
      vecmem::get_data(m_v1), vecmem::get_data(m_v2),
      vecmem::get_data(output_host), vecmem::get_data(output_device));

  // Compare the outputs.
  for (std::size_t i = 0; i < s_arraySize; ++i) {
    EXPECT_FLOAT_EQ(output_host[i], output_device[i]);
  }
}

/// Test for some operations with @c algebra::cmath::polar2
TEST_F(cuda_array_cmath, polar2) {

  // Set up the output vector(s).
  vecmem::vector<algebra::scalar> output_host(s_arraySize, &m_resource),
      output_device(s_arraySize, &m_resource);

  // Run the test on the host, and on the/a device.
  array_cmath_polar2_ops(vecmem::get_data(m_t1), vecmem::get_data(m_t2),
                         vecmem::get_data(m_t3), vecmem::get_data(m_v1),
                         vecmem::get_data(m_v2), vecmem::get_data(output_host),
                         vecmem::get_data(output_device));

  // Compare the outputs.
  for (std::size_t i = 0; i < s_arraySize; ++i) {
    EXPECT_FLOAT_EQ(output_host[i], output_device[i]);
  }
}
