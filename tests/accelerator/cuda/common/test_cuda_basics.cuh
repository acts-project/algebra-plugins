/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/algebra_qualifiers.hpp"

// Local include(s).
#include "execute_cuda_test.cuh"
#include "execute_host_test.hpp"
#include "test_base.hpp"
#include "test_cuda_basics_functors.hpp"

// VecMem include(s).
#include <vecmem/containers/vector.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <cstddef>

/// Test case class, to be specialised for the different plugins
template <typename T>
class test_cuda_basics : public testing::Test, public test_base<T> {

 public:
  /// Constructor, setting up the inputs for all of the tests
  test_cuda_basics() {

    // Initialise the input vectors with some dummy values.
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
  /// Function setting things up for a test
  virtual void SetUp() override {

    // Reset the values in the output vectors.
    for (std::size_t i = 0; i < s_arraySize; ++i) {
      m_output_host[i] = 0;
      m_output_device[i] = 0;
    }
  }

  /// Memory resource for all of the tests.
  vecmem::cuda::managed_memory_resource m_resource;

  /// Size for the tested arrays.
  static constexpr std::size_t s_arraySize = 5000;

  /// @name Inputs for the tests
  /// @{

  vecmem::vector<typename T::vector3> m_t1{s_arraySize, &m_resource};
  vecmem::vector<typename T::vector3> m_t2{s_arraySize, &m_resource};
  vecmem::vector<typename T::vector3> m_t3{s_arraySize, &m_resource};

  vecmem::vector<typename T::point2> m_p1{s_arraySize, &m_resource};
  vecmem::vector<typename T::point2> m_p2{s_arraySize, &m_resource};

  vecmem::vector<typename T::vector3> m_v1{s_arraySize, &m_resource};
  vecmem::vector<typename T::vector3> m_v2{s_arraySize, &m_resource};

  /// @}

  /// @name Outputs for the tests
  /// @{

  vecmem::vector<typename T::scalar> m_output_host{s_arraySize, &m_resource};
  vecmem::vector<typename T::scalar> m_output_device{s_arraySize, &m_resource};

  /// @}
};
TYPED_TEST_SUITE_P(test_cuda_basics);

/// Test for some basic 2D "vector operations"
TYPED_TEST_P(test_cuda_basics, vector_2d_ops) {

  // Run the test on the host, and on the/a device.
  execute_host_test<cuda::vector_2d_ops_functor<TypeParam> >(
      this->m_p1.size(), vecmem::get_data(this->m_p1),
      vecmem::get_data(this->m_p2), vecmem::get_data(this->m_output_host));
  execute_cuda_test<cuda::vector_2d_ops_functor<TypeParam> >(
      this->m_p1.size(), vecmem::get_data(this->m_p1),
      vecmem::get_data(this->m_p2), vecmem::get_data(this->m_output_device));

  // Compare the outputs.
  for (std::size_t i = 0; i < this->s_arraySize; ++i) {
    EXPECT_FLOAT_EQ(this->m_output_host[i], this->m_output_device[i]);
  }
}

/// Test for some basic 3D "vector operations"
TYPED_TEST_P(test_cuda_basics, vector_3d_ops) {

  // Run the test on the host, and on the/a device.
  execute_host_test<cuda::vector_3d_ops_functor<TypeParam> >(
      this->m_v1.size(), vecmem::get_data(this->m_v1),
      vecmem::get_data(this->m_v2), vecmem::get_data(this->m_output_host));
  execute_cuda_test<cuda::vector_3d_ops_functor<TypeParam> >(
      this->m_v1.size(), vecmem::get_data(this->m_v1),
      vecmem::get_data(this->m_v2), vecmem::get_data(this->m_output_device));

  // Compare the outputs.
  for (std::size_t i = 0; i < this->s_arraySize; ++i) {
    EXPECT_FLOAT_EQ(this->m_output_host[i], this->m_output_device[i]);
  }
}

/// Test for some operations with @c transform3
TYPED_TEST_P(test_cuda_basics, transform3) {

  // Run the test on the host, and on the/a device.
  execute_host_test<cuda::transform3_ops_functor<TypeParam> >(
      this->m_t1.size(), vecmem::get_data(this->m_t1),
      vecmem::get_data(this->m_t2), vecmem::get_data(this->m_t3),
      vecmem::get_data(this->m_v1), vecmem::get_data(this->m_v2),
      vecmem::get_data(this->m_output_host));
  execute_cuda_test<cuda::transform3_ops_functor<TypeParam> >(
      this->m_t1.size(), vecmem::get_data(this->m_t1),
      vecmem::get_data(this->m_t2), vecmem::get_data(this->m_t3),
      vecmem::get_data(this->m_v1), vecmem::get_data(this->m_v2),
      vecmem::get_data(this->m_output_device));

  // Compare the outputs.
  for (std::size_t i = 0; i < this->s_arraySize; ++i) {
    EXPECT_FLOAT_EQ(this->m_output_host[i], this->m_output_device[i]);
  }
}

/// Test for some operations with @c cartesian2
TYPED_TEST_P(test_cuda_basics, cartesian2) {

  // Run the test on the host, and on the/a device.
  execute_host_test<cuda::cartesian2_ops_functor<TypeParam> >(
      this->m_t1.size(), vecmem::get_data(this->m_t1),
      vecmem::get_data(this->m_t2), vecmem::get_data(this->m_t3),
      vecmem::get_data(this->m_v1), vecmem::get_data(this->m_v2),
      vecmem::get_data(this->m_output_host));
  execute_cuda_test<cuda::cartesian2_ops_functor<TypeParam> >(
      this->m_t1.size(), vecmem::get_data(this->m_t1),
      vecmem::get_data(this->m_t2), vecmem::get_data(this->m_t3),
      vecmem::get_data(this->m_v1), vecmem::get_data(this->m_v2),
      vecmem::get_data(this->m_output_device));

  // Compare the outputs.
  for (std::size_t i = 0; i < this->s_arraySize; ++i) {
    EXPECT_FLOAT_EQ(this->m_output_host[i], this->m_output_device[i]);
  }
}

/// Test for some operations with @c cylindrical2
TYPED_TEST_P(test_cuda_basics, cylindrical2) {

  // Run the test on the host, and on the/a device.
  execute_host_test<cuda::cylindrical2_ops_functor<TypeParam> >(
      this->m_t1.size(), vecmem::get_data(this->m_t1),
      vecmem::get_data(this->m_t2), vecmem::get_data(this->m_t3),
      vecmem::get_data(this->m_v1), vecmem::get_data(this->m_v2),
      vecmem::get_data(this->m_output_host));
  execute_cuda_test<cuda::cylindrical2_ops_functor<TypeParam> >(
      this->m_t1.size(), vecmem::get_data(this->m_t1),
      vecmem::get_data(this->m_t2), vecmem::get_data(this->m_t3),
      vecmem::get_data(this->m_v1), vecmem::get_data(this->m_v2),
      vecmem::get_data(this->m_output_device));

  // Compare the outputs.
  for (std::size_t i = 0; i < this->s_arraySize; ++i) {
    EXPECT_FLOAT_EQ(this->m_output_host[i], this->m_output_device[i]);
  }
}

/// Test for some operations with @c polar2
TYPED_TEST_P(test_cuda_basics, polar2) {

  // Run the test on the host, and on the/a device.
  execute_host_test<cuda::polar2_ops_functor<TypeParam> >(
      this->m_t1.size(), vecmem::get_data(this->m_t1),
      vecmem::get_data(this->m_t2), vecmem::get_data(this->m_t3),
      vecmem::get_data(this->m_v1), vecmem::get_data(this->m_v2),
      vecmem::get_data(this->m_output_host));
  execute_cuda_test<cuda::polar2_ops_functor<TypeParam> >(
      this->m_t1.size(), vecmem::get_data(this->m_t1),
      vecmem::get_data(this->m_t2), vecmem::get_data(this->m_t3),
      vecmem::get_data(this->m_v1), vecmem::get_data(this->m_v2),
      vecmem::get_data(this->m_output_device));

  // Compare the outputs.
  for (std::size_t i = 0; i < this->s_arraySize; ++i) {
    EXPECT_FLOAT_EQ(this->m_output_host[i], this->m_output_device[i]);
  }
}

REGISTER_TYPED_TEST_SUITE_P(test_cuda_basics, vector_2d_ops, vector_3d_ops,
                            transform3, cartesian2, cylindrical2, polar2);
