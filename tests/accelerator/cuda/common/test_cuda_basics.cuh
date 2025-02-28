/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Test include(s).
#include "cuda_error_check.hpp"
#include "execute_cuda_test.cuh"
#include "execute_host_test.hpp"
#include "test_basics_base.hpp"
#include "test_basics_functors.hpp"

// GoogleTest include(s).
#include <cuda_runtime_api.h>
#include <gtest/gtest.h>

/// Test case class, to be specialised for the different plugins
template <typename T>
class test_cuda_basics : public test_basics_base<T> {

 public:
  /// Constructor, setting up the inputs for all of the tests
  test_cuda_basics() : test_basics_base<T>() {}

 protected:
  /// Function setting things up for a test
  virtual void SetUp() override {

    // Create the vectors.
    CUDA_ERROR_CHECK(cudaMallocManaged(
        &this->m_t1, this->s_arraySize * sizeof(typename T::vector3)));
    CUDA_ERROR_CHECK(cudaMallocManaged(
        &this->m_t2, this->s_arraySize * sizeof(typename T::vector3)));
    CUDA_ERROR_CHECK(cudaMallocManaged(
        &this->m_t3, this->s_arraySize * sizeof(typename T::vector3)));

    CUDA_ERROR_CHECK(cudaMallocManaged(
        &this->m_p1, this->s_arraySize * sizeof(typename T::point2)));
    CUDA_ERROR_CHECK(cudaMallocManaged(
        &this->m_p2, this->s_arraySize * sizeof(typename T::point2)));

    CUDA_ERROR_CHECK(cudaMallocManaged(
        &this->m_v1, this->s_arraySize * sizeof(typename T::vector3)));
    CUDA_ERROR_CHECK(cudaMallocManaged(
        &this->m_v2, this->s_arraySize * sizeof(typename T::vector3)));

    CUDA_ERROR_CHECK(cudaMallocManaged(
        &this->m_m1,
        this->s_arraySize * sizeof(typename T::template matrix<6, 4>)));
    CUDA_ERROR_CHECK(cudaMallocManaged(
        &this->m_m2,
        this->s_arraySize * sizeof(typename T::template matrix<2, 2>)));

    CUDA_ERROR_CHECK(
        cudaMallocManaged(&this->m_output_device,
                          this->s_arraySize * sizeof(typename T::scalar)));

    this->m_output_host = static_cast<T::scalar *>(
        std::calloc(this->s_arraySize, sizeof(typename T::scalar)));

    test_basics_base<T>::fill_pointers();
  }

  /// Function tearing things down after the test
  virtual void TearDown() override {
    // Delete the vectors.
    CUDA_ERROR_CHECK(cudaFree(this->m_t1));
    CUDA_ERROR_CHECK(cudaFree(this->m_t2));
    CUDA_ERROR_CHECK(cudaFree(this->m_t3));
    CUDA_ERROR_CHECK(cudaFree(this->m_p1));
    CUDA_ERROR_CHECK(cudaFree(this->m_p2));
    CUDA_ERROR_CHECK(cudaFree(this->m_v1));
    CUDA_ERROR_CHECK(cudaFree(this->m_v2));
    CUDA_ERROR_CHECK(cudaFree(this->m_m1));
    CUDA_ERROR_CHECK(cudaFree(this->m_m2));
    CUDA_ERROR_CHECK(cudaFree(this->m_output_device));

    free(this->m_output_host);
  }
};
TYPED_TEST_SUITE_P(test_cuda_basics);

/// Test for some basic 2D "vector operations"
TYPED_TEST_P(test_cuda_basics, vector_2d_ops) {

  // Run the test on the host, and on the/a device.
  execute_host_test<vector_2d_ops_functor<TypeParam> >(
      this->m_p1->size(), this->m_p1, this->m_p2, this->m_output_host);
  execute_cuda_test<vector_2d_ops_functor<TypeParam> >(
      this->m_p1->size(), this->m_p1, this->m_p2, this->m_output_device);

  // Compare the outputs.
  this->compareOutputs();
}

/// Test for some basic 3D "vector operations"
TYPED_TEST_P(test_cuda_basics, vector_3d_ops) {

  // This test is just not numerically stable at float precision in optimized
  // mode for some reason. :-(
#ifdef NDEBUG
  if (typeid(typename TypeParam::scalar) == typeid(float)) {
    GTEST_SKIP();
  }
#endif  // NDEBUG

  // Run the test on the host, and on the/a device.
  execute_host_test<vector_3d_ops_functor<TypeParam> >(
      this->m_v1->size(), this->m_v1, this->m_v2, this->m_output_host);
  execute_cuda_test<vector_3d_ops_functor<TypeParam> >(
      this->m_v1->size(), this->m_v1, this->m_v2, this->m_output_device);

  // Compare the outputs.
  this->compareOutputs();
}

/// Test for handling matrices
TYPED_TEST_P(test_cuda_basics, matrix64_ops) {

  // Run the test on the host, and on the/a device.
  execute_host_test<matrix64_ops_functor<TypeParam> >(
      this->m_m1->size(), this->m_m1, this->m_output_host);
  execute_cuda_test<matrix64_ops_functor<TypeParam> >(
      this->m_m1->size(), this->m_m1, this->m_output_device);

  // Compare the outputs.
  this->compareOutputs();
}

/// Test for handling matrices
TYPED_TEST_P(test_cuda_basics, matrix22_ops) {

  // Run the test on the host, and on the/a device.
  execute_host_test<matrix22_ops_functor<TypeParam> >(
      this->m_m2->size(), this->m_m2, this->m_output_host);
  execute_cuda_test<matrix22_ops_functor<TypeParam> >(
      this->m_m2->size(), this->m_m2, this->m_output_device);

  // Compare the outputs.
  this->compareOutputs();
}

/// Test for some operations with @c transform3
TYPED_TEST_P(test_cuda_basics, transform3) {

  // Run the test on the host, and on the/a device.
  execute_host_test<transform3_ops_functor<TypeParam> >(
      this->m_t1->size(), this->m_t1, this->m_t2, this->m_t3, this->m_v1,
      this->m_v2, this->m_output_host);
  execute_cuda_test<transform3_ops_functor<TypeParam> >(
      this->m_t1->size(), this->m_t1, this->m_t2, this->m_t3, this->m_v1,
      this->m_v2, this->m_output_device);

  // Compare the outputs.
  this->compareOutputs();
}

REGISTER_TYPED_TEST_SUITE_P(test_cuda_basics, vector_2d_ops, vector_3d_ops,
                            matrix64_ops, matrix22_ops, transform3);
