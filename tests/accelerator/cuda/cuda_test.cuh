/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Test include(s).
#include "cuda_error_check.hpp"

template <class FUNCTOR, typename... ARGS>
__global__ void testKernel(std::size_t arraySizes, ARGS... args) {

  // Find the current index that we need to process.
  const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= arraySizes) {
    return;
  }

  // Execute the test functor for this index.
  FUNCTOR()(i, args...);
  return;
}

/// Execute a test functor on a device, on @c arraySizes threads
template <class FUNCTOR, class... ARGS>
void runOnDevice(std::size_t arraySizes, ARGS... args) {

  // Number of threads per execution block.
  int nThreadsPerBlock = 1024;

  // If the arrays are not even this large, then reduce the value to the
  // size of the arrays.
  if (arraySizes < nThreadsPerBlock) {
    nThreadsPerBlock = arraySizes;
  }

  // Launch the test on the device.
  const int nBlocks = ((arraySizes + nThreadsPerBlock - 1) / nThreadsPerBlock);
  testKernel<FUNCTOR><<<nBlocks, nThreadsPerBlock>>>(arraySizes, args...);

  // Check whether it succeeded to run.
  CUDA_ERROR_CHECK(cudaGetLastError());
  CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

/// Execute a test functor on the host
template <class FUNCTOR, class... ARGS>
void runOnHost(std::size_t arraySizes, ARGS... args) {

  // Instantiate the functor.
  auto functor = FUNCTOR();

  // Execute the functor on all elements of the array(s).
  for (std::size_t i = 0; i < arraySizes; ++i) {
    functor(i, args...);
  }
}
