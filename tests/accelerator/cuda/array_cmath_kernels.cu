/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Test include(s).
#include "array_cmath_kernels.hpp"
#include "cuda_test.cuh"
#include "test_device_basics.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

namespace {

/// Functor running @c test_device_basics::vector_2d_ops
struct array_cmath_2d_vector_ops_func {

  __host__ __device__ void operator()(
      std::size_t i, vecmem::data::vector_view<const algebra::array::point2> a,
      vecmem::data::vector_view<const algebra::array::point2> b,
      vecmem::data::vector_view<algebra::scalar>& output) {

    // Instantiate the tester object.
    test_device_basics<test_types<
        algebra::scalar, algebra::array::point2, algebra::array::point3,
        algebra::array::vector2, algebra::array::vector3,
        algebra::array::transform3, algebra::array::cartesian2,
        algebra::array::polar2, algebra::array::cylindrical2>>
        tester;

    // Create the VecMem vector(s).
    vecmem::device_vector<const algebra::array::point2> vec_a(a), vec_b(b);
    vecmem::device_vector<algebra::scalar> vec_output(output);

    // Perform the operation.
    vec_output[i] = tester.vector_2d_ops(vec_a[i], vec_b[i]);
  }
};

}  // namespace

void array_cmath_2d_vector_ops(
    vecmem::data::vector_view<const algebra::array::point2> a,
    vecmem::data::vector_view<const algebra::array::point2> b,
    vecmem::data::vector_view<algebra::scalar> output_host,
    vecmem::data::vector_view<algebra::scalar> output_device) {

  // Run the test on the host, and on the/a device.
  runOnHost<array_cmath_2d_vector_ops_func>(a.size(), a, b, output_host);
  runOnDevice<array_cmath_2d_vector_ops_func>(a.size(), a, b, output_device);
}
