/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Test include(s).
#include "array_cmath_kernels.hpp"
#include "execute_cuda_test.cuh"
#include "execute_host_test.hpp"
#include "test_device_basics.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

namespace {

/// Base class for all of the functors
struct array_cmath_func_base {

  /// Tester object
  test_device_basics<test_types<
      algebra::scalar, algebra::array::point2, algebra::array::point3,
      algebra::array::vector2, algebra::array::vector3,
      algebra::array::transform3, algebra::array::cartesian2,
      algebra::array::polar2, algebra::array::cylindrical2>>
      m_tester;
};

/// Functor running @c test_device_basics::vector_2d_ops
struct array_cmath_2d_vector_ops_func : public array_cmath_func_base {

  __host__ __device__ void operator()(
      std::size_t i, vecmem::data::vector_view<const algebra::array::point2> a,
      vecmem::data::vector_view<const algebra::array::point2> b,
      vecmem::data::vector_view<algebra::scalar>& output) {

    // Create the VecMem vector(s).
    vecmem::device_vector<const algebra::array::point2> vec_a(a), vec_b(b);
    vecmem::device_vector<algebra::scalar> vec_output(output);

    // Perform the operation.
    vec_output[i] = m_tester.vector_2d_ops(vec_a[i], vec_b[i]);
  }
};

/// Functor running @c test_device_basics::vector_3d_ops
struct array_cmath_3d_vector_ops_func : public array_cmath_func_base {

  __host__ __device__ void operator()(
      std::size_t i, vecmem::data::vector_view<const algebra::array::vector3> a,
      vecmem::data::vector_view<const algebra::array::vector3> b,
      vecmem::data::vector_view<algebra::scalar>& output) {

    // Create the VecMem vector(s).
    vecmem::device_vector<const algebra::array::vector3> vec_a(a), vec_b(b);
    vecmem::device_vector<algebra::scalar> vec_output(output);

    // Perform the operation.
    vec_output[i] = m_tester.vector_3d_ops(vec_a[i], vec_b[i]);
  }
};

/// Functor running @c test_device_basics::transform3_ops
struct array_cmath_transform3_ops_func : public array_cmath_func_base {

  __host__ __device__ void operator()(
      std::size_t i,
      vecmem::data::vector_view<const algebra::array::vector3> t1,
      vecmem::data::vector_view<const algebra::array::vector3> t2,
      vecmem::data::vector_view<const algebra::array::vector3> t3,
      vecmem::data::vector_view<const algebra::array::vector3> a,
      vecmem::data::vector_view<const algebra::array::vector3> b,
      vecmem::data::vector_view<algebra::scalar>& output) {

    // Create the VecMem vector(s).
    vecmem::device_vector<const algebra::array::vector3> vec_t1(t1), vec_t2(t2),
        vec_t3(t3), vec_a(a), vec_b(b);
    vecmem::device_vector<algebra::scalar> vec_output(output);

    // Perform the operation.
    vec_output[i] = m_tester.transform3_ops(vec_t1[i], vec_t2[i], vec_t3[i],
                                            vec_a[i], vec_b[i]);
  }
};

/// Functor running @c test_device_basics::cartesian2_ops
struct array_cmath_cartesian2_ops_func : public array_cmath_func_base {

  __host__ __device__ void operator()(
      std::size_t i,
      vecmem::data::vector_view<const algebra::array::vector3> t1,
      vecmem::data::vector_view<const algebra::array::vector3> t2,
      vecmem::data::vector_view<const algebra::array::vector3> t3,
      vecmem::data::vector_view<const algebra::array::vector3> a,
      vecmem::data::vector_view<const algebra::array::vector3> b,
      vecmem::data::vector_view<algebra::scalar>& output) {

    // Create the VecMem vector(s).
    vecmem::device_vector<const algebra::array::vector3> vec_t1(t1), vec_t2(t2),
        vec_t3(t3), vec_a(a), vec_b(b);
    vecmem::device_vector<algebra::scalar> vec_output(output);

    // Perform the operation.
    vec_output[i] = m_tester.cartesian2_ops(vec_t1[i], vec_t2[i], vec_t3[i],
                                            vec_a[i], vec_b[i]);
  }
};

/// Functor running @c test_device_basics::cylindrical2_ops
struct array_cmath_cylindrical2_ops_func : public array_cmath_func_base {

  __host__ __device__ void operator()(
      std::size_t i,
      vecmem::data::vector_view<const algebra::array::vector3> t1,
      vecmem::data::vector_view<const algebra::array::vector3> t2,
      vecmem::data::vector_view<const algebra::array::vector3> t3,
      vecmem::data::vector_view<const algebra::array::vector3> a,
      vecmem::data::vector_view<const algebra::array::vector3> b,
      vecmem::data::vector_view<algebra::scalar>& output) {

    // Create the VecMem vector(s).
    vecmem::device_vector<const algebra::array::vector3> vec_t1(t1), vec_t2(t2),
        vec_t3(t3), vec_a(a), vec_b(b);
    vecmem::device_vector<algebra::scalar> vec_output(output);

    // Perform the operation.
    vec_output[i] = m_tester.cylindrical2_ops(vec_t1[i], vec_t2[i], vec_t3[i],
                                              vec_a[i], vec_b[i]);
  }
};

/// Functor running @c test_device_basics::polar2_ops
struct array_cmath_polar2_ops_func : public array_cmath_func_base {

  __host__ __device__ void operator()(
      std::size_t i,
      vecmem::data::vector_view<const algebra::array::vector3> t1,
      vecmem::data::vector_view<const algebra::array::vector3> t2,
      vecmem::data::vector_view<const algebra::array::vector3> t3,
      vecmem::data::vector_view<const algebra::array::vector3> a,
      vecmem::data::vector_view<const algebra::array::vector3> b,
      vecmem::data::vector_view<algebra::scalar>& output) {

    // Create the VecMem vector(s).
    vecmem::device_vector<const algebra::array::vector3> vec_t1(t1), vec_t2(t2),
        vec_t3(t3), vec_a(a), vec_b(b);
    vecmem::device_vector<algebra::scalar> vec_output(output);

    // Perform the operation.
    vec_output[i] = m_tester.polar2_ops(vec_t1[i], vec_t2[i], vec_t3[i],
                                        vec_a[i], vec_b[i]);
  }
};

}  // namespace

void array_cmath_2d_vector_ops(
    vecmem::data::vector_view<const algebra::array::point2> a,
    vecmem::data::vector_view<const algebra::array::point2> b,
    vecmem::data::vector_view<algebra::scalar> output_host,
    vecmem::data::vector_view<algebra::scalar> output_device) {

  // Run the test on the host, and on the/a device.
  execute_host_test<array_cmath_2d_vector_ops_func>(a.size(), a, b,
                                                    output_host);
  execute_cuda_test<array_cmath_2d_vector_ops_func>(a.size(), a, b,
                                                    output_device);
}

void array_cmath_3d_vector_ops(
    vecmem::data::vector_view<const algebra::array::vector3> a,
    vecmem::data::vector_view<const algebra::array::vector3> b,
    vecmem::data::vector_view<algebra::scalar> output_host,
    vecmem::data::vector_view<algebra::scalar> output_device) {

  // Run the test on the host, and on the/a device.
  execute_host_test<array_cmath_3d_vector_ops_func>(a.size(), a, b,
                                                    output_host);
  execute_cuda_test<array_cmath_3d_vector_ops_func>(a.size(), a, b,
                                                    output_device);
}

void array_cmath_transform3_ops(
    vecmem::data::vector_view<const algebra::array::vector3> t1,
    vecmem::data::vector_view<const algebra::array::vector3> t2,
    vecmem::data::vector_view<const algebra::array::vector3> t3,
    vecmem::data::vector_view<const algebra::array::vector3> a,
    vecmem::data::vector_view<const algebra::array::vector3> b,
    vecmem::data::vector_view<algebra::scalar> output_host,
    vecmem::data::vector_view<algebra::scalar> output_device) {

  // Run the test on the host, and on the/a device.
  execute_host_test<array_cmath_transform3_ops_func>(a.size(), t1, t2, t3, a, b,
                                                     output_host);
  execute_cuda_test<array_cmath_transform3_ops_func>(a.size(), t1, t2, t3, a, b,
                                                     output_device);
}

void array_cmath_cartesian2_ops(
    vecmem::data::vector_view<const algebra::array::vector3> t1,
    vecmem::data::vector_view<const algebra::array::vector3> t2,
    vecmem::data::vector_view<const algebra::array::vector3> t3,
    vecmem::data::vector_view<const algebra::array::vector3> a,
    vecmem::data::vector_view<const algebra::array::vector3> b,
    vecmem::data::vector_view<algebra::scalar> output_host,
    vecmem::data::vector_view<algebra::scalar> output_device) {

  // Run the test on the host, and on the/a device.
  execute_host_test<array_cmath_cartesian2_ops_func>(a.size(), t1, t2, t3, a, b,
                                                     output_host);
  execute_cuda_test<array_cmath_cartesian2_ops_func>(a.size(), t1, t2, t3, a, b,
                                                     output_device);
}

void array_cmath_cylindrical2_ops(
    vecmem::data::vector_view<const algebra::array::vector3> t1,
    vecmem::data::vector_view<const algebra::array::vector3> t2,
    vecmem::data::vector_view<const algebra::array::vector3> t3,
    vecmem::data::vector_view<const algebra::array::vector3> a,
    vecmem::data::vector_view<const algebra::array::vector3> b,
    vecmem::data::vector_view<algebra::scalar> output_host,
    vecmem::data::vector_view<algebra::scalar> output_device) {

  // Run the test on the host, and on the/a device.
  execute_host_test<array_cmath_cylindrical2_ops_func>(a.size(), t1, t2, t3, a,
                                                       b, output_host);
  execute_cuda_test<array_cmath_cylindrical2_ops_func>(a.size(), t1, t2, t3, a,
                                                       b, output_device);
}

void array_cmath_polar2_ops(
    vecmem::data::vector_view<const algebra::array::vector3> t1,
    vecmem::data::vector_view<const algebra::array::vector3> t2,
    vecmem::data::vector_view<const algebra::array::vector3> t3,
    vecmem::data::vector_view<const algebra::array::vector3> a,
    vecmem::data::vector_view<const algebra::array::vector3> b,
    vecmem::data::vector_view<algebra::scalar> output_host,
    vecmem::data::vector_view<algebra::scalar> output_device) {

  // Run the test on the host, and on the/a device.
  execute_host_test<array_cmath_polar2_ops_func>(a.size(), t1, t2, t3, a, b,
                                                 output_host);
  execute_cuda_test<array_cmath_polar2_ops_func>(a.size(), t1, t2, t3, a, b,
                                                 output_device);
}
