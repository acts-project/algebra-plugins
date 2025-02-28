/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/qualifiers.hpp"

// Local include(s).
#include "test_device_basics.hpp"

/// Base class for all of the functors
template <typename T>
class functor_base {

 protected:
  /// Tester object
  test_device_basics<T> m_tester;
};

/// Functor running @c test_device_basics::vector_2d_ops
template <typename T>
class vector_2d_ops_functor : public functor_base<T> {

 public:
  ALGEBRA_HOST_DEVICE void operator()(std::size_t i,
                                      const typename T::point2* a,
                                      const typename T::point2* b,
                                      typename T::scalar* output) const {

    output[i] = this->m_tester.vector_2d_ops(a[i], b[i]);
  }
};

/// Functor running @c test_device_basics::vector_3d_ops
template <typename T>
class vector_3d_ops_functor : public functor_base<T> {

 public:
  ALGEBRA_HOST_DEVICE void operator()(std::size_t i,
                                      const typename T::vector3* a,
                                      const typename T::vector3* b,
                                      typename T::scalar* output) const {

    output[i] = this->m_tester.vector_3d_ops(a[i], b[i]);
  }
};

/// Functor running @c test_device_basics::matrix64_ops
template <typename T>
class matrix64_ops_functor : public functor_base<T> {

 public:
  ALGEBRA_HOST_DEVICE void operator()(
      std::size_t i, const typename T::template matrix<6, 4>* m,
      typename T::scalar* output) const {

    output[i] = this->m_tester.matrix64_ops(m[i]);
  }
};

/// Functor running @c test_device_basics::matrix22_ops
template <typename T>
class matrix22_ops_functor : public functor_base<T> {

 public:
  ALGEBRA_HOST_DEVICE void operator()(
      std::size_t i, const typename T::template matrix<2, 2>* m,
      typename T::scalar* output) const {

    output[i] = this->m_tester.matrix22_ops(m[i]);
  }
};

/// Functor running @c test_device_basics::transform3_ops
template <typename T>
class transform3_ops_functor : public functor_base<T> {

 public:
  ALGEBRA_HOST_DEVICE void operator()(std::size_t i,
                                      const typename T::vector3* t1,
                                      const typename T::vector3* t2,
                                      const typename T::vector3* t3,
                                      const typename T::vector3* a,
                                      const typename T::vector3* b,
                                      typename T::scalar* output) const {

    output[i] = this->m_tester.transform3_ops(t1[i], t2[i], t3[i], a[i], b[i]);
  }
};
