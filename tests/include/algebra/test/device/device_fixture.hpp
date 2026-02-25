/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "algebra/concepts.hpp"
#include "algebra/utils/approximately_equal.hpp"

// Test include(s).
#include "algebra/test/framework/test_base.hpp"

// VecMem include(s).
#include <vecmem/containers/vector.hpp>
#include <vecmem/memory/memory_resource.hpp>

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <cmath>
#include <memory>
#include <type_traits>

namespace algebra::test {

/// Data fixture for device test cases
template <algebra::concepts::algebra A, typename R>
  requires std::is_default_constructible_v<R>
class device_fixture : public testing::Test, public test_base<A> {

 public:
  /// Constructor, setting up the vecmem resource
  device_fixture(vecmem::memory_resource& mr) : m_resource(mr) {}

 protected:
  /// Set up the host and device results
  virtual void SetUp() override {
    m_output_host =
        std::make_unique<vecmem::vector<R>>(this->size(), &m_resource);
    m_output_device =
        std::make_unique<vecmem::vector<R>>(this->size(), &m_resource);

    for (std::size_t i = 0; i < this->size(); ++i) {
      m_output_host->at(i) = {};
      m_output_device->at(i) = {};
    }
  }

  /// Function resetting the output after the test
  virtual void TearDown() override {
    m_output_host.reset();
    m_output_device.reset();
  }

  /// @returns the number of results
  constexpr std::size_t size() const { return m_size; }

  /// @returns access to the memory resource of the device test suite
  constexpr vecmem::memory_resource& resource() { return m_resource; }

  /// Compare the outputs, after the data processing is finished.
  void compareOutputs() const {
    for (std::size_t i = 0u; i < this->size(); ++i) {
      // Use the inbuilt comparator for algebra types
      if constexpr (algebra::concepts::value<R> ||
                    algebra::concepts::scalar<R> ||
                    algebra::concepts::vector<R> ||
                    algebra::concepts::matrix<R>) {

        using value_t = algebra::traits::value_t<R>;
        constexpr value_t max_rel_error{this->m_epsilon};

        EXPECT_TRUE(algebra::approx_equal(
            m_output_host->at(i), m_output_device->at(i), max_rel_error))
            << "error at i = " << i << "\nHOST:\n"
            << m_output_host->at(i) << "\nDEVICE:\n"
            << m_output_device->at(i);

      } else {
        EXPECT_EQ(m_output_host->at(i), m_output_device->at(i))
            << "error at i = " << i << "\nHOST:\n"
            << m_output_host->at(i) << "\nDEVICE:\n"
            << m_output_device->at(i);
      }
    }
  }

  /// Memory resource provided by the derived class
  vecmem::memory_resource& m_resource;

  /// Size for the tested arrays.
  static constexpr std::size_t m_size = 5000u;

  /// @name Outputs for the tests
  /// @{
  std::unique_ptr<vecmem::vector<R>> m_output_host;
  std::unique_ptr<vecmem::vector<R>> m_output_device;
  /// @}

};  // class device_fixture

}  // namespace algebra::test
