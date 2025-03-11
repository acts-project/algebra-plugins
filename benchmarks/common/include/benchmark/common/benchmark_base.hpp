/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Benchmark include
#include <benchmark/benchmark.h>

// System include(s)
#include <ostream>
#include <string>

namespace algebra {

/// Base type for linear algebra benchmarks with google benchmark
struct benchmark_base {
  /// Local configuration type
  struct configuration {
    /// Size of data sample to be used in benchmark
    std::size_t m_samples{100u};
    // Sleep after building data sample
    bool m_sleep = false;
    // Number of seconds to sleep
    std::size_t m_n_sleep{1u};

    /// Setters
    /// @{
    configuration& n_samples(std::size_t n) {
      m_samples = n;
      return *this;
    }
    configuration& do_sleep(bool b) {
      m_sleep = b;
      return *this;
    }
    configuration& n_sleep(std::size_t n) {
      m_n_sleep = n;
      m_sleep = true;
      return *this;
    }
    /// @}

    /// Getters
    /// @{
    std::size_t n_samples() const { return m_samples; }
    constexpr bool do_sleep() const { return m_sleep; }
    constexpr std::size_t n_sleep() const { return m_n_sleep; }
    /// @}

    /// Print configuration
    friend std::ostream& operator<<(std::ostream& os,
                                    const benchmark_base::configuration& cfg) {
      os << " -> running:\t " << cfg.n_samples() << " samples" << std::endl;
      if (cfg.do_sleep()) {
        os << " -> cool down:\t " << cfg.n_sleep() << "s" << std::endl;
      }
      os << std::endl;
      return os;
    }
  };

  /// The benchmark configuration
  configuration m_cfg{};

  /// Default construction
  benchmark_base() = default;

  /// Construct from an externally provided configuration @param cfg
  explicit benchmark_base(configuration cfg) : m_cfg{cfg} {}

  /// @returns the benchmark configuration
  configuration& config() { return m_cfg; }

  /// Default destructor
  virtual ~benchmark_base() = default;

  /// @returns the benchmark name
  virtual constexpr std::string name() const = 0;

  /// Benchmark case
  virtual void operator()(::benchmark::State&) const = 0;
};

}  // namespace algebra
