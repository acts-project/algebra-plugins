/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "benchmark_base.hpp"

// System include(s)
#include <chrono>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

namespace algebra {

template <typename vector_t>
void fill_random_vec(std::vector<vector_t> &);

/// Benchmark for vector operations
template <typename vector_t>
struct vector_bm : public benchmark_base {

  /// Prefix for the benchmark name
  inline static const std::string name{"vector"};

  std::vector<vector_t> a, b, results;

  /// No default construction: Cannot prepare data
  vector_bm() = delete;

  /// Construct from an externally provided configuration @param cfg
  vector_bm(benchmark_base::configuration cfg) : benchmark_base{cfg} {

    const std::size_t n_data{this->m_cfg.n_samples()};

    a.reserve(n_data);
    b.reserve(n_data);

    fill_random_vec(a);
    fill_random_vec(b);
  }

  /// Clear state
  virtual ~vector_bm() {
    a.clear();
    b.clear();
  }
};

/// Benchmark elementwise addition of vectors
template <template <typename> class vector_t, typename scalar_t,
          typename unaryOP>
struct vector_unaryOP_bm : public vector_bm<vector_t<scalar_t>> {
  using base_type = vector_bm<vector_t<scalar_t>>;

  vector_unaryOP_bm() = delete;
  vector_unaryOP_bm(benchmark_base::configuration cfg) : base_type{cfg} {}
  std::string name() const override {
    return base_type::name + "_" + unaryOP::name;
  }

  void operator()(::benchmark::State &state) override {

    const std::size_t n_samples{this->m_cfg.n_samples()};

    // Run the benchmark
    for (auto _ : state) {
      for (std::size_t i{0}; i < n_samples; ++i) {
        ::benchmark::DoNotOptimize(unaryOP{}(this->a[i]));
      }
    }
  }
};

/// Benchmark elementwise addition of vectors
template <template <typename> class vector_t, typename scalar_t,
          typename binaryOP>
struct vector_binaryOP_bm : public vector_bm<vector_t<scalar_t>> {
  using base_type = vector_bm<vector_t<scalar_t>>;

  vector_binaryOP_bm() = delete;
  vector_binaryOP_bm(benchmark_base::configuration cfg) : base_type{cfg} {}
  std::string name() const override {
    return base_type::name + "_" + binaryOP::name;
  }

  void operator()(::benchmark::State &state) override {

    const std::size_t n_samples{this->m_cfg.n_samples()};

    // Run the benchmark
    for (auto _ : state) {
      for (std::size_t i{0}; i < n_samples; ++i) {
        ::benchmark::DoNotOptimize(binaryOP{}(this->a[i], this->b[i]));
      }
    }
  }
};

// Functions to be benchmarked
namespace bench_op {

struct add {
  inline static const std::string name{"add"};
  template <typename vector_t>
  auto operator()(const vector_t &a, const vector_t &b) const {
    return a + b;
  }
};
struct sub {
  inline static const std::string name{"sub"};
  template <typename vector_t>
  auto operator()(const vector_t &a, const vector_t &b) const {
    return a - b;
  }
};
struct dot {
  inline static const std::string name{"dot"};
  template <typename vector_t>
  auto operator()(const vector_t &a, const vector_t &b) const {
    return algebra::vector::dot(a, b);
  }
};
struct cross {
  inline static const std::string name{"cross"};
  template <typename vector_t>
  auto operator()(const vector_t &a, const vector_t &b) const {
    return algebra::vector::cross(a, b);
  }
};

// Macro for declaring vector unary ops
#define ALGEBRA_PLUGINS_BENCH_VECTOR(OP)       \
  struct OP {                                  \
    inline static const std::string name{#OP}; \
    template <typename vector_t>               \
    auto operator()(const vector_t &a) const { \
      return algebra::vector::OP(a);           \
    }                                          \
  };

ALGEBRA_PLUGINS_BENCH_VECTOR(phi)
ALGEBRA_PLUGINS_BENCH_VECTOR(theta)
ALGEBRA_PLUGINS_BENCH_VECTOR(eta)
ALGEBRA_PLUGINS_BENCH_VECTOR(perp)
ALGEBRA_PLUGINS_BENCH_VECTOR(norm)
ALGEBRA_PLUGINS_BENCH_VECTOR(normalize)

}  // namespace bench_op

}  // namespace algebra
