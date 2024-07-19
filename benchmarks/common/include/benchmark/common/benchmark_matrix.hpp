/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2024 CERN for the benefit of the ACTS project
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

template <typename matrix_t>
void fill_random_matrix(std::vector<matrix_t> &);

template <typename vector_t>
void fill_random_vec(std::vector<vector_t> &);

/// Benchmark for vector operations
template <typename matrix_t>
struct matrix_bm : public benchmark_base {

  /// Prefix for the benchmark name
  inline static const std::string name{"matrix"};

  std::vector<matrix_t> a, b;

  /// No default construction: Cannot prepare data
  matrix_bm() = delete;

  /// Construct from an externally provided configuration @param cfg
  matrix_bm(benchmark_base::configuration cfg) : benchmark_base{cfg} {

    const std::size_t n_data{this->m_cfg.n_samples()};

    a.reserve(n_data);
    b.reserve(n_data);

    fill_random_matrix(a);
    fill_random_matrix(b);
  }

  /// Clear state
  virtual ~matrix_bm() {
    a.clear();
    b.clear();
  }
};

/// Benchmark elementwise addition of vectors

template <typename matrix_t, typename unaryOP>
struct matrix_unaryOP_bm : public matrix_bm<matrix_t> {
  using base_type = matrix_bm<matrix_t>;

  matrix_unaryOP_bm() = delete;
  matrix_unaryOP_bm(benchmark_base::configuration cfg) : base_type{cfg} {}
  std::string name() const override {
    return base_type::name + "_" + unaryOP::name;
  }

  void operator()(::benchmark::State &state) override {

    using result_t = std::invoke_result_t<unaryOP, matrix_t>;

    const std::size_t n_samples{this->m_cfg.n_samples()};

    // Run the benchmark
    for (auto _ : state) {
      for (std::size_t i{0}; i < n_samples; ++i) {
        result_t result = unaryOP{}(this->a[i]);
        ::benchmark::DoNotOptimize(result);
      }
    }
  }
};

/// Benchmark elementwise addition of vectors
template <typename matrix_t, typename binaryOP>
struct matrix_binaryOP_bm : public matrix_bm<matrix_t> {
  using base_type = matrix_bm<matrix_t>;

  matrix_binaryOP_bm() = delete;
  matrix_binaryOP_bm(benchmark_base::configuration cfg) : base_type{cfg} {}
  std::string name() const override {
    return base_type::name + "_" + binaryOP::name;
  }

  void operator()(::benchmark::State &state) override {

    using result_t = std::invoke_result_t<binaryOP, matrix_t, matrix_t>;

    const std::size_t n_samples{this->m_cfg.n_samples()};

    // Run the benchmark
    for (auto _ : state) {
      for (std::size_t i{0}; i < n_samples; ++i) {
        result_t result = binaryOP{}(this->a[i], this->b[i]);
        ::benchmark::DoNotOptimize(result);
      }
    }
  }
};

/// Benchmark matrix vector multiplication
template <typename matrix_t, typename vector_t>
struct matrix_vector_bm : public matrix_bm<matrix_t> {
  using base_type = matrix_bm<matrix_t>;

  std::vector<vector_t> v;

  matrix_vector_bm() = delete;
  matrix_vector_bm(benchmark_base::configuration cfg) : base_type{cfg} {

    v.reserve(this->m_cfg.n_samples());

    fill_random_vec(v);
  }
  /// Clear state
  virtual ~matrix_vector_bm() { v.clear(); }
  std::string name() const override { return base_type::name + "_vector"; }

  void operator()(::benchmark::State &state) override {

    const std::size_t n_samples{this->m_cfg.n_samples()};

    // Run the benchmark
    for (auto _ : state) {
      for (std::size_t i{0}; i < n_samples; ++i) {
        vector_t result = this->a[i] * this->v[i];
        ::benchmark::DoNotOptimize(result);
      }
    }
  }
};

// Functions to be benchmarked
namespace bench_op {

struct add {
  inline static const std::string name{"add"};
  template <typename matrix_t>
  matrix_t operator()(const matrix_t &a, const matrix_t &b) const {
    return a + b;
  }
};
struct sub {
  inline static const std::string name{"sub"};
  template <typename matrix_t>
  matrix_t operator()(const matrix_t &a, const matrix_t &b) const {
    return a - b;
  }
};
struct mul {
  inline static const std::string name{"mul"};
  template <typename matrix_t>
  matrix_t operator()(const matrix_t &a, const matrix_t &b) const {
    return a * b;
  }
};
struct transpose {
  inline static const std::string name{"transpose"};
  template <typename matrix_t>
  auto operator()(const matrix_t &a) const {
    return algebra::matrix::transpose(a);
  }
};

}  // namespace bench_op

}  // namespace algebra
