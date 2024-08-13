/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "benchmark_base.hpp"
#include "register_benchmark.hpp"

// System include(s)
#include <chrono>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

namespace algebra {

template <concepts::vector vector_t>
void fill_random_vec(std::vector<vector_t> &);

/// Benchmark for vector operations
template <concepts::vector vector_t>
struct vector_bm : public benchmark_base {

  /// Prefix for the benchmark name
  inline static const std::string name{"vector"};

  std::vector<vector_t> a, b, results;

  /// No default construction: Cannot prepare data
  vector_bm() = delete;

  /// Construct from an externally provided configuration @param cfg
  explicit vector_bm(benchmark_base::configuration cfg) : benchmark_base{cfg} {

    const std::size_t n_data{this->m_cfg.n_samples()};

    a.reserve(n_data);
    b.reserve(n_data);

    fill_random_vec(a);
    fill_random_vec(b);
  }
  vector_bm(const vector_bm &bm) = default;
  vector_bm &operator=(vector_bm &other) = default;

  /// Clear state
  ~vector_bm() override {
    a.clear();
    b.clear();
  }
};

/// Benchmark elementwise addition of vectors
template <template <typename> class vector_t, concepts::scalar scalar_t,
          typename unaryOP>
requires std::invocable<unaryOP, vector_t<scalar_t>> struct vector_unaryOP_bm
    : public vector_bm<vector_t<scalar_t>> {
  using base_type = vector_bm<vector_t<scalar_t>>;

  vector_unaryOP_bm() = delete;
  explicit vector_unaryOP_bm(benchmark_base::configuration cfg)
      : base_type{cfg} {}
  vector_unaryOP_bm(const vector_unaryOP_bm &bm) = default;
  vector_unaryOP_bm &operator=(vector_unaryOP_bm &other) = default;

  std::string name() const override {
    return base_type::name + "_" + unaryOP::name;
  }

  void operator()(::benchmark::State &state) override {

    using result_t = std::invoke_result_t<unaryOP, vector_t<scalar_t>>;

    const std::size_t n_samples{this->m_cfg.n_samples()};

    // Run the benchmark
    for (auto _ : state) {
      for (std::size_t i{0}; i < n_samples; ++i) {
        result_t result = unaryOP{}(this->a[i]);
        ::benchmark::DoNotOptimize(const_cast<const result_t &>(result));
      }
    }
  }
};

/// Benchmark elementwise addition of vectors
template <template <typename> class vector_t, concepts::scalar scalar_t,
          typename binaryOP>
requires std::invocable<binaryOP, vector_t<scalar_t>,
                        vector_t<scalar_t>> struct vector_binaryOP_bm
    : public vector_bm<vector_t<scalar_t>> {
  using base_type = vector_bm<vector_t<scalar_t>>;

  vector_binaryOP_bm() = delete;
  explicit vector_binaryOP_bm(benchmark_base::configuration cfg)
      : base_type{cfg} {}
  vector_binaryOP_bm(const vector_binaryOP_bm &bm) = default;
  vector_binaryOP_bm &operator=(vector_binaryOP_bm &other) = default;

  std::string name() const override {
    return base_type::name + "_" + binaryOP::name;
  }

  void operator()(::benchmark::State &state) override {

    using result_t =
        std::invoke_result_t<binaryOP, vector_t<scalar_t>, vector_t<scalar_t>>;

    const std::size_t n_samples{this->m_cfg.n_samples()};

    // Run the benchmark
    for (auto _ : state) {
      for (std::size_t i{0}; i < n_samples; ++i) {
        result_t result = binaryOP{}(this->a[i], this->b[i]);
        ::benchmark::DoNotOptimize(const_cast<const result_t &>(result));
      }
    }
  }
};

// Functions to be benchmarked
namespace bench_op {

struct add {
  inline static const std::string name{"add"};
  template <concepts::vector vector_t>
  auto operator()(const vector_t &a, const vector_t &b) const {
    return a + b;
  }
};
struct sub {
  inline static const std::string name{"sub"};
  template <concepts::vector vector_t>
  auto operator()(const vector_t &a, const vector_t &b) const {
    return a - b;
  }
};
struct dot {
  inline static const std::string name{"dot"};
  template <concepts::vector vector_t>
  auto operator()(const vector_t &a, const vector_t &b) const {
    return algebra::vector::dot(a, b);
  }
};
struct cross {
  inline static const std::string name{"cross"};
  template <concepts::vector vector_t>
  auto operator()(const vector_t &a, const vector_t &b) const {
    return algebra::vector::cross(a, b);
  }
};

// Macro for declaring vector unary ops
#define ALGEBRA_PLUGINS_BENCH_VECTOR(OP)       \
  struct OP {                                  \
    inline static const std::string name{#OP}; \
    template <concepts::vector vector_t>       \
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

// Macro for registering all vector benchmarks
#define ALGEBRA_PLUGINS_REGISTER_VECTOR_BENCH(CFG)         \
  algebra::register_benchmark<add_f_t>(CFG, "_single");    \
  algebra::register_benchmark<add_d_t>(CFG, "_double");    \
  algebra::register_benchmark<sub_f_t>(CFG, "_single");    \
  algebra::register_benchmark<sub_d_t>(CFG, "_double");    \
  algebra::register_benchmark<dot_f_t>(CFG, "_single");    \
  algebra::register_benchmark<dot_d_t>(CFG, "_double");    \
  algebra::register_benchmark<cross_f_t>(CFG, "_single");  \
  algebra::register_benchmark<cross_d_t>(CFG, "_double");  \
  algebra::register_benchmark<normlz_f_t>(CFG, "_single"); \
  algebra::register_benchmark<normlz_d_t>(CFG, "_double"); \
                                                           \
  algebra::register_benchmark<phi_f_t>(CFG, "_single");    \
  algebra::register_benchmark<phi_d_t>(CFG, "_double");    \
  algebra::register_benchmark<theta_f_t>(CFG, "_single");  \
  algebra::register_benchmark<theta_d_t>(CFG, "_double");  \
  algebra::register_benchmark<perp_f_t>(CFG, "_single");   \
  algebra::register_benchmark<perp_d_t>(CFG, "_double");   \
  algebra::register_benchmark<norm_f_t>(CFG, "_single");   \
  algebra::register_benchmark<norm_d_t>(CFG, "_double");   \
  algebra::register_benchmark<eta_f_t>(CFG, "_single");    \
  algebra::register_benchmark<eta_d_t>(CFG, "_double");

}  // namespace algebra
