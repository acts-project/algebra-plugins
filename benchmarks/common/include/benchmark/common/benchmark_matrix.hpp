/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "benchmark_base.hpp"
#include "register_benchmark.hpp"

// System include(s)
#include <string_view>
#include <vector>

namespace algebra {

template <concepts::matrix matrix_t>
void fill_random_matrix(std::vector<matrix_t>&);

template <concepts::vector vector_t>
void fill_random_vec(std::vector<vector_t>&);

/// Benchmark for matrix operations
template <concepts::matrix matrix_t>
struct matrix_bm : public benchmark_base {

  /// Prefix for the benchmark name
  static constexpr std::string_view name{"matrix"};

  std::vector<matrix_t> a;
  std::vector<matrix_t> b;

  /// No default construction: Cannot prepare data
  matrix_bm() = delete;

  /// Construct from an externally provided configuration @param cfg
  explicit matrix_bm(benchmark_base::configuration cfg) : benchmark_base{cfg} {

    const std::size_t n_data{this->m_cfg.n_samples()};

    a.reserve(n_data);
    b.reserve(n_data);

    fill_random_matrix(a);
    fill_random_matrix(b);
  }

  matrix_bm(const matrix_bm& bm) = default;
  matrix_bm& operator=(matrix_bm& other) = default;

  /// Clear state
  ~matrix_bm() override {
    a.clear();
    b.clear();
  }
};

/// Benchmark operations on a single matrix (transpose, inverse etc)
template <concepts::matrix matrix_t, typename unaryOP>
requires std::invocable<unaryOP, matrix_t> struct matrix_unaryOP_bm
    : public matrix_bm<matrix_t> {
  using base_type = matrix_bm<matrix_t>;

  matrix_unaryOP_bm() = delete;
  explicit matrix_unaryOP_bm(benchmark_base::configuration cfg)
      : base_type{cfg} {}
  matrix_unaryOP_bm(const matrix_unaryOP_bm& bm) = default;
  matrix_unaryOP_bm& operator=(matrix_unaryOP_bm& other) = default;

  constexpr std::string name() const override {
    return std::string{base_type::name} + "_" + std::string{unaryOP::name};
  }

  inline void operator()(::benchmark::State& state) const override {

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

/// Benchmark elementwise addition/subtraction/multiplication of matrices
template <concepts::matrix matrix_t, typename binaryOP>
requires std::invocable<binaryOP, matrix_t, matrix_t> struct matrix_binaryOP_bm
    : public matrix_bm<matrix_t> {
  using base_type = matrix_bm<matrix_t>;

  matrix_binaryOP_bm() = delete;
  explicit matrix_binaryOP_bm(benchmark_base::configuration cfg)
      : base_type{cfg} {}
  matrix_binaryOP_bm(const matrix_binaryOP_bm& bm) = default;
  matrix_binaryOP_bm& operator=(matrix_binaryOP_bm& other) = default;

  constexpr std::string name() const override {
    return std::string{base_type::name} + "_" + std::string{binaryOP::name};
  }

  inline void operator()(::benchmark::State& state) const override {

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
template <concepts::matrix matrix_t, concepts::vector vector_t>
struct matrix_vector_bm : public matrix_bm<matrix_t> {
  using base_type = matrix_bm<matrix_t>;

  std::vector<vector_t> v;

  matrix_vector_bm() = delete;
  explicit matrix_vector_bm(benchmark_base::configuration cfg)
      : base_type{cfg} {

    v.reserve(this->m_cfg.n_samples());

    fill_random_vec(v);
  }
  matrix_vector_bm(const matrix_vector_bm& bm) = default;
  matrix_vector_bm& operator=(matrix_vector_bm& other) = default;

  /// Clear state
  ~matrix_vector_bm() override { v.clear(); }

  constexpr std::string name() const override {
    return std::string{base_type::name} + "_vector";
  }

  inline void operator()(::benchmark::State& state) const override {

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
  static constexpr std::string_view name{"add"};
  template <concepts::matrix matrix_t>
  constexpr matrix_t operator()(const matrix_t& a, const matrix_t& b) const {
    return a + b;
  }
};
struct sub {
  static constexpr std::string_view name{"sub"};
  template <concepts::matrix matrix_t>
  constexpr matrix_t operator()(const matrix_t& a, const matrix_t& b) const {
    return a - b;
  }
};
struct mul {
  static constexpr std::string_view name{"mul"};
  template <concepts::matrix matrix_t>
  constexpr matrix_t operator()(const matrix_t& a, const matrix_t& b) const {
    return a * b;
  }
};
struct transpose {
  static constexpr std::string_view name{"transpose"};
  template <concepts::matrix matrix_t>
  constexpr auto operator()(const matrix_t& a) const {
    return algebra::matrix::transpose(a);
  }
};
struct determinant {
  static constexpr std::string_view name{"determinant"};
  template <concepts::matrix matrix_t>
  constexpr auto operator()(const matrix_t& a) const {
    return algebra::matrix::determinant(a);
  }
};
struct invert {
  static constexpr std::string_view name{"invert"};
  template <concepts::matrix matrix_t>
  constexpr auto operator()(const matrix_t& a) const {
    return algebra::matrix::inverse(a);
  }
};

}  // namespace bench_op

// Macro for registering all vector benchmarks
#define ALGEBRA_PLUGINS_REGISTER_MATRIX_BENCH(CFG)                   \
  algebra::register_benchmark<mat44_transp_f_t>(CFG, "_4x4_single"); \
  algebra::register_benchmark<mat44_transp_d_t>(CFG, "_4x4_double"); \
  algebra::register_benchmark<mat66_transp_f_t>(CFG, "_6x6_single"); \
  algebra::register_benchmark<mat66_transp_d_t>(CFG, "_6x6_double"); \
  algebra::register_benchmark<mat88_transp_f_t>(CFG, "_8x8_single"); \
  algebra::register_benchmark<mat88_transp_d_t>(CFG, "_8x8_double"); \
                                                                     \
  algebra::register_benchmark<mat44_inv_f_t>(CFG, "_4x4_single");    \
  algebra::register_benchmark<mat44_inv_d_t>(CFG, "_4x4_double");    \
  algebra::register_benchmark<mat66_inv_f_t>(CFG, "_6x6_single");    \
  algebra::register_benchmark<mat66_inv_d_t>(CFG, "_6x6_double");    \
  algebra::register_benchmark<mat88_inv_f_t>(CFG, "_8x8_single");    \
  algebra::register_benchmark<mat88_inv_d_t>(CFG, "_8x8_double");    \
                                                                     \
  algebra::register_benchmark<mat44_det_f_t>(CFG, "_4x4_single");    \
  algebra::register_benchmark<mat44_det_d_t>(CFG, "_4x4_double");    \
  algebra::register_benchmark<mat66_det_f_t>(CFG, "_6x6_single");    \
  algebra::register_benchmark<mat66_det_d_t>(CFG, "_6x6_double");    \
  algebra::register_benchmark<mat88_det_f_t>(CFG, "_8x8_single");    \
  algebra::register_benchmark<mat88_det_d_t>(CFG, "_8x8_double");    \
                                                                     \
  algebra::register_benchmark<mat44_add_f_t>(CFG, "_4x4_single");    \
  algebra::register_benchmark<mat44_add_d_t>(CFG, "_4x4_double");    \
  algebra::register_benchmark<mat66_add_f_t>(CFG, "_6x6_single");    \
  algebra::register_benchmark<mat66_add_d_t>(CFG, "_6x6_double");    \
  algebra::register_benchmark<mat88_add_f_t>(CFG, "_8x8_single");    \
  algebra::register_benchmark<mat88_add_d_t>(CFG, "_8x8_double");    \
                                                                     \
  algebra::register_benchmark<mat44_mul_f_t>(CFG, "_4x4_single");    \
  algebra::register_benchmark<mat44_mul_d_t>(CFG, "_4x4_double");    \
  algebra::register_benchmark<mat66_mul_f_t>(CFG, "_6x6_single");    \
  algebra::register_benchmark<mat66_mul_d_t>(CFG, "_6x6_double");    \
  algebra::register_benchmark<mat88_mul_f_t>(CFG, "_8x8_single");    \
  algebra::register_benchmark<mat88_mul_d_t>(CFG, "_8x8_double");    \
                                                                     \
  algebra::register_benchmark<mat44_vec_f_t>(CFG, "_4x4_single");    \
  algebra::register_benchmark<mat44_vec_d_t>(CFG, "_4x4_double");    \
  algebra::register_benchmark<mat66_vec_f_t>(CFG, "_6x6_single");    \
  algebra::register_benchmark<mat66_vec_d_t>(CFG, "_6x6_double");    \
  algebra::register_benchmark<mat88_vec_f_t>(CFG, "_8x8_single");    \
  algebra::register_benchmark<mat88_vec_d_t>(CFG, "_8x8_double");

}  // namespace algebra
