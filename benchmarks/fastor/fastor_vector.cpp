/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "algebra/fastor_fastor.hpp"
#include "benchmark/common/benchmark_vector.hpp"
#include "benchmark/fastor/data_generator.hpp"

// Benchmark include
#include <benchmark/benchmark.h>

// System include(s)
#include <iostream>

using namespace algebra;

/// Run vector benchmarks
int main(int argc, char** argv) {

  //
  // Prepare benchmarks
  //
  algebra::benchmark_base::configuration cfg{};
  cfg.n_samples(100000);

  using phi_f_t = vector_unaryOP_bm<fastor::vector3, float, bench_op::phi>;
  using theta_f_t = vector_unaryOP_bm<fastor::vector3, float, bench_op::theta>;
  using perp_f_t = vector_unaryOP_bm<fastor::vector3, float, bench_op::perp>;
  using norm_f_t = vector_unaryOP_bm<fastor::vector3, float, bench_op::norm>;
  using eta_f_t = vector_unaryOP_bm<fastor::vector3, float, bench_op::eta>;

  using add_f_t = vector_binaryOP_bm<fastor::vector3, float, bench_op::add>;
  using sub_f_t = vector_binaryOP_bm<fastor::vector3, float, bench_op::sub>;
  using dot_f_t = vector_binaryOP_bm<fastor::vector3, float, bench_op::dot>;
  using cross_f_t = vector_binaryOP_bm<fastor::vector3, float, bench_op::cross>;
  using normlz_f_t =
      vector_unaryOP_bm<fastor::vector3, float, bench_op::normalize>;

  using phi_d_t = vector_unaryOP_bm<fastor::vector3, double, bench_op::phi>;
  using theta_d_t = vector_unaryOP_bm<fastor::vector3, double, bench_op::theta>;
  using perp_d_t = vector_unaryOP_bm<fastor::vector3, double, bench_op::perp>;
  using norm_d_t = vector_unaryOP_bm<fastor::vector3, double, bench_op::norm>;
  using eta_d_t = vector_unaryOP_bm<fastor::vector3, double, bench_op::eta>;

  using add_d_t = vector_binaryOP_bm<fastor::vector3, double, bench_op::add>;
  using sub_d_t = vector_binaryOP_bm<fastor::vector3, double, bench_op::sub>;
  using dot_d_t = vector_binaryOP_bm<fastor::vector3, double, bench_op::dot>;
  using cross_d_t =
      vector_binaryOP_bm<fastor::vector3, double, bench_op::cross>;
  using normlz_d_t =
      vector_unaryOP_bm<fastor::vector3, double, bench_op::normalize>;

  std::cout << "-------------------------------------------\n"
            << "Algebra-Plugins 'vector' benchmark (Fastor)\n"
            << "-------------------------------------------\n\n"
            << cfg;

  //
  // Register all benchmarks
  //
  ALGEBRA_PLUGINS_REGISTER_VECTOR_BENCH(cfg)

  ::benchmark::Initialize(&argc, argv);
  ::benchmark::RunSpecifiedBenchmarks();
  ::benchmark::Shutdown();
}
