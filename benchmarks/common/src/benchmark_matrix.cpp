/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
// clang-format off
#include "benchmark/common/benchmark_types.hpp"
#include "benchmark/common/benchmark_matrix.hpp"
// clang-format on

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

  std::cout << "-----------------------------------------------\n"
            << "Algebra-Plugins 'matrix' benchmark ("
            << algebra::benchmark::plugin_name << ")\n"
            << "-----------------------------------------------\n\n"
            << cfg;

//
// Define and register all benchmarks
//
#if ALGEBRA_BENCHMARK_ARRAY
  ALGEBRA_PLUGINS_DEFINE_MATRIX_BENCH(array)
#elif ALGEBRA_BENCHMARK_EIGEN
  ALGEBRA_PLUGINS_DEFINE_MATRIX_BENCH(eigen)
#elif ALGEBRA_BENCHMARK_FASTOR
  ALGEBRA_PLUGINS_DEFINE_MATRIX_BENCH(fastor)
#elif ALGEBRA_BENCHMARK_SMATRIX
  ALGEBRA_PLUGINS_DEFINE_MATRIX_BENCH(smatrix)
#elif ALGEBRA_BENCHMARK_VC_AOS
  ALGEBRA_PLUGINS_DEFINE_MATRIX_BENCH(vc_aos)
#endif

  ALGEBRA_PLUGINS_REGISTER_MATRIX_BENCH(cfg, cfg)

  ::benchmark::Initialize(&argc, argv);
  ::benchmark::RunSpecifiedBenchmarks();
  ::benchmark::Shutdown();
}
