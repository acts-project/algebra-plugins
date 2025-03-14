/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "algebra/array_cmath.hpp"
#include "benchmark/array/data_generator.hpp"
#include "benchmark/common/benchmark_transform3.hpp"
#include "benchmark/common/register_benchmark.hpp"

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

  using trf_f_t = transform3_bm<array::transform3<float>>;
  using trf_d_t = transform3_bm<array::transform3<double>>;

  std::cout << "---------------------------------------------------\n"
            << "Algebra-Plugins 'transform3' benchmark (std::array)\n"
            << "---------------------------------------------------\n\n"
            << cfg;

  //
  // Register all benchmarks
  //
  algebra::register_benchmark<trf_f_t>(cfg, "_single");
  algebra::register_benchmark<trf_d_t>(cfg, "_double");

  ::benchmark::Initialize(&argc, argv);
  ::benchmark::RunSpecifiedBenchmarks();
  ::benchmark::Shutdown();
}
