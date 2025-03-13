/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "algebra/eigen_eigen.hpp"
#include "benchmark/common/benchmark_getter.hpp"
#include "benchmark/eigen/data_generator.hpp"

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

  std::cout << "------------------------------------------\n"
            << "Algebra-Plugins 'getter' benchmark (Eigen)\n"
            << "------------------------------------------\n\n"
            << cfg;

  //
  // Register all benchmarks
  //
  ::benchmark::Initialize(&argc, argv);
  ::benchmark::RunSpecifiedBenchmarks();
  ::benchmark::Shutdown();
}
