/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "algebra/vc_soa.hpp"
#include "benchmark/common/benchmark_transform3.hpp"
#include "benchmark/common/register_benchmark.hpp"
#include "benchmark/vc_soa/data_generator.hpp"

// Benchmark include
#include <benchmark/benchmark.h>

// System include(s)
#include <iostream>

using namespace algebra;

/// Run vector benchmarks
int main(int argc, char** argv) {

  constexpr std::size_t n_samples{100000};

  //
  // Prepare benchmarks
  //
  algebra::benchmark_base::configuration cfg_s{};
  // Reduce the number of samples, since a single SoA struct contains multiple
  // vectors
  cfg_s.n_samples(n_samples / Vc::float_v::Size);

  // For double precision we need more samples (less vectors per SoA)
  algebra::benchmark_base::configuration cfg_d{cfg_s};
  cfg_d.n_samples(n_samples / Vc::double_v::Size);

  using trf_f_t = transform3_bm<vc_soa::transform3<float>>;
  using trf_d_t = transform3_bm<vc_soa::transform3<double>>;

  std::cout << "-------------------------------------------\n"
            << "Algebra-Plugins 'transform3' benchmark (Vc SoA)\n"
            << "-----------------------------------------------\n\n"
            << "(single)\n"
            << cfg_s << "(double)\n"
            << cfg_d;

  //
  // Register all benchmarks
  //
  algebra::register_benchmark<trf_f_t>(cfg_s, "_single");
  algebra::register_benchmark<trf_d_t>(cfg_d, "_double");

  ::benchmark::Initialize(&argc, argv);
  ::benchmark::RunSpecifiedBenchmarks();
  ::benchmark::Shutdown();
}
