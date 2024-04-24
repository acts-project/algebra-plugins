/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "algebra/vc_soa.hpp"
#include "benchmark/common/benchmark_getter.hpp"
#include "benchmark/common/register_benchmark.hpp"
#include "benchmark/vc_soa/data_generator.hpp"

// Benchmark include
#include <benchmark/benchmark.h>

// System include(s)
#include <string>

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

  using phi_f_t = vector_unaryOP_bm<vc_soa::vector3, float, bench_op::phi>;
  using theta_f_t = vector_unaryOP_bm<vc_soa::vector3, float, bench_op::theta>;
  using perp_f_t = vector_unaryOP_bm<vc_soa::vector3, float, bench_op::perp>;
  using norm_f_t = vector_unaryOP_bm<vc_soa::vector3, float, bench_op::norm>;
  using eta_f_t = vector_unaryOP_bm<vc_soa::vector3, float, bench_op::eta>;

  using phi_d_t = vector_unaryOP_bm<vc_soa::vector3, double, bench_op::phi>;
  using theta_d_t = vector_unaryOP_bm<vc_soa::vector3, double, bench_op::theta>;
  using perp_d_t = vector_unaryOP_bm<vc_soa::vector3, double, bench_op::perp>;
  using norm_d_t = vector_unaryOP_bm<vc_soa::vector3, double, bench_op::norm>;
  using eta_d_t = vector_unaryOP_bm<vc_soa::vector3, double, bench_op::eta>;

  std::cout << "-------------------------------------------\n"
            << "Algebra-Plugins 'getter' benchmark (Vc SoA)\n"
            << "-------------------------------------------\n\n"
            << "(single)\n"
            << cfg_s << "(double)\n"
            << cfg_d;

  //
  // Register all benchmarks
  //
  algebra::register_benchmark<phi_f_t>(cfg_s, "_single");
  algebra::register_benchmark<phi_d_t>(cfg_d, "_double");
  algebra::register_benchmark<theta_f_t>(cfg_s, "_single");
  algebra::register_benchmark<theta_d_t>(cfg_d, "_double");
  algebra::register_benchmark<perp_f_t>(cfg_s, "_single");
  algebra::register_benchmark<perp_d_t>(cfg_d, "_double");
  algebra::register_benchmark<norm_f_t>(cfg_s, "_single");
  algebra::register_benchmark<norm_d_t>(cfg_d, "_double");
  algebra::register_benchmark<eta_f_t>(cfg_s, "_single");
  algebra::register_benchmark<eta_d_t>(cfg_d, "_double");

  ::benchmark::Initialize(&argc, argv);
  ::benchmark::RunSpecifiedBenchmarks();
  ::benchmark::Shutdown();
}
