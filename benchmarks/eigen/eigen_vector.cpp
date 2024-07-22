/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "algebra/eigen_eigen.hpp"
#include "benchmark/common/benchmark_vector.hpp"
#include "benchmark/common/register_benchmark.hpp"
#include "benchmark/eigen/data_generator.hpp"

// Benchmark include
#include <benchmark/benchmark.h>

// System include(s)
#include <string>

using namespace algebra;

/// Run vector benchmarks
int main(int argc, char** argv) {

  //
  // Prepare benchmarks
  //
  algebra::benchmark_base::configuration cfg{};
  cfg.n_samples(100000);

  using phi_f_t = vector_unaryOP_bm<eigen::vector3, float, bench_op::phi>;
  using theta_f_t = vector_unaryOP_bm<eigen::vector3, float, bench_op::theta>;
  using perp_f_t = vector_unaryOP_bm<eigen::vector3, float, bench_op::perp>;
  using norm_f_t = vector_unaryOP_bm<eigen::vector3, float, bench_op::norm>;
  using eta_f_t = vector_unaryOP_bm<eigen::vector3, float, bench_op::eta>;

  using add_f_t = vector_binaryOP_bm<eigen::vector3, float, bench_op::add>;
  using sub_f_t = vector_binaryOP_bm<eigen::vector3, float, bench_op::sub>;
  using dot_f_t = vector_binaryOP_bm<eigen::vector3, float, bench_op::dot>;
  using cross_f_t = vector_binaryOP_bm<eigen::vector3, float, bench_op::cross>;
  using normlz_f_t =
      vector_unaryOP_bm<eigen::vector3, float, bench_op::normalize>;

  using phi_d_t = vector_unaryOP_bm<eigen::vector3, double, bench_op::phi>;
  using theta_d_t = vector_unaryOP_bm<eigen::vector3, double, bench_op::theta>;
  using perp_d_t = vector_unaryOP_bm<eigen::vector3, double, bench_op::perp>;
  using norm_d_t = vector_unaryOP_bm<eigen::vector3, double, bench_op::norm>;
  using eta_d_t = vector_unaryOP_bm<eigen::vector3, double, bench_op::eta>;

  using add_d_t = vector_binaryOP_bm<eigen::vector3, double, bench_op::add>;
  using sub_d_t = vector_binaryOP_bm<eigen::vector3, double, bench_op::sub>;
  using dot_d_t = vector_binaryOP_bm<eigen::vector3, double, bench_op::dot>;
  using cross_d_t = vector_binaryOP_bm<eigen::vector3, double, bench_op::cross>;
  using normlz_d_t =
      vector_unaryOP_bm<eigen::vector3, double, bench_op::normalize>;

  std::cout << "------------------------------------------\n"
            << "Algebra-Plugins 'vector' benchmark (Eigen)\n"
            << "------------------------------------------\n\n"
            << cfg;

  //
  // Register all benchmarks
  //
  algebra::register_benchmark<phi_f_t>(cfg, "_single");
  algebra::register_benchmark<phi_d_t>(cfg, "_double");
  algebra::register_benchmark<theta_f_t>(cfg, "_single");
  algebra::register_benchmark<theta_d_t>(cfg, "_double");
  algebra::register_benchmark<perp_f_t>(cfg, "_single");
  algebra::register_benchmark<perp_d_t>(cfg, "_double");
  algebra::register_benchmark<norm_f_t>(cfg, "_single");
  algebra::register_benchmark<norm_d_t>(cfg, "_double");
  algebra::register_benchmark<eta_f_t>(cfg, "_single");
  algebra::register_benchmark<eta_d_t>(cfg, "_double");

  algebra::register_benchmark<add_f_t>(cfg, "_single");
  algebra::register_benchmark<add_d_t>(cfg, "_double");
  algebra::register_benchmark<sub_f_t>(cfg, "_single");
  algebra::register_benchmark<sub_d_t>(cfg, "_double");
  algebra::register_benchmark<dot_f_t>(cfg, "_single");
  algebra::register_benchmark<dot_d_t>(cfg, "_double");
  algebra::register_benchmark<cross_f_t>(cfg, "_single");
  algebra::register_benchmark<cross_d_t>(cfg, "_double");
  algebra::register_benchmark<normlz_f_t>(cfg, "_single");
  algebra::register_benchmark<normlz_d_t>(cfg, "_double");

  ::benchmark::Initialize(&argc, argv);
  ::benchmark::RunSpecifiedBenchmarks();
  ::benchmark::Shutdown();
}
