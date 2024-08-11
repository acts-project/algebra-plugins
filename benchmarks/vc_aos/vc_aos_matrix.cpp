/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "algebra/vc_aos.hpp"
#include "benchmark/common/benchmark_matrix.hpp"
#include "benchmark/vc_aos/data_generator.hpp"

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

  using mat44_transp_f_t =
      matrix_unaryOP_bm<vc_aos::matrix_type<float, 4, 4>, bench_op::transpose>;
  using mat44_transp_d_t =
      matrix_unaryOP_bm<vc_aos::matrix_type<double, 4, 4>, bench_op::transpose>;
  using mat66_transp_f_t =
      matrix_unaryOP_bm<vc_aos::matrix_type<float, 6, 6>, bench_op::transpose>;
  using mat66_transp_d_t =
      matrix_unaryOP_bm<vc_aos::matrix_type<double, 6, 6>, bench_op::transpose>;
  using mat88_transp_f_t =
      matrix_unaryOP_bm<vc_aos::matrix_type<float, 8, 8>, bench_op::transpose>;
  using mat88_transp_d_t =
      matrix_unaryOP_bm<vc_aos::matrix_type<double, 8, 8>, bench_op::transpose>;

  using mat44_add_f_t =
      matrix_binaryOP_bm<vc_aos::matrix_type<float, 4, 4>, bench_op::add>;
  using mat44_add_d_t =
      matrix_binaryOP_bm<vc_aos::matrix_type<double, 4, 4>, bench_op::add>;
  using mat66_add_f_t =
      matrix_binaryOP_bm<vc_aos::matrix_type<float, 6, 6>, bench_op::add>;
  using mat66_add_d_t =
      matrix_binaryOP_bm<vc_aos::matrix_type<double, 6, 6>, bench_op::add>;
  using mat88_add_f_t =
      matrix_binaryOP_bm<vc_aos::matrix_type<float, 8, 8>, bench_op::add>;
  using mat88_add_d_t =
      matrix_binaryOP_bm<vc_aos::matrix_type<double, 8, 8>, bench_op::add>;

  using mat44_mul_f_t =
      matrix_binaryOP_bm<vc_aos::matrix_type<float, 4, 4>, bench_op::mul>;
  using mat44_mul_d_t =
      matrix_binaryOP_bm<vc_aos::matrix_type<double, 4, 4>, bench_op::mul>;
  using mat66_mul_f_t =
      matrix_binaryOP_bm<vc_aos::matrix_type<float, 6, 6>, bench_op::mul>;
  using mat66_mul_d_t =
      matrix_binaryOP_bm<vc_aos::matrix_type<double, 6, 6>, bench_op::mul>;
  using mat88_mul_f_t =
      matrix_binaryOP_bm<vc_aos::matrix_type<float, 8, 8>, bench_op::mul>;
  using mat88_mul_d_t =
      matrix_binaryOP_bm<vc_aos::matrix_type<double, 8, 8>, bench_op::mul>;

  using mat44_vec_f_t = matrix_vector_bm<vc_aos::matrix_type<float, 4, 4>,
                                         vc_aos::vector_type<float, 4>>;
  using mat44_vec_d_t = matrix_vector_bm<vc_aos::matrix_type<double, 4, 4>,
                                         vc_aos::vector_type<double, 4>>;
  using mat66_vec_f_t = matrix_vector_bm<vc_aos::matrix_type<float, 6, 6>,
                                         vc_aos::vector_type<float, 6>>;
  using mat66_vec_d_t = matrix_vector_bm<vc_aos::matrix_type<double, 6, 6>,
                                         vc_aos::vector_type<double, 6>>;
  using mat88_vec_f_t = matrix_vector_bm<vc_aos::matrix_type<float, 8, 8>,
                                         vc_aos::vector_type<float, 8>>;
  using mat88_vec_d_t = matrix_vector_bm<vc_aos::matrix_type<double, 8, 8>,
                                         vc_aos::vector_type<double, 8>>;

  std::cout << "------------------------------------------\n"
            << "Algebra-Plugins 'matrix' benchmark (Vc AoS)\n"
            << "-------------------------------------------\n\n"
            << cfg;

  //
  // Register all benchmarks
  //
  ALGEBRA_PLUGINS_REGISTER_MATRIX_BENCH(cfg)

  ::benchmark::Initialize(&argc, argv);
  ::benchmark::RunSpecifiedBenchmarks();
  ::benchmark::Shutdown();
}