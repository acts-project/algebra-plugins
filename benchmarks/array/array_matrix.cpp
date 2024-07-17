/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "algebra/array_cmath.hpp"
#include "benchmark/array/data_generator.hpp"
#include "benchmark/common/benchmark_matrix.hpp"
#include "benchmark/common/register_benchmark.hpp"

// Benchmark include
#include <benchmark/benchmark.h>

using namespace algebra;

/// Run vector benchmarks
int main(int argc, char** argv) {

  //
  // Prepare benchmarks
  //
  algebra::benchmark_base::configuration cfg{};
  cfg.n_samples(100000);

  using mat44_add_f_t =
      matrix_binaryOP_bm<array::matrix_type<float, 4, 4>, bench_op::add>;
  using mat44_add_d_t =
      matrix_binaryOP_bm<array::matrix_type<double, 4, 4>, bench_op::add>;
  using mat66_add_f_t =
      matrix_binaryOP_bm<array::matrix_type<float, 6, 6>, bench_op::add>;
  using mat66_add_d_t =
      matrix_binaryOP_bm<array::matrix_type<double, 6, 6>, bench_op::add>;
  using mat88_add_f_t =
      matrix_binaryOP_bm<array::matrix_type<float, 8, 8>, bench_op::add>;
  using mat88_add_d_t =
      matrix_binaryOP_bm<array::matrix_type<double, 8, 8>, bench_op::add>;

  using mat44_mul_f_t =
      matrix_binaryOP_bm<array::matrix_type<float, 4, 4>, bench_op::mul>;
  using mat44_mul_d_t =
      matrix_binaryOP_bm<array::matrix_type<double, 4, 4>, bench_op::mul>;
  using mat66_mul_f_t =
      matrix_binaryOP_bm<array::matrix_type<float, 6, 6>, bench_op::mul>;
  using mat66_mul_d_t =
      matrix_binaryOP_bm<array::matrix_type<double, 6, 6>, bench_op::mul>;
  using mat88_mul_f_t =
      matrix_binaryOP_bm<array::matrix_type<float, 8, 8>, bench_op::mul>;
  using mat88_mul_d_t =
      matrix_binaryOP_bm<array::matrix_type<double, 8, 8>, bench_op::mul>;

  using mat44_vec_f_t = matrix_vector_bm<array::matrix_type<float, 4, 4>,
                                         array::vector_type<float, 4>>;
  using mat44_vec_d_t = matrix_vector_bm<array::matrix_type<double, 4, 4>,
                                         array::vector_type<double, 4>>;
  using mat66_vec_f_t = matrix_vector_bm<array::matrix_type<float, 6, 6>,
                                         array::vector_type<float, 6>>;
  using mat66_vec_d_t = matrix_vector_bm<array::matrix_type<double, 6, 6>,
                                         array::vector_type<double, 6>>;
  using mat88_vec_f_t = matrix_vector_bm<array::matrix_type<float, 8, 8>,
                                         array::vector_type<float, 8>>;
  using mat88_vec_d_t = matrix_vector_bm<array::matrix_type<double, 8, 8>,
                                         array::vector_type<double, 8>>;

  std::cout << "-----------------------------------------------\n"
            << "Algebra-Plugins 'matrix' benchmark (std::array)\n"
            << "-----------------------------------------------\n\n"
            << cfg;

  //
  // Register all benchmarks
  //
  algebra::register_benchmark<mat44_add_f_t>(cfg, "_4x4_add_single");
  algebra::register_benchmark<mat44_add_d_t>(cfg, "_4x4_add_double");
  algebra::register_benchmark<mat66_add_f_t>(cfg, "_6x6_add_single");
  algebra::register_benchmark<mat66_add_d_t>(cfg, "_6x6_add_double");
  algebra::register_benchmark<mat88_add_f_t>(cfg, "_8x8_add_single");
  algebra::register_benchmark<mat88_add_d_t>(cfg, "_8x8_add_double");

  algebra::register_benchmark<mat44_mul_f_t>(cfg, "_4x4_mul_single");
  algebra::register_benchmark<mat44_mul_d_t>(cfg, "_4x4_mul_double");
  algebra::register_benchmark<mat66_mul_f_t>(cfg, "_6x6_mul_single");
  algebra::register_benchmark<mat66_mul_d_t>(cfg, "_6x6_mul_double");
  algebra::register_benchmark<mat88_mul_f_t>(cfg, "_8x8_mul_single");
  algebra::register_benchmark<mat88_mul_d_t>(cfg, "_8x8_mul_double");

  algebra::register_benchmark<mat44_vec_f_t>(cfg, "_4x4_vec_single");
  algebra::register_benchmark<mat44_vec_d_t>(cfg, "_4x4_vec_double");
  algebra::register_benchmark<mat66_vec_f_t>(cfg, "_6x6_vec_single");
  algebra::register_benchmark<mat66_vec_d_t>(cfg, "_6x6_vec_double");
  algebra::register_benchmark<mat88_vec_f_t>(cfg, "_8x8_vec_single");
  algebra::register_benchmark<mat88_vec_d_t>(cfg, "_8x8_vec_double");

  ::benchmark::Initialize(&argc, argv);
  ::benchmark::RunSpecifiedBenchmarks();
  ::benchmark::Shutdown();
}
