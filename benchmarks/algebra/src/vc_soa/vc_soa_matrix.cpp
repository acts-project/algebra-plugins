/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
// clang-format off
#include "algebra/vc_soa.hpp"
#include "algebra/common/benchmark_matrix.hpp"
// clang-format on
#include "algebra/common/register_benchmark.hpp"
#include "algebra/vc_soa/data_generator.hpp"

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

  using mat44_transp_f_t = matrix_unaryOP_bm<vc_soa::matrix_type<float, 4, 4>,
                                             bench_op::matrix::transpose>;
  using mat44_transp_d_t = matrix_unaryOP_bm<vc_soa::matrix_type<double, 4, 4>,
                                             bench_op::matrix::transpose>;
  using mat66_transp_f_t = matrix_unaryOP_bm<vc_soa::matrix_type<float, 6, 6>,
                                             bench_op::matrix::transpose>;
  using mat66_transp_d_t = matrix_unaryOP_bm<vc_soa::matrix_type<double, 6, 6>,
                                             bench_op::matrix::transpose>;
  using mat88_transp_f_t = matrix_unaryOP_bm<vc_soa::matrix_type<float, 8, 8>,
                                             bench_op::matrix::transpose>;
  using mat88_transp_d_t = matrix_unaryOP_bm<vc_soa::matrix_type<double, 8, 8>,
                                             bench_op::matrix::transpose>;

  using mat44_add_f_t = matrix_binaryOP_bm<vc_soa::matrix_type<float, 4, 4>,
                                           bench_op::matrix::add>;
  using mat44_add_d_t = matrix_binaryOP_bm<vc_soa::matrix_type<double, 4, 4>,
                                           bench_op::matrix::add>;
  using mat66_add_f_t = matrix_binaryOP_bm<vc_soa::matrix_type<float, 6, 6>,
                                           bench_op::matrix::add>;
  using mat66_add_d_t = matrix_binaryOP_bm<vc_soa::matrix_type<double, 6, 6>,
                                           bench_op::matrix::add>;
  using mat88_add_f_t = matrix_binaryOP_bm<vc_soa::matrix_type<float, 8, 8>,
                                           bench_op::matrix::add>;
  using mat88_add_d_t = matrix_binaryOP_bm<vc_soa::matrix_type<double, 8, 8>,
                                           bench_op::matrix::add>;

  using mat44_mul_f_t = matrix_binaryOP_bm<vc_soa::matrix_type<float, 4, 4>,
                                           bench_op::matrix::mul>;
  using mat44_mul_d_t = matrix_binaryOP_bm<vc_soa::matrix_type<double, 4, 4>,
                                           bench_op::matrix::mul>;
  using mat66_mul_f_t = matrix_binaryOP_bm<vc_soa::matrix_type<float, 6, 6>,
                                           bench_op::matrix::mul>;
  using mat66_mul_d_t = matrix_binaryOP_bm<vc_soa::matrix_type<double, 6, 6>,
                                           bench_op::matrix::mul>;
  using mat88_mul_f_t = matrix_binaryOP_bm<vc_soa::matrix_type<float, 8, 8>,
                                           bench_op::matrix::mul>;
  using mat88_mul_d_t = matrix_binaryOP_bm<vc_soa::matrix_type<double, 8, 8>,
                                           bench_op::matrix::mul>;

  using mat44_vec_f_t = matrix_vector_bm<vc_soa::matrix_type<float, 4, 4>,
                                         vc_soa::vector_type<float, 4>>;
  using mat44_vec_d_t = matrix_vector_bm<vc_soa::matrix_type<double, 4, 4>,
                                         vc_soa::vector_type<double, 4>>;
  using mat66_vec_f_t = matrix_vector_bm<vc_soa::matrix_type<float, 6, 6>,
                                         vc_soa::vector_type<float, 6>>;
  using mat66_vec_d_t = matrix_vector_bm<vc_soa::matrix_type<double, 6, 6>,
                                         vc_soa::vector_type<double, 6>>;
  using mat88_vec_f_t = matrix_vector_bm<vc_soa::matrix_type<float, 8, 8>,
                                         vc_soa::vector_type<float, 8>>;
  using mat88_vec_d_t = matrix_vector_bm<vc_soa::matrix_type<double, 8, 8>,
                                         vc_soa::vector_type<double, 8>>;

  std::cout << "-------------------------------------------\n"
            << "Algebra-Plugins 'matrix' benchmark (Vc SoA)\n"
            << "-------------------------------------------\n\n"
            << "(single)\n"
            << cfg_s << "(double)\n"
            << cfg_d;

  //
  // Register all benchmarks
  //
  algebra::register_benchmark<mat44_transp_f_t>(cfg_s, "_4x4_single");
  algebra::register_benchmark<mat44_transp_d_t>(cfg_d, "_4x4_double");
  algebra::register_benchmark<mat66_transp_f_t>(cfg_s, "_6x6_single");
  algebra::register_benchmark<mat66_transp_d_t>(cfg_d, "_6x6_double");
  algebra::register_benchmark<mat88_transp_f_t>(cfg_s, "_8x8_single");
  algebra::register_benchmark<mat88_transp_d_t>(cfg_d, "_8x8_double");

  algebra::register_benchmark<mat44_add_f_t>(cfg_s, "_4x4_single");
  algebra::register_benchmark<mat44_add_d_t>(cfg_d, "_4x4_double");
  algebra::register_benchmark<mat66_add_f_t>(cfg_s, "_6x6_single");
  algebra::register_benchmark<mat66_add_d_t>(cfg_d, "_6x6_double");
  algebra::register_benchmark<mat88_add_f_t>(cfg_s, "_8x8_single");
  algebra::register_benchmark<mat88_add_d_t>(cfg_d, "_8x8_double");

  algebra::register_benchmark<mat44_mul_f_t>(cfg_s, "_4x4_single");
  algebra::register_benchmark<mat44_mul_d_t>(cfg_d, "_4x4_double");
  algebra::register_benchmark<mat66_mul_f_t>(cfg_s, "_6x6_single");
  algebra::register_benchmark<mat66_mul_d_t>(cfg_d, "_6x6_double");
  algebra::register_benchmark<mat88_mul_f_t>(cfg_s, "_8x8_single");
  algebra::register_benchmark<mat88_mul_d_t>(cfg_d, "_8x8_double");

  algebra::register_benchmark<mat44_vec_f_t>(cfg_s, "_4x4_single");
  algebra::register_benchmark<mat44_vec_d_t>(cfg_d, "_4x4_double");
  algebra::register_benchmark<mat66_vec_f_t>(cfg_s, "_6x6_single");
  algebra::register_benchmark<mat66_vec_d_t>(cfg_d, "_6x6_double");
  algebra::register_benchmark<mat88_vec_f_t>(cfg_s, "_8x8_single");
  algebra::register_benchmark<mat88_vec_d_t>(cfg_d, "_8x8_double");

  ::benchmark::Initialize(&argc, argv);
  ::benchmark::RunSpecifiedBenchmarks();
  ::benchmark::Shutdown();
}
