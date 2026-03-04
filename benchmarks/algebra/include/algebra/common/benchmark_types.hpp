/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Algebra plugin selection
#if ALGEBRA_BENCHMARK_ARRAY
#include "algebra/array.hpp"
#include "algebra/array/data_generator.hpp"
#endif
#if ALGEBRA_BENCHMARK_EIGEN
#include "algebra/eigen.hpp"
#include "algebra/eigen/data_generator.hpp"
#endif
#if ALGEBRA_BENCHMARK_FASTOR
#include "algebra/fastor.hpp"
#include "algebra/fastor/data_generator.hpp"
#endif
#if ALGEBRA_BENCHMARK_SMATRIX
#include "algebra/smatrix.hpp"
#include "algebra/smatrix/data_generator.hpp"
#endif
#if ALGEBRA_BENCHMARK_VC_AOS
#include "algebra/vc_aos.hpp"
#include "algebra/vc_aos/data_generator.hpp"
#endif

// Project include(s).
#include "algebra/common/benchmark_vector.hpp"

// Benchmark include
#include <benchmark/benchmark.h>

// System include(s)
#include <iostream>
#include <string>

namespace algebra::benchmark {

#define ALGEBRA_DEFINE_BENCHMARK_NAME(ALGEBRA) \
    static std::string plugin_name(#ALGEBRA);

// Select algebra-plugin to compile the test with
#if ALGEBRA_BENCHMARK_ARRAY
ALGEBRA_DEFINE_BENCHMARK_NAME(array)
#elif ALGEBRA_BENCHMARK_EIGEN
ALGEBRA_DEFINE_BENCHMARK_NAME(eigen)
#elif ALGEBRA_BENCHMARK_FASTOR
ALGEBRA_DEFINE_BENCHMARK_NAME(fastor)
#elif ALGEBRA_BENCHMARK_SMATRIX
ALGEBRA_DEFINE_BENCHMARK_NAME(smatrix)
#elif ALGEBRA_BENCHMARK_VC_AOS
ALGEBRA_DEFINE_BENCHMARK_NAME(vc_aos)
#elif ALGEBRA_BENCHMARK_VC_SOA
ALGEBRA_DEFINE_BENCHMARK_NAME(vc_soa)
#endif

}  // namespace algebra::benchmark
