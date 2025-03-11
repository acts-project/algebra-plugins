/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "benchmark_base.hpp"

// Benchmark include
#include <benchmark/benchmark.h>

// System include(s)
#include <string>

namespace algebra {

template <typename benchmark_t>
requires std::derived_from<benchmark_t, benchmark_base> inline void
register_benchmark(const benchmark_base::configuration& cfg,
                   const std::string& suffix) {
  benchmark_t bench{cfg};

  ::benchmark::RegisterBenchmark((bench.name() + suffix).c_str(), bench)
      ->UseRealTime()
      ->MeasureProcessCPUTime()
      ->ThreadPerCpu();
}

}  // namespace algebra
