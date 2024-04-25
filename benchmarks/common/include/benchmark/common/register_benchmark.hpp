/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Benchmark include
#include <benchmark/benchmark.h>

// System include(s)
#include <string>

namespace algebra {

template <typename benchmark_t, typename config_t>
void register_benchmark(const config_t& cfg, const std::string& suffix) {
  benchmark_t bench{cfg};

  ::benchmark::RegisterBenchmark((bench.name() + suffix).c_str(), bench)
      ->UseRealTime()
      ->MeasureProcessCPUTime()
      ->ThreadPerCpu();
}

}  // namespace algebra
