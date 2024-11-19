/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "benchmark_vector.hpp"
#include "register_benchmark.hpp"

namespace algebra::bench_op {

// Macro for declaring the predefined materials (with Density effect data)
#define ALGEBRA_PLUGINS_BENCH_GETTER(GETTER_NAME)       \
  struct GETTER_NAME {                                  \
    inline static const std::string name{#GETTER_NAME}; \
    template <typename vector_t>                        \
    auto operator()(const vector_t &a) const {          \
      return algebra::getter::GETTER_NAME(a);           \
    }                                                   \
  };

// Functions to be benchmarked
ALGEBRA_PLUGINS_BENCH_GETTER(phi)
ALGEBRA_PLUGINS_BENCH_GETTER(theta)
ALGEBRA_PLUGINS_BENCH_GETTER(perp)
ALGEBRA_PLUGINS_BENCH_GETTER(norm)
ALGEBRA_PLUGINS_BENCH_GETTER(eta)

// Macro for registering all getter benchmarks
#define ALGEBRA_PLUGINS_REGISTER_GETTER_BENCH(CFG)        \
  algebra::register_benchmark<phi_f_t>(CFG, "_single");   \
  algebra::register_benchmark<phi_d_t>(CFG, "_double");   \
  algebra::register_benchmark<theta_f_t>(CFG, "_single"); \
  algebra::register_benchmark<theta_d_t>(CFG, "_double"); \
  algebra::register_benchmark<perp_f_t>(CFG, "_single");  \
  algebra::register_benchmark<perp_d_t>(CFG, "_double");  \
  algebra::register_benchmark<norm_f_t>(CFG, "_single");  \
  algebra::register_benchmark<norm_d_t>(CFG, "_double");  \
  algebra::register_benchmark<eta_f_t>(CFG, "_single");   \
  algebra::register_benchmark<eta_d_t>(CFG, "_double");

}  // namespace algebra::bench_op
