/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "benchmark_vector.hpp"

// System include(s)
#include <chrono>
#include <iostream>
#include <string_view>
#include <thread>
#include <vector>

namespace algebra {

template <concepts::transform3D transform3_t>
void fill_random_trf(std::vector<transform3_t>&);

/// Benchmark for vector operations
template <concepts::transform3D transform3_t>
struct transform3_bm : public vector_bm<typename transform3_t::vector3> {
 private:
  using base_type = vector_bm<typename transform3_t::vector3>;

 public:
  /// Prefix for the benchmark name
  static constexpr std::string_view bm_name{"transform3"};

  std::vector<transform3_t> trfs;

  /// No default construction: Cannot prepare data
  transform3_bm() = delete;
  /// Construct from an externally provided configuration @param cfg
  explicit transform3_bm(benchmark_base::configuration cfg) : base_type{cfg} {

    trfs.reserve(this->m_cfg.n_samples());

    fill_random_trf(trfs);
  }
  transform3_bm(const transform3_bm& bm) = default;
  transform3_bm& operator=(transform3_bm& other) = default;

  /// Clear state
  ~transform3_bm() override { trfs.clear(); }

  constexpr std::string name() const override {
    return std::string{base_type::name} + "_" + std::string{bm_name};
  }

  /// Benchmark case
  inline void operator()(::benchmark::State& state) const override {

    using vector_t = typename transform3_t::vector3;
    using point_t = typename transform3_t::point3;

    const std::size_t n_samples{this->m_cfg.n_samples()};

    // Run the benchmark
    for (auto _ : state) {
      for (std::size_t i{0}; i < n_samples; ++i) {

        point_t result1 = this->trfs[i].point_to_global(this->a[i]);
        point_t result2 = this->trfs[i].point_to_local(this->a[i]);
        vector_t result3 = this->trfs[i].vector_to_global(this->a[i]);
        vector_t result4 = this->trfs[i].vector_to_local(this->a[i]);

        ::benchmark::DoNotOptimize(result1);
        ::benchmark::DoNotOptimize(result2);
        ::benchmark::DoNotOptimize(result3);
        ::benchmark::DoNotOptimize(result4);
      }
    }
  }
};

}  // namespace algebra
