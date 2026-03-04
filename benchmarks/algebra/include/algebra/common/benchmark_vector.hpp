/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "algebra/common/benchmark_base.hpp"
#include "algebra/common/register_benchmark.hpp"
#include "algebra/concepts.hpp"

// System include(s)
#include <chrono>
#include <iostream>
#include <string_view>
#include <thread>
#include <vector>

namespace algebra {

template <concepts::vector vector_t>
void fill_random_vec(std::vector<vector_t> &);

/// Benchmark for vector operations
template <concepts::vector vector_t>
struct vector_bm : public benchmark_base {

    /// Prefix for the benchmark name
    static constexpr std::string_view name{"vector"};

    std::vector<vector_t> a;
    std::vector<vector_t> b;
    std::vector<vector_t> results;

    /// No default construction: Cannot prepare data
    vector_bm() = delete;

    /// Construct from an externally provided configuration @param cfg
    explicit vector_bm(benchmark_base::configuration cfg)
        : benchmark_base{cfg} {

        const std::size_t n_data{this->m_cfg.n_samples()};

        a.reserve(n_data);
        b.reserve(n_data);

        fill_random_vec(a);
        fill_random_vec(b);
    }
    vector_bm(const vector_bm &bm) = default;
    vector_bm &operator=(vector_bm &other) = default;

    /// Clear state
    ~vector_bm() override {
        a.clear();
        b.clear();
    }
};

/// Benchmark elementwise addition of vectors
template <template <typename> class vector_t, concepts::scalar scalar_t,
          typename unaryOP>
    requires std::invocable<unaryOP, vector_t<scalar_t>>
struct vector_unaryOP_bm : public vector_bm<vector_t<scalar_t>> {
    using base_type = vector_bm<vector_t<scalar_t>>;

    vector_unaryOP_bm() = delete;
    explicit vector_unaryOP_bm(benchmark_base::configuration cfg)
        : base_type{cfg} {}
    vector_unaryOP_bm(const vector_unaryOP_bm &bm) = default;
    vector_unaryOP_bm &operator=(vector_unaryOP_bm &other) = default;

    constexpr std::string name() const override {
        return std::string{base_type::name} + "_" + std::string{unaryOP::name};
    }

    inline void operator()(::benchmark::State &state) const override {

        using result_t = std::invoke_result_t<unaryOP, vector_t<scalar_t>>;

        const std::size_t n_samples{this->m_cfg.n_samples()};

        // Run the benchmark
        for (auto _ : state) {
            for (std::size_t i{0}; i < n_samples; ++i) {
                result_t result = unaryOP{}(this->a[i]);
                ::benchmark::DoNotOptimize(result);
            }
        }
    }
};

/// Benchmark elementwise addition of vectors
template <template <typename> class vector_t, concepts::scalar scalar_t,
          typename binaryOP>
    requires std::invocable<binaryOP, vector_t<scalar_t>, vector_t<scalar_t>>
struct vector_binaryOP_bm : public vector_bm<vector_t<scalar_t>> {
    using base_type = vector_bm<vector_t<scalar_t>>;

    vector_binaryOP_bm() = delete;
    explicit vector_binaryOP_bm(benchmark_base::configuration cfg)
        : base_type{cfg} {}
    vector_binaryOP_bm(const vector_binaryOP_bm &bm) = default;
    vector_binaryOP_bm &operator=(vector_binaryOP_bm &other) = default;

    constexpr std::string name() const override {
        return std::string{base_type::name} + "_" + std::string{binaryOP::name};
    }

    inline void operator()(::benchmark::State &state) const override {

        using result_t = std::invoke_result_t<binaryOP, vector_t<scalar_t>,
                                              vector_t<scalar_t>>;

        const std::size_t n_samples{this->m_cfg.n_samples()};

        // Run the benchmark
        for (auto _ : state) {
            for (std::size_t i{0}; i < n_samples; ++i) {
                result_t result = binaryOP{}(this->a[i], this->b[i]);
                ::benchmark::DoNotOptimize(result);
            }
        }
    }
};

// Functions to be benchmarked
namespace bench_op::vector {

struct add {
    static constexpr std::string_view name{"add"};
    template <concepts::vector vector_t>
    constexpr vector_t operator()(const vector_t &a, const vector_t &b) const {
        return a + b;
    }
};
struct sub {
    static constexpr std::string_view name{"sub"};
    template <concepts::vector vector_t>
    constexpr vector_t operator()(const vector_t &a, const vector_t &b) const {
        return a - b;
    }
};
struct dot {
    static constexpr std::string_view name{"dot"};
    template <concepts::vector vector_t>
    constexpr algebra::traits::scalar_t<vector_t> operator()(
        const vector_t &a, const vector_t &b) const {
        return algebra::vector::dot(a, b);
    }
};
struct cross {
    static constexpr std::string_view name{"cross"};
    template <concepts::vector vector_t>
    constexpr vector_t operator()(const vector_t &a, const vector_t &b) const {
        return algebra::vector::cross(a, b);
    }
};

// Macro for declaring vector unary ops
#define ALGEBRA_PLUGINS_BENCH_VECTOR(OP, RES)               \
    struct OP {                                             \
        static constexpr std::string_view name{#OP};        \
        template <concepts::vector vector_t>                \
        constexpr RES operator()(const vector_t &a) const { \
            return algebra::vector::OP(a);                  \
        }                                                   \
    };

ALGEBRA_PLUGINS_BENCH_VECTOR(phi, algebra::traits::scalar_t<vector_t>)
ALGEBRA_PLUGINS_BENCH_VECTOR(theta, algebra::traits::scalar_t<vector_t>)
ALGEBRA_PLUGINS_BENCH_VECTOR(eta, algebra::traits::scalar_t<vector_t>)
ALGEBRA_PLUGINS_BENCH_VECTOR(perp, algebra::traits::scalar_t<vector_t>)
ALGEBRA_PLUGINS_BENCH_VECTOR(norm, algebra::traits::scalar_t<vector_t>)
ALGEBRA_PLUGINS_BENCH_VECTOR(normalize, vector_t)

}  // namespace bench_op::vector

// Macro for defining all vector benchmark types
#define ALGEBRA_PLUGINS_DEFINE_VECTOR_BENCH(PLUGIN)                           \
    using phi_f_t =                                                           \
        vector_unaryOP_bm<PLUGIN::vector3, float, bench_op::vector::phi>;     \
    using theta_f_t =                                                         \
        vector_unaryOP_bm<PLUGIN::vector3, float, bench_op::vector::theta>;   \
    using perp_f_t =                                                          \
        vector_unaryOP_bm<PLUGIN::vector3, float, bench_op::vector::perp>;    \
    using norm_f_t =                                                          \
        vector_unaryOP_bm<PLUGIN::vector3, float, bench_op::vector::norm>;    \
    using eta_f_t =                                                           \
        vector_unaryOP_bm<PLUGIN::vector3, float, bench_op::vector::eta>;     \
                                                                              \
    using add_f_t =                                                           \
        vector_binaryOP_bm<PLUGIN::vector3, float, bench_op::vector::add>;    \
    using sub_f_t =                                                           \
        vector_binaryOP_bm<PLUGIN::vector3, float, bench_op::vector::sub>;    \
    using dot_f_t =                                                           \
        vector_binaryOP_bm<PLUGIN::vector3, float, bench_op::vector::dot>;    \
    using cross_f_t =                                                         \
        vector_binaryOP_bm<PLUGIN::vector3, float, bench_op::vector::cross>;  \
    using normlz_f_t = vector_unaryOP_bm<PLUGIN::vector3, float,              \
                                         bench_op::vector::normalize>;        \
                                                                              \
    using phi_d_t =                                                           \
        vector_unaryOP_bm<PLUGIN::vector3, double, bench_op::vector::phi>;    \
    using theta_d_t =                                                         \
        vector_unaryOP_bm<PLUGIN::vector3, double, bench_op::vector::theta>;  \
    using perp_d_t =                                                          \
        vector_unaryOP_bm<PLUGIN::vector3, double, bench_op::vector::perp>;   \
    using norm_d_t =                                                          \
        vector_unaryOP_bm<PLUGIN::vector3, double, bench_op::vector::norm>;   \
    using eta_d_t =                                                           \
        vector_unaryOP_bm<PLUGIN::vector3, double, bench_op::vector::eta>;    \
                                                                              \
    using add_d_t =                                                           \
        vector_binaryOP_bm<PLUGIN::vector3, double, bench_op::vector::add>;   \
    using sub_d_t =                                                           \
        vector_binaryOP_bm<PLUGIN::vector3, double, bench_op::vector::sub>;   \
    using dot_d_t =                                                           \
        vector_binaryOP_bm<PLUGIN::vector3, double, bench_op::vector::dot>;   \
    using cross_d_t =                                                         \
        vector_binaryOP_bm<PLUGIN::vector3, double, bench_op::vector::cross>; \
    using normlz_d_t = vector_unaryOP_bm<PLUGIN::vector3, double,             \
                                         bench_op::vector::normalize>;

// Macro for registering all vector benchmarks
#define ALGEBRA_PLUGINS_REGISTER_VECTOR_BENCH(CFGS, CFGD)     \
    algebra::register_benchmark<add_f_t>(CFGS, "_single");    \
    algebra::register_benchmark<add_d_t>(CFGD, "_double");    \
    algebra::register_benchmark<sub_f_t>(CFGS, "_single");    \
    algebra::register_benchmark<sub_d_t>(CFGD, "_double");    \
    algebra::register_benchmark<dot_f_t>(CFGS, "_single");    \
    algebra::register_benchmark<dot_d_t>(CFGD, "_double");    \
    algebra::register_benchmark<cross_f_t>(CFGS, "_single");  \
    algebra::register_benchmark<cross_d_t>(CFGD, "_double");  \
    algebra::register_benchmark<normlz_f_t>(CFGS, "_single"); \
    algebra::register_benchmark<normlz_d_t>(CFGD, "_double"); \
                                                              \
    algebra::register_benchmark<phi_f_t>(CFGS, "_single");    \
    algebra::register_benchmark<phi_d_t>(CFGD, "_double");    \
    algebra::register_benchmark<theta_f_t>(CFGS, "_single");  \
    algebra::register_benchmark<theta_d_t>(CFGD, "_double");  \
    algebra::register_benchmark<perp_f_t>(CFGS, "_single");   \
    algebra::register_benchmark<perp_d_t>(CFGD, "_double");   \
    algebra::register_benchmark<norm_f_t>(CFGS, "_single");   \
    algebra::register_benchmark<norm_d_t>(CFGD, "_double");   \
    algebra::register_benchmark<eta_f_t>(CFGS, "_single");    \
    algebra::register_benchmark<eta_d_t>(CFGD, "_double");

}  // namespace algebra
