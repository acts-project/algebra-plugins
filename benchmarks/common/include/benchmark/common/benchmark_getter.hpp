/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "benchmark_vector.hpp"

namespace algebra::bench_op {

// Functions to be benchmarked

struct phi {
  inline static const std::string name{"phi"};
  template <typename vector_t>
  auto operator()(const vector_t &a) const {
    return algebra::getter::phi(a);
  }
};
struct theta {
  inline static const std::string name{"theta"};
  template <typename vector_t>
  auto operator()(const vector_t &a) const {
    return algebra::getter::theta(a);
  }
};
struct perp {
  inline static const std::string name{"perp"};
  template <typename vector_t>
  auto operator()(const vector_t &a) const {
    return algebra::getter::perp(a);
  }
};
struct norm {
  inline static const std::string name{"norm"};
  template <typename vector_t>
  auto operator()(const vector_t &a) const {
    return algebra::getter::norm(a);
  }
};
struct eta {
  inline static const std::string name{"eta"};
  template <typename vector_t>
  auto operator()(const vector_t &a) const {
    return algebra::getter::eta(a);
  }
};

}  // namespace algebra::bench_op