# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the library.
algebra_add_library(algebra_common_math common_math
   # Math
   "include/algebra/math/common.hpp"
   # Utils
   "include/algebra/utils/column_wise_operator.hpp"
   "include/algebra/utils/utils.hpp")
target_link_libraries(algebra_common_math
   INTERFACE algebra::common)
