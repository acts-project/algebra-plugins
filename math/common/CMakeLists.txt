# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the library.
algebra_add_library(algebra_common_math common_math
   "include/algebra/coordinates/cartesian2.hpp"
   "include/algebra/coordinates/coordinate_base.hpp"
   "include/algebra/coordinates/coordinates.hpp"
   "include/algebra/coordinates/cylindrical2.hpp"
   "include/algebra/coordinates/line2.hpp"
   "include/algebra/coordinates/polar2.hpp"
   "include/algebra/math/common.hpp")
target_link_libraries(algebra_common_math
   INTERFACE algebra::common)
