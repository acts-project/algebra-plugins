# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the library.
algebra_add_library( algebra_common_math common_math
   # coordinate Include(s)
   "include/algebra/coordinate/cartesian2.hpp"
   "include/algebra/coordinate/coordinates.hpp"
   "include/algebra/coordinate/cylindrical2.hpp"
   "include/algebra/coordinate/polar2.hpp"
   "include/algebra/coordinate/radial.hpp"
   # math Include(s)
   "include/algebra/math/common.hpp" )
target_link_libraries( algebra_common_math
   INTERFACE algebra::common )
