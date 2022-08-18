# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the library.
algebra_add_library(algebra_common_math common_math
   # Coordinates
   "include/algebra/coordinates/cartesian2.hpp"
   "include/algebra/coordinates/cartesian3.hpp"
   "include/algebra/coordinates/coordinate_base.hpp"
   "include/algebra/coordinates/coordinates.hpp"
   "include/algebra/coordinates/cylindrical2.hpp"
   "include/algebra/coordinates/cylindrical3.hpp"
   "include/algebra/coordinates/line2.hpp"
   "include/algebra/coordinates/polar2.hpp"
   # Math
   "include/algebra/math/common.hpp"
   # Tracks
   "include/algebra/tracks/bound_track_parameters.hpp"
   "include/algebra/tracks/detail/track_helper.hpp"
   "include/algebra/tracks/free_track_parameters.hpp"
   "include/algebra/tracks/tracks.hpp"
   # Utils
   "include/algebra/utils/column_wise_operator.hpp")
target_link_libraries(algebra_common_math
   INTERFACE algebra::common)
