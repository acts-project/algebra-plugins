# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Let the user know what's happening.
message( STATUS "Building 'algebra::array' plugin" )

# Set up the library.
algebra_add_library( algebra_array array
   "include/algebra/definitions/std_array.hpp" )
target_link_libraries( algebra_array INTERFACE algebra::common )
