# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the library.
algebra_add_library( algebra_std_array std_array
   "include/algebra/storage/std_array.hpp" )
target_link_libraries( algebra_std_array
   INTERFACE algebra::common )
