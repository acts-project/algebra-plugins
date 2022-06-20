# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the library.
algebra_add_library( algebra_array_storage array_storage
   "include/algebra/storage/array.hpp" )
target_link_libraries( algebra_array_storage
   INTERFACE algebra::common )
