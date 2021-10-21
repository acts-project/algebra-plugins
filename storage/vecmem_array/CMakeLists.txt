# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the library.
algebra_add_library( algebra_vecmem_array vecmem_array
   "include/algebra/storage/vecmem_array.hpp" )
target_link_libraries( algebra_vecmem_array
   INTERFACE algebra::common vecmem::core )
