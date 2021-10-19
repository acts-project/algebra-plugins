# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Let the user know what's happening.
message( STATUS "Building 'algebra::vecmem_array' plugin" )

# Set up the library.
algebra_add_library( algebra_vecmem_array vecmem_array
   "include/algebra/definitions/vecmem_array.hpp" )
target_compile_definitions( algebra_vecmem_array INTERFACE
   ALGEBRA_PLUGIN_USE_VECMEM=TRUE )
target_link_libraries( algebra_vecmem_array INTERFACE
   algebra::common vecmem::core )
