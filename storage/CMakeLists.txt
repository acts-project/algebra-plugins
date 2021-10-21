# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up all enabled libraries.
if( ALGEBRA_PLUGIN_INCLUDE_ARRAY )
   add_subdirectory( std_array )
endif()
if( ALGEBRA_PLUGIN_INCLUDE_VECMEM )
   add_subdirectory( vecmem_array )
endif()
