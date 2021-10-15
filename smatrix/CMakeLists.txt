# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Let the user know what's happening.
message( STATUS "Building 'algebra::smatrix' plugin" )

# Find Smatrix.
find_package( ROOT COMPONENTS Smatrix REQUIRED)

# Set up the library.
algebra_add_library( algebra_smatrix smatrix
   "include/algebra/definitions/smatrix.hpp" )
target_link_libraries( algebra_smatrix INTERFACE
   algebra::common ROOT::Smatrix )
if( ALGEBRA_PLUGIN_INCLUDE_VC )
   target_link_libraries( algebra_smatrix INTERFACE algebra::vc_array )
endif()
if( ALGEBRA_PLUGIN_INCLUDE_VECMEM )
   target_link_libraries( algebra_smatrix INTERFACE algebra::vecmem_array )
endif()
