# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up all enabled libraries.
add_subdirectory( common )
add_subdirectory( cmath )
add_subdirectory( generic )
if( ALGEBRA_PLUGINS_INCLUDE_EIGEN )
   add_subdirectory( eigen )
endif()
if( ALGEBRA_PLUGINS_INCLUDE_FASTOR )
   add_subdirectory( fastor )
endif()
if( ALGEBRA_PLUGINS_INCLUDE_SMATRIX )
   add_subdirectory( smatrix )
endif()
if( ALGEBRA_PLUGINS_INCLUDE_VC )
   add_subdirectory( vc_aos )
   add_subdirectory( vc_soa )
endif()
