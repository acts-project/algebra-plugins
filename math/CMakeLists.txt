# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up all enabled libraries.
add_subdirectory( cmath )
if( ALGEBRA_PLUGIN_INCLUDE_EIGEN )
   add_subdirectory( eigen )
endif()
if( ALGEBRA_PLUGIN_INCLUDE_SMATRIX )
   add_subdirectory( smatrix )
endif()
if( ALGEBRA_PLUGIN_INCLUDE_VC )
   add_subdirectory( vc )
endif()
