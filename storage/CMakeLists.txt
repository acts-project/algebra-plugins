# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up all enabled libraries.
add_subdirectory( array )
if( ALGEBRA_PLUGINS_INCLUDE_EIGEN )
   add_subdirectory( eigen )
endif()
if( ALGEBRA_PLUGINS_INCLUDE_SMATRIX )
   add_subdirectory( smatrix )
endif()
if( ALGEBRA_PLUGINS_INCLUDE_VC )
   add_subdirectory( vc )
endif()
if( ALGEBRA_PLUGINS_INCLUDE_FASTOR )
   add_subdirectory( fastor )
endif()
if( ALGEBRA_PLUGINS_INCLUDE_VECMEM )
   add_subdirectory( vecmem )
endif()
