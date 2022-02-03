# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Include the helper function(s).
include( algebra-plugins-functions )

# Turn on the correct setting for the __cplusplus macro with MSVC.
if( "${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC" )
   algebra_add_flag( CMAKE_CXX_FLAGS "/Zc:__cplusplus" )
   algebra_add_flag( CMAKE_CUDA_FLAGS "-Xcompiler /Zc:__cplusplus" )
endif()

# Turn on a number of warnings for the "known compilers".
if( ( "${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU" ) OR
    ( "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" ) )

   # Basic flags for all build modes.
   algebra_add_flag( CMAKE_CXX_FLAGS "-Wall" )
   algebra_add_flag( CMAKE_CXX_FLAGS "-Wextra" )
   algebra_add_flag( CMAKE_CXX_FLAGS "-Wshadow" )
   algebra_add_flag( CMAKE_CXX_FLAGS "-Wunused-local-typedefs" )

   # More rigorous tests for the Debug builds.
   algebra_add_flag( CMAKE_CXX_FLAGS_DEBUG "-Werror" )
   algebra_add_flag( CMAKE_CXX_FLAGS_DEBUG "-pedantic" )

elseif( "${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC" )

   # Basic flags for all build modes.
   string( REGEX REPLACE "/W[0-9]" "" CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS}" )
   algebra_add_flag( CMAKE_CXX_FLAGS "/W4" )

   # More rigorous tests for the Debug builds.
   algebra_add_flag( CMAKE_CXX_FLAGS_DEBUG "/WX" )

endif()
