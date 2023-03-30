# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Include the helper function(s).
include( algebra-plugins-functions )

# Basic flags for all build modes.
foreach( mode RELEASE RELWITHDEBINFO MINSIZEREL DEBUG )
   algebra_add_flag( CMAKE_SYCL_FLAGS_${mode} "-Wall" )
   algebra_add_flag( CMAKE_SYCL_FLAGS_${mode} "-Wextra" )
   algebra_add_flag( CMAKE_SYCL_FLAGS_${mode} "-Wno-unknown-cuda-version" )
   algebra_add_flag( CMAKE_SYCL_FLAGS_${mode} "-Wshadow" )
   algebra_add_flag( CMAKE_SYCL_FLAGS_${mode} "-Wunused-local-typedefs" )
endforeach()
if( NOT WIN32 )
   foreach( mode RELEASE RELWITHDEBINFO MINSIZEREL DEBUG )
      algebra_add_flag( CMAKE_SYCL_FLAGS_${mode} "-pedantic" )
   endforeach()
endif()

# Avoid issues coming from MSVC<->DPC++ argument differences.
if( "${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC" )
   foreach( mode RELEASE RELWITHDEBINFO MINSIZEREL DEBUG )
      algebra_add_flag( CMAKE_SYCL_FLAGS_${mode}
         "-Wno-unused-command-line-argument" )
   endforeach()
endif()

# Fail on warnings, if asked for that behaviour.
if( ALGEBRA_PLUGINS_FAIL_ON_WARNINGS )
   foreach( mode RELEASE RELWITHDEBINFO MINSIZEREL DEBUG )
      algebra_add_flag( CMAKE_SYCL_FLAGS_${mode} "-Werror" )
   endforeach()
endif()
