# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required( VERSION 3.14 )
include( FetchContent )

# Tell the user what's happening.
message( STATUS "Building Vc as part of the Algebra Plugins project" )

# Declare where to get Vc from.
set( ALGEBRA_PLUGINS_VC_SOURCE
   "GIT_REPOSITORY;https://github.com/VcDevel/Vc.git;GIT_TAG;b84dcd0"
   CACHE STRING "Source for Vc, when built as part of this project" )
mark_as_advanced( ALGEBRA_PLUGINS_VC_SOURCE )
FetchContent_Declare( Vc ${ALGEBRA_PLUGINS_VC_SOURCE} )

# Turn off the unit tests for Vc.
if( DEFINED CACHE{BUILD_TESTING} )
   set( _buildTestingValue ${BUILD_TESTING} )
endif()
set( BUILD_TESTING FALSE CACHE INTERNAL "Forceful setting of BUILD_TESTING" )

# Get it into the current directory.
FetchContent_MakeAvailable( Vc )

# Reset the BUILD_TESTING variable.
if( DEFINED _buildTestingValue )
   set( BUILD_TESTING ${_buildTestingValue} CACHE BOOL "Turn tests on/off"
      FORCE )
   unset( _buildTestingValue )
else()
   unset( BUILD_TESTING CACHE )
endif()

# Set up an alias for the Vc target, with the same name that it will have
# when "finding it".
add_library( Vc::Vc ALIAS Vc )
