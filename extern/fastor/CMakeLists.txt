# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2022-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required( VERSION 3.24 )
include( FetchContent )

# Declare where to get Fastor from, in case it is not found with find_package.
# We need to use this alternative syntax for FetchContent_Declare because the
# latest release for Fastor does not have a CMakeLists.txt file.
set( ALGEBRA_PLUGINS_FASTOR_SOURCE
   "fastor;GIT_REPOSITORY;https://github.com/romeric/Fastor.git;GIT_TAG;e96e63fe8a05898d732d0b9faba3667d9a25a38f"
   CACHE STRING "Source for Fastor, when built as part of this project" )
mark_as_advanced( ALGEBRA_PLUGINS_FASTOR_SOURCE )
FetchContent_Declare( Fastor ${ALGEBRA_PLUGINS_FASTOR_SOURCE}
   FIND_PACKAGE_ARGS 0.6.3 )

# Turn off the unit tests for Fastor.
if( DEFINED CACHE{BUILD_TESTING} )
   set( _buildTestingValue ${BUILD_TESTING} )
endif()
set( BUILD_TESTING FALSE CACHE INTERNAL "Forceful setting of BUILD_TESTING" )

# Set it up.
FetchContent_MakeAvailable( Fastor )

# Reset the BUILD_TESTING variable.
if( DEFINED _buildTestingValue )
   set( BUILD_TESTING ${_buildTestingValue} CACHE BOOL "Turn tests on/off"
      FORCE )
   unset( _buildTestingValue )
else()
   unset( BUILD_TESTING CACHE )
endif()

# Treat the Fastor headers as "system headers", to avoid getting warnings from
# them.
if( TARGET Fastor )
   get_target_property( _incDirs Fastor INTERFACE_INCLUDE_DIRECTORIES )
   target_include_directories( Fastor
      SYSTEM INTERFACE ${_incDirs} )
   unset( _incDirs )
endif()

# Set up an alias for the Fastor target, with the same name that it will have
# when "finding it".
if( NOT TARGET Fastor::Fastor )
   add_library( Fastor::Fastor ALIAS Fastor )
endif()
