# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required( VERSION 3.24 )
include( FetchContent )

# Declare where to get Eigen3 from, in case it is not found with find_package.
set( ALGEBRA_PLUGINS_EIGEN_SOURCE
   "URL;https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2;URL_MD5;132dde48fe2b563211675626d29f1707"
   CACHE STRING "Source for Eigen, when built as part of this project" )
mark_as_advanced( ALGEBRA_PLUGINS_EIGEN_SOURCE )
FetchContent_Declare( Eigen3 ${ALGEBRA_PLUGINS_EIGEN_SOURCE}
   FIND_PACKAGE_ARGS )

# Turn off the unit tests for Eigen3.
if( DEFINED CACHE{BUILD_TESTING} )
   set( _buildTestingValue ${BUILD_TESTING} )
endif()
set( BUILD_TESTING FALSE CACHE INTERNAL "Forceful setting of BUILD_TESTING" )

# Set it up.
FetchContent_MakeAvailable( Eigen3 )

# Reset the BUILD_TESTING variable.
if( DEFINED _buildTestingValue )
   set( BUILD_TESTING ${_buildTestingValue} CACHE BOOL "Turn tests on/off"
      FORCE )
   unset( _buildTestingValue )
else()
   unset( BUILD_TESTING CACHE )
endif()

# Treat the Eigen headers as "system headers", to avoid getting warnings from
# them. But only bother about this when building Eigen as part of this project.
if( TARGET eigen )
   get_target_property( _incDirs eigen INTERFACE_INCLUDE_DIRECTORIES )
   target_include_directories( eigen
      SYSTEM INTERFACE ${_incDirs} )
   unset( _incDirs )
endif()
