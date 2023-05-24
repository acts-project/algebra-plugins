# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required( VERSION 3.14 )
include( FetchContent )

# Suppress developer warnings for Vc by default.
# Unless the user explicitly requested otherwise.
if( NOT DEFINED CMAKE_SUPPRESS_DEVELOPER_WARNINGS )
   set( CMAKE_SUPPRESS_DEVELOPER_WARNINGS TRUE )
endif()

# Declare where to get Vc from, in case it is not found with find_package.
set( ALGEBRA_PLUGINS_VC_SOURCE
   "URL;https://github.com/VcDevel/Vc/archive/b84dcd0a65d8dc5de6a2bd4d367882b3748f812c.tar.gz;URL_MD5;ae1e8317352df53113ee28ba3f0d490a"
   CACHE STRING "Source for Vc, when built as part of this project" )
mark_as_advanced( ALGEBRA_PLUGINS_VC_SOURCE )
FetchContent_Declare( Vc ${ALGEBRA_PLUGINS_VC_SOURCE}
   FIND_PACKAGE_ARGS 1.4.2 )

# Turn off the unit tests for Vc.
if( DEFINED CACHE{BUILD_TESTING} )
   set( _buildTestingValue ${BUILD_TESTING} )
endif()
set( BUILD_TESTING FALSE CACHE INTERNAL "Forceful setting of BUILD_TESTING" )

# Set it up.
FetchContent_MakeAvailable( Vc )

# Reset the BUILD_TESTING variable.
if( DEFINED _buildTestingValue )
   set( BUILD_TESTING ${_buildTestingValue} CACHE BOOL "Turn tests on/off"
      FORCE )
   unset( _buildTestingValue )
else()
   unset( BUILD_TESTING CACHE )
endif()

# Treat the Vc headers as "system headers", to avoid getting warnings from them.
if( TARGET Vc )
   target_include_directories( Vc
      SYSTEM INTERFACE $<TARGET_PROPERTY:Vc,INCLUDE_DIRECTORIES> )
endif()

# Disable a warning for GCC which the Vc build turns on explicitly for that
# compiler, just to then go and trigger it... :-/
if( ( "${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU" ) AND ( TARGET Vc ) )
   target_compile_options( Vc PRIVATE "-Wno-old-style-cast" )
endif()

# Set up an alias for the Vc target, with the same name that it will have
# when "finding it".
if( NOT TARGET Vc::Vc )
   add_library( Vc::Vc ALIAS Vc )
endif()
