# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required( VERSION 3.24 )
include( FetchContent )

# Declare where to get GoogleTest from.
set( ALGEBRA_PLUGINS_GOOGLETEST_SOURCE
   "URL;https://github.com/google/googletest/archive/release-1.11.0.tar.gz;URL_MD5;e8a8df240b6938bb6384155d4c37d937"
   CACHE STRING "Source for GoogleTest, when built as part of this project" )
mark_as_advanced( ALGEBRA_PLUGINS_GOOGLETEST_SOURCE )
FetchContent_Declare( GTest ${ALGEBRA_PLUGINS_GOOGLETEST_SOURCE}
   FIND_PACKAGE_ARGS )

# Options used in the build of GoogleTest.
set( BUILD_GMOCK FALSE CACHE BOOL "Turn off the build of GMock" )
set( INSTALL_GTEST FALSE CACHE BOOL "Turn off the installation of GoogleTest" )
if( WIN32 )
   set( gtest_force_shared_crt TRUE CACHE BOOL
      "Use shared (DLL) run-time library, even with static libraries" )
endif()

# Silence some warnings with modern versions of CMake on macOS.
set( CMAKE_MACOSX_RPATH TRUE )

# Set it up.
FetchContent_MakeAvailable( GTest )

# Set up aliases for the GTest targets with the same name that they have
# when we find GTest pre-installed.
if( NOT TARGET GTest::gtest )
   add_library( GTest::gtest ALIAS gtest )
endif()
if( NOT TARGET GTest::gtest_main )
   add_library( GTest::gtest_main ALIAS gtest_main )
endif()
