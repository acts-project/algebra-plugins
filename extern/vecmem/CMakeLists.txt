# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required( VERSION 3.14 )
include( FetchContent )

# Tell the user what's happening.
message( STATUS "Building VecMem as part of the Algebra Plugins project" )

# Declare where to get VecMem from.
set( ALGEBRA_PLUGINS_VECMEM_SOURCE
   "URL;https://github.com/acts-project/vecmem/archive/refs/tags/v0.27.0.tar.gz;URL_MD5;cd1520efbd46d5d09ac5727452939ede"
   CACHE STRING "Source for VecMem, when built as part of this project" )
mark_as_advanced( ALGEBRA_PLUGINS_VECMEM_SOURCE )
FetchContent_Declare( VecMem ${ALGEBRA_PLUGINS_VECMEM_SOURCE} )

# Options used in the build of VecMem.
set( VECMEM_BUILD_TESTING FALSE CACHE BOOL
   "Turn off the build of the VecMem unit tests" )

# Get it into the current directory.
FetchContent_MakeAvailable( VecMem )
