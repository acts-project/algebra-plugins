# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required( VERSION 3.24 )
include( FetchContent )

# Declare where to get VecMem from, in case it is not found with find_package.
set( ALGEBRA_PLUGINS_VECMEM_SOURCE
   "URL;https://github.com/acts-project/vecmem/archive/refs/tags/v0.24.0.tar.gz;URL_MD5;f6936ea57921d7e5110606c7a9e9a371"
   CACHE STRING "Source for VecMem, when built as part of this project" )
mark_as_advanced( ALGEBRA_PLUGINS_VECMEM_SOURCE )
FetchContent_Declare( vecmem ${ALGEBRA_PLUGINS_VECMEM_SOURCE}
   FIND_PACKAGE_ARGS )

# Options used in the build of VecMem.
set( VECMEM_BUILD_TESTING FALSE CACHE BOOL
   "Turn off the build of the VecMem unit tests" )

# Set it up.
FetchContent_MakeAvailable( vecmem )
