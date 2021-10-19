# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Let the user know what's happening.
message( STATUS "Building 'algebra::vc_array' plugin" )

# Set up the library.
algebra_add_library( algebra_vc_array vc_array
   "include/algebra/definitions/vc_array.hpp" )
target_compile_definitions( algebra_vc_array INTERFACE
   ALGEBRA_PLUGIN_INCLUDE_VC=TRUE )
target_link_libraries( algebra_vc_array INTERFACE
   algebra::common Vc )
