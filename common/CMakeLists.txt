# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the library.
algebra_add_library( algebra_common common
   "include/algebra/common/definitions/array.hpp"
   "include/algebra/common/simd/simd_array_wrapper.hpp"
   "include/algebra/common/simd/simd_types.hpp"
   "include/algebra/common/types.hpp" )
target_compile_definitions( algebra_common INTERFACE
   ALGEBRA_PLUGIN_CUSTOM_SCALARTYPE=${ALGEBRA_PLUGIN_CUSTOM_SCALARTYPE} )
