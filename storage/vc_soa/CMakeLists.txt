# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the library.
algebra_add_library( algebra_vc_soa_storage vc_soa_storage
   "include/algebra/storage/vc_soa.hpp" )
target_link_libraries( algebra_vc_soa_storage
   INTERFACE algebra::common algebra::common_storage Vc::Vc )
algebra_test_public_headers( algebra_vc_soa_storage
   "algebra/storage/vc_soa.hpp" )
