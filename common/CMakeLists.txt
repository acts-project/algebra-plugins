# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the library.
algebra_add_library( algebra_common common
   "include/algebra/assert.hpp"
   "include/algebra/concepts.hpp"
   "include/algebra/qualifiers.hpp"
   "include/algebra/type_traits.hpp"
)

algebra_test_public_headers( algebra_common
   "algebra/assert.hpp"
   "algebra/concepts.hpp"
   "algebra/qualifiers.hpp"
   "algebra/type_traits.hpp"
)
