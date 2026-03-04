# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the common library.
algebra_add_library( algebra_common common
   "include/algebra/array_operators.hpp"
   "include/algebra/assert.hpp"
   "include/algebra/boolean.hpp"
   "include/algebra/concepts.hpp"
   "include/algebra/math.hpp"
   "include/algebra/matrix_getter.hpp"
   "include/algebra/matrix.hpp"
   "include/algebra/qualifiers.hpp"
   "include/algebra/type_traits.hpp"
   "include/algebra/vector.hpp"
)

algebra_test_public_headers( algebra_common
   "algebra/array_operators.hpp"
   "algebra/assert.hpp"
   "algebra/boolean.hpp"
   "algebra/concepts.hpp"
   "algebra/math.hpp"
   "algebra/matrix_getter.hpp"
   "algebra/matrix.hpp"
   "algebra/qualifiers.hpp"
   "algebra/type_traits.hpp"
   "algebra/vector.hpp"
)
