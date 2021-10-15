# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Helper function for setting up the algebra plugins libraries.
#
# Usage: algebra_add_library( algebra_array array "header1.hpp"... )
#
function( algebra_add_library fullname basename )

   # Create the library.
   add_library( ${fullname} INTERFACE ${ARG_UNPARSED_ARGUMENTS} )

   # Set up how clients should find its headers.
   target_include_directories( ${fullname} INTERFACE
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}> )

   # Make sure that the library is available as "algebra::${basename}" in every
   # situation.
   set_target_properties( ${fullname} PROPERTIES EXPORT_NAME ${basename} )
   add_library( algebra::${basename} ALIAS ${fullname} )

   # Set up the installation of the library and its headers.
   install( TARGETS ${fullname}
      EXPORT algebra-plugins-exports
      LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" )
   install( DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/"
      DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" )

endfunction( algebra_add_library )

# Helper function for setting up the algebra plugin tests.
#
# Usage: algebra_add_test( array source1.cpp source2.cpp
#                          LINK_LIBRARIES algebra::array )
#
function( algebra_add_test name )

   # Parse the function's options.
   cmake_parse_arguments( ARG "" "" "LINK_LIBRARIES" ${ARGN} )

   # Create the test executable.
   set( test_exe_name "algebra_test_${name}" )
   add_executable( ${test_exe_name} ${ARG_UNPARSED_ARGUMENTS} )
   if( ARG_LINK_LIBRARIES )
      target_link_libraries( ${test_exe_name} PRIVATE ${ARG_LINK_LIBRARIES} )
   endif()

   # Run the executable as the test.
   add_test( NAME ${test_exe_name}
      COMMAND $<TARGET_FILE:${test_exe_name}> )

endfunction( algebra_add_test )
