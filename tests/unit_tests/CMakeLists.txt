message(STATUS "Testing: 'algebra' unit tests")

mark_as_advanced(
    BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
    gmock_build_tests gtest_build_samples gtest_build_tests
    gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
)

include(GoogleTest)

macro(add_algebra_test TESTNAME FILES PLUGIN_LIBRARY)
    add_executable(${TESTNAME} ${FILES})
    target_link_libraries(${TESTNAME} PRIVATE ${PLUGIN_EXTRA_LIBRARIES})
    target_link_libraries(${TESTNAME} PRIVATE gtest gmock gtest_main)
    target_link_libraries(${TESTNAME} PRIVATE algebra::tests_common)
    target_link_libraries(${TESTNAME} PRIVATE  ${PLUGIN_LIBRARY} ${UNIT_TEST_EXTRA_LIBRARIES})
    set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
    add_test(unit_test_${TESTNAME} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${TESTNAME})
    set_tests_properties(unit_test_${TESTNAME} PROPERTIES 
        ENVIRONMENT ALGEBRA_PLUGIN_TEST_DATA_DIR=${ALGEBRA_PLUGIN_SOURCE_DIR}/data/)
endmacro()

set(all_unit_tests "plugin")

if(ALGEBRA_PLUGIN_INCLUDE_ARRAY)
    add_subdirectory(array)
endif()
if(ALGEBRA_PLUGIN_INCLUDE_EIGEN)
    add_subdirectory(eigen)
endif()
if(ALGEBRA_PLUGIN_INCLUDE_SMATRIX)
    add_subdirectory(smatrix)
endif()
if(ALGEBRA_PLUGIN_INCLUDE_VC)
    add_subdirectory(vc)
endif()
if(ALGEBRA_PLUGIN_USE_VECMEM)
    add_subdirectory(vecmem)
endif()
