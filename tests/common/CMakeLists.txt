add_library(algebra_tests_common INTERFACE)

target_include_directories(algebra_tests_common
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

if(DETRAY_CUSTOM_SCALARTYPE)
  target_compile_definitions(
    algebra_tests_common
    INTERFACE -DALGEBRA_PLUGIN_CUSTOM_SCALARTYPE=${ALGEBRA_PLUGIN_CUSTOM_SCALARTYPE})
endif()

add_library(algebra::tests_common ALIAS algebra_tests_common)
