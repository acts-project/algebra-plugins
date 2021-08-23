message(STATUS "Building 'algebra::array' plugin")

add_library(algebra_array INTERFACE)

target_include_directories(algebra_array
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    ${ALGEBRA_PLUGIN_SOURCE_DIR}/core/common/include/algebra)

install(
  DIRECTORY include/plugins/algebra
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(ALGEBRA_PLUGIN_CUSTOM_SCALARTYPE)
  target_compile_definitions(
    algebra_array
    INTERFACE -DALGEBRA_PLUGIN_CUSTOM_SCALARTYPE=${ALGEBRA_PLUGIN_CUSTOM_SCALARTYPE})
endif()

if(ALGEBRA_PLUGIN_USE_VECMEM)
  target_compile_options(algebra_array INTERFACE -DALGEBRA_PLUGIN_USE_VECMEM)
  target_link_libraries(algebra_array INTERFACE vecmem::core)
endif()

add_library(algebra::array ALIAS algebra_array)
