message(STATUS "Building 'algebra::vecmem_array' plugin")

add_library(vecmem_array INTERFACE)

target_include_directories(vecmem_array
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    ${ALGEBRA_PLUGIN_SOURCE_DIR}/common/include/algebra)

install(
  DIRECTORY include/plugins/algebra
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(ALGEBRA_PLUGIN_CUSTOM_SCALARTYPE)
  target_compile_definitions(
    vecmem_array
    INTERFACE -DALGEBRA_PLUGIN_CUSTOM_SCALARTYPE=${ALGEBRA_PLUGIN_CUSTOM_SCALARTYPE})
endif()

target_compile_options(vecmem_array INTERFACE -DALGEBRA_PLUGIN_USE_VECMEM)
target_link_libraries(vecmem_array INTERFACE vecmem::core)

add_library(algebra::vecmem_array ALIAS vecmem_array)
