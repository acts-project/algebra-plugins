message(STATUS "Building 'algebra::vc_array' plugin" ${})

add_library(vc_array INTERFACE)

target_include_directories(vc_array
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    ${ALGEBRA_PLUGIN_SOURCE_DIR}/common/include/algebra)

install(
  DIRECTORY include/plugins/algebra
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(ALGEBRA_PLUGIN_CUSTOM_SCALARTYPE)
  target_compile_definitions(
    vc_array
    INTERFACE -DALGEBRA_PLUGIN_CUSTOM_SCALARTYPE=${ALGEBRA_PLUGIN_CUSTOM_SCALARTYPE})
endif()

target_compile_options(vc_array INTERFACE ${Vc_ARCHITECTURE_FLAGS}  -DALGEBRA_PLUGIN_INCLUDE_VC)

target_link_libraries(Vc)

add_library(algebra::vc_array ALIAS vc_array)
