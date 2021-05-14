message(STATUS "Building 'algebra::array' plugin")

add_library(algebra_array INTERFACE)

target_include_directories(algebra_array
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

install(
  DIRECTORY include/plugins
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(DETRAY_CUSTOM_SCALARTYPE)
  target_compile_definitions(
    algebra_array
    INTERFACE -DALGEBRA_CUSTOM_SCALARTYPE=${ALGEBRA_CUSTOM_SCALARTYPE})
endif()

add_library(algebra::array ALIAS algebra_array)
