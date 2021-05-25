message(STATUS "Building 'algebra::smatrix' plugin")

add_library(algebra_smatrix INTERFACE)

target_include_directories(algebra_smatrix
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

install(
  DIRECTORY include/plugins
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(ALGEBRA_CUSTOM_SCALARTYPE)
  target_compile_definitions(
    algebra_smatrix
    INTERFACE -DALGEBRA_CUSTOM_SCALARTYPE=${ALGEBRA_CUSTOM_SCALARTYPE})
endif()

add_library(algebra::smatrix ALIAS algebra_smatrix)
