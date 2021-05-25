

if (ALGEBRA_INCLUDE_ARRAY_PLUGIN)
    add_subdirectory(array)
endif()

if (ALGEBRA_INCLUDE_EIGEN_PLUGIN)
    add_subdirectory(eigen)
endif()

if (ALGEBRA_INCLUDE_SMATRIX_PLUGIN)
    add_subdirectory(smatrix)
endif()
