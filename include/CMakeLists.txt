

if (INCLUDE_ARRAY_PLUGIN)
    add_subdirectory(array)
endif()

if (INCLUDE_EIGEN_PLUGIN)
    add_subdirectory(eigen)
endif()

if (INCLUDE_SMATRIX_PLUGIN)
    add_subdirectory(smatrix)
endif()
