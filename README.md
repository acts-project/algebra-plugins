# algebra-plugins

This repository provides different algebra plugins with minimal functionality
for the R&D projects [detray](https://github.com/acts-project/detray) and
[traccc](https://github.com/acts-project/traccc).

To build it standalone, run e.g.

```
mkdir build
cd build
cmake -DALGEBRA_PLUGIN_INCLUDE_ARRAY=ON ..
make -j5
```

Available options:

- `ALGEBRA_PLUGIN_CUSTOM_SCALARTYPE`: Scalar value type
  (`"float"` or `"double"`, `"double"` by default)
- `ALGEBRA_PLUGIN_INCLUDE_<XXX>`: Boolean to turn on/off the build of one of
  the following plugins:
  * `ARRAY`: Plugin using
    [std::array](https://en.cppreference.com/w/cpp/container/array)
    (`ON` by default)
  * `EIGEN`: Plugin using [Eigen](https://eigen.tuxfamily.org)
    (`OFF` by default)
  * `SMATRIX`: Plugin using
    [Smatrix](https://root.cern/doc/master/group__SMatrixGroup.html)
    (`OFF` by default)
  * `VC`: Plugin using [Vc](https://github.com/VcDevel/Vc)
    (`OFF` by default)
  * `VECMEM`: Plugin using [VecMem](https://github.com/acts-project/vecmem)
    (`OFF` by default)
- `ALGEBRA_PLUGIN_SETUP_<XXX>`: Boolean to turn on/off the explicit "setup" of
  the externals (`GOOGLETEST`, `EIGEN3`, `VC` and `VECMEM`)
- `ALGEBRA_PLUGIN_USE_SYSTEM_<XXX>`: Boolean configuring how to set up a given
  external
  * `ON`: The external is searched for "on the system" using
    [find_package](https://cmake.org/cmake/help/latest/command/find_package.html);
  * `OFF`: The package is set up for build as part of this project, using
    [FetchContent](https://cmake.org/cmake/help/latest/module/FetchContent.html).
- `ALGEBRA_PLUGIN_BUILD_TESTING`: Turn the build/setup of the unit tests on/off
  (`ON` by default)
