enable_testing()

foreach(etest ${all_unit_tests})
    add_algebra_test(eigen_algebra_${etest}
                     eigen_algebra_${etest}.cpp algebra::eigen)
endforeach(etest)
