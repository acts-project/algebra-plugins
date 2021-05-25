enable_testing()

foreach(etest ${all_unit_tests})
    add_algebra_test(array_algebra_${etest}
                     array_algebra_${etest}.cpp 
                     algebra::array)
endforeach(etest)

