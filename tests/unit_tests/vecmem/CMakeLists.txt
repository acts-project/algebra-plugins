enable_testing()

foreach(etest ${all_unit_tests})
    add_algebra_test(vecmem_array_algebra_${etest}
                     vecmem_array_algebra_${etest}.cpp 
                     algebra::vecmem_array)
endforeach(etest)

