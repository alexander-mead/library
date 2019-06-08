#!/bin/bash

tests+=('cosmology_functions')
tests+=('logical_operations')
tests+=('array_operations')

compiler='gfortran'

normal='-Warray-bounds -ffree-line-length-none -fmax-errors=4 -ffpe-trap=invalid,zero,overflow -fimplicit-none -std=gnu'   
slow='-Wall -fcheck=all -fbounds-check -fbacktrace -Og'
library='-L/usr/local/lib -lfftw3 -lfftw3f -lfftw3l -I/Users/Mead/Physics/library /Users/Mead/Physics/library/meadlib -L/Users/Mead/Physics/library'

for i in "${tests[@]}"; do
    i_test=${i}_test
    echo 'Compiling:' $i_test
    $compiler $i_test/$i_test.f90 $normal $slow $library -o $i_test/$i_test.e
    ./$i_test/$i_test.e
done
