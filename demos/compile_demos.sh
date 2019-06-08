#!/bin/bash

compiler='gfortran'

normal='-Warray-bounds -ffree-line-length-none -fmax-errors=4 -ffpe-trap=invalid,zero,overflow -fimplicit-none -std=gnu'   
slow='-Wall -fcheck=all -fbounds-check -fbacktrace -Og'
library='-L/usr/local/lib -lfftw3 -lfftw3f -lfftw3l -I/Users/Mead/Physics/library /Users/Mead/Physics/library/meadlib -L/Users/Mead/Physics/library'

demos+=('GRF_demo')
demos+=('HMx_demo')
demos+=('accept_reject_demo')
demos+=('array_operations_demo')
demos+=('calculus_demo')
demos+=('calculus_table_demo')
demos+=('constants_demo')
demos+=('cosmology_functions_demo')
demos+=('dice_demo')
demos+=('fft_demo')
demos+=('field_operations_demo')
demos+=('file_info_demo')
demos+=('fitting_demo')
demos+=('fix_polynomial_demo')
demos+=('gadget_demo')
demos+=('generate_randoms_demo')
demos+=('interpolate_demo')
demos+=('logical_operations_demo')
demos+=('multidark_demo')
demos+=('numerology_demo')
demos+=('ode_solvers_demo')
demos+=('orbits_demo')
demos+=('physics_demo')
demos+=('power_demo')
demos+=('random_integers_demo')
demos+=('random_numbers_demo')
demos+=('shot_noise_demo')
demos+=('simulations_demo')
demos+=('solve_equations_demo')
demos+=('sorting_demo')
demos+=('special_functions_demo')
demos+=('statistics_demo')
demos+=('string_operations_demo')
demos+=('table_integer_demo')
demos+=('vectors_demo')

for i in "${demos[@]}"; do
    echo 'Compiling:' $i
    $compiler $i/$i.f90 $normal $slow $library -o $i/$i.e
done
