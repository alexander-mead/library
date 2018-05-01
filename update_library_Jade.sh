#!/bin/bash

libs+=('precision')
libs+=('constants')
libs+=('fix_polynomial')
libs+=('table_integer')
libs+=('array_operations')
libs+=('ode_solvers')
libs+=('file_info')
libs+=('interpolate')
libs+=('logical_operations')
libs+=('random_numbers')
libs+=('solve_equations')
libs+=('sorting')
libs+=('special_functions')
libs+=('statistics')
libs+=('string_operations')
libs+=('vectors')
libs+=('numerology')
libs+=('calculus')
libs+=('calculus_table')
libs+=('fitting')
libs+=('fft')
libs+=('gadget')
libs+=('owls')
libs+=('field_operations')
libs+=('cosmology_functions')

scp /Users/Mead/Physics/library/compile_library.sh amead@jade.phas.ubc.ca:~/library/.
for lib in "${libs[@]}"; do
    scp -r /Users/Mead/Physics/library/$lib/$lib.f90 amead@jade.phas.ubc.ca:~/library/$lib/.
done
