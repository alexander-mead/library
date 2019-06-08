#!/bin/bash

#Set the compiler
compiler=gfortran

#Set the code to compile
code=random_distributions_test

#Location of my libraries
mead=/Users/Mead/Physics/library

#Set the precision
precision='-fdefault-real-8'

$compiler $mead/constants.f90 $mead/special_functions.f90 $mead/random_numbers.f90 $code.f90 $precision -o $code.e

rm *.mod
