#!/bin/bash

#Initial white space
#echo

#Set the compiler and check it is either gfortran or ifort
compiler=$1
if [ -z compiler ]; then
    echo 'Please specify a compiler (gfortran or ifort)'
    exit 1
fi
if [ "$compiler" = "gfortran" ]; then
    echo 'Compiler: ' $compiler
elif [ "$compiler" = "ifort" ]; then
    echo 'Compiler: ' $compiler
else
    echo 'Please specify a compiler (gfortran or ifort)'
    exit 1
fi

#Print compiler version to screen
$compiler --version

#Check the fast-slow option is either 1 or 0
fast_slow=$2
if [ -z $fast_slow ]; then
    echo 'Please specify fast (1) or slow (0) compiling'
    exit 1
fi
if [ $fast_slow -eq 1 ]; then
    echo 'Compiling optimised libraries'
elif [ $fast_slow -eq 0 ]; then
    echo 'Compiling libraries slowly for bug checking'
else
    echo 'Please specify fast (1) or slow (0) compiling'
    exit 1
fi

#Name for eventual library
meadlib='meadlib'

#Remove all the old module files and the old library
rm *.mod
rm $meadlib

#Make an array of all the libraries
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

#Numerical recipes
#nr='/Users/Mead/Physics/numerical_recipes/recipes_compiled'

#Normal or slow gfortran compile options
if [ "$compiler" = "gfortran" ]; then    
    normal='-Warray-bounds -ffree-line-length-none -fmax-errors=4 -ffpe-trap=invalid,zero,overflow -fimplicit-none'   
    slow='-Wall -fcheck=all -fbounds-check -fbacktrace -Og' #Could add -g to generate 'debug symbols' and get .dSYM directory (?)
fi

#Normal or slow ifort compile options
if [ "$compiler" = "ifort" ]; then
    normal=''
    slow=''
fi

#Write out compile flags
echo 'Normal compiler flags:' $normal
if [ $fast_slow -eq 0 ]; then
    echo 'Slow compiler flags:' $slow
fi

#Loop over and compile all libraries
for i in "${libs[@]}"; do
    echo 'Compiling:' $i
    if [ $fast_slow -eq 1 ]; then
	#$compiler $i/$i.f90 $normal -lfftw3 -O3 -c
	$compiler $i.f90 $normal -O3 -c
    else
	#$compiler $i/$i.f90 $normal $slow -lfftw3 -c
	$compiler $i.f90 $normal $slow -c
    fi
done

#Make the big library from all the .o files
ar rc $meadlib *.o

#Remove all the .o files
rm *.o
