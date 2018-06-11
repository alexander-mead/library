# Makefile for meadlib

# Set the Fortran compiler
FC = gfortran

# Default Fortran compile flags
FFLAGS = -Warray-bounds -fmax-errors=4 -ffpe-trap=invalid,zero,overflow -fimplicit-none -std=gnu -ffree-line-length-none

# Debug flags
DEBUG_FLAGS = -Wall -fcheck=all -fbounds-check -fbacktrace -Og

# All the object files that the library is composed of 
OBJS =  precision.o \
	constants.o \
	fix_polynomial.o \
	array_operations.o \
	calculus.o \
	file_info.o \
	camb_stuff.o \
	random_numbers.o \
	table_integer.o \
	interpolate.o \
	special_functions.o \
	calculus_table.o \
	cosmology_functions.o \
	fft.o \
	statistics.o \
	fitting.o \
	gadget.o \
	logical_operations.o \
	numerology.o \
	ode_solvers.o \
	vectors.o \
	orbits.o \
	owls.o \
	solve_equations.o \
	sorting.o \
	string_operations.o \
	field_operations.o \
	simulations

# Default compile option
all: FFLAGS += -O3
all: meadlib

# Debug mode
debug: FFLAGS += $(DEBUG_FLAGS)
debug: meadlib

# Make the library
meadlib: $(OBJS).o
	$(FC) --version
	ar rc meadlib $(OBJS).o

array_operations.o: fix_polynomial.mod array_operations.f90
	$(FC) $(FFLAGS) -c array_operations.f90

calculus.o: calculus.f90
	$(FC) $(FFLAGS) -c calculus.f90

calculus_table.o: calculus_table.f90
	$(FC) $(FFLAGS) -c calculus_table.f90

camb_stuff.o: file_info.mod camb_stuff.f90
	$(FC) $(FFLAGS) -c camb_stuff.f90

constants.o: constants.f90
	$(FC) $(FFLAGS) -c constants.f90

cosmology_functions.o: cosmology_functions.f90
	$(FC) $(FFLAGS) -c cosmology_functions.f90

fft.o: fft.f90
	$(FC) $(FFLAGS) -c fft.f90

field_operations.o: field_operations.f90
	$(FC) $(FFLAGS) -c field_operations.f90

file_info.o: file_info.f90
	$(FC) $(FFLAGS) -c file_info.f90

fitting.o: fitting.f90
	$(FC) $(FFLAGS) -c fitting.f90

fix_polynomial.o: fix_polynomial.f90
	$(FC) $(FFLAGS) -c fix_polynomial.f90

gadget.o: gadget.f90
	$(FC) $(FFLAGS) -c gadget.f90

interpolate.o: interpolate.f90
	$(FC) $(FFLAGS) -c interpolate.f90

logical_operations.o: logical_operations.f90
	$(FC) $(FFLAGS) -c logical_operations.f90

numerology.o: numerology.f90
	$(FC) $(FFLAGS) -c numerology.f90

ode_solvers.o: ode_solvers.f90
	$(FC) $(FFLAGS) -c ode_solvers.f90

orbits.o: orbits.f90
	$(FC) $(FFLAGS) -c orbits.f90

owls.o: owls.f90
	$(FC) $(FFLAGS) -c owls.f90

precision.o: precision.f90
	$(FC) $(FFLAGS) -c precision.f90

random_numbers.o: random_numbers.f90
	$(FC) $(FFLAGS) -c random_numbers.f90

simulations.o: simulations.f90
	$(FC) $(FFLAGS) -c simulations.f90

solve_equations.o: solve_equations.f90
	$(FC) $(FFLAGS) -c solve_equations.f90

sorting.o: sorting.f90
	$(FC) $(FFLAGS) -c sorting.f90

special_functions.o: special_functions.f90
	$(FC) $(FFLAGS) -c special_functions.f90

statistics.o: statistics.f90
	$(FC) $(FFLAGS) -c statistics.f90

string_operations.o: string_operations.f90
	$(FC) $(FFLAGS) -c string_operations.f90

table_integer.o: table_integer.f90
	$(FC) $(FFLAGS) -c table_integer.f90

vectors.o: vectors.f90
	$(FC) $(FFLAGS) -c vectors.f90

# Clean up
clean:
	rm meadlib
	rm *.mod
	rm *.o
