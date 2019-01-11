# Makefile for meadlib

# Set the Fortran compiler
FC = gfortran

# Standard Fortran compile flags
FFLAGS = \
	-Warray-bounds \
	-fmax-errors=4 \
	-ffpe-trap=invalid,zero,overflow \
	-fimplicit-none \
	-std=gnu \
	-ffree-line-length-none \
	-I/usr/local/include \
	-I/usr/include

# Debug flags
DEBUG_FLAGS = \
	-Wall \
	-fcheck=all \
	-fbounds-check \
	-fbacktrace \
	-Og

# All the object files that the library is composed of 
OBJS =  precision.o \
	constants.o \
	physics.o \
	logical_operations.o \
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
	numerology.o \
	ode_solvers.o \
	vectors.o \
	orbits.o \
	owls.o \
	solve_equations.o \
	sorting.o \
	string_operations.o \
	field_operations.o \
	simulations.o \
	cosmic_emu_stuff.o \
	multidark_stuff.o \
	Limber.o \
	HMx

# Default compile option
all: FFLAGS += -O3 
all: meadlib

# Debug mode
debug: FFLAGS += $(DEBUG_FLAGS)
debug: meadlib

# Make the library
meadlib: $(OBJS).o
	@echo
	@$(FC) --version
	@echo 'compiling library'
	@ar rc meadlib $(OBJS).o
	@echo 'done'
	@echo

# Rule to create the object file
%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $< 

# Clean up
clean:
	rm meadlib
	rm *.mod
	rm *.o
