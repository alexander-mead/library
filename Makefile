# Makefile for meadlib

# Set the Fortran compiler
FC = gfortran

# Library name
LIB = meadlib

# Debug library name
LIB_DEBUG = $(LIB)_debug

# Standard Fortran compile flags
FFLAGS = \
	-fcheck=all \
	-fmax-errors=4 \
	-ffpe-trap=invalid,zero,overflow \
	-fimplicit-none \
	-O3 \
	-fdefault-real-8 \
	-fdefault-double-8 \
	-std=gnu \
	-lgfortran \
	-ffree-line-length-none \
	-I/usr/local/include \
	-I/usr/include# \ 
	#-Warray-bounds \
	#-std=f2008 \
	#- pedantic 

# Debug flags
DEBUG_FLAGS = \
	-Wall \
	-fcheck=all \
	-fbounds-check \
	-fbacktrace \
	-Og

# Source-code directory
SRC_DIR = src

# Build directory
BUILD_DIR = build

# All the object files that the library is composed of 
_OBJ = \
	precision.o \
	constants.o \
	physics.o \
	fix_polynomial.o \
	array_operations.o \
	logical_operations.o \
	random_numbers.o \
	calculus.o \
	file_info.o \
	camb_stuff.o \
	table_integer.o \
	interpolate.o \
	special_functions.o \
	calculus_table.o \
	string_operations.o \
	cosmology_functions.o \
	fft.o \
	statistics.o \
	fitting.o \
	gadget.o \
	ode_solvers.o \
	vectors.o \
	orbits.o \
	solve_equations.o \
	sorting.o \
	field_operations.o \
	simulations.o \
	cosmic_emu_stuff.o \
	multidark_stuff.o \
	limber.o \
	hmx.o \
	owls.o \
	owls_extras

# Add prefix of build directory to objects
OBJ = $(addprefix $(BUILD_DIR)/,$(_OBJ))

# Default compile option
#all: FFLAGS += -O3 
all: $(LIB)

# Debug mode
debug: FFLAGS += $(DEBUG_FLAGS)
debug: $(LIB_DEBUG)

# Make the library
meadlib: $(OBJ).o
	@echo
	@$(FC) --version
	@echo 'compiling library'
	@ar rc $(LIB) $(OBJ).o
	@echo 'done'
	@echo

# Make the library
meadlib_debug: $(OBJ).o
	@echo
	@$(FC) --version
	@echo 'compiling debugging library'
	@ar rc $(LIB_DEBUG) $(OBJ).o
	@echo 'done'
	@echo

# Rule to create the object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -c -o $@ $< -J$(BUILD_DIR) $(FFLAGS)

# Clean up
clean:
	rm -f $(LIB)
	rm -f $(LIB_DEBUG)
	rm -f $(BUILD_DIR)/*.mod
	rm -f $(BUILD_DIR)/*.o
