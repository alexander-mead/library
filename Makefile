# Makefile for my library

# Set the Fortran compiler
FC = gfortran

# Library name
LIB = libmead.a

# Debug library name
LIB_DEBUG = libmead_debug.a

# Standard Fortran compile flags
FFLAGS = \
	-fimplicit-none \
	-std=gnu \
	-lgfortran \
	-fmax-errors=4 \
	-fdefault-real-8 \
	-fdefault-double-8 \
	-ffree-line-length-none \
	-I/usr/local/include \
	-I/usr/include

# Additional flags for standard compilation
FLAGS_ALL = \
	-fcheck=all \
	-ffpe-trap=invalid,zero,overflow \
	-O3

# Additional flags for debug compilation
FLAGS_DEBUG = \
	-Wall \
	-fcheck=all \
	-ffpe-trap=invalid,zero,overflow \
	-fbounds-check \
	-fbacktrace \
	-Og \
	#-pedantic \
	#-Warray-bounds \
	#-std=f2008

# Source-code directory
SRC_DIR = src

# Build directory
BUILD_DIR = build

# Debug build directory
BUILD_DIR_DEBUG = build_debug

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

# Default compile option
all: FFLAGS += $(FLAGS_ALL)
all: $(LIB)

# Debug mode
debug: FFLAGS += $(FLAGS_DEBUG)
debug: $(LIB_DEBUG)

# Add prefix of build directory to objects
OBJ = $(addprefix $(BUILD_DIR)/,$(_OBJ))
OBJ_DEBUG = $(addprefix $(BUILD_DIR_DEBUG)/,$(_OBJ))

# Make the library
$(LIB): $(OBJ).o
	@echo
	@$(FC) --version
	@echo 'compiling library'
	@ar rc $(LIB) $(OBJ).o
	@echo 'done'
	@echo

# Make the debugging library
$(LIB_DEBUG): $(OBJ_DEBUG).o
	@echo
	@$(FC) --version
	@echo 'compiling debugging library'
	@ar rc $(LIB_DEBUG) $(OBJ_DEBUG).o
	@echo 'done'
	@echo

# Rule to create the object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -c -o $@ $< -J$(BUILD_DIR)

# Rule to create the object files
$(BUILD_DIR_DEBUG)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -c -o $@ $< -J$(BUILD_DIR_DEBUG)

# Clean up
clean:
	rm -f $(LIB)
	rm -f $(BUILD_DIR)/*.mod
	rm -f $(BUILD_DIR)/*.o
	rm -f $(LIB_DEBUG)
	rm -f $(BUILD_DIR_DEBUG)/*.mod
	rm -f $(BUILD_DIR_DEBUG)/*.o
