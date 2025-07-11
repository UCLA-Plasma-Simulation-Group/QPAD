##### compiler configuration #####
# This template configuration file is based on the GNU compiler series.

# ------------------------------------------------------------------------------
# Compiler configuration
# ------------------------------------------------------------------------------
# Specify the MPI-Fortran compiler wrapper. Frequently used wrapper includes
# "mpif90" for GNU family. The name could vary on different system so one should
# seek the user guide for the wrapper name.
FC = mpif90

# Configure the compiler options
FC_OPTS = -g -O3 -fdefault-real-8 -fdefault-double-8 -fopenmp -ffree-form

# Enable the generation of run-time checks. The argument shall be a comma-delimited 
# list of the following keywords:
# 'all': Enable all run-time test of -fcheck.
# 'array-temps': Warns at run time when for passing an actual argument a temporary array had to be generated.
# 'bits': Enable generation of run-time checks for invalid arguments to the bit manipulation intrinsics.
# 'bounds': Enable generation of run-time checks for array subscripts and against the declared minimum and maximum values.
# 'do': Enable generation of run-time checks for invalid modification of loop iteration variables.
# 'mem': Enable generation of run-time checks for memory allocation.
# 'pointer': Enable generation of run-time checks for pointers and allocatables.
# 'recursion': Enable generation of run-time checks for recursively called subroutines and functions which are not marked as recursive.
# FC_OPTS += -fcheck=all

# Specify a list of floating point exception traps to enable. The comma-separated
# list of the following exceptions: 
# ‘invalid’: invalid floating point operation, such as SQRT(-1.0).
# ‘zero’: division by zero.
# ‘overflow’: overflow in a floating point operation.
# ‘underflow’: underflow in a floating point operation.
# ‘inexact’: loss of precision during operation.
# ‘denormal’: operation performed on a denormal value.
# FC_OPTS += -ffpe-trap=invalid,zero,overflow

# Specify warning options. -Wall enables commonly used warning options including:
# -Waliasing, -Wampersand, -Wconversion, -Wsurprising, -Wc-binding-type, -Wintrinsics-std, 
# -Wtabs, -Wintrinsic-shadow, -Wline-truncation, -Wtarget-lifetime, -Winteger-division,
# -Wreal-q-constant, -Wunused and -Wundefined-do-loop.
# FC_OPTS += -Wall 

# C compiler configuration
CC = gcc
CC_OPTS = -O3 -std=c99

# Preprocessor configuration
FPP = gcc -C -E -cpp

# ------------------------------------------------------------------------------
# Configure the modules that will NOT be installed
# ------------------------------------------------------------------------------

# DISABLE_TEMPLATE = true
# DISABLE_POPAS = true

# ------------------------------------------------------------------------------
# Configure the linker
# ------------------------------------------------------------------------------
LINKER = mpif90
LINKER_OPTS = -fopenmp -O3
LDF = -L$(HYPRE_LIB) -L$(JSON_LIB) -L$(HDF5_LIB) -lHYPRE -ljsonfortran -lhdf5_fortran -lhdf5

# ------------------------------------------------------------------------------
# Specify the path of each libraries
# ------------------------------------------------------------------------------
HYPRE_LIB = /path_to_hypre_library
JSON_LIB = /path_to_json_fortran_library
JSON_INC = /path_to_json_fortran_head_files
HDF5_LIB = /path_to_hdf5_library
HDF5_INC = /path_to_hdf5_head_files


INCF = -I$(HDF5_INC) -I$(JSON_INC)
