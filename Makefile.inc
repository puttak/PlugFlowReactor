# Path to blas/lapack library if used to compile Cantera.
BLAS_LAPACK   := -lopenblas

# Path to Cantera install where include/lib directories are found.
# If include/lib are not in same path, manually edit the next variables.
ROOT_CANTERA  := ../CanteraPFR/external/Nix
LIB_CANTERA   := -L$(ROOT_CANTERA)/lib
INC_CANTERA   := -I$(ROOT_CANTERA)/include

# Path to Sundials install where include/lib directories are found.
# If include/lib are not in same path, manually edit the next variables.
ROOT_SUNDIALS := $(ROOT_CANTERA)
LIB_SUNDIALS  := -L$(ROOT_SUNDIALS)/lib
INC_SUNDIALS  := -I$(ROOT_SUNDIALS)/include

# Compiler/archiver and options.
CXX           := g++
CXXFLAGS      := -g -O3 -std=c++14
AR            := ar
ARFLAGS       := sq

# Python distribution.
PYTHON        := python
