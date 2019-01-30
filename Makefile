PYTHON        := python
OPTIONS       := -g -O3 -std=c++14
BLAS_LAPACK   := -lopenblas
ROOT_SELF     := PlugFlowReactor
ROOT_CANTERA  := ../CanteraPFR/external/Nix
ROOT_SUNDIALS := $(ROOT_CANTERA)
LIB_SUNDIALS  := -L$(ROOT_SUNDIALS)/lib -lsundials_ida -lsundials_nvecserial
INC_SUNDIALS  := -I$(ROOT_SUNDIALS)/include
LIB_CANTERA   := -L$(ROOT_CANTERA)/lib -lcantera_shared
INC_CANTERA   := -I$(ROOT_CANTERA)/include
LIB_SELF      :=
INC_SELF      := -I$(ROOT_SELF)/include
LIBRARIES     := $(LIB_CANTERA) $(LIB_SUNDIALS) $(LIB_SELF) $(BLAS_LAPACK)
INCLUDES      := $(INC_CANTERA) $(INC_SUNDIALS) $(INC_SELF)

LIBRARIES += -pthread
CPPFLAGS  := $(OPTIONS) $(INCLUDES)
LIBFLAGS  := -shared -fPIC $(CPPFLAGS)

LIB_PFR   := $(ROOT_SELF)/libPlugFlowReactor_shared.so
SRC_LIB   := $(ROOT_SELF)/src/cPlugFlowReactor.c
SRC_TEST  := $(ROOT_SELF)/src/test.cpp

CXX       := g++
AR        := ar
ARFLAGS   := sq
RM        := rm -rf

all: shared_library
	$(PYTHON) setup.py install

test: $(SRC_TEST)
	$(CXX) $(CPPFLAGS) $(SRC_TEST) $(LIBRARIES) -o test.exe

shared_library: $(SRC_LIB)
	$(CXX) $(LIBFLAGS) $(SRC_LIB) $(LIBRARIES) -o $(LIB_PFR)

clean:
	$(RM) *.exe *.so
