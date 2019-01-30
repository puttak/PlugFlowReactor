include Makefile.inc

ROOT_SELF     := PlugFlowReactor
LIB_SELF      :=
INC_SELF      := -I$(ROOT_SELF)/include

LIB_CANTERA   += -lcantera_shared
LIB_SUNDIALS  += -lsundials_ida -lsundials_nvecserial

LIBRARIES     := $(LIB_CANTERA) $(LIB_SUNDIALS) $(LIB_SELF) $(BLAS_LAPACK)
INCLUDES      := $(INC_CANTERA) $(INC_SUNDIALS) $(INC_SELF)

LIBRARIES += -pthread
CPPFLAGS  := $(OPTIONS) $(INCLUDES)
LIBFLAGS  := -shared -fPIC $(CPPFLAGS)

LIB_PFR   := $(ROOT_SELF)/libPlugFlowReactor_shared.so
SRC_LIB   := $(ROOT_SELF)/src/cPlugFlowReactor.c
SRC_TEST  := $(ROOT_SELF)/src/test.cpp

RM        := rm -rf

all: shared_library
	$(PYTHON) setup.py install

shared_library: $(SRC_LIB)
	$(CXX) $(LIBFLAGS) $(SRC_LIB) $(LIBRARIES) -o $(LIB_PFR)

test: $(SRC_TEST)
	$(CXX) $(CPPFLAGS) $(SRC_TEST) $(LIBRARIES) -o test.exe


clean:
	$(RM) *.exe *.so
