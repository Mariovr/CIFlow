# Makefile for the CIFlow project tests
#

# Rudimentary platform support.

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S), Linux)
	ifeq ($(VSC_INSTITUTE_LOCAL), gent)
	CXX=mpiicpc
	MPIOPT=-DMPICH_IGNORE_CXX_SEEK -I/usr/lib/openmpi/include
	LIBS=-pthread -L$(EBROOTHDF5)/lib -L$(EBROOTBOOST)/lib -lhdf5 -larpack -lboost_unit_test_framework
else
	#CXX=mpiicpc
	CXX = icpc
	#CXX = g++
	MPIOPT=-I/usr/lib/openmpi/include -pthread 
	LIBS=  -L/usr/lib -lhdf5 -lmpi -lmpi_cxx -lboost_unit_test_framework -lm -larpack
endif
endif

ifeq ($(UNAME_S), Darwin)
	CXX = g++
	LIBS= -lm -L/opt/local/lib -lhdf5 -larpack -lboost_unit_test_framework-mt
	MPIOPT=-pthread 
endif

SRC = $(patsubst src/%.cpp, %.cpp, $(wildcard src/*.cpp))
BIN = $(patsubst %.cpp, bin/%.x, $(SRC))
INCLUDE_PATH = -I../include -I/opt/local/include 
LIB_PATH = -L../lib 
CXX_FLAGS = -Wall -Wextra -std=c++11 -g -march=native $(INCLUDE_PATH) $(MPIOPT)
CXX_FLAGS += -D_DEBUG

ifneq ($(CXX) , g++) 
	CXX_FLAGS += -D__INTEL_COMPILER__ -D__OMP__ -openmp 
	LIBS += -mkl -openmp 
else
	LIBS += -llapack -lblas #-fopenmp
	#CXX_FLAGS += -std=libc++ -fopenmp
endif

LINKER_FLAGS = $(LIB_PATH) -lciflow $(LIBS)
RM = rm -f

MKDIR_P = mkdir -p
RMDIR = rm -rf 

all: directories $(BIN)

directories :
		${RMDIR} bin
	        ${MKDIR_P} bin

$(BIN) : bin/%.x : src/%.cpp
	$(CXX) $(CXX_FLAGS) -o $@ $^ $(LINKER_FLAGS)

mp : 
	$(CXX) $(CXX_FLAGS) -o bin/Testmp.x src/Testmp.cpp $(LINKER_FLAGS)

density : 
	$(CXX) $(CXX_FLAGS) -o bin/TestDensity.x src/TestDensity.cpp $(LINKER_FLAGS)

transform : 
	$(CXX) $(CXX_FLAGS) -o bin/TestTransformations.x src/TestTransformations.cpp $(LINKER_FLAGS)

ham: 
	$(CXX) $(CXX_FLAGS) -o bin/TestHam.x src/TestHam.cpp $(LINKER_FLAGS)

sub:
	$(CXX) $(CXX_FLAGS) -o bin/TestSubotnik.x src/TestSubotnik.cpp $(LINKER_FLAGS)

sim:
	$(CXX) $(CXX_FLAGS) -o bin/TestSimulatedAnnealing.x src/TestSimulatedAnnealing.cpp $(LINKER_FLAGS)

scf:
	$(CXX) $(CXX_FLAGS) -o bin/TestRHF.x src/TestRHF.cpp $(LINKER_FLAGS)

er:
	$(CXX) $(CXX_FLAGS) -o bin/TestER.x src/TestER.cpp $(LINKER_FLAGS)

sub_diis:
	$(CXX) $(CXX_FLAGS) -o bin/TestSubotnikDIIS.x src/TestSubotnikDIIS.cpp $(LINKER_FLAGS)

redbcs:
	$(CXX) $(CXX_FLAGS) -o bin/TestRedBCS.x src/TestRedBCS.cpp $(LINKER_FLAGS)

clean :
	$(RM) bin/*.x
