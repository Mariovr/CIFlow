# Makefile for the CIFlow project
# 

# Rudimentary platform support.

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S), Linux)
	ifeq ($(VSC_INSTITUTE_LOCAL), gent)
		CXX=mpiicpc
		MPIOPT=-DMPICH_IGNORE_CXX_SEEK -I/usr/lib/openmpi/include
		LIBS= -lm -pthread -L$(EBROOTHDF5)/lib -lhdf5 -larpack -L$(EBROOTBOOST)/lib 
	else
		CXX=icpc
		#CXX = g++
		MPIOPT=-I/usr/lib/openmpi/include -pthread 
		LIBS=  -lm -lhdf5 -L/usr/lib -lmpi -lmpi_cxx -larpack
	endif
endif

ifeq ($(UNAME_S), Darwin)
	CXX = g++
	LIBS= -lm -L/opt/local/lib/ -lhdf5 -larpack
	MPIOPT=-pthread 
endif

SRC = $(patsubst src/%.cpp, %.cpp, $(wildcard src/*.cpp))
OBJ_SRC := $(filter-out Main.cpp, $(SRC))
OBJ = $(patsubst %.cpp, obj/%.o, $(OBJ_SRC))
INCLUDE_PATH = -I./include -I/opt/local/include 
LIB_PATH = -L./lib 
CXX_FLAGS = -Wall -Wextra -std=c++11 -march=native $(INCLUDE_PATH) $(MPIOPT) -O3 -fpic
#CXX_FLAGS += -g -ggdb3 #-D_DEBUG

ifneq ($(CXX) , g++) 
	CXX_FLAGS += -D__INTEL_COMPILER__ -D__OMP__ -openmp 
	LIBS += -mkl -openmp 
else
	LIBS += -llapack -lblas #-fopenmp
	#CXX_FLAGS += -fopenmp
endif

LINKER_FLAGS = $(LIB_PATH) -lciflow $(LIBS)
DOCGEN=doxygen
DOCFILE=Doxyfile

all : library program doc tests

.PHONY : library
library : $(OBJ)
	@echo "Building library..."
	@if [ ! -d "lib" ]; then \
		mkdir -p lib; \
	fi
	@ar rcs lib/libciflow.a $^

.PHONY : program
program : library
	@echo "Building program..."
	@if [ ! -d "bin" ]; then \
		mkdir -p bin; \
	fi
	@$(CXX) $(CXX_FLAGS) -o bin/ciflow.x src/Main.cpp -v $(LINKER_FLAGS)

$(OBJ) : obj/%.o : src/%.cpp
	@if [ ! -d "obj" ]; then \
		mkdir -p obj; \
	fi
	$(CXX) $(CXX_FLAGS) -c $< -o $@

.PHONY : doc
doc :
	@echo "Building documentation..."
	@ cd doc && $(DOCGEN) $(DOCFILE)

.PHONY : tests
tests :
	@echo "Building tests..."
	@cd tests && ./run.sh

# If you have PSI4.
.PHONY: mointegrals
mointegrals:
	cd ./mointegrals && $(MAKE)

.PHONY : clean
clean:
	$(RM) obj/*.o
	$(RM) lib/*.a
	$(RM) bin/*.x
