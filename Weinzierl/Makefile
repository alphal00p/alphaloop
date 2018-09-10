# Compilers:                                                                                                                                              
UNAME_S := $(shell uname -s)
PWD := $(shell pwd)

ifeq ($(UNAME_S),Linux)
  CXX=g++
else ifeq ($(UNAME_S),Darwin)
  CXX=g++
endif

# Cuba
CUBA_METHOD = 1 #1:Vegas 2:Suave 3:Divone 4:Cuhre

AND_CUBA_LIB = /home/andrea/Programs/Cuba-4.2/
AND_CUBA_INC = /home/andrea/Programs/Cuba-4.2/

VAL_CUBA_LIB = /Users/valentin/Documents/HEP_softs/Cuba-4.2/lib
VAL_CUBA_INC = /Users/valentin/Documents/HEP_softs/Cuba-4.2/include

# DCD
DCD_INC = $(PWD)/src

export LIB_FLAG = -L$(AND_CUBA_LIB) -L$(VAL_CUBA_LIB)
export INC_FLAG = -I$(AND_CUBA_INC) -I$(VAL_CUBA_INC) -I$(DCD_INC)

# Library
LFLAGS = -lcuba -lboost_program_options

# Compiler flags:
CXXFLAGS = -std=c++11 -O3

SOURCE = src/DCD.o src/RC4vector.o
PY_LIB = src/DCD_interface.o
INTEGRAND = integrand/box1L.o integrand/Integrand.o
INTEGRAND_INTERFACE = integrand/integrand_interface.o src/DCD.o src/RC4vector.o integrand/box1L.o

EXEC = DCD_interface.so

all: $(EXEC)

#####################################
# Executable:
DCDtest: DCD.cpp $(SOURCE)
	$(CXX) $(CXXFLAGS) -D _DCD_MAIN_ $< -o $@

cuba_integrate: cuba.c $(INTEGRAND) $(SOURCE)
	$(CXX) $(CXXFLAGS) $(INTEGRAND) $(SOURCE) $(INC_FLAG) $(LIB_FLAG) -o $@ $< $(LFLAGS)

#####################################
# Shared Library: C Wrap for python
DCD_interface.so: $(PY_LIB) $(SOURCE)
	$(CXX) $(CXXFLAGS) $(PY_LIB) $(SOURCE) -fPIC -shared -o $@ 

integrand_interface.so: $(INTEGRAND_INTERFACE)
	$(CXX) $(CXXFLAGS) $(INTEGRAND_INTERFACE) -fPIC -shared -o $@ 


#####################################
# Create Links:

src/%.o: src/%.cpp src/%.h
	$(MAKE) -C ./src/ 

integrand/%.o: integrand/%.cpp integrand/Integrand.h integrand/integrand_interface.h
	$(MAKE) -C ./integrand/

#####################################
# Clean:
clean :
	$(MAKE) -C ./src/ clean
	$(MAKE) -C ./integrand/ clean
	-rm -f ${SRC}*.o
