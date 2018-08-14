# Compilers:                                                                                                                                              
UNAME_S := $(shell uname -s)
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

LIB_FLAG = -L$(AND_CUBA_LIB) -L$(VAL_CUBA_LIB)
INC_FLAG = -I$(AND_CUBA_INC) -I$(VAL_CUBA_INC)

# Scattering Angle
SCATTERING_ANGLE=M_PI/2.

# Library
LFLAGS = -lcuba

# Compiler flags:
CXXFLAGS = -std=c++11 -O3

EXEC = DCD_interface.so

all: $(EXEC)

#####################################
# Executable:
DCDtest: DCD.cpp 
	$(CXX) $(CXXFLAGS) -D _DCD_MAIN_ $^ -o $@

cuba_integrate: cuba_integrand.o DCD.o cuba_integrate.o
	$(CXX) $(CXXFLAGS) $(INC_FLAG) $(LIB_FLAG) -o $@ $^ $(LFLAGS)

#####################################
# Shared Library:
DCD_interface.so: DCD_interface.o DCD.o 
	$(CXX) $(CXXFLAGS) -fPIC -shared $^ -o $@ 

#####################################
# Create Links:

#cuba integrator
cuba_integrate.o : cuba.c
	$(CXX) $(CXXFLAGS) -D CUBA_WHAT=$(CUBA_METHOD) -c $^ -o $@ $(INC_FLAG)
#function to be integrated
cuba_integrand.o : cuba_integrand.cpp
	$(CXX) $(CXXFLAGS) -D ROTATION_ANGLE=$(SCATTERING_ANGLE) -c $^ -o $@ $(INC_FLAG)
#path deformation for direct integration
DCD.o: DCD.cpp DCD.h
	$(CXX) $(CXXFLAGS) -fPIC -D _PySTUFF_ -c $< -o $@ 
#C Wrap for python
DCD_interface.o: DCD_interface.cpp DCD_interface.h
	$(CXX) $(CXXFLAGS) -fPIC -c $< -o $@ 

clean :
	-rm -f ${SRC}*.o
