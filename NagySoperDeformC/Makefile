CXX=g++ -std=c++11 -O3

#Cuba
CUBA_METHOD=1 #1:Vegas 2:Suave 3:Divone 4:Cuhre
MY_LIB=/Users/valentin/Documents/HEP_softs/Cuba-4.2/lib
MY_INC=/Users/valentin/Documents/HEP_softs/Cuba-4.2/include

#Scattering Angle
SCATTERING_ANGLE=M_PI/2.

# Library
LIB=-lcuba

# Compiler flags:
CXXFLAGS=

EXEC=cuba_integrate

all: $(EXEC)

#####################################
# Link together:

#Create values out of trees over finitefield
cuba_integrate: cuba_integrand.o DI_PathDef.o cuba_integrate.o
	$(CXX) $(CXXFLAGS) -D CUBA_WHAT=$CUBA_METHOD -I$(MY_INC) -L$(MY_LIB) -o $@  $(LIB) $^

#####################################
# Create Links:

#cuba integrator
cuba_integrate.o : cuba.c
	$(CXX) $(CXXFLAGS) -c $^ -o $@ -I$(MY_INC)
#function to be integrated
cuba_integrand.o : cuba_integrand.cpp
	$(CXX) $(CXXFLAGS) -D ROTATION_ANGLE=$(SCATTERING_ANGLE) -c $^ -o $@ -I$(MY_INC)
#path deformation for direct integration
DI_PathDef.o : DI_PathDef.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@ 

clean :
	-rm -f ${SRC}*.o
