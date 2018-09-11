# Compilers:                                                                                                                                              
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
  CXX=g++
else ifeq ($(UNAME_S),Darwin)
  CXX=g++
endif

# Compiler flags:
CXXFLAGS = -std=c++11 -O3 -fPIC 

EXEC = box1L.o integrand_interface.o

all: $(EXEC)

integrand_interface.o : integrand_interface.cpp ../src/DCD.h integrand_interface.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INC_FLAG)

%.o : %.cpp ../src/DCD.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INC_FLAG)
clean :
	-rm -f ${SRC}*.o
