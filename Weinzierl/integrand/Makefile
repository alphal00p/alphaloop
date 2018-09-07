# Compilers:                                                                                                                                              
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
  CXX=g++
else ifeq ($(UNAME_S),Darwin)
  CXX=g++
endif

# Compiler flags:
CXXFLAGS = -std=c++11 -O3

EXEC = box1L.o Integrand.o integrand_interface.o

all: $(EXEC)

%.o : %.cpp Integrand.h ../src/DCD.h integrand_interface.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INC_FLAG)
clean :
	-rm -f ${SRC}*.o
