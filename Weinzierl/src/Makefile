# Compilers:                                                                                                                                              
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
  CXX=g++
else ifeq ($(UNAME_S),Darwin)
  CXX=g++
endif

# Compiler flags:
CXXFLAGS = -std=c++11 -O3

SOURCE = DCD_interface.o DCD.o

all: $(SOURCE)

%.o : %.cpp %.h
	$(CXX) $(CXXFLAGS) -fPIC -c $< -o $@
clean :
	-rm -f ${SRC}*.o
