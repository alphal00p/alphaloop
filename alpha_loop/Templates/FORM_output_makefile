include ../Source/make_opts

ifeq ($(origin GCC),undefined)
GCC=gcc
endif
LIBBPATH=../lib
SRC=$(wildcard numerator*.c)
OBJ=$(patsubst %.c,%.o, $(SRC))
LIBRARY=libFORM_numerators.so

all: $(LIBBPATH)/$(LIBRARY)

%.o: %.c
	$(GCC) -c -O3 -fcx-fortran-rules -fcx-limited-range -lm $< -o $@

$(LIBBPATH)/$(LIBRARY): $(OBJ)
	$(GCC) --shared -fPIC -O3 -fcx-fortran-rules -fcx-limited-range -lm -o $@ $^

clean:
	rm -f $(LIBBPATH)/$(LIBRARY)
