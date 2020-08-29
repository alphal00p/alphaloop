include ../Source/make_opts

ifeq ($(origin GCC),undefined)
GCC=gcc
endif
LIBBPATH=../lib
SRC=$(wildcard numerator*.c)
OBJ=$(patsubst %.c,%.o, $(SRC))
INTEGRANDSRC=$(wildcard integrand*.c)
INTEGRANDOBJ=$(patsubst %.c,%.o, $(INTEGRANDSRC))
LIBRARY=libFORM_numerators.so
INTEGRANDLIBRARY=libFORM_integrands.so

all: $(LIBBPATH)/$(LIBRARY) $(LIBBPATH)/$(INTEGRANDLIBRARY)
numerator: $(LIBBPATH)/$(LIBRARY)
integrand: $(LIBBPATH)/$(INTEGRANDLIBRARY)

OPTIMIZATION_LVL = 3
CFLAGS=-Wconversion -O$(OPTIMIZATION_LVL) -fcx-fortran-rules -fcx-limited-range -lm -lquadmath

OPTIMIZATION_LVL_f128 = 0
CFLAGSf128=-Wconversion -O$(OPTIMIZATION_LVL_f128) -fcx-fortran-rules -fcx-limited-range -lm -lquadmath

%_f64.o: %_f64.c
	$(GCC) -c $(CFLAGS) $< -o $@

%_f128.o: %_f128.c
	$(GCC) -c $(CFLAGSf128) $< -o $@

$(LIBBPATH)/$(LIBRARY): $(OBJ)
	$(GCC) --shared -fPIC $(CFLAGS) -o $@ $^

$(LIBBPATH)/$(INTEGRANDLIBRARY): $(INTEGRANDOBJ)
	$(GCC) --shared -fPIC $(CFLAGS) -o $@ $^

clean:
	rm -f $(LIBBPATH)/$(LIBRARY)
	rm -f $(LIBBPATH)/$(INTEGRANDLIBRARY)
