# ----------------------------------------------------------------------------
#
# Makefile for model library
#
# ----------------------------------------------------------------------------

# Check for ../make_opts
ifeq ($(wildcard ../make_opts), ../make_opts)
  include ../make_opts
else
  FFLAGS+= -ffixed-line-length-132
  FC=gfortran
endif

include makeinc.inc

LIBDIR=../../lib/
LIBRARY=libmodel.$(libext)

ALL: $(LIBDIR)$(LIBRARY) $(LIBDIR)$(patsubst %.a,%_dynamic.so,$(LIBRARY))

helas_couplings: helas_couplings.o $(LIBRARY)
	$(FC) $(FFLAGS) -o $@ $^

testprog: testprog.o $(LIBRARY)
	$(FC) $(FFLAGS) -o $@ $^

$(LIBRARY): $(MODEL)
	ar cru $(LIBRARY) $(MODEL)
	ranlib $(LIBRARY)

$(LIBDIR)$(LIBRARY): $(MODEL)
	$(call CREATELIB, $@, $^)

$(LIBDIR)$(patsubst %.a,%_dynamic.so,$(LIBRARY)): $(MODEL)
	$(call CREATEDYNLIB, $@, $^)

clean: 
	$(RM) *.o $(LIBDIR)$(LIBRARY) $(LIBDIR)$(patsubst %.a,%_dynamic.so,$(LIBRARY))

