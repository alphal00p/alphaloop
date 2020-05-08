ALL: test_pol_sum

FFLAGS= -O -w -fbounds-check -fPIC -ffixed-line-length-132

test_pol_sum: spinor_polarization_sum.f ../../Templates/aloha_functions.f
	gfortran $(FFLAGS) -o test_pol_sum spinor_polarization_sum.f ../../Templates/aloha_functions.f
