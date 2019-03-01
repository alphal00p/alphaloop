F2PY = f2py
FC = gfortran

all: compute_chain gamma_chain.so

gamma_chain.so: compute_chain.f90  
	touch __init__.py
	$(F2PY) -m gamma_chain -c compute_chain.f90 --fcompiler=$(FC)

compute_chain:	compute_chain.f90
	$(FC) -o compute_chain compute_chain.f90 -O5

clean:
	rm -rf *.o compute_chain gamma_chain.so *.dSYM *.pyc
