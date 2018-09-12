F2PY ?= f2py
FC = gfortran
ONELOOP_LIB_PATH ?= ../../../HEPTools/oneloop
ONELOOP_INC_PATH ?= ../../../HEPTools/oneloop

all: pyAVH_OneLOop_hook.so

pyAVH_OneLOop_hook.so: AVHOneLOopHook.f90
	touch __init__.py
	$(F2PY) -m pyAVH_OneLOop_hook -c AVHOneLOopHook.f90 --fcompiler=$(FC) -L$(ONELOOP_LIB_PATH) -lavh_olo -I$(ONELOOP_INC_PATH)
