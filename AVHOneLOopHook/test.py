#!/usr/bin/env python 
import pyAVH_OneLOop_hook
print dir(pyAVH_OneLOop_hook)
result = pyAVH_OneLOop_hook.compute_one_loop_box(100.**2,200.**2,300.**2,400.**2,812438.8277451771e0,-406042.8596833937e0,0.,0.,0.,0.)
print "Single pole obtained: %.16e %s %.16e*i"%(result[1].real, '+' if result[1].real > 0. else '-', abs(result[1].imag))
print "Double pole obtained: %.16e %s %.16e*i"%(result[2].real, '+' if result[2].real > 0. else '-', abs(result[2].imag))
print "The poles should be zero."
print "The finite part oof the result obained is: %.16e %s %.16e*i"%(result[0].real, '+' if result[0].real > 0. else '-', abs(result[0].imag))
print "The finite part of the result should be  : %.16e - %.16e*i"%( -5.0195667135415558e-11, 9.2157738446831588e-11) 
