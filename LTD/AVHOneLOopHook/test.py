#!/usr/bin/env python 
import math
import pyAVH_OneLOop_hook
result = pyAVH_OneLOop_hook.compute_one_loop_box(100.**2,200.**2,300.**2,400.**2,812438.8277451771e0,-406042.8596833937e0,0.,0.,0.,0.)
print "BOX: Single pole obtained: %.16e %s %.16e*i"%(result[1].real, '+' if result[1].real > 0. else '-', abs(result[1].imag))
print "BOX: Double pole obtained: %.16e %s %.16e*i"%(result[2].real, '+' if result[2].real > 0. else '-', abs(result[2].imag))
print "BOX: The poles should be zero."
print "BOX: The finite part of the result obained is: %.16e %s %.16e*i"%(result[0].real, '+' if result[0].real > 0. else '-', abs(result[0].imag))
print "BOX: The finite part of the result should be : %.16e - %.16e*i"%( -5.0195667135415558e-11, 9.2157738446831588e-11) 
result = pyAVH_OneLOop_hook.compute_one_loop_triangle(1,1./4.,1./16.,0.,0.,0.)
print(result)
#result = pyAVH_OneLOop_hook.compute_one_loop_triangle(1.,1.,1.,0.,0.,0.)
# Below is an exceptional point that leads to unstable evaluations...
#result = pyAVH_OneLOop_hook.compute_one_loop_triangle(1.,4.,9.,0.,0.,0.)
print "TRIANGLE: Single pole obtained: %.16e %s %.16e*i"%(result[1].real, '+' if result[1].real > 0. else '-', abs(result[1].imag))
print "TRIANGLE: Double pole obtained: %.16e %s %.16e*i"%(result[2].real, '+' if result[2].real > 0. else '-', abs(result[2].imag))
print "TRIANGLE: The poles should be zero."
print "TRIANGLE: The finite part of the result obained is: %.16e %s %.16e*i"%(result[0].real, '+' if result[0].real > 0. else '-', abs(result[0].imag))
print "TRIANGLE: The finite part of the result should be : %.16e + %.16e*i"%( 7.8701840491118205e+00 , 0.0000000000000000e+00 )

print(pyAVH_OneLOop_hook.compute_one_loop_triangle(-1.,1.,1.,0.,0.,0.))

print("----------------")
#result = pyAVH_OneLOop_hook.compute_one_loop_box(-100.,0.**2,0.**2,0.**2, -100.,-100.,0.,0.,0.,0.)
result = pyAVH_OneLOop_hook.compute_one_loop_triangle(-49.,-148.,-367.,0.0,0.0,0.0)
result = pyAVH_OneLOop_hook.compute_one_loop_box(
    100.0**2,
    300.0**2,
    500.0**2,
    700.0**2,
    10240000.000000000,
    -2826958.6951293740,
    800.**2,200.**2,400.**2,600.**2)
print result
result = pyAVH_OneLOop_hook.compute_one_loop_box(
    100.0**2,
    300.0**2,
    500.0**2,
    700.0**2,
    10240000.000000000,
    -5950922.6967353411,
    800.**2,200.**2,400.**2,600.**2)
print result
result = pyAVH_OneLOop_hook.compute_one_loop_box(
    -10000.000000000000,
    -90000.000000000000,
    -4320649.8702493673,
    -4560649.8702493701,
    -10240000.000000000,
    1880272.8264859747,
    800.**2,200.**2,400.**2,600.**2)
print result 
print('==========')
result = pyAVH_OneLOop_hook.compute_one_loop_box(
    -4869.270805985198,
    -4029.186884648,
    4319.195273527799,
    29764.057148062195,
    8849.472192272002,
    4533.972680586998,
     9.82998**2,9.82998**2,9.82998**2,9.82998**2)
result[0]*=(1.j/(16.*(math.pi**2)))
print result
print '%.16e + %.16e'%(result[0].real, result[0].imag)
print('==========')
result = pyAVH_OneLOop_hook.compute_one_loop_box(
    0.7,
    -0.12,
    -0.37,
    -2.2,
    -0.32,
    -0.93,
     0.,0.,0.,0.)
result[0]*=(1.j/(16.*(math.pi**2)))
print result
print '%.16e + %.16e'%(result[0].real, result[0].imag)
print('==========')

result = pyAVH_OneLOop_hook.compute_one_loop_box(
    3.74,
    -0.33,
    -1.01,
    2.86,
    5.03,
    -0.68,
     0.,0.,0.,0.)
result[0]*=(1.j/(16.*(math.pi**2)))
print result
print '%.16e + %.16e'%(result[0].real, result[0].imag)
print('==========')
result = pyAVH_OneLOop_hook.compute_one_loop_box(
    0.5875,
    -0.2089,
    0.1171,
    -2.4673,
    -1.0874,
    -0.9688,
     0.,0.,0.,0.)
result[0]*=(1.j/(16.*(math.pi**2)))
print result
print '%.16e + %.16e'%(result[0].real, result[0].imag)
print('==========')
result = pyAVH_OneLOop_hook.compute_one_loop_box(
    0.84,
    -2.75,
    0.56,
    -2.61,
    -3.31,
    -3.45,
     0.,0.,0.,0.)
result[0]*=(1.j/(16.*(math.pi**2)))
print result
print '%.16e + %.16e'%(result[0].real, result[0].imag)

print('==========')
result = pyAVH_OneLOop_hook.compute_one_loop_box(
    -0.29,
    -0.12,
    -0.37,
    -2.2,
    -0.77,
    -0.93,
     0.,0.,0.,0.)
result[0]*=(1.j/(16.*(math.pi**2)))
print result
print '%.16e + %.16e'%(result[0].real, result[0].imag)

print('==========')
result = pyAVH_OneLOop_hook.compute_one_loop_triangle(
    -0.29,
    -0.12,
    -0.77,
     0.,0.,0.)
result[0]*=(1.j/(16.*(math.pi**2)))
print result
print '%.16e + %.16e'%(result[0].real, result[0].imag)
print('==========')
result = pyAVH_OneLOop_hook.compute_one_loop_triangle(
    0.7,
    -0.12,
    -0.32,
     0.,0.,0.)
result[0]*=(1.j/(16.*(math.pi**2)))
print result
print '%.16e + %.16e'%(result[0].real, result[0].imag)
print('==========')
result = pyAVH_OneLOop_hook.compute_one_loop_triangle(
    0.72308431,
    -0.10124341,
    -0.55763794,
     0.,0.,0.)
result[0]*=(1.j/(16.*(math.pi**2)))
print result
print '%.16e + %.16e'%(result[0].real, result[0].imag)


print('==========')
result = pyAVH_OneLOop_hook.compute_one_loop_box(
    0.79462243,
    -0.08114367,
    -0.38902201,
    -1.95377605,
    -0.55985718,
    -0.99810302,
     0., 0.,0.,0.)
result[0]*=(1.j/(16.*(math.pi**2)))
print result
print '%.16e + %.16e'%(result[0].real, result[0].imag)


print('==========')
result = pyAVH_OneLOop_hook.compute_one_loop_box(
    0.7,
    -0.12,
    -0.37,
    -1.57,
    -0.32,
    -0.93,
     0., 0.,0.,0.)
result[0]*=(1.j/(16.*(math.pi**2)))
print result
print '%.16e + %.16e'%(result[0].real, result[0].imag)


print('==========')
result = pyAVH_OneLOop_hook.compute_one_loop_box(
    0.84,
    -2.75,
    0.56,
    -2.61,
    -3.31,
    -3.45,
     0., 0.,0.,0.)
result[0]*=(1.j/(16.*(math.pi**2)))
print result
print '%.16e + %.16e'%(result[0].real, result[0].imag)

print('==========')
result = pyAVH_OneLOop_hook.compute_one_loop_box(
    0.5875,
    -0.2089,
    0.1171,
    -0.7573,
    -1.0874,
    -0.9688,
     0., 0.,0.,0.)
result[0]*=(1.j/(16.*(math.pi**2)))
print result
print '%.16e + %.16e'%(result[0].real, result[0].imag)


print('==========')
result = pyAVH_OneLOop_hook.compute_one_loop_box(
    0.01,
    0.01,
    0.01,
    0.01,
    4.0,
    -1.1897703974737204,
     0., 0.,0.,0.)
result[0]*=(1.j/(16.*(math.pi**2)))
print result
print '%.16e + %.16e'%(result[0].real, result[0].imag)


print('==========')
result = pyAVH_OneLOop_hook.compute_one_loop_triangle(
    -49.0,
    -148.0,
     -367.0,
    0.,0.,0.)
result[0]*=(1.j/(16.*(math.pi**2)))
print result
print '%.16e + %.16e'%(result[0].real, result[0].imag)

print('==========')
result = pyAVH_OneLOop_hook.compute_one_loop_triangle(
     11.0,
     23.0,
     70.0,
    0.,0.,0.)
result[0]*=(1.j/(16.*(math.pi**2)))
print result
print '%.16e + %.16e'%(result[0].real, result[0].imag)
