#!/usr/bin/env python 
import pyAVH_OneLOop_hook
print "Offshell box computation"
print "------------------------"
result = pyAVH_OneLOop_hook.compute_one_loop_box(100.**2,200.**2,300.**2,400.**2,812438.8277451771e0,-406042.8596833937e0,0.,0.,0.,0.)
print "Single pole obtained: %.16e %s %.16e*i"%(result[1].real, '+' if result[1].real > 0. else '-', abs(result[1].imag))
print "Double pole obtained: %.16e %s %.16e*i"%(result[2].real, '+' if result[2].real > 0. else '-', abs(result[2].imag))
print "The poles should be zero."
print "The finite part oof the result obained is: %.16e %s %.16e*i"%(result[0].real, '+' if result[0].imag > 0. else '-', abs(result[0].imag))
print "The finite part of the result should be  : %.16e - %.16e*i"%( -5.0195667135415558e-11, 9.2157738446831588e-11)


class vector(list):
    def square(self):
        return self[0]**2-self[1]**2-self[2]**2-self[3]**2
    def __add__(self, other):
        return vector(self[i]+other[i] for i in range(4))

# Testing the onshell subtracted box
P_point = [
    vector([-4.8678217388064405e+02,    -0.0000000000000000e+00,   -0.0000000000000000e+00,  -4.8678217388064405e+02]),
    vector([-4.8678217388064405e+02,    -0.0000000000000000e+00,   -0.0000000000000000e+00,   4.8678217388064405e+02]),
    vector([ 4.8678217388064405e+02,    1.9365322696179936e+02 ,   1.1431607376733305e+02 ,  -4.3172577844468481e+02]),
    vector([ 4.8678217388064405e+02,   -1.9365322696179936e+02 ,  -1.1431607376733305e+02 ,   4.3172577844468481e+02]),
]

box_res = list(pyAVH_OneLOop_hook.compute_one_loop_box(
    P_point[0].square(),P_point[1].square(),P_point[2].square(),P_point[3].square(),
    (P_point[0]+P_point[1]).square(),(P_point[1]+P_point[2]).square(),0.,0.,0.,0.))

s = (P_point[2]+P_point[3]).square()
t = (P_point[2]+P_point[1]).square()
triangles_res = [
    [(1./s)*p for p in pyAVH_OneLOop_hook.compute_one_loop_triangle(
        P_point[1].square(),P_point[2].square(),(P_point[0]+P_point[3]).square(),0.,0.,0.)],
    [(1./s)*p for p in pyAVH_OneLOop_hook.compute_one_loop_triangle(
        P_point[0].square(),P_point[3].square(),(P_point[1]+P_point[2]).square(),0.,0.,0.)],
    [(1./t)*p for p in pyAVH_OneLOop_hook.compute_one_loop_triangle(
        P_point[0].square(),P_point[1].square(),(P_point[2]+P_point[3]).square(),0.,0.,0.)],
    [(1./t)*p for p in pyAVH_OneLOop_hook.compute_one_loop_triangle(
        P_point[2].square(),P_point[3].square(),(P_point[0]+P_point[1]).square(),0.,0.,0.)],
]
print box_res
for tr in triangles_res:
    print tr
result = [(box_res[i]-sum(tr[i] for tr in triangles_res)) for i in range(3)]

print "Subbtracted onshell box computation"
print "-----------------------------------"
print "Single pole obtained: %.16e %s %.16e*i"%(result[1].real, '+' if result[1].real > 0. else '-', abs(result[1].imag))
print "Double pole obtained: %.16e %s %.16e*i"%(result[2].real, '+' if result[2].real > 0. else '-', abs(result[2].imag))
print "The poles should be zero."
print "The finite part of the result obained is: %.16e %s %.16e*i"%(result[0].real, '+' if result[0].imag > 0. else '-', abs(result[0].imag))
print "The finite part of the result should be  : %.16e - %.16e*i"%( 1.6242404963212e-10, 3.55266416995851e-10)
