from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp
import math
import vectors

# load c++ library
Pentabox_physical = IntegralLibrary('Pentabox_physical/Pentabox_physical_pylink.so')

# choose integrator
Pentabox_physical.use_Vegas(flags=2, epsrel=1e-3, epsabs=1e-25,
nstart=100000, nincrease=100000, maxeval=500000000,
real_complex_together=False)

p1 = vectors.LorentzVector(
    [0.149500000000000E+01,    0.000000000000000E+00,    0.000000000000000E+00,    0.149165176901313E+01]
) 
p2 = vectors.LorentzVector(
    [0.150500000000000E+01,    0.000000000000000E+00,    0.000000000000000E+00,   -0.149165176901313E+01]
)
p3 = vectors.LorentzVector(
    [-0.126041949101381e+01,    -0.452362952912639e+00,    -0.101350243653045e+01,   0.516563513332600e+00]
)
p4 = vectors.LorentzVector(
    [-0.105098730574850e+01,   0.489324061520790e-01,   0.928212188578101e+00,    -0.283905035967510e+00]
)
loop_mass = 0.35
# mass: 80.419 outgoing W-
p5 = -p4-p3-p2-p1

sqrt = math.sqrt
print(sqrt(p1.square()))
print(sqrt(p2.square()))
print(sqrt(p3.square()))
print(sqrt(p4.square()))
print(sqrt(p5.square()))
print(p1+p2+p3+p4+p5)

# integrate
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = Pentabox_physical(real_parameters=[

 p1.square(),
 p2.square(),
 p3.square(),
 p4.square(),
# p5.square(),

 p1.dot(p2),
 p1.dot(p3),
 p1.dot(p4),
# p1.dot(p5),

 p2.dot(p3),
 p2.dot(p4),
# p2.dot(p5),

 p3.dot(p4),
# p3.dot(p5),
 
# p4.dot(p5),

],
complex_parameters = [complex(loop_mass**2,0.0),],
number_of_presamples=10**6,deformation_parameters_maximum=0.3
)

# convert complex numbers from c++ to sympy notation
str_integral_with_prefactor = str_integral_with_prefactor.replace(',','+I*')
str_prefactor = str_prefactor.replace(',','+I*')
str_integral_without_prefactor = str_integral_without_prefactor.replace(',','+I*')

# convert result to sympy expressions
integral_with_prefactor = sp.sympify(str_integral_with_prefactor.replace('+/-','*value+error*'))
integral_with_prefactor_err = sp.sympify(str_integral_with_prefactor.replace('+/-','*value+error*'))
prefactor = sp.sympify(str_prefactor)
integral_without_prefactor = sp.sympify(str_integral_without_prefactor.replace('+/-','*value+error*'))
integral_without_prefactor_err = sp.sympify(str_integral_without_prefactor.replace('+/-','*value+error*'))

# numerical result
conversion_factor = ((2.0*math.pi)**4)**2 / (math.pi**2)**2
print('Numerical Result')
print('Verify lack of poles:')
print('eps^-1:', integral_with_prefactor.coeff('eps',-1).coeff('value') / conversion_factor, 
      '+/- (', integral_with_prefactor_err.coeff('eps',-1).coeff('error') / conversion_factor, 
      ')')
print('Finite:')
print('eps^0:', integral_with_prefactor.coeff('eps',0).coeff('value') / conversion_factor, 
      '+/- (', integral_with_prefactor_err.coeff('eps',0).coeff('error') / conversion_factor, 
      ')')
#print('eps^1:', integral_with_prefactor.coeff('eps',1).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',1).coeff('error'), ')')
#print('eps^2:', integral_with_prefactor.coeff('eps',2).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',2).coeff('error'), ')')
#print('eps^3:', integral_with_prefactor.coeff('eps',3).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',3).coeff('error'), ')')
#print('eps^4:', integral_with_prefactor.coeff('eps',4).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',4).coeff('error'), ')')
