import math
import mpmath
import numpy
from scipy.special import zeta

#############################################################################################################
# List analytical expressions
#############################################################################################################

def ladder_phi(x,y,l):
    ladder_lambda = lambda x,y: numpy.sqrt((1.-x-y)**2-4*x*y+0j)
    ladder_rho = lambda x,y: 2./(1.-x-y+ladder_lambda(x,y))

    bb = -1./(math.factorial(l)*ladder_lambda(x,y))
    summand = 0.
    for j in xrange(l,2*l+1):
        cc = (-1.)**j*math.factorial(j)*numpy.log(y/x+0j)**(2*l-j)
        dd = math.factorial(j-l)*math.factorial(2*l-j)
        ee = mpmath.polylog(j,-1./(x*ladder_rho(x,y)))-mpmath.polylog(j,-y*ladder_rho(x,y))
        summand += cc*ee/dd
    return bb*summand

def analytic_four_point_ladder(s1,s2,s3,s4,s,t,l):
    # leg labelling according to Weinzierl https://arxiv.org/abs/1211.0509v3
    # p2 --> |------| <-- p3
    #        |-Diag-|
    # p1 --> |------| <-- p4
    # s1 = p1^2, s2 = p2^2, s3 = p3^2, s4 = p4^2, s = (p1+p2)^2, t = (p2+p3)^2
    # l = nr of loops
    factor = 1./t*(1j/(16.*math.pi**2*s))**l
    X = s1*s3/(s*t)
    Y = s2*s4/(s*t)
    return factor*ladder_phi(X,Y,l)

def analytic_three_point_ladder(s1,s2,s3,l):
    # leg labelling according to Weinzierl https://arxiv.org/abs/1211.0509v3
    #        |------| <-- p1
    # p3 --> |-Diag-|
    #        |------| <-- p2
    # s1 = p1^2, s2 = p2^2, s3 = p3^2
    # l = nr of loops
    factor = (1j/(16*math.pi**2*s3))**l
    x = s1/s3
    y = s2/s3
    return factor*ladder_phi(x,y,l)

def analytic_two_point_ladder(s,l):
    # leg labelling according to Weinzierl https://arxiv.org/abs/1211.0509v3
    # p --> |-Diag-| --> p
    # s = p^2
    # l = nr of loops
    result = (1j/(16.*math.pi**2*s))**l
    result *= s
    result *= math.factorial(2*l)/math.factorial(l)**2
    result *= zeta(2*l-1)
    return result

