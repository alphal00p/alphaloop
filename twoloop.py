import os
import sys
import inspect
currentdir = os.path.dirname(os.path.abspath(
    inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, currentdir + "/rust_backend")
sys.path.insert(0, currentdir + "/../..")

import integrand
import deformation
import random

import loop_momenta_generator
import nloop_integrands

import time
from math import sqrt

import numdifftools as nd
import numpy as np
import numpy.linalg as linalg

import madgraph.integrator.vectors as vectors

channel = 0
mu_sq = -1e-5
qs = [[0., 0., 0., 0.], [-495., 0., 0., 495.],
      [5., 26.97011261, 42.00346171, 992.50208264], [505., 0., 0., 495.]]

qs_l = [vectors.LorentzVector(q) for q in qs]
external_momenta = [qs_l[i] - qs_l[i - 1] for i in range(4)]
masses = [0, 0, 0, 0]

# fixed like this in C++
e_cm_sq = (external_momenta[0] + external_momenta[1]).square()

# defaults to Weinzierl mapping
parameterization = deformation.Parameterization(
    e_cm_sq=e_cm_sq, region=0, channel=channel, qs_py=qs)

deformation = deformation.Deformation(
    e_cm_sq=e_cm_sq, mu_sq=mu_sq, region=0, qs_py=qs, masses=masses)

doublebox = integrand.Integrand(
    "box2L_direct_integration", channel=channel, region=0, mu_sq=mu_sq)

doublebox.set_externals([[x for x in e] for e in external_momenta])

t = time.time()
N = 10

random.seed(123)


def deform(loop_momenta):
    k_mapped = vectors.LorentzVector(loop_momenta[:4])
    l_mapped = vectors.LorentzVector(loop_momenta[4:])

    # for the shift of 1/l^2
    a = l_mapped - k_mapped

    # cycle 12 with l -> k - l
    C12_qs = [vectors.LorentzVector(), a, a + external_momenta[1] + external_momenta[2],
              a + external_momenta[1] + external_momenta[2] - external_momenta[0]]

    deformation.set_qs([[x for x in e] for e in C12_qs])
    mapped_r, jac_real, jac_imag = deformation.deform([x for x in k_mapped])
    C12_k = [imag for real, imag in mapped_r]  # just take the imaginary part

    # cycle 23 with k = k + 2l
    C23_qs = [vectors.LorentzVector(), -k_mapped, -k_mapped + external_momenta[1],
              -k_mapped + external_momenta[1] + external_momenta[2]]

    deformation.set_qs([[x for x in e] for e in C23_qs])
    mapped_r, jac_real, jac_imag = deformation.deform([x for x in l_mapped])
    C23_k = [imag for real, imag in mapped_r]

    # cycle 13 without shift
    C13_qs = [vectors.LorentzVector(), l_mapped, l_mapped + external_momenta[1],
              l_mapped + external_momenta[1] + external_momenta[2],
              external_momenta[1] + external_momenta[2],
              -external_momenta[0]]

    deformation.set_qs([[x for x in e] for e in C13_qs])
    mapped_r, jac_real, jac_imag = deformation.deform([x for x in k_mapped])
    C13_k = [imag for real, imag in mapped_r]

    k_1 = vectors.LorentzVector(C12_k) + vectors.LorentzVector(C13_k)

    # TODO: what about all the shifts?
    # should we use C12_l instead?
    k_2 = vectors.LorentzVector(C12_k) + vectors.LorentzVector(C23_k)

    lambda_overall = 1  # maximum lambda value
    ks = [k_1, k_2]
    lmom = [k_mapped, l_mapped]

    # TODO: the paper suggests the loop is over the propagators
    # but we only have two deformed ks instead of 7 and their propagator is simply 1/k^2
    for j in range(2):
        xj = (ks[j].dot(lmom[j]) / ks[j].square())**2
        yj = (
            (lmom[j]).square() - masses[j]**2) / ks[j].square()

        if 2 * xj < yj:
            lambda_overall = min(lambda_overall, sqrt(yj / 4.))
        elif yj < 0:
            lambda_overall = min(lambda_overall, sqrt(xj - yj/2.))
        else:
            lambda_overall = min(lambda_overall, sqrt(xj - yj/4.))

    k_1 = k_mapped + lambda_overall * k_1 * 1j
    k_2 = k_mapped + lambda_overall * k_2 * 1j

    return [x for x in k_1] + [x for x in k_2]


for _ in range(N):
    # first loop momentum
    k = [random.random(), random.random(),
         random.random(), random.random()]

    # second loop momentum
    l = [random.random(), random.random(),
         random.random(), random.random()]

    # TODO: does the per-channel treatment still make sense?
    # should we do channels over the two-loop integral?
    parameterization.set_mode("log")
    k_mapped, jac_k = parameterization.map(k)
    l_mapped, jac_l = parameterization.map(l)

    ks = deform(k_mapped + l_mapped)

    # compute numerical Jacobian
    def wrapped_function(loop_momenta):
        return np.r_[deform(loop_momenta)]

    local_point = k_mapped + l_mapped
    numerical_jacobian, info = nd.Jacobian(
        wrapped_function, full_output=True)(local_point)

    # And now compute the determinant
    numerical_jacobian_weight = linalg.det(numerical_jacobian)

    if np.max(info.error_estimate) > 1.e-3:
        print(
            "Large error of %f (for which det(jac)=%f) encountered in the numerical evaluation of the Jacobian for the inputs: %s" %
            (np.max(info.error_estimate), numerical_jacobian_weight, str(local_point)))

    # integrate the two-loop
    out_real, out_imag = doublebox.evaluate([ks[:4], ks[4:]])

    result = complex(out_real) * numerical_jacobian_weight * jac_k * jac_l

    print('Out', result)

print('Time', time.time() - t)
