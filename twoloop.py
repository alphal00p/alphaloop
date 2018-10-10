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
from math import sqrt, pi

import numdifftools as nd
import numpy as np
import numpy.linalg as linalg
from mpmath import polylog, ln

import madgraph.integrator.vectors as vectors
import madgraph.integrator.phase_space_generators as phase_space_generators
from madgraph.integrator.vegas3_integrator import Vegas3Integrator
import madgraph.integrator.integrands as integrands
from madgraph.various.cluster import MultiCore


class DoubleBox(integrands.VirtualIntegrand):
    def __init__(self):
        n_loops = 2
        dimensions = integrands.DimensionList(sum([[
            integrands.ContinuousDimension(
                'l%d_E' % i_loop, lower_bound=0.0, upper_bound=1.0),
            integrands.ContinuousDimension(
                'l%d_x' % i_loop, lower_bound=0.0, upper_bound=1.0),
            integrands.ContinuousDimension(
                'l%d_y' % i_loop, lower_bound=0.0, upper_bound=1.0),
            integrands.ContinuousDimension(
                'l%d_z' % i_loop, lower_bound=0.0, upper_bound=1.0),
        ] for i_loop in range(1, n_loops+1)], []))

        super(DoubleBox, self).__init__(dimensions=dimensions)

        self.dimension_name_to_position = {
            d.name: i for i, d in enumerate(dimensions)}

        # set up a configuration
        self.channel = 0
        self.mu_sq = -1e-5

        self.parameterization = None
        self.deformation = None
        self.internal_masses = [0., 0., 0., 0., 0., 0., 0.]

        #self.integrand = DummyIntegrand()
        self.integrand = integrand.Integrand(
            "box2L_direct_integration", channel=self.channel, region=0, mu_sq=self.mu_sq)

    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        # compute the real phase for now
        return self.evaluate([x for x in continuous_inputs[:4]], [x for x in continuous_inputs[4:]]).real

    def deform(self, loop_momenta):
        k_mapped = vectors.LorentzVector(loop_momenta[:4])
        l_mapped = vectors.LorentzVector(loop_momenta[4:])

        # compute the external momenta for the C12(k) cycle and express the qs in terms of that.
        # This is not a shift. We also set q[0] = 0 at the 1/k^2 line
        # k, l, k - l, k - l - p2, k - l - p2 - p3, k - p2 - p3, k + p1
        # TODO: check
        C12_qs = [vectors.LorentzVector(), k_mapped - l_mapped,
                  self.external_momenta[1] + self.external_momenta[2],
                  -self.external_momenta[0]]

        self.deformation.set_qs([[x for x in e] for e in C12_qs])
        mapped_r, jac_real, jac_imag = self.deformation.deform(
            [x for x in k_mapped])
        # just take the imaginary part
        C12_k = [imag for real, imag in mapped_r]

        # now do the same for C12(l)_qs
        # TODO: check
        C12_qs = [vectors.LorentzVector(), -k_mapped + l_mapped + self.external_momenta[1] + self.external_momenta[2],
                  -k_mapped + l_mapped - self.external_momenta[0],
                  -k_mapped + l_mapped]

        self.deformation.set_qs([[x for x in e] for e in C12_qs])
        mapped_r, jac_real, jac_imag = self.deformation.deform(
            [x for x in l_mapped])
        # just take the imaginary part
        C12_l = [imag for real, imag in mapped_r]

        # cycle C23(l), going with the momentum direction
        # TODO: check
        C23_qs = [vectors.LorentzVector(), -k_mapped + self.external_momenta[1] + self.external_momenta[2],
                  -k_mapped + self.external_momenta[1], -k_mapped]

        self.deformation.set_qs([[x for x in e] for e in C23_qs])
        mapped_r, jac_real, jac_imag = self.deformation.deform(
            [x for x in l_mapped])
        C23_l = [imag for real, imag in mapped_r]

        # cycle C13(k)
        C13_qs = [vectors.LorentzVector(), l_mapped, l_mapped + self.external_momenta[1],
                  l_mapped +
                  self.external_momenta[1] + self.external_momenta[2],
                  self.external_momenta[1] + self.external_momenta[2],
                  -self.external_momenta[0]]

        self.deformation.set_qs([[x for x in e] for e in C13_qs])
        mapped_r, jac_real, jac_imag = self.deformation.deform(
            [x for x in k_mapped])
        C13_k = [imag for real, imag in mapped_r]

        k_1 = vectors.LorentzVector(C12_k) + vectors.LorentzVector(C13_k)
        k_2 = vectors.LorentzVector(C12_l) + vectors.LorentzVector(C23_l)

        lambda_overall = 1  # maximum lambda value

        # propagators with substituted momenta split in real and imag
        # k, l, k - l, k - l - p2, k - l - p2 - p3, k - p2 - p3, k + p1
        props = [(k_mapped, k_1), (l_mapped, k_2), (k_mapped - l_mapped, k_1 - k_2),
                 (k_mapped - l_mapped - self.external_momenta[1], k_1 - k_2),
                 (k_mapped - l_mapped -
                  self.external_momenta[1] - self.external_momenta[2], k_1 - k_2),
                 (k_mapped - self.external_momenta[1] -
                  self.external_momenta[2], k_1),
                 (k_mapped + self.external_momenta[0], k_1)]

        for (qt, kt), m in zip(props, self.internal_masses):
            xj = (kt.dot(qt) / kt.square())**2
            yj = (
                qt.square() - m**2) / kt.square()

            if 2 * xj < yj:
                lambda_overall = min(lambda_overall, sqrt(yj / 4.))
            elif yj < 0:
                lambda_overall = min(lambda_overall, sqrt(xj - yj/2.))
            else:
                lambda_overall = min(lambda_overall, sqrt(xj - yj/4.))

        k_1 = k_mapped + lambda_overall * k_1 * 1j
        k_2 = k_mapped + lambda_overall * k_2 * 1j

        return [x for x in k_1] + [x for x in k_2]

    def evaluate(self, k, l):
        # TODO: does the per-channel treatment still make sense?
        # should we do channels over the two-loop integral?
        self.parameterization.set_mode("log")
        k_mapped, jac_k = self.parameterization.map(k)
        l_mapped, jac_l = self.parameterization.map(l)
        numerical_jacobian_weight = 1.
        ks = k_mapped + l_mapped

        if True:
            ks = self.deform(k_mapped + l_mapped)

            # compute numerical Jacobian
            def wrapped_function(loop_momenta):
                return np.r_[self.deform(loop_momenta)]

            local_point = k_mapped + l_mapped

            # note: method=forward reduces the number of samples from 244 to 122
            numerical_jacobian = nd.Jacobian(
                wrapped_function, full_output=False, method='central')(local_point)

            # And now compute the determinant
            numerical_jacobian_weight = linalg.det(numerical_jacobian)

            # NOTE: This happens a lot!
            # if np.max(info.error_estimate) > 1.e-3:
            #    print(
            #        "Large error of %f (for which det(jac)=%s) encountered in the numerical evaluation of the Jacobian for the inputs: %s" %
            #        (np.max(info.error_estimate), numerical_jacobian_weight, str(local_point)))

        # integrate the two-loop
        out_real, out_imag = self.integrand.evaluate([ks[:4], ks[4:]])

        result = complex(out_real, out_imag) * \
            numerical_jacobian_weight * jac_k * jac_l
        return result

    def integrate(self, sqrt_s, masses):
        phase_space_generator = phase_space_generators.FlatInvertiblePhasespace(
            masses[:2], masses[2:],
            [sqrt_s/2., sqrt_s/2.],
            beam_types=(1, 1)
        )

        # Specifying None to get a random PS point
        random_PS_point, _, _, _ = phase_space_generator.get_PS_point(None)
        random_PS_point = random_PS_point.to_dict()

        for i in random_PS_point:
            if i > 2:
                continue
            random_PS_point[i] = -random_PS_point[i]

        print(random_PS_point)

        self.external_momenta = [random_PS_point[x] for x in range(1, 5)]

        print("Analytic result: %s" %
              self.doublebox_analytic(*self.external_momenta))

        # defaults to Weinzierl mapping
        # set dummy qs
        self.parameterization = deformation.Parameterization(
            e_cm_sq=sqrt_s**2, region=0, channel=self.channel, qs_py=[[v for v in q] for q in self.external_momenta])

        self.deformation = deformation.Deformation(
            e_cm_sq=sqrt_s**2, mu_sq=self.mu_sq, region=0, qs_py=[[v for v in q] for q in self.external_momenta], masses=self.internal_masses)

        self.integrand.set_externals([[x for x in e]
                                      for e in self.external_momenta])

        # TODO: set sensible options
        integrator = Vegas3Integrator(self, verbosity=1, cluster=MultiCore(
            4), survey_n_points=10000, survey_n_iterations=10, refine_n_points=20000, refine_n_iterations=5)

        amplitude, error = integrator.integrate()
        print("Result: %s +/ %s" % (amplitude, error))

    def sample_points(self):
        """Sample random points"""

        qs = [[0., 0., 0., 0.], [-495., 0., 0., 495.],
              [5., 26.97011261, 42.00346171, 992.50208264], [505., 0., 0., 495.]]

        qs_l = [vectors.LorentzVector(q) for q in qs]
        external_momenta = [qs_l[i] - qs_l[i - 1] for i in range(4)]
        masses = [0, 0, 0, 0]

        # fixed like this in C++
        e_cm_sq = (external_momenta[0] + external_momenta[1]).square()
        self.external_momenta = external_momenta

        # defaults to Weinzierl mapping
        self.parameterization = deformation.Parameterization(
            e_cm_sq=e_cm_sq, region=0, channel=self.channel, qs_py=qs)

        self.deformation = deformation.Deformation(
            e_cm_sq=e_cm_sq, mu_sq=self.mu_sq, region=0, qs_py=qs, masses=masses)

        self.integrand.set_externals([[x for x in e]
                                      for e in external_momenta])

        t = time.time()
        N = 30

        random.seed(123)

        for _ in range(N):
            # first loop momentum
            k = [random.random(), random.random(),
                 random.random(), random.random()]

            # second loop momentum
            l = [random.random(), random.random(),
                 random.random(), random.random()]

            print('Out', self.evaluate(k, l))

        print('Time', time.time() - t)

    def doublebox_analytic(self, k1, k2, k3, k4):
        # from http://www.higgs.de/~davyd/preprints/ud1.pdf
        def lmbda(x, y):
            return sqrt((1.-x-y)**2 - 4. * x * y)

        def rho(x, y):
            return 2. / (1. - x - y + lmbda(x, y))

        s = (k1 + k2).square()
        t = (k2 + k3).square()

        x = k1.square() * k3.square() / float(s * t)
        y = k2.square() * k4.square() / float(s * t)

        l = lmbda(x, y)
        r = rho(x, y)
        rx = r * x
        ry = r * y

        ln_rx = ln(rx) if x >= 0 else ln(-rx) + 1j * pi * np.sign(s*t)
        ln_ry = ln(ry) if y >= 0 else ln(-ry) + 1j * pi * np.sign(s*t)

        return - t * (pi**2*1j/(s*t))**2 / l * (6. * (polylog(4, -rx) + polylog(4, -ry)) + 3 * ln(y/x)*(polylog(3, -rx) - polylog(3, -ry)) + 0.5*ln(y/x)**2*(polylog(2, -rx)+polylog(2, -ry))
                                                + 0.25*ln_rx**2*ln_ry**2 + pi**2/2.*ln_rx*ln_ry + pi**2/12*ln(y/x)**2 + 7 * pi**4/60.)


class DummyIntegrand():
    def set_externals(self, *args):
        pass

    def evaluate(self, l_moms):
        euclidian_product = sum(
            l_moms[0][i] ** 2 for i in range(4)) + sum(l_moms[1][i] ** 2 for i in range(4))
        Mass_scale = 1.
        regulator = 50.
        d = 6
        #
        # with the above parameter, the analytical result then reads for:
        #  Mass_scale = 1.   result --> 2.6654
        #  Mass_scale = 100. result --> 2.6087*10^16
        #  mscale=1, reg=500 result --> 1.23639
        integrand = (
            1. / ((euclidian_product / (Mass_scale ** 2)) ** d + regulator))
        return (integrand, 0.)


random.seed(123)
d = DoubleBox()
# d.sample_points()
d.integrate(1000., [100., 200., 300., 400.])
