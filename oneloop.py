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
import logging
import math

import numdifftools as nd
import numpy as np
import numpy.linalg as linalg
from mpmath import polylog, ln

import madgraph.integrator.vectors as vectors
import madgraph.integrator.phase_space_generators as phase_space_generators
from madgraph.integrator.vegas3_integrator import Vegas3Integrator
from madgraph.integrator.pyCubaIntegrator import pyCubaIntegrator
import madgraph.integrator.integrands as integrands
from madgraph.various.cluster import MultiCore

logger = logging.getLogger('doublebox')


class Box(integrands.VirtualIntegrand):
    def __init__(self):
        n_loops = 1
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

        super(Box, self).__init__(dimensions=dimensions)

        self.dimension_name_to_position = {
            d.name: i for i, d in enumerate(dimensions)}

        # set up a configuration
        self.channel = 0
        self.mu_sq = -1e6

        self.parameterization = None
        self.deformation = None
        self.internal_masses = [0., 0., 0., 0., 0., 0., 0.]

        self.integrand = integrand.Integrand(
            "box1L_direct_integration", channel=self.channel, region=0, mu_sq=self.mu_sq)

    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        # compute the real phase for now
        return self.evaluate([x for x in continuous_inputs]).real

    def evaluate(self, k):
        if True:
            # TODO: also sample -k_mapped in external region
            r1 = self.evaluate_region(k, 1) # external
            r2 = self.evaluate_region(k, 2) # internal
            return r1 + r2
        else:
            return self.evaluate_region(k, 0)

    def evaluate_region(self, k, region=0):
        # TODO: does the per-channel treatment still make sense?
        # should we do channels over the two-loop integral?
        if region == 1:
            # external region
            self.parameterization.set_mode("weinzierl")
        else:
            # for now we don't do channels
            self.parameterization.set_mode("log")
        self.parameterization.set_region(region)
        k_mapped, jac_k = self.parameterization.map(k)

        self.deformation.set_region(region)
        (ksv, jac_real, jac_imag) = self.deformation.deform(
            [x for x in k_mapped])

        ks = [complex(*x) for x in ksv]
        jac = complex(jac_real, jac_imag)

        # integrate the two-loop
        self.integrand.set_region(region)
        out_real, out_imag = self.integrand.evaluate([ks])

        result = complex(out_real, out_imag) * \
            jac * jac_k

        # for Cuba, since it evaluates the boundaries
        if math.isnan(result.real) or math.isnan(result.imag):
            return complex(0, 0)

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

        # compute qs
        qs = [vectors.LorentzVector()]
        for i, x in enumerate(self.external_momenta[:-1]):
            qs.append(qs[i] + x)

        qs = [[x for x in y] for y in qs]

        # defaults to Weinzierl mapping
        self.parameterization = deformation.Parameterization(
            e_cm_sq=sqrt_s**2, region=0, channel=self.channel, qs_py=qs)

        self.deformation = deformation.Deformation(
            e_cm_sq=sqrt_s**2, mu_sq=self.mu_sq, region=0, qs_py=qs, masses=self.internal_masses)

        self.deformation.set_external_momenta([[x for x in e]
                                               for e in self.external_momenta])

        self.integrand.set_externals([[x for x in e]
                                      for e in self.external_momenta])

        # TODO: set sensible options
        integrator = Vegas3Integrator(self, verbosity=2, cluster=MultiCore(
            4), survey_n_points=1000000, survey_n_iterations=10, refine_n_points=200000, refine_n_iterations=5)

        #integrator = pyCubaIntegrator(self, algorithm= 'Vegas', verbosity=2, target_accuracy=1.0e-3, max_eval=100000, n_start=1000, n_increase=1000, nb_core=4)

        amplitude, error = integrator.integrate()
        print("Result: %s +/ %s" % (amplitude, error))


random.seed(2)
d = Box()
# d.sample_points()
d.integrate(1000., [100., 200., 300., 400.])
