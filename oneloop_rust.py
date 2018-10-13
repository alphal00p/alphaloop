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
from integrator import Integrator
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
        self.mu_sq = -1e8

        self.parameterization = None
        self.deformation = None
        self.internal_masses = [0., 0., 0., 0., 0., 0., 0.]

        self.integrand = integrand.Integrand(
            "box1L_direct_integration", channel=0, region=0, mu_sq=self.mu_sq)

    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        # compute the real phase for now
        return self.aggregrator.evaluate([x for x in continuous_inputs])[0]

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

        self.aggregrator = Integrator("box1L_direct_integration",
                                      do_regions=True, do_multichanneling=True, mu_sq=-1e8, e_cm_sq=sqrt_s**2, ext_py=[[x for x in e]
                                                                                                                       for e in self.external_momenta])
        # TODO: set sensible options
        integrator = Vegas3Integrator(self, verbosity=2, survey_n_points=100000, cluster=MultiCore(4),
                                      survey_n_iterations=10, refine_n_points=200000, refine_n_iterations=5)

        #integrator = pyCubaIntegrator(self, algorithm= 'Vegas', verbosity=2, target_accuracy=1.0e-3, max_eval=100000, n_start=1000, n_increase=1000, nb_core=4)

        amplitude, error = integrator.integrate()
        print("Result: %s +/- %s" % (amplitude, error))


random.seed(2)
d = Box()
# d.sample_points()
d.integrate(1000., [100., 200., 300., 400.])
