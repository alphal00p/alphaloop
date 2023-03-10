#####################################################
#                                                   #
#  Source file of the pyNLoop GGVV MG5aMC plugin.   #
#  Use only with consent of its authors.            #
#                                                   #
#  author: Valentin Hirschi, Ben Ruij               #
#                                                   #
#####################################################

import math
import cmath
import shutil
import os
import logging

import madgraph.integrator.vectors as vectors
import madgraph.integrator.integrands as integrands
import madgraph.various.misc as misc
import loop_momenta_generator
from madgraph import InvalidCmd, MadGraph5Error
import itertools
import numpy as np
from pyNLoop import plugin_path
import ctypes

import utils

logger = logging.getLogger('pyNLoop.Integrand')

pjoin = os.path.join

class IntegrandError(MadGraph5Error):
    """ Error for the LoopMomentaGenerator class suite."""
    pass

try:
    import pySecDec
    from pySecDec.loop_integral import loop_package as pySecDec_loop_package
except ImportError:
    logger.warning("Import of pySecDec fails; you will not be able integrate loops depending on this integrator.")

class NLoopIntegrand(integrands.VirtualIntegrand):
    """ Class for implementing an N-loop integrand."""

    # Specify here which integrator this integrand supports
    _supported_integrators = ['Vegas3', 'pySecDec']
    
    # Of course only relevant for integrands suited for the pySecDec integrator
    _pySecDecOutputName    = None

    def __init__(self,
                 dimensions=None,
                 n_loops=1,
                 external_momenta=vectors.LorentzVectorDict(),
                 phase_computed='Real',
                 # Data-structure for specifying a topology to be determined
                 topology=None,
                 **opts):

        # If dimensions are not specified then generate a naive direct space in an hypercube
        if dimensions is None:
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

        super(NLoopIntegrand, self).__init__(dimensions=dimensions)

        self.dimension_name_to_position = {d.name: i for i, d in enumerate(dimensions)}

        self.n_loops = 1
        self.external_momenta = external_momenta
        self.phase_computed = phase_computed
        self.topology = topology
        self.integral_scaling = 1.0 # scaling of the integral

    def assign_kinematic_configuration(self, PS_point):
        """ Function allowing to redefine a new kinematic configuration, re-instantiating what must be re-instantiated."""
        raise NotImplementedError("Daughter integrand classes must define the function 'assign_kinematic_configuration' for"+
                                  " it to be available, as the actions that must be taken here are particular to each.")

    def get_integrated_counterterms(self):
        """ In case local numerical subtraction is employed, then the integrand of the
        correpsonding integrated counterterms must be returned as well."""
        
        # By default, an empty list is returned.
        return []

    def nice_string(self):
        """ Short string representation of this integrand. We could eventually add more
        attributes to the string returned. But for now simply return the class name."""
        return self.__class__.__name__

    def output(self, output_folder, force=False, **opts):
        """ Possibly output some low-level code representation of the integrand to make
        its evaluation faster."""

        self.output_folder = output_folder
        if os.path.exists(output_folder):
            if force:
                return False
            if raw_input("Do you want to recycle existing output '%s'?\n[Y/N] > " %
                                  output_folder).upper() not in ['N', 'NO']:
                logger.info("Recycling existing output '%s'." % output_folder)
                return False
            else:
                logger.warning(
                    "Folder '%s' already exists, it will be overwritten." % output_folder)
                shutil.rmtree(output_folder)
        os.mkdir(output_folder)
        return True

    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        """ Actual evaluation of the call."""

        raise NotImplementedError

    # Handling of the analytic result. By default, set it non-available
    def setup_analytic_computation(self, pyNLoop_command, *args, **opts):
        """ Performs tasks possibly necessary (like loading a library) for computing the analytic result"""
        pass
    def is_analytic_result_available(self, *args, **opts):
        """ Checks if an analytical result exists for this integrand."""
        return False
    def get_analytic_result(self, PS_point, *args, **opts):
        """ Return the real and imaginary part of the analytical result."""
        raise NotImplementedError('Analytic result apparently not available for integrand %s.'%self.nice_string())

class IntegratedCounterterm(object):
    """ Container class from which integrated counterterms inherit."""
    
    def __init__(self, *args, **opts):
        if 'integrated_CT_name' in opts:
            self.integrated_CT_name = opts.pop('integrated_CT_name')
        else:
            self.integrated_CT_name = 'generic_integrated_CT'
        return opts
    
    def nice_string(self):
        """ Short string representation of this integrand. Add as a suffix the name of 
        this integrated counterterm to the name as there can be several integrated CT per loop."""
        return '%s@%s'%(self.__class__.__name__, self.integrated_CT_name)

class box1L(NLoopIntegrand):

    # We plan on being able to use pySecDec integrator only for this case
    _supported_integrators = ['pySecDec']
    
    _pySecDecOutputName    = 'pySecDecLoopPackage'

    def __init__(self,
                 n_loops=1,
                 external_momenta=vectors.LorentzVectorDict(),
                 phase_computed='Real',
                 # Data-structure for specifying a topology to be determined
                 topology=None,
                 **opts):

        # Create the Feynman parameters dimensions. As we intend to use the pySecDec integrator
        # this has little relevance anyway.
        dimensions = integrands.DimensionList([
            integrands.ContinuousDimension(
                'x%d' % (i_prop+1), lower_bound=0.0, upper_bound=1.0)
            for i_prop in range(len(external_momenta))])

        # I am copying here all options just to make it explicit that could be used here as well.
        super(box1L, self).__init__(
            dimensions=dimensions,
            n_loops=n_loops,
            external_momenta=external_momenta,
            phase_computed=phase_computed,
            # Data-structure for specifying a topology to be determined
            topology=topology,
            **opts
        )

        # Now generate the pySecDec loop integrand object
        self.pySecDec_loop_integrand = pySecDec.loop_integral.LoopIntegralFromGraph(
            internal_lines=[['m', [1, 2]], [
                0, [2, 3]], [0, [3, 4]], [0, [4, 1]]],
            external_lines=[['p1', 1], ['p2', 2], ['p3', 3], ['p4', 4]],
            replacement_rules=[
                ('p1*p1', 's1'),
                ('p2*p2', 0),
                ('p3*p3', 0),
                ('p4*p4', 0),
                ('p3*p2', 't/2'),
                ('p1*p2', 's/2-s1/2'),
                ('p1*p4', 't/2-s1/2'),
                ('p2*p4', 's1/2-t/2-s/2'),
                ('p3*p4', 's/2'),
                ('m**2', 'msq')
            ]
        )
        self.integrand_parameters = ['s', 't', 's1', 'msq']

        self.integrand_parameters_values = [
            2.*self.external_momenta[3].dot(self.external_momenta[4]),
            2.*self.external_momenta[3].dot(self.external_momenta[2]),
            self.external_momenta[1].dot(self.external_momenta[1]),
            # Let's hardcode a mass for now until we figure out how we want to specify the
            # topology
            1.0
        ]

        # Normally we should use the above, but as a quick hack to exactly match the
        # result from the canonical pySecDec example, we will use the inputs below
        self.integrand_parameters_values = [4.0, -0.75, 1.25, 1.0]

    def output(self, output_folder, verbosity=0, **opts):
        """ Possibly output some low-level code representation of the integrand to make
        its evaluation faster."""

        if not super(box1L, self).output(output_folder, **opts):
            return False

        logger.info("Generating pySecDec output for integrand '%s' ..." %
                    self.nice_string())
        with misc.chdir(self.output_folder) as cwd:
            with misc.Silence(active=(verbosity == 0)):
                pySecDec_loop_package(
                    name='pySecDecLoopPackage',
                    loop_integral=self.pySecDec_loop_integrand,
                    real_parameters=self.integrand_parameters,
                    # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
                    requested_order=0,
                    # the optimization level to use in FORM (can be 0, 1, 2, 3)
                    form_optimization_level=2,
                    # the WorkSpace parameter for FORM
                    form_work_space='100M',
                    # the method to be used for the sector decomposition
                    # valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
                    decomposition_method='iterative',
                    # if you choose ``geometric[_ku]`` and 'normaliz' is not in your
                    # $PATH, you can set the path to the 'normaliz' command-line
                    # executable here
                    # normaliz_executable='/path/to/normaliz',
                )

        return True

    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        """ Actual evaluation of the call."""

        raise NotImplementedError


class box1L_offshell_massless(NLoopIntegrand):

    # We plan on being able to use pySecDec integrator only for this case
    _supported_integrators = ['pySecDec']

    _pySecDecOutputName    = 'pySecDecLoopPackage'

    def __init__(self,
                 n_loops=1,
                 external_momenta=vectors.LorentzVectorDict(),
                 phase_computed='All',
                 # Data-structure for specifying a topology to be determined
                 topology=None,
                 **opts):

        # Create the Feynman parameters dimensions. As we intend to use the pySecDec integrator
        # this has little relevance anyway.
        dimensions = integrands.DimensionList([
            integrands.ContinuousDimension(
                'x%d' % (i_prop+1), lower_bound=0.0, upper_bound=1.0)
            for i_prop in range(len(external_momenta))])

        # I am copying here all options just to make it explicit that could be used here as well.
        super(box1L_offshell_massless, self).__init__(
            dimensions=dimensions,
            n_loops=n_loops,
            external_momenta=external_momenta,
            phase_computed=phase_computed,
            # Data-structure for specifying a topology to be determined
            topology=topology,
            **opts
        )

        # Now generate the pySecDec loop integrand object
        self.pySecDec_loop_integrand = pySecDec.loop_integral.LoopIntegralFromGraph(
            internal_lines=[[0, [1, 2]], [0, [2, 3]],
                            [0, [3, 4]], [0, [4, 1]]],
            external_lines=[['p1', 1], ['p2', 2], ['p3', 3], ['p4', 4]],
            replacement_rules=[
                        ('p1*p1', 's1'),
                        ('p2*p2', 's2'),
                        ('p3*p3', 's3'),
                        ('p4*p4', 's4'),
                        ('p1*p2', 's/2-s1/2-s2/2'), # s direct
                        ('p1*p3', 's2/2+s4/2-t/2-s/2'), # u -> sum masses - s - t
                        ('p1*p4', 't/2-s1/2-s4/2'), # t direct
                        ('p2*p3', 't/2-s2/2-s3/2'), # t direct
                        ('p2*p4', 's1/2+s3/2-t/2-s/2'), # u -> sum masses - s - t
                        ('p3*p4', 's/2-s3/2-s4/2'), # s direct
            ],
#            regulator='eps',
#            regulator_power=0,
#            dimensionality='4-2*eps',
        )
        self.integrand_parameters =  ['s','t','s1','s2','s3','s4']

        self.integrand_parameters_values = [
            (self.external_momenta[1] + self.external_momenta[2]).square(), # s
            (self.external_momenta[2] + self.external_momenta[3]).square(), # t
            self.external_momenta[1].dot(self.external_momenta[1]), # s1
            self.external_momenta[2].dot(self.external_momenta[2]),
            self.external_momenta[3].dot(self.external_momenta[3]),
            self.external_momenta[4].dot(self.external_momenta[4]),
        ]

    def output(self, output_folder, verbosity=0, **opts):
        """ Possibly output some low-level code representation of the integrand to make
        its evaluation faster."""

        if not super(box1L_offshell_massless, self).output(output_folder, **opts):
            return False

        logger.info("Generating pySecDec output for integrand '%s' ..." %
                    self.nice_string())
        with misc.chdir(self.output_folder) as cwd:
            with misc.Silence(active=(verbosity == 0)):
                pySecDec_loop_package(
                    name='pySecDecLoopPackage',
                    loop_integral=self.pySecDec_loop_integrand,
                    real_parameters=self.integrand_parameters,
                    # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
                    requested_order=0,
                    # the optimization level to use in FORM (can be 0, 1, 2, 3)
                    form_optimization_level=2,
                    # the WorkSpace parameter for FORM
                    form_work_space='100M',
                    # the method to be used for the sector decomposition
                    # valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
                    decomposition_method='iterative',
                    # if you choose ``geometric[_ku]`` and 'normaliz' is not in your
                    # $PATH, you can set the path to the 'normaliz' command-line
                    # executable here
                    # normaliz_executable='/path/to/normaliz',
                )

        return True

    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        """ Actual evaluation of the call."""

        raise NotImplementedError

class box1L_onshell_massless(NLoopIntegrand):

    # We plan on being able to use pySecDec integrator only for this case
    _supported_integrators = ['pySecDec']

    _pySecDecOutputName    = 'pySecDecLoopPackage'

    def __init__(self,
                 n_loops=1,
                 external_momenta=vectors.LorentzVectorDict(),
                 phase_computed='All',
                 # Data-structure for specifying a topology to be determined
                 topology=None,
                 **opts):

        # Create the Feynman parameters dimensions. As we intend to use the pySecDec integrator
        # this has little relevance anyway.
        dimensions = integrands.DimensionList([
            integrands.ContinuousDimension(
                'x%d' % (i_prop+1), lower_bound=0.0, upper_bound=1.0)
            for i_prop in range(len(external_momenta))])

        # I am copying here all options just to make it explicit that could be used here as well.
        super(box1L_onshell_massless, self).__init__(
            dimensions=dimensions,
            n_loops=n_loops,
            external_momenta=external_momenta,
            phase_computed=phase_computed,
            # Data-structure for specifying a topology to be determined
            topology=topology,
            **opts
        )

        # Now generate the pySecDec loop integrand object
        self.pySecDec_loop_integrand = pySecDec.loop_integral.LoopIntegralFromGraph(
            internal_lines=[[0, [1, 2]], [0, [2, 3]],
                            [0, [3, 4]], [0, [4, 1]]],
            external_lines=[['p1', 1], ['p2', 2], ['p3', 3], ['p4', 4]],
            replacement_rules=[
                ('p1*p1', 0),
                ('p2*p2', 0),
                ('p3*p3', 0),
                ('p4*p4', 0),
                ('p3*p2', 't/2'),
                ('p1*p2', 's/2'),
                ('p1*p4', 't/2'),
                ('p2*p4', '-t/2-s/2'),
                ('p3*p4', 's/2')
            ],
#            regulator='eps',
#            regulator_power=0,
#            dimensionality='4-2*eps',
        )
        self.integrand_parameters = ['s', 't']

        self.integrand_parameters_values = [
            2.*self.external_momenta[3].dot(self.external_momenta[4]),
            2.*self.external_momenta[3].dot(self.external_momenta[2]),
            self.external_momenta[1].dot(self.external_momenta[1])
        ]

        # Normally we should use the above, but as a quick hack to exactly match the
        # result from the canonical pySecDec example, we will use the inputs below
        self.integrand_parameters_values = [4.0, -0.75]

    def output(self, output_folder, verbosity=0, **opts):
        """ Possibly output some low-level code representation of the integrand to make
        its evaluation faster."""

        if not super(box1L_onshell_massless, self).output(output_folder, **opts):
            return False

        logger.info("Generating pySecDec output for integrand '%s' ..." %
                    self.nice_string())
        with misc.chdir(self.output_folder) as cwd:
            with misc.Silence(active=(verbosity == 0)):
                pySecDec_loop_package(
                    name='pySecDecLoopPackage',
                    loop_integral=self.pySecDec_loop_integrand,
                    real_parameters=self.integrand_parameters,
                    # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
                    requested_order=0,
                    # the optimization level to use in FORM (can be 0, 1, 2, 3)
                    form_optimization_level=2,
                    # the WorkSpace parameter for FORM
                    form_work_space='100M',
                    # the method to be used for the sector decomposition
                    # valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
                    decomposition_method='iterative',
                    # if you choose ``geometric[_ku]`` and 'normaliz' is not in your
                    # $PATH, you can set the path to the 'normaliz' command-line
                    # executable here
                    # normaliz_executable='/path/to/normaliz',
                )

        return True

    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        """ Actual evaluation of the call."""

        raise NotImplementedError


class DummyNLoopIntegrand(NLoopIntegrand):

    # We plan on being able to use pySecDec integrator only for this case
    _supported_integrators = ['Vegas3']

    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        """ Actual evaluation of the call, with a dummy function"""

        # Let's put a dummy example for now
        l1_E = continuous_inputs[self.dimension_name_to_position['l1_E']]
        l1_x = continuous_inputs[self.dimension_name_to_position['l1_x']]
        l1_y = continuous_inputs[self.dimension_name_to_position['l1_y']]
        l1_z = continuous_inputs[self.dimension_name_to_position['l1_z']]

        res = math.sin(l1_E/l1_x) + math.cos(l1_y*l1_z) + math.tan(l1_x*l1_z)

        # misc.sprint(res)
        return res

###########################################################################################
#
# First example of a local + integrated CT ala Babis
#
###########################################################################################

class box1L_subtracted(NLoopIntegrand):
    # We plan on being able to use pySecDec integrator only for this case
    _supported_integrators = ['pySecDec']

    _pySecDecOutputName    = 'pySecDecLoopPackage_subtracted'

    def __init__(self,
                 n_loops=1,
                 external_momenta=vectors.LorentzVectorDict(),
                 phase_computed='Real',
                 # Data-structure for specifying a topology to be determined
                 topology=None,
                 **opts):

        # Create the Feynman parameters dimensions. As we intend to use the pySecDec integrator
        # this has little relevance anyway.
        dimensions = integrands.DimensionList([
            integrands.ContinuousDimension(
                'x%d' % (i_prop+1), lower_bound=0.0, upper_bound=1.0)
            for i_prop in range(len(external_momenta))])

        # I am copying here all options just to make it explicit that could be used here as well.
        super(box1L_subtracted, self).__init__(
            dimensions=dimensions,
            n_loops=n_loops,
            external_momenta=external_momenta,
            phase_computed=phase_computed,
            # Data-structure for specifying a topology to be determined
            topology=topology,
            **opts
        )

        self.define_loop_integrand()

    def define_loop_integrand(self):
        """ Define the pySecDecObjects identifying the loop integrand."""

        
        self.pySecDec_loop_integrand = pySecDec.loop_integral.LoopIntegralFromPropagators(

            loop_momenta        =   ['k1'],
            external_momenta    =   ['p1', 'p2', 'p3', 'p4'],
            Lorentz_indices     =   ['mu', 'mu1', 'mu2', 'mu3',
                                     'mu4', 'mu5', 'mu6', 'mu7', 'mu8', 'mu9'],
            propagators         =   ['k1**2', '(k1+p2)**2', '(k1+p2+p3)**2',
                                     '(k1-p1)**2', '-(k1(mu7)*p2(mu7))/s', '1-(k1(mu8)*p2(mu8))/s', 's/2-k1**2'],
            powerlist           =   [1, 1, 1, 1, 1, 1, 1],
            numerator           =   '1' + '-(%s)'%self.get_local_counterterms(),
            replacement_rules   =   [  ('p1*p1', 0),
                                       ('p2*p2', 0),
                                       ('p3*p3', 0),
                                       ('p4*p4', 0),
                                       ('p1*p2', 's/2'),
                                       ('p2*p3', 't/2'),
                                       ('p1*p3', '-t/2-s/2')
                                    ],
            regulator           =   'eps',
            regulator_power     =   0,
            dimensionality      =   '4-2*eps',
        )

        """
        self.pySecDec_loop_integrand = pySecDec.loop_integral.LoopIntegralFromGraph(
            internal_lines=[[0, [1, 2]], [0, [2, 3]],
                            [0, [3, 4]], [0, [4, 1]]],
            external_lines=[['p1', 1], ['p2', 2], ['p3', 3], ['p4', 4]],
            replacement_rules=[
                ('p1*p1', 0),
                ('p2*p2', 0),
                ('p3*p3', 0),
                ('p4*p4', 0),
                ('p3*p2', 't/2'),
                ('p1*p2', 's/2'),
                ('p1*p4', 't/2'),
                ('p2*p4', '-t/2-s/2'),
                ('p3*p4', 's/2')
            ],
            regulator='eps',
            regulator_power=0,
            dimensionality='4-2*eps',
        )
"""
        self.integrand_parameters = ['s', 't']

        self.integrand_parameters_values = [
            2.*self.external_momenta[3].dot(self.external_momenta[4]),
            2.*self.external_momenta[3].dot(self.external_momenta[2]),
            self.external_momenta[1].dot(self.external_momenta[1]),
            # Let's hardcode a mass for now until we figure out how we want to specify the
            # topology
            0.0
        ]

        # Normally we should use the above, but as a quick hack to exactly match the
        # result from the canonical pySecDec example, we will use the inputs below
        self.integrand_parameters_values = [4.0, -0.75]

    def get_local_counterterms(self):
        """ Get counterterms to the current loop integrand."""

        As = ['k1(mu)*k1(mu)',
              '(k1(mu)+p2(mu))*(k1(mu)+p2(mu))',
              '(k1(mu)+p2(mu)+p3(mu))*(k1(mu)+p2(mu)+p3(mu))',
              '(k1(mu)-p1(mu))*(k1(mu)-p1(mu))'
              ]

        softs = '(%s+%s)/t + (%s+%s)/s' % (As[0], As[2], As[1], As[3])

        # we choose p2 as the reference vector
        # p1.p2 = s/2
        def x1f(x1): return '-2*k1({0})*p2({0})/s'.format(x1)

        AsC = [x % (x1f('mu1'), x1f('mu2')) for x in
               ('%s*%s*p1(mu)*p1(mu)',
                '(%s*p1(mu)+p2(mu))*(%s*p1(mu)+p2(mu))',
                '(%s*p1(mu)+p2(mu)+p3(mu))*(%s*p1(mu)+p2(mu)+p3(mu))',
                '(%s*p1(mu)-p1(mu))*(%s*p1(mu)-p1(mu))')
               ]

        softcol = '(%s+%s)/t + (%s+%s)/s' % (AsC[0], AsC[2], AsC[1], AsC[3])

        # we set mu^2=s/2 and make sure we have a common denominator
        cols = '-(s/2*(%s)*%s*%s)/s/t' % (
            softcol, '(k1(mu3)+p2(mu3)+p3(mu3))*(k1(mu3)+p2(mu3)+p3(mu3))', '(k1(mu4)-p1(mu4))*(k1(mu4)-p1(mu4))')

        return '((%s)*(%s)*(%s)*(%s)-(%s))' % (softs, '-(k1(mu5)*p2(mu5))/s', '1+(k1(mu5)*p2(mu5))/s', 's/2-k1(mu9)*k1(mu9)', cols)

    def get_integrated_counterterms(self):
        """ Return instances of n_loop integrands corresponding to the integrated counterpart
        of the local counterterms used here."""
        return []
        # For now simply pass along the same information as this loop integrand, but eventually
        # this will be refined, especially the topology.
        return [
            box1L_integrated_CT(
                n_loops=self.n_loops,
                external_momenta=self.external_momenta,
                phase_computed=self.phase_computed,
                # Data-structure for specifying a topology to be determined
                topology=self.topology,
                # Additional name for identifying this particular integrated CT
                integrated_CT_name = 'allIntegratedCT'
            ),
        ]

    def output(self, output_folder, verbosity=0, **opts):
        """ Possibly output some low-level code representation of the integrand to make
        its evaluation faster."""

        super(box1L_subtracted, self).output(output_folder, **opts)
        
        if os.path.exists(pjoin(output_folder, self._pySecDecOutputName)):
            return False

        logger.info("Generating pySecDec output for integrand '%s' ..." %self.nice_string())
        logger.info("in folder '%s' ..." %pjoin(output_folder, self._pySecDecOutputName))
        
        with misc.chdir(self.output_folder) as cwd:
            with misc.Silence(active=(verbosity == 0)):
                pySecDec_loop_package(
                    name=self._pySecDecOutputName,
                    loop_integral=self.pySecDec_loop_integrand,
                    real_parameters=self.integrand_parameters,
                    # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
                    requested_order=0,
                    # the optimization level to use in FORM (can be 0, 1, 2, 3)
                    form_optimization_level=2,
                    # the WorkSpace parameter for FORM
                    form_work_space='100M',
                    # the method to be used for the sector decomposition
                    # valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
                    decomposition_method='iterative',
                    # if you choose ``geometric[_ku]`` and 'normaliz' is not in your
                    # $PATH, you can set the path to the 'normaliz' command-line
                    # executable here
                    # normaliz_executable='/path/to/normaliz',
                )

        return True

    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        """ Actual evaluation of the call."""

        raise NotImplementedError
    
class box1L_integrated_CT(IntegratedCounterterm, box1L_subtracted):

    _pySecDecOutputName = 'pySecDecLoopPackage_integratedCT'
    
    def __init__(self, *args, **opts):
        opts = IntegratedCounterterm.__init__(self, *args, **opts)
        box1L_subtracted.__init__(self, *args, **opts)

        if 'integral_scaling' in opts:
            self.integral_scaling = opts.pop('integral_scaling')

        # Rename the pySecDec output folder so as to avoid conflict when outputting
        # multiple integrated CT.
        self._pySecDecOutputName += '__%s'%self.integrated_CT_name
            
    def define_loop_integrand(self):

        # The definition below for now simply corresponds to the original box1L
        # It is just a token for performing structural tests for now.

        self.pySecDec_loop_integrand = pySecDec.loop_integral.LoopIntegralFromPropagators(

            loop_momenta        =   ['k1'],
            external_momenta    =   ['p1', 'p2', 'p3', 'p4'],
            Lorentz_indices     =   ['mu', 'mu1', 'mu2', 'mu3',
                                     'mu4', 'mu5', 'mu6', 'mu7', 'mu8', 'mu9'],
            propagators         =   ['TODO'],
            powerlist           =   ['TODO'],
            numerator           =   'TODO',
            replacement_rules   =   [  ('p1*p1', 0),
                                       ('p2*p2', 0),
                                       ('p3*p3', 0),
                                       ('p4*p4', 0),
                                       ('p1*p2', 's/2'),
                                       ('p2*p3', 't/2'),
                                       ('p1*p3', '-t/2-s/2')
                                    ],
            regulator           =   'eps',
            regulator_power     =   0,
            dimensionality      =   '4-2*eps',

        )

        """
        self.pySecDec_loop_integrand = pySecDec.loop_integral.LoopIntegralFromGraph(
            internal_lines=[[0, [1, 2]], [0, [2, 3]],
                            [0, [3, 4]], [0, [4, 1]]],
            external_lines=[['p1', 1], ['p2', 2], ['p3', 3], ['p4', 4]],
            replacement_rules=[
                ('p1*p1', 0),
                ('p2*p2', 0),
                ('p3*p3', 0),
                ('p4*p4', 0),
                ('p3*p2', 't/2'),
                ('p1*p2', 's/2'),
                ('p1*p4', 't/2'),
                ('p2*p4', '-t/2-s/2'),
                ('p3*p4', 's/2')
            ],
            regulator='eps',
            regulator_power=0,
            dimensionality='4-2*eps',
        )
"""
        self.integrand_parameters = ['s', 't']

        self.integrand_parameters_values = [
            2.*self.external_momenta[3].dot(self.external_momenta[4]),
            2.*self.external_momenta[3].dot(self.external_momenta[2]),
            self.external_momenta[1].dot(self.external_momenta[1]),
            # Let's hardcode a mass for now until we figure out how we want to specify the
            # topology
            0.0
        ]

        # Normally we should use the above, but as a quick hack to exactly match the
        # result from the canonical pySecDec example, we will use the inputs below
        self.integrand_parameters_values = [4.0, -0.75]
        
###########################################################################################
#
# First example of a local + integrated VH Attempt
#
###########################################################################################
class box1L_subtracted_VH(box1L_subtracted):
    # We plan on being able to use pySecDec integrator only for this case
    _supported_integrators = ['pySecDec']

    _pySecDecOutputName    = 'pySecDecLoopPackage_subtracted'

    def __init__(self,
                 n_loops=1,
                 external_momenta=vectors.LorentzVectorDict(),
                 phase_computed='Real',
                 # Data-structure for specifying a topology to be determined
                 topology=None,
                 **opts):

        # I am copying here all options just to make it explicit that could be used here as well.
        super(box1L_subtracted_VH, self).__init__(
            n_loops=n_loops,
            external_momenta=external_momenta,
            phase_computed=phase_computed,
            # Data-structure for specifying a topology to be determined
            topology=topology,
            **opts
        )

        # scale the external momenta to a secdec-friendly range
        max_external_momenta = max(abs(ki) for k in self.external_momenta.to_list() for ki in k)
        self.external_momenta = vectors.LorentzVectorDict({k: v / max_external_momenta for k,v in self.external_momenta.items()})
        self.integral_scaling = (1 / max_external_momenta )**4

        self.define_loop_integrand()

    def define_loop_integrand(self):
        """ Define the pySecDecObjects identifying the loop integrand."""

        
        self.pySecDec_loop_integrand = pySecDec.loop_integral.LoopIntegralFromPropagators(

            loop_momenta        =   ['k1'],
            external_momenta    =   ['p1', 'p2', 'p3', 'p4'],
            Lorentz_indices     =   ['mu',],
            propagators         =   ['k1**2', '(k1+p2)**2', '(k1+p2+p3)**2','(k1-p1)**2'],
            powerlist           =   [1, 1, 1, 1],
            numerator           =   '1' + '-(%s)'%self.get_local_counterterms(),
            replacement_rules   =   [  ('p1*p1', 0),
                                       ('p2*p2', 0),
                                       ('p3*p3', 0),
                                       ('p4*p4', 0),
                                       ('p1*p2', 's/2'),
                                       ('p2*p3', 't/2'),
                                       ('p1*p3', '-t/2-s/2'),
                                    ],
            regulator           =   'eps',
            regulator_power     =   0,
            dimensionality      =   '4-2*eps',
        )

        self.integrand_parameters = ['s', 't']

        self.integrand_parameters_values = [
            2.*self.external_momenta[3].dot(self.external_momenta[4]),
            2.*self.external_momenta[3].dot(self.external_momenta[2]),
        ]

    def get_local_counterterms(self):
        """ Get counterterms to the current loop integrand."""

        As = ['k1(mu)*k1(mu)',
              '(k1(mu)+p2(mu))*(k1(mu)+p2(mu))',
              '(k1(mu)+p2(mu)+p3(mu))*(k1(mu)+p2(mu)+p3(mu))',
              '(k1(mu)-p1(mu))*(k1(mu)-p1(mu))'
              ]

        softs = '(%s+%s)/t + (%s+%s)/s' % (As[0], As[2], As[1], As[3])

        return softs

    def get_integrated_counterterms(self):
        """ Return instances of n_loop integrands corresponding to the integrated counterpart
        of the local counterterms used here."""
        integrated_CT_options = {
                'n_loops'                   : self.n_loops,
                'external_momenta'          : self.external_momenta,
                'integral_scaling'          : self.integral_scaling,
                'phase_computed'            : self.phase_computed,
                # Data-structure for specifying a topology to be determined
                'topology'                  : self.topology,
                # Additional name for identifying this particular integrated CT
                'integrated_CT_name'        : 'TBD',
                'soft_propagator'           : 'TBD'
        }
        integrated_CT_options_A1 = dict(integrated_CT_options)
        integrated_CT_options_A1['integrated_CT_name'] = 'softA1'
        integrated_CT_options_A1['soft_propagator']    = 'A1'
        integrated_CT_options_A2 = dict(integrated_CT_options)
        integrated_CT_options_A2['integrated_CT_name'] = 'softA2'
        integrated_CT_options_A2['soft_propagator']    = 'A2'
        integrated_CT_options_A3 = dict(integrated_CT_options)
        integrated_CT_options_A3['integrated_CT_name'] = 'softA3'
        integrated_CT_options_A3['soft_propagator']    = 'A3'
        integrated_CT_options_A4 = dict(integrated_CT_options)
        integrated_CT_options_A4['integrated_CT_name'] = 'softA4'
        integrated_CT_options_A4['soft_propagator']    = 'A4'
        return [
            box1L_integrated_CT_VH(**integrated_CT_options_A1), 
            box1L_integrated_CT_VH(**integrated_CT_options_A2),
            box1L_integrated_CT_VH(**integrated_CT_options_A3),
            box1L_integrated_CT_VH(**integrated_CT_options_A4),           
        ]

class box1L_integrated_CT_VH(box1L_integrated_CT):

    def __init__(self, *args, **opts):
        if 'soft_propagator' not in opts and opts['soft_propagator'] not in ['A1','A2','A3','A4']:
            raise BaseException("The option 'soft_propagator' must be specified to the class box1L_integrated_CT_VH.")
        self.soft_propagator = opts.pop('soft_propagator')
        super(box1L_integrated_CT_VH, self).__init__(*args, **opts)

    def define_loop_integrand(self):

        propagators = None
        numerator = None
        if self.soft_propagator == 'A1':
            propagators = ['(k1+p2)**2','(k1+p2+p3)**2','(k1-p1)**2']
            numerator   = '1/t'
        elif self.soft_propagator == 'A2':
            propagators = ['k1**2','(k1+p2+p3)**2','(k1-p1)**2']
            numerator   = '1/s'
        elif self.soft_propagator == 'A3':
            propagators = ['k1**2','(k1+p2)**2','(k1-p1)**2']
            numerator   = '1/t'
        elif self.soft_propagator == 'A4':
            propagators = ['k1**2','(k1+p2)**2', '(k1+p2+p3)**2']
            numerator   = '1/s'

        # The definition below for now simply corresponds to the original box1L
        # It is just a token for performing structural tests for now.
        self.pySecDec_loop_integrand = pySecDec.loop_integral.LoopIntegralFromPropagators(

            loop_momenta        =   ['k1'],
            external_momenta    =   ['p1', 'p2', 'p3', 'p4'],
            Lorentz_indices     =   ['mu',],
            propagators         =   propagators,
            powerlist           =   [1,1,1],
            numerator           =   numerator,
            replacement_rules   =   [  ('p1*p1', 0),
                                       ('p2*p2', 0),
                                       ('p3*p3', 0),
                                       ('p4*p4', 0),
                                       ('p1*p2', 's/2'),
                                       ('p2*p3', 't/2'),
                                       ('p1*p3', '-t/2-s/2'),
                                    ],
            regulator           =   'eps',
            regulator_power     =   0,
            dimensionality      =   '4-2*eps',
        )        

        self.integrand_parameters = ['s', 't']

        self.integrand_parameters_values = [
            2.*self.external_momenta[3].dot(self.external_momenta[4]),
            2.*self.external_momenta[3].dot(self.external_momenta[2]),
        ]

###########################################################################################
#
# First example of a one-loop diagram integrated directly in momentum space
#
###########################################################################################
class box1L_direct_integration(NLoopIntegrand):

    # We plan on being able to use pySecDec integrator only for this case
    _supported_integrators = ['Vegas3','Cuba']

    NORMALIZATION_FACTOR =  (-1.j / math.pi ** 2)

    def __init__(self,
                 n_loops            = 1,
                 external_momenta   = vectors.LorentzVectorDict(),
                 phase_computed     = 'Real',
                 # Data-structure for specifying a topology to be determined
                 topology           = None,
                 loop_momenta_generator_class = None,
                 loop_momenta_generator_options = {},
                 cpp_integrand = False,
                 channel = None,
                 **opts):

        # Offer the possibility of specifying the class implementing the loop momentum deformation
        if loop_momenta_generator_class is None:
            loop_momenta_generator_class = loop_momenta_generator.OneLoopMomentumGenerator_WeinzierlCPP

        # Create the dimensions for integrating directly in momentum space.
        # The rescaling from the unit hypercube to infinity will be performed directly
        # in the loop momentum generator.
        component_names = {0:'E', 1:'x', 2:'y', 3:'z'}
        dimensions = integrands.DimensionList([
            integrands.ContinuousDimension(
                'k1_%s' % component_names[i], lower_bound=0.0, upper_bound=1.0)
            for i in range(0,4)])

        # I am copying here all options just to make it explicit that could be used here as well.
        super(box1L_direct_integration, self).__init__(
            dimensions          = dimensions,
            n_loops             = n_loops,
            external_momenta    = external_momenta,
            phase_computed      = phase_computed,
            topology            = topology,
            **opts
        )

        # For now, to avoid having any momentum back to back, boost them out of the center of mass rest frame.
#        boost_vector = (self.external_momenta[1]+2.0*self.external_momenta[2]).boostVector()
#        for p in self.external_momenta.values():
#            p.boost(-boost_vector)

        #
        self.cpp_integrand = cpp_integrand

        self.channel = channel
        if self.channel:
            logger.info("Integrand %s, on channel %d"%(self.__class__.__name__, self.channel))
        # Add the channel to the loop momentum generator options
        loop_momenta_generator_options['channel'] = self.channel


        self.loop_momenta_generator_class = loop_momenta_generator_class
        self.loop_momenta_generator_options = loop_momenta_generator_options
        self.assign_kinematic_configuration(self.external_momenta, self.channel)

    def assign_kinematic_configuration(self, PS_point, channel):
        """ Function allowing to redefine a new kinematic configuration, re-instantiating what must be re-instantiated."""
        self.channel = channel
        self.external_momenta = PS_point
        self.loop_momenta_generator_options['channel'] = channel
        self.loop_momentum_generator = self.loop_momenta_generator_class( self.topology,
                                                        self.external_momenta, **self.loop_momenta_generator_options)

        self.integrand_cpp_interface = None
        if self.cpp_integrand:
            self.integrand_cpp_interface = IntegrandCPPinterface(self.topology, self.loop_momentum_generator.q_is)
            self.integrand_cpp_interface.set_option('CHANNEL_ID', 0 if self.channel is None else self.channel)
            self.integrand_cpp_interface.set_option('UVSQ', -1.e3j * self.loop_momentum_generator.sqrt_S)
            self.integrand_cpp_interface.set_option('S12', (self.external_momenta[1] + self.external_momenta[2]).square())
            self.integrand_cpp_interface.set_option('S23', (self.external_momenta[2] + self.external_momenta[3]).square())

        self.define_loop_integrand(self.topology)

    def setup_analytic_computation(self, pyNLoop_command, *args, **opts):
        """ Performs tasks possibly necessary (like loading a library) for computing the analytic result"""

        self.avh_oneloop_hook = utils.AVHOneLOopHook(heptools_install_dir=
            pyNLoop_command.options['heptools_install_dir'] if 'heptools_install_dir' in pyNLoop_command.options else None,
            f2py_compiler=pyNLoop_command.options['f2py_compiler'] if pyNLoop_command.options['f2py_compiler'] else 'f2py')

        if not self.avh_oneloop_hook.is_available():
            logger.warning('AVH OneLOop library could not be properly loaded. '+
                                          'Analytic results for integrand %s will not be available.'%self.nice_string())
            self.avh_oneloop_hook = None

    def is_analytic_result_available(self, *args, **opts):
        """ Checks if an analytical result exists for this integrand."""
        return (self.avh_oneloop_hook is not None)

    def get_analytic_result(self, PS_point, *args, **opts):
        """ Return the real and imaginary part of the analytical result."""
        if not self.is_analytic_result_available():
            raise IntegrandError('Analytic computation for integrand %s is not available.'%self.nice_string())
        with misc.Silence():
            return self.avh_oneloop_hook.compute_one_loop_box(PS_point, loop_propagator_masses = self.loop_propagator_masses)

    def define_loop_integrand(self, topology):
        """Possibly generates new instance attributes book-keeping some representation of the loop integrand."""

        # Eventually we should extract the loop propagator masses from the topology supplied but for now we will
        # simply hardcode it
        self.loop_propagator_masses = [0.,]*4

    def output(self, output_folder, verbosity=0, **opts):
        """ Possibly output some low-level code representation of the integrand to make
        its evaluation faster."""

        # We will use a native python implementation for now
        opts['force'] = True
        super(box1L_direct_integration, self).output(output_folder, **opts)

        return False

    def get_integrated_counterterms(self):
        """ In case local numerical subtraction is employed, then the integrand of the
        correpsonding integrated counterterms must be returned as well."""

        # An empty list with no instances of IntegratedCounterterm for now.
        return []

    def get_local_counterterms(self, l_mom):
        """ Get counterterms to the current loop integrand."""

        return 0.

    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        """ Actual evaluation of the loop integrand."""

        # In the plot_deformation command we must be able to call the integrand while overriding the phase-choice.
        user_phase_choice = None
        if 'phase' in opts:
            user_phase_choice = opts.pop('phase')

        if 'input_already_in_infinite_hyperbox' not in opts:
            # Let's put a dummy example for now
            k1_E = continuous_inputs[self.dimension_name_to_position['k1_E']]
            k1_x = continuous_inputs[self.dimension_name_to_position['k1_x']]
            k1_y = continuous_inputs[self.dimension_name_to_position['k1_y']]
            k1_z = continuous_inputs[self.dimension_name_to_position['k1_z']]

            l_moms, jacobian_weight = self.loop_momentum_generator.generate_loop_momenta((k1_E,k1_x,k1_y,k1_z))
            l_mom = l_moms[0]
        else:
            # In this case 'continuous_inputs' will be the list of loop momenta
            l_moms           = continuous_inputs
            l_mom            = l_moms[0]
            jacobian_weight  = opts['jacobian']

        chosen_integrand = 'box'
        if chosen_integrand=='validation':
            # Here are some dummy function to try in order to test the phase-space volume
            euclidian_product = sum(l_mom[i] ** 2 for i in range(4))
            # Number 1:
            #    Warning, this function still has pole in the complex plane, so it may run into issue if the path happens
            #    to wander around it (and Vegas3 may push it there). For d>2, we have
            #
            #    analytical_result_d_gt_2 = math.pi**3 * (1./d) * (Mass_scale**4) * (regulator**( (2./d)-1 )) * (1./math.sin(2.*math.pi/d))
            #
            #    for d = 2, we have:
            #    analytical_result_d_eq_2 = -(1./2) * math.pi**2 * math.log(M_regulator)
            #
            Mass_scale = 1.
            regulator = 50.
            d = 3
            #
            # with the above parameter, the analytical result then reads for:
            #  Mass_scale = 100. result --> 3.2394732409401247722*10^8
            #  Mass_scale = 1.   result --> 3.2394732409401247722
            integrand = (1. / ((euclidian_product / (Mass_scale ** 2)) ** d + regulator))

        if chosen_integrand=='box':

            if self.integrand_cpp_interface:
                integrand = self.integrand_cpp_interface.evaluate(l_mom)
            else:
                numerator = 1.

                denoms = [((l_mom - q_i).square() - self.loop_propagator_masses[i_prop] ** 2) for i_prop, q_i in
                                                                                enumerate(self.loop_momentum_generator.q_is)]
                denominator = 1.
                for d in denoms:
                    denominator *= d
                if abs(denominator)==0.:
                    logger.critical('Exactly on-shell denominator encountered with the following inputs:'+
                                    '\n%s\nand resulting deformed loop momentum:\n%s\nSkipping this point.'%(
                                        str(continuous_inputs), str(l_mom) ))
                    return 0.
                integrand = (numerator / denominator)
                integrand -= self.get_local_counterterms(l_mom)

                # Normalize the loop integral properly
                integrand *= self.NORMALIZATION_FACTOR

                ############################################
                # Multi-channeling
                ############################################
                def compute_channel_weight(channel_id):
                    alpha = 2
                    return ( np.absolute((l_mom - self.loop_momentum_generator.q_is[channel_id]).square()) *
                            np.absolute((l_mom - self.loop_momentum_generator.q_is[channel_id + 1]).square()) )**alpha
                if self.channel is not None and self.channel > 0:
                    assert(self.channel + 1 < len(denoms))
                    channels = [1 / compute_channel_weight(i)  for i in range(0,len(denoms) - 1)]
                    MC_factor = channels[self.channel - 1] / sum(channels)
                    integrand *= MC_factor
                ############################################

        # Return a dummy function for now
        if user_phase_choice is None:
            chosen_phase = self.phase_computed
        else:
            chosen_phase = user_phase_choice

        if chosen_phase == 'Real':
            #misc.sprint("Returning real part: %e"%(( (-1.j/math.pi**2) * jacobian_weight * integrand_box).real))
            res = ( jacobian_weight * integrand).real
        elif chosen_phase == 'Imaginary':
            #misc.sprint("Returning complex part: %e"%(( (-1.j/math.pi**2) * jacobian_weight * integrand_box).imag))
            res = ( jacobian_weight * integrand).imag
        elif chosen_phase == 'All':
            #misc.sprint("Returning complex part: %e"%(( (-1.j/math.pi**2) * jacobian_weight * integrand_box).imag))
            res = ( jacobian_weight * integrand)
        else:
            raise IntegrandError("Unsupported phase computed option specified: %s"%self.phase_computed)

        if math.isnan(res.real) or math.isnan(res.imag):
            if 'input_already_in_infinite_hyperbox' not in opts:
                logger.warning('Integrand produced NaN result. Returning zero instead.'+
                                               ' Input random variables from the integrator: %s'%str(continuous_inputs))
            else:
                logger.warning('Integrand produced NaN result. Returning zero instead.'+
                                                 ' Complex momentum passed to the integrand: %s'%str(continuous_inputs))
            #raise IntegrandError('Integrand produced NaN result. Returning zero instead.'+
            #                                   ' Input random variables from the integrator: %s'%str(continuous_inputs))
            res = 0.

        return res

class box1L_direct_integration_subtracted(box1L_direct_integration):
    """ Implementation of the box1L with external momenta massless and onshell, using subtraction."""

    def setup_analytic_computation(self, pyNLoop_command, *args, **opts):
        """ Performs tasks possibly necessary (like loading a library) for computing the analytic result"""

        self.avh_oneloop_hook = utils.AVHOneLOopHook(heptools_install_dir=
            pyNLoop_command.options['heptools_install_dir'] if 'heptools_install_dir' in pyNLoop_command.options else None,
            f2py_compiler=pyNLoop_command.options['f2py_compiler'] if pyNLoop_command.options['f2py_compiler'] else 'f2py')

        if not self.avh_oneloop_hook.is_available():
            logger.warning('AVH OneLOop library could not be properly loaded. '+
                                          'Analytic results for integrand %s will not be available.'%self.nice_string())
            self.avh_oneloop_hook = None

    def is_analytic_result_available(self, *args, **opts):
        """ Checks if an analytical result exists for this integrand."""
        return ((self.avh_oneloop_hook is not None) and all(m==0. for m in self.loop_propagator_masses))

    def get_analytic_result(self, PS_point, *args, **opts):
        """ Return the real and imaginary part of the analytical result."""

        if not self.is_analytic_result_available():
            raise IntegrandError('Analytic computation for integrand %s is not available.'%self.nice_string())

        # First compute the box
        box_res_re, box_res_im = None, None
        with misc.Silence():
            box_res_re, box_res_im = self.avh_oneloop_hook.compute_one_loop_box(PS_point, loop_propagator_masses = (0.,0.,0.,0.))

        # Now subtract the pinched boxes, aka triangles
        s = (PS_point[3] + PS_point[4]).square()
        t = (PS_point[3] + PS_point[2]).square()

        with misc.Silence():
            tri_res_re, tri_res_im = self.avh_oneloop_hook.compute_one_loop_triangle( vectors.LorentzVectorDict(
                       { 1:PS_point[2], 2:PS_point[3], 3:(PS_point[1]+PS_point[4]) } ), loop_propagator_masses = (0.,0.,0.))
            box_res_re -= tri_res_re * (1. / s)
            box_res_im -= tri_res_im * (1. / s)
            tri_res_re, tri_res_im = self.avh_oneloop_hook.compute_one_loop_triangle( vectors.LorentzVectorDict(
                       { 1:PS_point[1], 2:PS_point[4], 3:(PS_point[2]+PS_point[3]) } ), loop_propagator_masses = (0.,0.,0.))
            box_res_re -= tri_res_re * (1. / s)
            box_res_im -= tri_res_im * (1. / s)
            tri_res_re, tri_res_im = self.avh_oneloop_hook.compute_one_loop_triangle( vectors.LorentzVectorDict(
                       { 1:PS_point[1], 2:PS_point[2], 3:(PS_point[3]+PS_point[4]) } ), loop_propagator_masses = (0.,0.,0.))
            box_res_re -= tri_res_re * (1. / t)
            box_res_im -= tri_res_im * (1. / t)
            tri_res_re, tri_res_im = self.avh_oneloop_hook.compute_one_loop_triangle( vectors.LorentzVectorDict(
                       { 1:PS_point[3], 2:PS_point[4], 3:(PS_point[1]+PS_point[2]) } ), loop_propagator_masses = (0.,0.,0.))
            box_res_re -= tri_res_re * (1. / t)
            box_res_im -= tri_res_im * (1. / t)

        return box_res_re, box_res_im

    def get_local_counterterms(self, l_mom):
        """ Get counterterms to the current loop integrand."""

        s = 2.*self.external_momenta[3].dot(self.external_momenta[4])
        t = 2.*self.external_momenta[3].dot(self.external_momenta[2])
        denoms = [((l_mom - q_i).square() - self.loop_propagator_masses[i_prop] ** 2) for i_prop, q_i in
                                                                           enumerate(self.loop_momentum_generator.q_is)]
        denoms_product = 1.
        for d in denoms:
            denoms_product *= d

        return ( ( denoms[0]+denoms[2] ) / t + ( denoms[1]+denoms[3] ) / s ) / denoms_product


class box1L_direct_integration_one_offshell_subtracted(box1L_direct_integration_subtracted):

    def setup_analytic_computation(self, pyNLoop_command, *args, **opts):
        """ Performs tasks possibly necessary (like loading a library) for computing the analytic result"""
        pass
    def is_analytic_result_available(self, *args, **opts):
        """ Checks if an analytical result exists for this integrand."""
        return False
    def get_analytic_result(self, PS_point, *args, **opts):
        """ Return the real and imaginary part of the analytical result."""
        raise NotImplementedError('Analytic result apparently not available for integrand %s.'%self.nice_string())

    def __init__(self, *args, **opts):
        """ Here, we only need to make sure that the CPP integrand is used."""
        super(box1L_direct_integration_one_offshell_subtracted, self).__init__(*args, **opts)

        if not self.integrand_cpp_interface:
            raise IntegrandError("The class %s requires the CPP integrand option to be enabled."%
                                                                        self.__class__.__name__)

class IntegrandCPPinterface(object):

    _debug_cpp = False
    _CPP_Weinzierl_src = pjoin(plugin_path,'Weinzierl')

    # List valid hyperparameters and their ID
    _valid_parameters = {
        'INTEGRAND_ID' : 1,
        'CHANNEL_ID' : 2,
        'UVSQ': 3,
        'S12': 4,
        'S23': 5
    }

    _topo_map = {
        'box1L_direct_integration': 1,
        'box1L_direct_integration_subtracted': 3,
        'box1L_direct_integration_one_offshell_subtracted': 4,
        'box1L_direct_integration_subtracted_uv_int': 5,
    }

    def compile_CPP_integrand_library(path):
        """Compiles the C++ integrand library if necessary."""

        logger.info('Now compiling shared library integrand_interface.so in path %s'%path)
        misc.compile(arg=['integrand_interface.so'], cwd=path)

        if not os.path.isfile(pjoin(path, 'integrand_interface.so')):
            raise IntegrandError("Could no compile C++ integrand source code in %s with command 'make'."%path)

    compile_CPP_integrand_library(_CPP_Weinzierl_src)

    _hook = ctypes.CDLL(pjoin(_CPP_Weinzierl_src,'integrand_interface.so'))
    # append Q and reset Q
    _hook.append_Q.argtypes = (ctypes.POINTER(ctypes.c_double),ctypes.c_int)
    _hook.append_Q.restype  = (ctypes.c_int)
    _hook.reset_Q.argtypes = ()
    _hook.reset_Q.restype  = (ctypes.c_void_p)
    # set optional factors
    _hook.set_factor_int.argtypes = (ctypes.c_int, ctypes.c_int)
    _hook.set_factor_int.restype = (ctypes.c_int)
    _hook.set_factor_complex.argtypes = (ctypes.c_int, ctypes.c_double, ctypes.c_double)
    _hook.set_factor_complex.restype = (ctypes.c_int)
    # evaluate
    _hook.evaluate.argtypes = (ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double))
    _hook.evaluate.restype  = (ctypes.c_int)


    def __init__(self, topology, qs):
        if topology not in self._topo_map:
            raise IntegrandError("Integrand {} not implemented in C++".format(topology))

        self.set_option('INTEGRAND_ID', self._topo_map[topology])
        self.reset_qs()
        self.set_qs(qs)

    def get_hook(self):
        """Build, start and return the C++ interface."""
        return self._hook

    #
    # Wrapper around all functions exposed in the C++ library
    #
    # ==================================================================================================================
    def reset_qs(self):
        self._hook.reset_Q()

    def register_qs(self, q):
        if self._debug_cpp: logger.debug(self.__class__.__name__+': In register_qs with q=%s'%str(q))
        dim = len(q)
        array_type = ctypes.c_double * dim
        return self._hook.append_Q(array_type(*q), dim)

    def set_option(self, option_name, option_values):
        if option_name not in self._valid_parameters:
            raise IntegrandError("Integrand option '%s' not recognized."%option_name)
        option_ID = self._valid_parameters[option_name]

        if isinstance(option_values, int):
            return_code = self._hook.set_factor_int(option_ID, option_values)
        elif isinstance(option_values, complex):
            return_code = self._hook.set_factor_complex(option_ID, option_values.real, option_values.imag)
        elif isinstance(option_values, float):
            return_code = self._hook.set_factor_complex(option_ID, option_values, 0.)
        else:
            raise IntegrandError("Unsupported type of option : '%s' set to '%s'."%(option_name,str(option_values)))
        if return_code != 0:
            raise IntegrandError("Error when setting option '%s' to '%s'."%(option_name,str(option_values)))

    def set_qs(self, qs):
        if self._debug_cpp: logger.debug(self.__class__.__name__+': In set_q_is with q_is=%s'%str(qs))
        for q_i in qs:
            if self.register_qs(q_i) != 0:
                raise IntegrandError('Error registering Qs')

    def evaluate(self, loop_momenta):
        real_loop_momenta = [l.real for l in loop_momenta]
        imag_loop_momenta = [l.imag for l in loop_momenta]

        dim = len(loop_momenta)
        array_type = ctypes.c_double * dim

        factor_real = ctypes.c_double(0.)
        factor_imag = ctypes.c_double(0.)
        return_code = self._hook.evaluate(array_type(*real_loop_momenta), array_type(*imag_loop_momenta), ctypes.byref(factor_real), ctypes.byref(factor_imag))
        if return_code != 0:
            raise IntegrandError("Error during evaluation of C++ integrand")
        return complex(factor_real.value, factor_imag.value)
