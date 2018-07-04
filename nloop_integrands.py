#####################################################
#                                                   #
#  Source file of the pyNLoop GGVV MG5aMC plugin.   #
#  Use only with consent of its authors.            #
#                                                   #
#  author: Valentin Hirschi, Ben Ruij               #
#                                                   #
#####################################################

import math
import shutil
import os
import logging

import madgraph.integrator.vectors as vectors
import madgraph.integrator.integrands as integrands
import madgraph.various.misc as misc
import loop_momenta_generator

logger = logging.getLogger('pyNLoop.Integrand')

pjoin = os.path.join

import pySecDec
from pySecDec.loop_integral import loop_package as pySecDec_loop_package

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

        super(NLoopIntegrand, self).__init__(dimensions=dimensions, **opts)

        self.dimension_name_to_position = {d.name: i for i, d in enumerate(dimensions)}

        self.n_loops = 1
        self.external_momenta = external_momenta
        self.phase_computed = phase_computed
        self.topology = topology

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
        

class HardCodedOffshellScalarBox(NLoopIntegrand):

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
        super(HardCodedOffshellScalarBox, self).__init__(
            dimensions=dimensions,
            n_loops=n_loops,
            external_momenta=external_momenta,
            phase_computed=phase_computed,
            # Data-structure for specifying a topology to be determined
            topology=topology,
            **opts
        )

    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        """ Actual evaluation of the call."""

        raise NotImplementedError

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


class box1L_onshell_massless(NLoopIntegrand):

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
                                       ('mur**2', 's/2'),
                                    ],
            regulator           =   'eps',
            regulator_power     =   0,
            dimensionality      =   '4-2*eps',
        )

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

        return softs

    def get_integrated_counterterms(self):
        """ Return instances of n_loop integrands corresponding to the integrated counterpart
        of the local counterterms used here."""
        integrated_CT_options = {
                'n_loops'                   : self.n_loops,
                'external_momenta'          : self.external_momenta,
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
                                       ('mur**2', 's/2'),
                                    ],
            regulator           =   'eps',
            regulator_power     =   0,
            dimensionality      =   '4-2*eps',
        )        

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
# First example of a one-loop diagram integrated directly in momentum space
#
###########################################################################################

class box1L_direct_integration(NLoopIntegrand):

    # We plan on being able to use pySecDec integrator only for this case
    _supported_integrators = ['Vegas3']

    def __init__(self,
                 n_loops            = 1,
                 external_momenta   = vectors.LorentzVectorDict(),
                 phase_computed     = 'Real',
                 # Data-structure for specifying a topology to be determined
                 topology           = None,
                 **opts):

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

        self.loop_momentum_generator = loop_momenta_generator.OneLoopMomentumGenerator(topology, self.external_momenta)
        self.define_loop_integrand(topology)

    def define_loop_integrand(self, topology):
        """Possibly generates new instance attributes book-keeping some representation of the loop integrand."""

        # Eventually we should extract the loop propagator masses from the topology supplied but for now we will
        # simply hardcode it
        self.loop_propagator_masses = [0.,]*4

    def output(self, output_folder, verbosity=0, **opts):
        """ Possibly output some low-level code representation of the integrand to make
        its evaluation faster."""

        # We will use a native python implementation for now
        pass
        #super(box1L_direct_integration, self).output(output_folder, **opts)

        return False

    def get_integrated_counterterms(self):
        """ In case local numerical subtraction is employed, then the integrand of the
        correpsonding integrated counterterms must be returned as well."""

        # An empty list with no instances of IntegratedCounterterm for now.
        return []

    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        """ Actual evaluation of the loop integrand."""

        # Let's put a dummy example for now
        k1_E = continuous_inputs[self.dimension_name_to_position['k1_E']]
        k1_x = continuous_inputs[self.dimension_name_to_position['k1_x']]
        k1_y = continuous_inputs[self.dimension_name_to_position['k1_y']]
        k1_z = continuous_inputs[self.dimension_name_to_position['k1_z']]

        l_moms, jacobian_weight = self.loop_momentum_generator.generate_loop_momenta((k1_E,k1_x,k1_y,k1_z))

        l_mom = l_moms[0]

        # Return a dummy function for now
        if self.phase_computed == 'Real':
            return jacobian_weight*(1./l_mom.dot(l_mom)).real
        else:
            return jacobian_weight*(1. / l_mom.dot(l_mom)).imag