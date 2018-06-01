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

logger = logging.getLogger('pyNLoop.Integrand')

import pySecDec
from pySecDec.loop_integral import loop_package as pySecDec_loop_package
    
class NLoopIntegrand(integrands.VirtualIntegrand):
    """ Class for implementing an N-loop integrand."""

    # Specify here which integrator this integrand supports
    _supported_integrators = ['Vegas3', 'pySecDec']

    def __init__(self,
                 dimensions         =   None,
                 n_loops            =   1,
                 external_momenta   =   vectors.LorentzVectorDict(),
                 phase_computed     =   'Real', 
                 # Data-structure for specifying a topology to be determined 
                 topology           =   None,
                 **opts):

        # If dimensions are not specified then generate a naive direct space in an hypercube
        if dimensions is None:
            dimensions = integrands.DimensionList(sum([[
                integrands.ContinuousDimension('l%d_E'%i_loop,lower_bound=0.0, upper_bound=1.0),
                integrands.ContinuousDimension('l%d_x'%i_loop,lower_bound=0.0, upper_bound=1.0),
                integrands.ContinuousDimension('l%d_y'%i_loop,lower_bound=0.0, upper_bound=1.0),
                integrands.ContinuousDimension('l%d_z'%i_loop,lower_bound=0.0, upper_bound=1.0),
             ] for i_loop in range(1,n_loops+1)],[]))

        super(NLoopIntegrand, self).__init__(dimensions=dimensions, **opts)

        self.dimension_name_to_position = {d.name : i for i,d in enumerate(dimensions)}

        self.n_loops                    = 1
        self.external_momenta           = external_momenta
        self.phase_computed             = phase_computed
        self.topology                   = topology

    def nice_string(self):
        """ Short string representation of this integrand. We could eventually add more
        attributes to the string returned. But for now simply return the class name."""
        return self.__class__.__name__

    def output(self, output_folder, force=False, **opts):
        """ Possibly output some low-level code representation of the integrand to make
        its evaluation faster."""

        self.output_folder = output_folder
        if os.path.exists(output_folder):
            if force or raw_input("Do you want to recycle existing output '%s'?\n[Y/N] > "%
                                                  output_folder).upper() not in ['N','NO']:
                logger.info("Recycling existing output '%s'."%output_folder)
                return False
            else:
                logger.warning("Folder '%s' already exists, it will be overwritten."%output_folder)
                shutil.rmtree(output_folder)
        os.mkdir(output_folder)
        return True
    
    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        """ Actual evaluation of the call."""

        raise NotImplementedError
    
class HardCodedOffshellScalarBox(NLoopIntegrand):

    # We plan on being able to use pySecDec integrator only for this case
    _supported_integrators = ['pySecDec']

    def __init__(self, 
                 n_loops            =   1,
                 external_momenta   =   vectors.LorentzVectorDict(),
                 phase_computed     =   'Real', 
                 # Data-structure for specifying a topology to be determined 
                 topology           =   None,
                 **opts):
        
        # Create the Feynman parameters dimensions. As we intend to use the pySecDec integrator
        # this has little relevance anyway.
        dimensions = integrands.DimensionList([
            integrands.ContinuousDimension('x%d'%(i_prop+1),lower_bound=0.0, upper_bound=1.0) 
                                              for i_prop in range(len(external_momenta)) ])

        # I am copying here all options just to make it explicit that could be used here as well.
        super(HardCodedOffshellScalarBox, self).__init__(
            dimensions          =   dimensions,
            n_loops             =   n_loops,
            external_momenta    =   external_momenta,
            phase_computed      =   phase_computed, 
            # Data-structure for specifying a topology to be determined 
            topology            =   topology,
            **opts
        )

    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        """ Actual evaluation of the call."""

        raise NotImplementedError

class box1L_Babis(NLoopIntegrand):
    # We plan on being able to use pySecDec integrator only for this case
    _supported_integrators = ['pySecDec']

    def __init__(self, 
                 n_loops            =   1,
                 external_momenta   =   vectors.LorentzVectorDict(),
                 phase_computed     =   'Real', 
                 # Data-structure for specifying a topology to be determined 
                 topology           =   None,
                 **opts):
        
        # Create the Feynman parameters dimensions. As we intend to use the pySecDec integrator
        # this has little relevance anyway.
        dimensions = integrands.DimensionList([
            integrands.ContinuousDimension('x%d'%(i_prop+1),lower_bound=0.0, upper_bound=1.0) 
                                              for i_prop in range(len(external_momenta)) ])

        # I am copying here all options just to make it explicit that could be used here as well.
        super(box1L_Babis, self).__init__(
            dimensions          =   dimensions,
            n_loops             =   n_loops,
            external_momenta    =   external_momenta,
            phase_computed      =   phase_computed, 
            # Data-structure for specifying a topology to be determined 
            topology            =   topology,
            **opts
        )

        self.pySecDec_loop_integrand = pySecDec.loop_integral.LoopIntegralFromPropagators(
        
            loop_momenta        = ['k1'],
            external_momenta    = ['p1','p2','p3','p4'],
            Lorentz_indices     = ['mu'],
            propagators         = ['k1**2', '(k1+p2)**2', '(k1+p2+p3)**2', '(k1-p1)**2' ],
            powerlist           = [1,1,1,1],
            numerator           = '1-'+self.get_counterterms(),
            replacement_rules   = [ ('p1*p1', 0),
                                    ('p2*p2', 0),
                                    ('p3*p3', 0),
                                    ('p4*p4', 0),
                                    ('p1*p2', 's/2'),
                                    ('p2*p3', 't/2'),
                                    ('p1*p3', '-t/2-s/2')
                                ],
            regulator           = 'eps',
            regulator_power     = 0,
            dimensionality      = '4-2*eps',
        )
        
        self.integrand_parameters        = ['s','t']
        
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
        self.integrand_parameters_values = [4.0, -0.75]
        
    def get_counterterms(self):
        """ Get counterterms to the current loop integrand."""
        
        As = [  'k1(mu)*k1(mu)', 
                '(k1(mu)+p2(mu))*(k1(mu)+p2(mu))',
                '(k1(mu)+p2(mu)+p3(mu))*(k1(mu)+p2(mu)+p3(mu))',
                '(k1(mu)-p1(mu))*(k1(mu)-p1(mu))'
             ]

        softs = '(%s+%s)/t + (%s+%s)/s'%(As[0],As[2],As[1],As[3])

        #misc.sprint('(%s)'%softs)
        #stop
        return '(%s)'%softs
        
    def output(self, output_folder, verbosity=0, **opts):
        """ Possibly output some low-level code representation of the integrand to make
        its evaluation faster."""
        
        if not super(box1L_Babis, self).output(output_folder, **opts):
            return False

        logger.info("Generating pySecDec output for integrand '%s' ..."%self.nice_string())
        with misc.chdir(self.output_folder) as cwd:
            with misc.Silence(active=(verbosity==0)):
                pySecDec_loop_package(
                    name            = 'pySecDecLoopPackage',
                    loop_integral   = self.pySecDec_loop_integrand,
                    real_parameters = self.integrand_parameters,
                    # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
                    requested_order = 0,
                    # the optimization level to use in FORM (can be 0, 1, 2, 3)
                    form_optimization_level = 2,
                    # the WorkSpace parameter for FORM
                    form_work_space = '100M',
                    # the method to be used for the sector decomposition
                    # valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
                    decomposition_method = 'iterative',
                    # if you choose ``geometric[_ku]`` and 'normaliz' is not in your
                    # $PATH, you can set the path to the 'normaliz' command-line
                    # executable here
                    #normaliz_executable='/path/to/normaliz',
                )
        
        return True

    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        """ Actual evaluation of the call."""

        raise NotImplementedError

class box1L(NLoopIntegrand):

    # We plan on being able to use pySecDec integrator only for this case
    _supported_integrators = ['pySecDec']

    def __init__(self, 
                 n_loops            =   1,
                 external_momenta   =   vectors.LorentzVectorDict(),
                 phase_computed     =   'Real', 
                 # Data-structure for specifying a topology to be determined 
                 topology           =   None,
                 **opts):
        
        # Create the Feynman parameters dimensions. As we intend to use the pySecDec integrator
        # this has little relevance anyway.
        dimensions = integrands.DimensionList([
            integrands.ContinuousDimension('x%d'%(i_prop+1),lower_bound=0.0, upper_bound=1.0) 
                                              for i_prop in range(len(external_momenta)) ])

        # I am copying here all options just to make it explicit that could be used here as well.
        super(box1L, self).__init__(
            dimensions          =   dimensions,
            n_loops             =   n_loops,
            external_momenta    =   external_momenta,
            phase_computed      =   phase_computed, 
            # Data-structure for specifying a topology to be determined 
            topology            =   topology,
            **opts
        )
        
        # Now generate the pySecDec loop integrand object
        self.pySecDec_loop_integrand = pySecDec.loop_integral.LoopIntegralFromGraph(
            internal_lines = [['m',[1,2]],[0,[2,3]],[0,[3,4]],[0,[4,1]]],
            external_lines = [['p1',1],['p2',2],['p3',3],['p4',4]],
            replacement_rules = [
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
        self.integrand_parameters        = ['s','t','s1','msq']
        
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

        logger.info("Generating pySecDec output for integrand '%s' ..."%self.nice_string())
        with misc.chdir(self.output_folder) as cwd:
            with misc.Silence(active=(verbosity==0)):
                pySecDec_loop_package(
                    name            = 'pySecDecLoopPackage',
                    loop_integral   = self.pySecDec_loop_integrand,
                    real_parameters = self.integrand_parameters,
                    # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
                    requested_order = 0,
                    # the optimization level to use in FORM (can be 0, 1, 2, 3)
                    form_optimization_level = 2,
                    # the WorkSpace parameter for FORM
                    form_work_space = '100M',
                    # the method to be used for the sector decomposition
                    # valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
                    decomposition_method = 'iterative',
                    # if you choose ``geometric[_ku]`` and 'normaliz' is not in your
                    # $PATH, you can set the path to the 'normaliz' command-line
                    # executable here
                    #normaliz_executable='/path/to/normaliz',
                )
        
        return True

    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        """ Actual evaluation of the call."""

        raise NotImplementedError

class box1L_onshell_massless(NLoopIntegrand):

    # We plan on being able to use pySecDec integrator only for this case
    _supported_integrators = ['pySecDec']

    def __init__(self, 
                 n_loops            =   1,
                 external_momenta   =   vectors.LorentzVectorDict(),
                 phase_computed     =   'Real', 
                 # Data-structure for specifying a topology to be determined 
                 topology           =   None,
                 **opts):
        
        # Create the Feynman parameters dimensions. As we intend to use the pySecDec integrator
        # this has little relevance anyway.
        dimensions = integrands.DimensionList([
            integrands.ContinuousDimension('x%d'%(i_prop+1),lower_bound=0.0, upper_bound=1.0) 
                                              for i_prop in range(len(external_momenta)) ])

        # I am copying here all options just to make it explicit that could be used here as well.
        super(box1L_onshell_massless, self).__init__(
            dimensions          =   dimensions,
            n_loops             =   n_loops,
            external_momenta    =   external_momenta,
            phase_computed      =   phase_computed, 
            # Data-structure for specifying a topology to be determined 
            topology            =   topology,
            **opts
        )
        
        # Now generate the pySecDec loop integrand object
        self.pySecDec_loop_integrand = pySecDec.loop_integral.LoopIntegralFromGraph(
            internal_lines = [[0,[1,2]],[0,[2,3]],[0,[3,4]],[0,[4,1]]],
            external_lines = [['p1',1],['p2',2],['p3',3],['p4',4]],
            replacement_rules = [
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
            regulator           = 'eps',
            regulator_power     = 0,
            dimensionality      = '4-2*eps',
        )
        self.integrand_parameters        = ['s','t']
        
        self.integrand_parameters_values = [
            2.*self.external_momenta[3].dot(self.external_momenta[4]),
            2.*self.external_momenta[3].dot(self.external_momenta[2]),
            self.external_momenta[1].dot(self.external_momenta[1])
        ]
        
        # Normally we should use the above, but as a quick hack to exactly match the
        # result from the canonical pySecDec example, we will use the inputs below
        self.integrand_parameters_values = [4.0, -0.75 ]

    def output(self, output_folder, verbosity=0, **opts):
        """ Possibly output some low-level code representation of the integrand to make
        its evaluation faster."""
        
        if not super(box1L_onshell_massless, self).output(output_folder, **opts):
            return False

        logger.info("Generating pySecDec output for integrand '%s' ..."%self.nice_string())
        with misc.chdir(self.output_folder) as cwd:
            with misc.Silence(active=(verbosity==0)):
                pySecDec_loop_package(
                    name            = 'pySecDecLoopPackage',
                    loop_integral   = self.pySecDec_loop_integrand,
                    real_parameters = self.integrand_parameters,
                    # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
                    requested_order = 0,
                    # the optimization level to use in FORM (can be 0, 1, 2, 3)
                    form_optimization_level = 2,
                    # the WorkSpace parameter for FORM
                    form_work_space = '100M',
                    # the method to be used for the sector decomposition
                    # valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
                    decomposition_method = 'iterative',
                    # if you choose ``geometric[_ku]`` and 'normaliz' is not in your
                    # $PATH, you can set the path to the 'normaliz' command-line
                    # executable here
                    #normaliz_executable='/path/to/normaliz',
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

        res =  math.sin(l1_E/l1_x) + math.cos(l1_y*l1_z) + math.tan(l1_x*l1_z)
        
        #misc.sprint(res)
        return res
