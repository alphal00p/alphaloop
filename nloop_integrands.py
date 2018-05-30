#####################################################
#                                                   #
#  Source file of the pyNLoop GGVV MG5aMC plugin.   #
#  Use only with consent of its authors.            #
#                                                   #
#  author: Valentin Hirschi, Ben Ruij               #
#                                                   #
#####################################################

import math

import madgraph.integrator.vectors as vectors
import madgraph.integrator.integrands as integrands
import madgraph.various.misc as misc

class NLoopIntegrand(integrands.VirtualIntegrand):
    """ Class for implementing an N-loop integrand."""

    def __init__(self, 
                 n_loops            =   1,
                 external_momenta   =   vectors.LorentzVectorDict(),
                 phase_computed     =   'Real', 
                 # Data-structure for specifying a topology to be determined 
                 topology           =   None,
                 **opts):

        # Generate the dimension automatically in direct space, assuming a rescaling
        # to unit hypercube for now.

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
        self.topology                   = topology

    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        """ Actual evaluation of the call."""

        raise NotImplementedError
    
class HardCodedOffshellScalarBox(NLoopIntegrand):

    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        """ Actual evaluation of the call."""

        raise NotImplementedError
    
class DummyNLoopIntegrand(NLoopIntegrand):
    
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
