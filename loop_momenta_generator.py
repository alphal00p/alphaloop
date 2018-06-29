#####################################################
#                                                   #
#  Source file of the pyNLoop GGVV MG5aMC plugin.   #
#  Use only with consent of its authors.            #
#                                                   #
#  author: Valentin Hirschi, Ben Ruij               #
#                                                   #
#####################################################

import logging
import madgraph.integrator.vectors as vectors

from madgraph import InvalidCmd, MadGraph5Error

logger = logging.getLogger('pyNLoop.LoopMomentaGenerator')

class LoopMomentaGeneratorError(MadGraph5Error):
    """ Error for the LoopMomentaGenerator class suite."""
    pass

class LoopMomentaGenerator(object):
    """ Class for recursively using OneLoopMomentumGenerator for generaring several complex-valued loop momenta
    from random variables in the unit hypercube, following a deformation that ensures that all propagators are
    evaluated in the physical region."""

    def __init__(self, topology, **opts):
        """ Instantiate the class, specifying various aspects of the one-loop topology for which the deformation
        must be generated."""


class OneLoopMomentumGenerator(LoopMomentaGenerator):
    """ Class for generating a complex-valued one-loop momentum from random variables in the unit hypercube,
    following a deformation that ensures that all propagators are evaluated in the physical region."""

    def __init__(self, topology, external_momenta, **opts):
        """ Instantiate the class, specifying various aspects of the one-loop topology for which the deformation
        must be generated."""

        self.external_momenta = external_momenta
        # Eventually this information must be extracted from the topology, for now it is hard-coded.
        self.loop_propagator_masses = [0.,]*4

    def generate_loop_momenta(self, random_variables):
        """ From the random variables passed in argument, this the deform one-loop four-momentum in the form
         of a LorentzVector."""

        if len(random_variables)!=4 or any((r>1. or r<0.) for r in random_variables):
            raise LoopMomentaGeneratorError(
                "The function 'generate_loop_momenta' class '%s' requires exactly 4"%self.__class__.__name__+
                                                                                 " input random variables in [0., 1.].")
        rv = random_variables
        # Return some dummy momentum for now for test purposes
        return vectors.LorentzVector([
            rv[0]*(1.+2.j),
            rv[1]*(2.+4.j),
            rv[2]*(3.+6.j),
            rv[3]*(4.+8.j)]
        )

    def test_deformation(self):
        """ Validation function that tests that the deformation yields a complex part of each denominator with the right
        sign when approaching the point where this denominator becomes onshell."""
        pass