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
from math import sqrt

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
        self.loop_propagator_masses = [0., ]*4

    def generate_loop_momenta(self, random_variables):
        """ From the random variables passed in argument, this the deform one-loop four-momentum in the form
         of a LorentzVector."""

        if len(random_variables) != 4 or any((r > 1. or r < 0.) for r in random_variables):
            raise LoopMomentaGeneratorError(
                "The function 'generate_loop_momenta' class '%s' requires exactly 4" % self.__class__.__name__ +
                " input random variables in [0., 1.].")
        rv = random_variables

        Pplus = self.findP(self.external_momenta)
        Pmin = self.findP(self.external_momenta, plus=False)
        muP = sqrt((Pmin - Pplus).square())  # characteristic scale

        sqrtS = 1  # TODO: determine the centre-of-mass energy
        M3 = 0.035*max(muP, sqrtS)

        cplus, cmin = 1, 1
        # TODO: k_i is the ith propagator, so rv - q_i?
        for qi, mi in zip(self.external_momenta, self.loop_propagator_masses):
            # note the sign reversal
            cplus *= self.h_delta(-1, rv - qi, mi, M3)
            cmin *= self.h_delta(1, rv - qi, mi, M3)

        kplus = rv - Pplus
        kmin = rv - Pmin
        kext = vectors.LorentzVector([-cplus*kplus[0] - cmin*kmin[0],
                                      cplus*kplus[1] + cmin*kmin[1],
                                      cplus*kplus[2] + cmin*kmin[2],
                                      cplus*kplus[3] + cmin*kmin[3]])

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

    def findP(self, qs, plus=True):
        """Find a vector P such that all external momenta `qs` are in
        the forward lightcone of P if `plus`, or are in the backwards lightcone
        if not `plus`.
        """
        vecs = [v for v in qs]
        while len(vecs) > 1:
            # filter all vectors that are in the forward/backward light-cone of another vector
            newvec = []
            for i, v in enumerate(vecs):
                for v1 in vecs[i + 1:]:
                    if (v-v1).square() < 0 and ((plus and v[0] > v1[0]) or (not plus and v[0] < v1[0])):
                        break
                else:
                    # this vector is space-like relative to all others
                    newvec.append(v)

            # find the pair with the smallest space-like seperation
            space_sep = [(i, j, -(newvec[i]-newvec[j]).square()) for i in range(newvec)
                         for j in range(i+1, newvec)]
            smallest = min(space_sep, key=lambda x: x[2])

            def Z(x, y):
                n = 1.0 / y.rho2()
                if plus:
                    return 0.5 * vectors.LorentzVector([x[0] + n*y.square() - n*y[0]**2,
                                                        x[1] - n*y[0]*y[1], x[2] - n*y[0]*y[2], x[3] - n*y[0]*y[3]])
                else:
                    return 0.5 * vectors.LorentzVector([x[0] - n*y.square() + n*y[0]**2,
                                                        x[1] + n*y[0]*y[1], x[2] + n*y[0]*y[2], x[3] + n*y[0]*y[3]])

            # replace first vector and drop the second
            newvec[smallest[0]] = Z(
                newvec[smallest[0]] + newvec[smallest[1]], newvec[smallest[0]] - newvec[smallest[1]])
            del newvec[smallest[1]]
            vecs = newvec
        return vecs[0]

    def h_delta(self, sign, k, m, M):
        """The three helper functions h_delta-, h_delta+, and h_delta0, indicated
           by `sign`.
        """
        if sign == 0:
            v = (abs(k[0]) - sqrt(k.rho2() + m))**2
        else:
            v = (sign * k[0] - sqrt(k.rho2() + m))**2
        return v / (v + M)

    def h_theta(self, t, M):
        if t < 0:
            return 0.
        else:
            return t / (t + M)

    def g(self, k, gamma, M):
        gamma * M / (k[0] * k[0] + k.rho2() + M)
