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
import madgraph.various.misc as misc
from math import sqrt, cos, sin, exp
from madgraph import InvalidCmd, MadGraph5Error

# Suppress harmless lapack warning
import warnings
warnings.filterwarnings(action="ignore", module="scipy",
                        message="^internal gelsd")


class LoopMomentaGeneratorError(MadGraph5Error):
    """ Error for the LoopMomentaGenerator class suite."""
    pass


try:
    import numdifftools as nd
    import numpy as np
    import numpy.linalg as linalg
except ImportError:
    raise LoopMomentaGeneratorError("The loop momenta generator requires the numdifftools python package for" +
                                    " the numerical computation of derivatives.")

logger = logging.getLogger('pyNLoop.LoopMomentaGenerator')


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

    contour_hyper_parameters = {
        'self.M1_factor':   0.035,
        'self.M2_factor':   0.7,
        'self.M3_factor':   0.035,
        'self.gamma1':   0.7,
        'self.gamma2':   0.008,
        'self.Esoft_factor':   0.003,
    }

    def __init__(self, topology, external_momenta, **opts):
        """ Instantiate the class, specifying various aspects of the one-loop topology for which the deformation
        must be generated."""

        self.external_momenta = external_momenta.to_list()

        self.define_global_quantities_for_contour_deformation()

        # Eventually this information must be extracted from the topology, for now it is hard-coded.
        self.loop_propagator_masses = [0., ]*4

    def define_global_quantities_for_contour_deformation(self):
        """Define some global quantities independent of the loop momenta and useful for computing the deformed contour."""

        self.sqrt_S = sqrt(sum(v for i, v in enumerate(
            self.external_momenta) if i <= 1).square())

        # Define the q_i's defined from the external momenta as:
        self.q_i = [self.external_momenta[0], ]
        for i, p_i in enumerate(self.external_momenta[1:]):
            self.q_i.append(self.q_i[i]+p_i)

        self.P_plus = self.findP(plus=True)
        self.P_minus = self.findP(plus=False)

        # Characteristic scale of the process
        self.mu_P = sqrt((self.P_minus - self.P_plus).square())

        self.M1 = self.contour_hyper_parameters['self.M1_factor'] * self.sqrt_S
        self.M2 = self.contour_hyper_parameters['self.M2_factor'] * max(
            self.mu_P, self.sqrt_S)
        self.M3 = self.contour_hyper_parameters['self.M3_factor'] * max(
            self.mu_P, self.sqrt_S)
        self.gamma1 = self.contour_hyper_parameters['self.gamma1']
        self.gamma2 = self.contour_hyper_parameters['self.gamma2']
        self.Esoft = self.contour_hyper_parameters['self.Esoft_factor'] * self.sqrt_S

        self.soft_vectors = [self.Esoft * vectors.LorentzVector([1, 0, 0, 0]),
                             self.Esoft * vectors.LorentzVector([0, 1, 0, 0]),
                             self.Esoft * vectors.LorentzVector([0, 0, 1, 0]),
                             self.Esoft * vectors.LorentzVector([0, 0, 0, 1])]

    def map_to_infinite_hyperbox(self, random_variables):
        """ Maps a set of four random variables in the unit hyperbox to an infinite dimensional cube that corresponds
         to the original loop momentum integration space."""

        # Analytically compute this trivial jacobian
        # dk / dr = -1 / ( 2 * ( r - 1/2. )**2 )
        jacobian = 1.
        for rv in random_variables:
            jacobian *= ((1. / rv**2) + (1 / ((rv - 1.)**2)))

        return [vectors.LorentzVector([((1./(1.-rv)) - 1./rv) for rv in random_variables[i:i+4]])
                for i in range(0, len(random_variables), 4)], jacobian

    def generate_loop_momenta(self, random_variables):
        """ From the random variables passed in argument, this returns the one-loop four-momentum in the form
         of a list of a single LorentzVector, together with the jacobian of the deformation and conformal map."""

        if len(random_variables) != 4 or any((r > 1. or r < 0.) for r in random_variables):
            raise LoopMomentaGeneratorError(
                "The function 'generate_loop_momenta' class '%s' requires exactly 4" % self.__class__.__name__ +
                " input random variables in [0., 1.].")

        k_loops, remapping_weight = self.map_to_infinite_hyperbox(
            random_variables)

        deformed_k_loops, defomation_jacobian = self.apply_deformation(k_loops)

#        misc.sprint("Returning k_deformed= \n    %s\nwith jacobian: %f"%(
#            '\n    '.join('%s'%ki for ki in deformed_k_loops[0]),
#            remapping_weight*defomation_jacobian))
        return deformed_k_loops, remapping_weight*defomation_jacobian

    def apply_deformation(self, loop_momenta):
        """ This function delegates the deformation of the starting loop momenta passed in argument, and returns it as a Lorentz
        4-vector, along with the corresponding jacobian computed numerically."""

        # First use numdifftool to compute the jacobian matrix
        def wrapped_function(loop_momentum):
            return np.r_[self.deform_loop_momenta([vectors.LorentzVector(loop_momentum), ])[0]]

        local_point = list(loop_momenta[0])
        jacobian, info = nd.Jacobian(
            wrapped_function, full_output=True)(local_point)

        # And now compute the determinant
        jacobian_weight = abs(linalg.det(jacobian))

        if np.max(info.error_estimate) > 1.e-3:
            logger.warning(
                "Large error of %f (for which det(jac)=%f) encountered in the numerical evaluation of the Jacobian for the inputs: %s" %
                (np.max(info.error_estimate), jacobian_weight, str(loop_momenta[0])))

        deformed_k_loops = self.deform_loop_momenta(loop_momenta)

        return deformed_k_loops, jacobian_weight

    def deform_loop_momenta(self, loop_momenta):
        """ This function deforms the starting loop momentum passed in argument and returns it as a Lorentz
        4-vector. In the implementation of this class, we consider a single loop."""

        # A single loop for now
        assert len(
            loop_momenta) == 1, "This class %s can only handle a single loop momentum" % self.__class__.__name__
        k_loop = loop_momenta[0]

        c_plus, c_minus = 1, 1
        # k_i is the ith propagator, so k_loop - q_i
        for qi, mi in zip(self.q_i, self.loop_propagator_masses):
            # note the sign reversal
            c_plus *= self.h_delta(-1, k_loop - qi, mi * mi, self.M3 * self.M3)
            c_minus *= self.h_delta(1, k_loop - qi, mi * mi, self.M3 * self.M3)

        k_plus = k_loop - self.P_plus
        k_minus = k_loop - self.P_minus
        k_ext = vectors.LorentzVector([c_plus*k_plus[0] + c_minus*k_minus[0],
                                       - c_plus*k_plus[1] - c_minus*k_minus[1],
                                       - c_plus*k_plus[2] - c_minus*k_minus[2],
                                       - c_plus*k_plus[3] - c_minus*k_minus[3]])

        k_centre = 0.5*(k_plus + k_minus)

        # Warning, it may be that the square root needs an absolute value
        # v = [[0 if i == j else 0.5*(qi+qj - (mi-mj)/sqrt((qi-qj).square())*(qi-qj))
        #      for i, (qi, mi) in enumerate(zip(self.q_i, self.loop_propagator_masses))]
        #     for j, (qj, mj) in enumerate(zip(self.q_i, self.loop_propagator_masses))]

        def d(i, l, q, m):
            if l == i and m[l] == 0:
                return 1
            if (q[i] - q[l]).square() == 0 and q[i][0] < q[l][0] and m[l] == 0:
                return self.h_delta(1, k_loop - q[l],  m[l] * m[l], self.M1**2)
            if (q[i] - q[l]).square() == 0 and q[i][0] > q[l][0] and m[l] == 0:
                return self.h_delta(-1, k_loop - q[l], m[l] * m[l], self.M1**2)
            return max(self.h_delta(0, k_loop - q[l], m[l] * m[l], self.M1**2),
                       self.h_theta(-2 * (k_loop - q[l]).dot(k_loop - q[i]), self.M1**2))

        def d2(i, j, l, q, m):
            return self.h_theta((q[i]-q[j]).square() - (m[i]+m[j])**2, self.M1**2) * \
                max(self.h_delta(0, k_loop - q[l], m[l] * m[l], self.M1**2),
                    self.h_theta(-2 * (k_loop - q[l]).dot(k_loop - v[i, j]), self.M1**2))

        # construct the coefficients
        n = len(self.q_i)
        f = self.g(k_centre, self.gamma1, self.M2 * self.M2)
        c = [f, ]*n
        c2 = [[f, ]*n, ]*n
        da_plus = [1, ]*len(self.soft_vectors)
        da_min = [1, ]*len(self.soft_vectors)
        ca = [f, ]*len(self.soft_vectors)

        k_int = 0
        for i in range(n):
            for l in range(n):  # note: in paper l starts at 1
                c[i] *= d(i, l, self.q_i, self.loop_propagator_masses)

                # Deformation taking care of massive hyperbola
                # for j in range(i + 1, n):
                #    c2[i][j] *= d2(i, j, l, self.q_i, self.loop_propagator_masses)

            k_int += - c[i] * (k_loop - self.q_i[i])  # massless part

            # Deformation taking care of two hyperbolae
            # for j in range(i + 1, n):
            #    k_int += - c2[i][j] * (k_loop - v[i][j])

            # Deformation for more than two hyperbolae
            # for a, ka in enumerate(self.soft_vectors):
            #    da_plus[a] = max(self.h_delta(0, k_loop - self.q_i[l], self.loop_propagator_masses[l]**2, self.gamma2 * self.M1**2),
            #                     self.h_theta(ka.dot(2*(k_loop - self.q_i[l])), self.gamma2 * self.M1**2))
            #    da_min[a] = max(self.h_delta(0, k_loop - self.q_i[l], self.loop_propagator_masses[l]**2, self.gamma2 * self.M1**2),
            #                    self.h_theta(-ka.dot(2*(k_loop - self.q_i[l])), self.gamma2 * self.M1**2))

        # Add the soft part
        # for a, ka in enumerate(self.soft_vectors):
        #    c[a] *= da_plus[a] - da_min[a]
        #    k_int += c[a] * ka

        # compute lambda
        k0 = k_int + k_ext
        lambda_cycle = 1  # maximum lambda value

        for j in range(n):
            xj = (k0.dot(k_loop - self.q_i[j]) / k0.square())**2
            yj = (
                (k_loop - self.q_i[j]).square() - self.loop_propagator_masses[j]**2) / k0.square()**2

            if 2 * xj < yj:
                lambda_cycle = min(lambda_cycle, sqrt(yj / 4))
            elif yj < 0:
                lambda_cycle = min(lambda_cycle, sqrt(xj - yj/2))
            else:
                lambda_cycle = min(lambda_cycle, sqrt(xj - yj/4))

        N = sum(c)  # + sum(sum(x) for x in c2) + sum(abs(x) for x in ca)
        lambda_cycle = min(lambda_cycle, 1/(4 * N))

        # TODO: bound lambda by imaginary part of UV scale

        deformed_k_loop = k_loop + 1j * lambda_cycle * (k_int + k_ext)

        #misc.sprint('Started with k_loop:',k_loop)
        #misc.sprint('Deformed k_loop:', deformed_k_loop)

        return [deformed_k_loop, ]

    def test_deformation(self):
        """ Validation function that tests that the deformation yields a complex part of each denominator with the right
        sign when approaching the point where this denominator becomes onshell."""
        pass

    def findP(self, plus=True):
        """Find a vector P such that all external momenta `qs` are in
        the forward lightcone of P if `plus`, or are in the backwards lightcone
        if not `plus`.
        """

        # helper function
        def Z(x, y):
            # WARNING: watchout for configurations with input momenta that are back-to-back as then rho2() is zero
            n = 1.0 / y.rho2()
            if plus:
                return 0.5 * vectors.LorentzVector([x[0] - n * y.square() + n * y[0] ** 2,
                                                    x[1] + n * y[0] *
                                                    y[1], x[2] + n *
                                                    y[0] * y[2],
                                                    x[3] + n * y[0] * y[3]])
            else:
                return 0.5 * vectors.LorentzVector([x[0] + n * y.square() - n * y[0] ** 2,
                                                    x[1] - n * y[0] *
                                                    y[1], x[2] - n *
                                                    y[0] * y[2],
                                                    x[3] - n * y[0] * y[3]])

        vecs = list(self.q_i)
        while True:
            # filter all vectors that are in the forward/backward light-cone of another vector
            newvec = []
            for i, v in enumerate(vecs):
                for v1 in vecs[i + 1:]:
                    if (v-v1).square() > 0 and ((plus and v[0] > v1[0]) or (not plus and v[0] < v1[0])):
                        break
                else:
                    # this vector is space-like relative to all others
                    newvec.append(v)

            vecs = newvec
            if len(vecs) <= 1:
                break

            # find the pair with the smallest space-like seperation
            space_sep = [(i, i+j+1, -(v1-v2).square())
                         for i, v1 in enumerate(vecs) for j, v2 in enumerate(vecs[i+1:])]
            smallest = min(space_sep, key=lambda x: x[2])

            # replace first vector and drop the second
            vecs[smallest[0]] = Z(
                vecs[smallest[0]] + vecs[smallest[1]], vecs[smallest[0]] - vecs[smallest[1]])
            del vecs[smallest[1]]
            if len(vecs) <= 1:
                break

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
        return gamma * M / (k[0] * k[0] + k.rho2() + M)


class OneLoopMomentumGenerator_NoDeformation(OneLoopMomentumGenerator):
    """ One loop momentum generator which only applies the conformal map but no deformation. """

    def define_global_quantities_for_contour_defomration(self):
        """Nothing to do in this case. Move along."""
        pass

    def generate_loop_momenta(self, random_variables):
        """ From the random variables passed in argument, this returns the one-loop four-momentum in the form
         of a list of a single LorentzVector, together with the jacobian of the conformal map."""

        if len(random_variables) != 4 or any((r > 1. or r < 0.) for r in random_variables):
            raise LoopMomentaGeneratorError(
                "The function 'generate_loop_momenta' class '%s' requires exactly 4" % self.__class__.__name__ +
                " input random variables in [0., 1.].")

        k_loops, remapping_weight = self.map_to_infinite_hyperbox(
            random_variables)

        #misc.sprint("Returning k= \n    %s\nwith jacobian: %f"%('\n    '.join('%s'%ki for ki in k_loops[0]),remapping_weight))
        return k_loops, remapping_weight


class OneLoopMomentumGenerator_SimpleDeformation(OneLoopMomentumGenerator):
    """ One loop momentum generator which only applies the conformal map but no deformation. """

    def define_global_quantities_for_contour_defomration(self):
        """Nothing to do in this case. Move along."""
        pass

    def apply_deformation(self, loop_momenta):
        """ This function delegates the deformation of the starting loop momenta passed in argument, and returns it as a Lorentz
        4-vector, along with the corresponding jacobian computed numerically.
        This is a carbon copy of the function in the mother class, overloaded here only to allow further debug statements."""

        # First use numdifftool to compute the jacobian matrix
        def wrapped_function(loop_momentum):
            return np.r_[self.deform_loop_momenta([vectors.LorentzVector(loop_momentum), ])[0]]

        local_point = list(loop_momenta[0])
        jacobian, info = nd.Jacobian(
            wrapped_function, full_output=True)(local_point)

        # And now compute the determinant
        jacobian_weight = abs(linalg.det(jacobian))

        if np.max(info.error_estimate) > 1.e-3:
            logger.warning(
                "Large error of %f (for which det(jac)=%f) encountered in the numerical evaluation of the Jacobian for the inputs: %s" %
                (np.max(info.error_estimate), jacobian_weight, str(loop_momenta[0])))

        deformed_k_loops = self.deform_loop_momenta(loop_momenta)

        # misc.sprint("jacobian_weight=%f"%jacobian_weight)
        return deformed_k_loops, jacobian_weight

    def deform_loop_momenta(self, loop_momenta):
        """ This function deforms the starting loop momentum passed in argument and returns it as a Lorentz
        4-vector. In the implementation of this class, we consider a single loop."""

        # A single loop for now
        assert len(
            loop_momenta) == 1, "This class %s can only handle a single loop momentum" % self.__class__.__name__
        k_loop = loop_momenta[0]

        euclidian_product = sum(k_loop[i]**2 for i in range(4))

        scaling_factor = 0.1e-1

        deformed_k_loop = k_loop+1j * scaling_factor * vectors.LorentzVector([
            k_loop[0] * (cos(k_loop[0]**2/euclidian_product)),
            k_loop[1] * (k_loop[1]**2/euclidian_product),
            k_loop[2] * (sin(k_loop[2]**2/euclidian_product)),
            k_loop[3] * ((k_loop[3]**2/euclidian_product)**2),
        ])

#        deformed_k_loop = k_loop+1j * scaling_factor * vectors.LorentzVector([
#            k_loop[0] * ( 1. ),
#            k_loop[1] * ( 2. ),
#            k_loop[2] * ( 3. ),
#            k_loop[3] * ( 4. ),
#        ])

#        deformed_k_loop = k_loop

        return [deformed_k_loop, ]
