#!/usr/bin/env python3

import sys
import os
from multiprocessing import Value 
import time
from pprint import pformat

import progressbar
import yaml
import shutil
import random
from scipy import optimize

import alpha_loop.integrator.phase_space_generators as PS
import madgraph.core.base_objects as base_objects
import madgraph.various.cluster as cluster
import madgraph.various.misc as misc
import re
import math

import LTD.vectors as vectors
from LTD.vectors import LorentzVector, LorentzVectorList

import alpha_loop.integrator.integrands as integrands
import alpha_loop.integrator.vegas3_integrator as vegas3

#import matplotlib.pyplot as plt
import numpy as np

import logging

pjoin = os.path.join

import copy
import math
import random

FINAL=True
INITIAL=False

DUMMY=99

logger = logging.getLogger("alpha_loop.sampler")

class HFunction(object):
    r"""
    Implements sampling from the noramlised distibution:

            PDF(t) = frac { t^\sigma } { 1 + t^2\sigma } frac { 2 \sigma } { \pi } Cos[ frac {\pi} {2 \sigma} ]
 
    The resulting CDF is expressed for generic \sigma in terms of a 2_F_1, which is not very practical, so we 
    implement it here only for \sigma=2 and \sigma=3
    """

    max_iteration = 500

    def __init__(self, sigma, debug=0):

        self.debug = debug

        self.sigma = sigma

        self.normalisation = ( (2*self.sigma) / math.pi ) * math.cos( math.pi / (2 * self.sigma) )

        self.PDF = (
            lambda t : (t**self.sigma / (1 + t**(2*self.sigma))) * self.normalisation
        )

        if self.sigma == 2:
            self.CDF = lambda t, x : (
                (1 / (2*math.pi))*(
                      math.log( 1 - ((2.*math.sqrt(2)*t)/(t*(t+math.sqrt(2))+1)) ) 
                    - 2.*math.atan(1-math.sqrt(2)*t) 
                    + 2.*math.atan(math.sqrt(2)*t+1)
                ) - x
            )
            self.CDF_prime = lambda t,x : (
                (2.*math.sqrt(2)*(t**2)) / (math.pi*((t**4)+1))
            )
            self.CDF_double_prime = lambda t,x : (
                - ( ( 4. * math.sqrt(2) * t * ( t**4 - 1 ) ) / ( math.pi * (t**4 + 1)**2 ) )
            )
        elif self.sigma == 3:
            self.CDF = lambda t,x : (
                (1 / (4*math.pi))*(
                    - 2*math.sqrt(3)*math.log(t**2+1)
                    + math.sqrt(3)*math.log(t**4-t**2+1)
                    - 6*math.atan(math.sqrt(3)-2*t)
                    - 6*math.atan(2*t+math.sqrt(3))
                    + 4*math.pi
                ) - x
            )
            self.CDF_prime = lambda t,x : (
                (3.*math.sqrt(3)*(t**3)) / (math.pi*(((t**6)+1)**2))
            )
            self.CDF_double_prime = lambda t,x : (
                - ( ( 9. * math.sqrt(3) * (t**2) * ( t**6 - 1 ) ) / ( math.pi * (t**6 + 1)**2 ) )
            )
        else:
            raise BaseException("Currently the H function is only implemented for sigma in [2,3].")

    def __call__(self, x):

        # Find the inverse point of the CDF
        bracket = [0.,1000000.0]
        if self.CDF(bracket[0],x) > 0.:
            logger.critical("ERROR: For x=%.16f, lower value %.16f for Netwons method invalid as it yields a positive value %.16f"%(x,bracket[0],self.CDF(bracket[0],x)))
            raise Exception("Incorrect seeds for Newtons method.")
        
        count=0
        while self.CDF(bracket[1],x) < 0. and count<10:
            bracket[1] *= 100
            count += 1        
        if self.CDF(bracket[1],x) < 0.:
            logger.critical("ERROR: For x=%.16f, upper value %.16f for Netwons method invalid as it yields a negative value %.16f"%(x,bracket[1],self.CDF(bracket[1],x)))
            raise Exception("Incorrect seeds for Newtons method.")

        root_result = optimize.root_scalar(self.CDF,
#            method="newton",
#            x0=0.01,x1=100.0,
            bracket=bracket,
            args=tuple([x,]),
            fprime=self.CDF_prime,
            fprime2=self.CDF_double_prime,
            maxiter=self.max_iteration)

    
        if self.debug or not root_result.converged:
            logger.debug("--------------------------------------------------------------")
            logger.debug("Numerical solver called with x=%.16e"%x)
            logger.debug("CDF(bracket[0]),CDF(bracket[1])=%.16e,%.16e"%(self.CDF(bracket[0],x),self.CDF(bracket[1],x)))
            logger.debug("CDF_prime(bracket[0]),CDF_prime(bracket[1])=%.16e,%.16e"%(self.CDF_prime(bracket[0],x),self.CDF_prime(bracket[1],x)))
            logger.debug("CDF_double_prime(bracket[0]),CDF_double_prime(bracket[1])=%.16e,%.16e"%(self.CDF_double_prime(bracket[0],x),self.CDF_double_prime(bracket[1],x)))
            logger.debug("Numerical solution found root=%.16e"%root_result.root)
            logger.debug("Number of iterations required=%d"%root_result.iterations)
            logger.debug("Number of function calls required=%d"%root_result.function_calls)
            logger.debug("Did converge? = %s"%root_result.converged)
            logger.debug("Algorithm termination cause: %s"%root_result.flag)
            logger.debug("PDF of for that t solution: %.16e"%self.PDF(root_result.root))
            logger.debug("--------------------------------------------------------------")
        
        return root_result.root, 1./self.PDF(root_result.root)

class TestHFuncIntegrand(integrands.VirtualIntegrand):
    """An integrand for this phase-space volume test."""

    def __init__(self, h_function, debug=0):
    
        super(TestHFuncIntegrand, self).__init__( integrands.DimensionList([ integrands.ContinuousDimension('x_1',lower_bound=0.0, upper_bound=1.0), ]) )

        self.h_function = h_function
        self.debug = debug

        self.n_evals = Value('i', 0)

        #if type(self.phase_space_generator).__name__ == 'SingleChannelPhasespace':
        #    self.my_random_path = self.phase_space_generator.generate_random_path()
    
    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        #if type(self.phase_space_generator).__name__ == 'SingleChannelPhasespace':
        #    PS_point, wgt, x1, x2 = self.phase_space_generator.get_PS_point(continuous_inputs,self.my_random_path)
        #else:
        self.n_evals.value += 1

        x = list(continuous_inputs)[0]
        t, wgt = self.h_function(x)


        final_res = (1./(1+t**2))*(1./(math.pi/2.)) * wgt

        if self.debug: logger.debug("Final wgt returned to integrator: %s"%str(final_res))
        if self.debug > 2: time.sleep(self.debug)

        return final_res

class DefaultALIntegrand(integrands.VirtualIntegrand):
    """An integrand for this phase-space volume test."""

    def __init__(self, rust_worker, generator, debug=0, phase='real'):

        self.generator = generator
        self.dimensions = generator.dimensions
        self.debug = debug
    
        super(DefaultALIntegrand, self).__init__( self.dimensions )

        self.rust_worker = rust_worker
        self.phase = phase
        self.n_evals = Value('i', 0)

        #if type(self.phase_space_generator).__name__ == 'SingleChannelPhasespace':
        #    self.my_random_path = self.phase_space_generator.generate_random_path()
    
        self.first_final_res = None

    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        #if type(self.phase_space_generator).__name__ == 'SingleChannelPhasespace':
        #    PS_point, wgt, x1, x2 = self.phase_space_generator.get_PS_point(continuous_inputs,self.my_random_path)
        #else:
        self.n_evals.value += 1

        xs, wgt = self.generator(continuous_inputs)
        re, im = self.rust_worker.evaluate_integrand( xs )

        if self.phase=='real':
            final_res = re*wgt
        else:
            final_res = im*wgt
        
        if self.debug: logger.debug("Final wgt returned to integrator: %.16e"%final_res)
        if hasattr(self.generator,"test_inverse_jacobian") and self.generator.test_inverse_jacobian:
            if final_res == 0.0:
                self.generator.i_count = 0
            else:
                if self.first_final_res is None:
                    self.first_final_res = final_res
                else:
                    if self.debug: logger.debug("Diff w.r.t previous final wgt   : %.16e"%(self.first_final_res-final_res))
        if self.debug > 2: time.sleep(self.debug/10.0)

        return final_res

class CustomGenerator(object):
    pass

class generator_aL(CustomGenerator):

    def __init__(self, dimensions, rust_worker, SG_info, model, h_function, hyperparameters, debug=0, **opts):

        self.rust_worker = rust_worker
        self.model = model
        self.SG_info = SG_info
        self.hyperparameters = hyperparameters
        self.dimensions = dimensions
        self.debug = debug

    def __call__(self, random_variables, **opts):
        """ 
        Generates a point using a sampling reflecting the topolgy of SG_QG0.
        It will return the point directly in x-space as well as the final weight to combine
        the result with.
        """
        return list(random_variables), 1.0

class generator_spherical(CustomGenerator):

    def __init__(self, dimensions, rust_worker, SG_info, model, h_function, hyperparameters, debug=0, **opts):

        self.debug = debug
        self.rust_worker = rust_worker
        self.model = model
        self.SG_info = SG_info
        self.hyperparameters = hyperparameters
        self.dimensions = dimensions
        self.h_function = h_function

        E_cm = math.sqrt(max(sum(LorentzVector(vec) for vec in hyperparameters['CrossSection']['incoming_momenta']).square(),0.0))
        self.E_cm = E_cm
        if debug: logger.debug("E_cm detected=%.16e"%self.E_cm)

        # This generator is actually only used for its 'self.generator.dim_name_to_position' attribute
        self.generator = PS.FlatInvertiblePhasespace(
            [0.]*2, [0.]*(SG_info['topo']['n_loops']+1), 
            beam_Es =(E_cm/2.,E_cm/2.), beam_types=(0,0),
            dimensions = self.dimensions,
        )

    def __call__(self, random_variables):
        """ 
        Generates a point using a sampling reflecting the topolgy of SG_QG0.
        It will return the point directly in x-space as well as the final weight to combine
        the result with.
        """

#        self.debug=True
        x_t = random_variables[self.generator.dim_name_to_position['t']]
        if self.debug: logger.debug('x_t=%s'%x_t)

        # We must now promote the point to full 3^L using the causal flow.
        # First transform t in [0,1] to t in [0,\inf]
#        rescaling_t, wgt_t = self.h_function(x_t)
        rescaling_t = self.E_cm * ( 1. / ( 1. - x_t ) - 1. / x_t )
        # Not clear where the factor 1./2. come from.
        wgt_t = (1./2.) * self.E_cm * ( 1. / x_t**2 + 1. / ( ( 1.-x_t ) * ( 1.-x_t ) ) )

        if self.debug: logger.debug('t, wgt_t=%s, %s'%(rescaling_t, wgt_t))


        theta = random_variables[self.generator.dim_name_to_position['x_1']]*math.pi
        phi = random_variables[self.generator.dim_name_to_position['x_2']]*2*math.pi

        rescaled_PS_point_LMB = [ [
            rescaling_t*math.sin(theta)*math.cos(phi) , 
            rescaling_t*math.sin(theta)*math.sin(phi) , 
            rescaling_t*math.cos(theta)
        ], ]

        wgt_param = rescaling_t**2 * math.sin(theta) * (2*math.pi) * math.pi

        if self.debug: logger.debug("before: %s"%rescaled_PS_point_LMB[0]) 

        # Now we inverse parameterise this point using alphaLoop (non-multichanneled internal param)
        aL_xs = []
        inv_aL_jac = 1.0
        for i_v, v in enumerate(rescaled_PS_point_LMB):
            kx, ky, kz, inv_jac = self.rust_worker.inv_parameterize(v, i_v, self.E_cm**2)
            inv_aL_jac *= inv_jac
            aL_xs.extend([kx, ky, kz])
        
#        if self.debug: logger.debug("inv parameterise: %s, %.16e"%(str([kx, ky, kz]), inv_aL_jac**-1))
#        kkx, kky, kkz, kjac = self.rust_worker.parameterize([kx, ky, kz], 0, self.E_cm**2)
#        if self.debug: logger.debug("parameterise: %s, %.16e"%(str([kkx, kky, kkz]), kjac))

        # The final jacobian must then be our param. jac together with that of t divided by the one from alphaloop.
        final_jacobian = wgt_param * wgt_t * inv_aL_jac

        if self.debug: logger.debug('xs=%s'%aL_xs)
        if self.debug: logger.debug('jacobian=%s'%final_jacobian)

        return aL_xs, final_jacobian


class generator_flat(CustomGenerator):

    def __init__(self, dimensions, rust_worker, SG_info, model, h_function, hyperparameters, debug=0, **opts):

        self.debug = debug
        self.rust_worker = rust_worker
        self.model = model
        self.SG_info = SG_info
        self.hyperparameters = hyperparameters
        self.dimensions = dimensions
        self.h_function = h_function

        E_cm = math.sqrt(max(sum(LorentzVector(vec) for vec in hyperparameters['CrossSection']['incoming_momenta']).square(),0.0))
        self.E_cm = E_cm
        if debug: logger.debug("E_cm detected=%.16e"%self.E_cm)

        self.generator = PS.FlatInvertiblePhasespace(
            [0.]*2, [0.]*(SG_info['topo']['n_loops']+1), 
            beam_Es =(E_cm/2.,E_cm/2.), beam_types=(0,0),
            dimensions = self.dimensions,
        )

        # Variables for testing that the inverse jacobian works well.
        # Set test_inverse_jacobian to True in order to debug the inverse jacobian.
        self.test_inverse_jacobian = False
        if self.test_inverse_jacobian:
            logger.debug("WARNING: Running in debug mode 'test_inverse_jacobian'.")
            time.sleep(1.0)
        self.rv = None
        self.i_count = 0

    def __call__(self, random_variables):
        """ 
        Generates a point using a sampling reflecting the topolgy of SG_QG0.
        It will return the point directly in x-space as well as the final weight to combine
        the result with.
        """
        #random_variables=[4.27617846e-01, 5.55958133e-01, 1.57593910e-01, 1.13340867e-01, 2.74746432e-01, 2.16284116e-01, 1.99368173e-04, 8.08726361e-02, 1.40842072e-01]

        if self.test_inverse_jacobian:
            if self.i_count==0:
                self.rv = copy.deepcopy(random_variables)
            else:
                random_variables = self.rv
            self.i_count += 1

        x_t = random_variables[self.generator.dim_name_to_position['t']]

        if self.debug: logger.debug('x_t=%s'%x_t)

        # We must now promote the point to full 3^L using the causal flow.
        # First transform t in [0,1] to t in [0,\inf]
        rescaling_t, wgt_t = self.h_function(x_t)

        if self.test_inverse_jacobian:
            rescaling_t /= float(self.i_count)

        if self.debug: logger.debug('t, wgt_t=%s, %s'%(rescaling_t, wgt_t))

        PS_point, wgt_param, x_1, x_2 = self.generator.get_PS_point(random_variables)
        if self.debug: logger.debug("orig wgt_param: %s"%wgt_param)

        # Correct for the 1/(2E) of each cut propagator that alphaLoop includes but which is already included in the normal PS parameterisation
        for v in PS_point[2:]:
            wgt_param *= 2*v[0]
        if self.debug: logger.debug("after wgt_param: %s"%wgt_param)

        if self.debug: logger.debug("PS point from param:\n%s"%str(LorentzVectorList(PS_point)))

        # Now rotate this point to the lMB of the topology.
        # In this case, since we generate flat, we don't care and simply take the first :-1 vectors
        # (and remove the initial external states)
        # and strip off the energies

#       START for ./run_ddx_NLO_SG_QG2.sh        
#        PS_point[2]=-PS_point[2]
#       END for ./run_ddx_NLO_SG_QG2.sh        

#       START for run_jjj_NLO_SG_QG36.sh
        PS_point[2]=-PS_point[2]
#       END for run_jjj_NLO_SG_QG36.sh

        PS_point_LMB = [ list(v[1:]) for v in PS_point[2:-1] ]
        if self.debug: logger.debug("PS point in LMB:\n%s"%pformat(PS_point_LMB))

        # Apply the rescaling
        rescaled_PS_point_LMB = [ [ rescaling_t*ki for ki in k ] for k in PS_point_LMB]

        # Exclude normalising func when testing the inverse jacobian.
        if not self.test_inverse_jacobian:
            normalising_func = self.h_function.PDF(rescaling_t)
        else:
            normalising_func = 1.0
        if self.debug: logger.debug('normalising_func=%s'%str(normalising_func))

        # Now we inverse parameterise this point using alphaLoop (non-multichanneled internal param)
        aL_xs = []
        # The rescaling jacobian comes both the rescaling change of variable and solving the delta function
        # TODO adjust so as to align with the particular LMB chosen for this topology.
        # TODO massive case not supported, issue a crash and warning in this case.
        dependent_momentum = -sum(vectors.Vector(v) for v in PS_point_LMB)

#       START for ./run_ddx_NLO_SG_QG2.sh
#        dependent_momentum = vectors.Vector(PS_point_LMB[0])-vectors.Vector(PS_point_LMB[1])
#       END for ./run_ddx_NLO_SG_QG2.sh        

#       START for run_jjj_NLO_SG_QG36.sh
        dependent_momentum = -vectors.Vector(PS_point_LMB[0])+vectors.Vector(PS_point_LMB[1])+vectors.Vector(PS_point_LMB[2])
#       END for run_jjj_NLO_SG_QG36.sh

        k_norms = [ math.sqrt(sum([ ki**2 for ki in k ])) for k in (PS_point_LMB+[list(dependent_momentum),]) ]

        if self.debug: logger.debug('h(t)=%s'%str(self.h_function.PDF(rescaling_t)))
        if self.debug: logger.debug('1/t=%s'%str(1./rescaling_t))
        if self.debug: logger.debug('h(1/t)=%s'%(self.h_function.PDF(1./rescaling_t)))
        if self.debug: logger.debug('rescaling_t**(3*len(PS_point_LMB))=%s'%str(rescaling_t**(3*len(PS_point_LMB))))
        if self.debug: logger.debug('sum(k_norms)=%s'%str(sum(k_norms)))
        if self.debug: logger.debug('1./(rescaling_t*sum(k_norms))=%s'%str(1./(rescaling_t*sum(k_norms))))

        inv_aL_jac = ( 1.0 / self.h_function.PDF(1./rescaling_t) ) * (rescaling_t**(len(PS_point_LMB)*3)) * ( rescaling_t * sum(k_norms) )

        for i_v, v in enumerate(rescaled_PS_point_LMB):
            kx, ky, kz, inv_jac = self.rust_worker.inv_parameterize(v, i_v, self.E_cm**2)
            inv_aL_jac *= inv_jac
            aL_xs.extend([kx, ky, kz])

        if self.debug: logger.debug('inv_aL_jac=%s'%str(inv_aL_jac))

        # The final jacobian must then be our param. jac together with that of t divided by the one from alphaloop.
        final_jacobian = wgt_param * wgt_t * normalising_func * inv_aL_jac

        if self.debug: logger.debug('xs=%s'%str(aL_xs))
        if self.debug: logger.debug('jacobian=%s'%str(final_jacobian))

        return aL_xs, final_jacobian

class generator_epem_a_ddx_SG_QG0(CustomGenerator):

    def __init__(self, dimensions, rust_worker, SG_info, model, h_function, hyperparameters, debug=0, **opts):

        self.rust_worker = rust_worker
        self.model = model
        self.SG_info = SG_info
        self.dimensions = dimensions
        self.h_function = h_function
        self.debug = debug

        self.topology = (
                    # s-channels first:
                    base_objects.VertexList([]),
                    # t-channels then:
                    base_objects.VertexList([
                        # The dummy vertex below is just to close the diagram and connect
                        # with the number -3 which is to be understood here as the initial state #2.
                        base_objects.Vertex({
                            'id': DUMMY, # Irrelevant
                            'legs': base_objects.LegList([
                                base_objects.Leg({
                                    'id': 11,
                                    'number': -2,
                                    'state': FINAL,
                                }),
                                base_objects.Leg({
                                    'id': 23,
                                    'number': -1,
                                    'state': FINAL,
                                }),
                                base_objects.Leg({
                                    'id': -11,
                                    'number': -3,
                                    'state': INITIAL,
                                })
                            ])
                        }),
                    ])
                )

        E_cm = math.sqrt(max(sum(LorentzVector(vec) for vec in hyperparameters['CrossSection']['incoming_momenta']).square(),0.0))


        self.generator = PS.SingleChannelPhasespace([0.]*2, [0.]*3,
                beam_Es =(E_cm/2.,E_cm/2.), beam_types=(0,0), 
                model=model, topology=self.topology, path = [[0,],[]],
                dimensions = self.dimensions
            )
        
        logger.debug("Considering the following topology:")
        logger.debug("-"*10)
        logger.debug(self.generator.get_topology_string(self.generator.topology, path_to_print=self.generator.path))
        logger.debug("-"*10)

    def __call__(self, random_variables):
        """ 
        Generates a point using a sampling reflecting the topolgy of SG_QG0.
        It will return the point directly in x-space as well as the final weight to combine
        the result with.
        """

        return [1.0,2.0,2.0], 1.0

class generator_epem_a_jjjj_SG_QG36(CustomGenerator):

    def __init__(self, dimensions, rust_worker, SG_info, model, h_function, hyperparameters, debug=0, **opts):

        self.rust_worker = rust_worker
        self.model = model
        self.SG_info = SG_info
        self.dimensions = dimensions
        self.h_function = h_function
        self.debug = debug

        # The process is e+(1) e-(2) > a > d(3) d~(4) g(5) g(6)
        # The SG_QG36 topoology assigns externals as follows:
        # LMB vec #1 = pq3 = -d(3)
        # LMB vec #2 = pq6 = +d~(4)
        # LMB vec #3 = pq10 = +g(5)
        # dependnet momentum = pq9 = -(-#1 +#2 +#3 +Q) 

        # Topology is 
        # g(5) g(6) -> -1
        # -1 d~(4) -> -2
        # -2 d(3) -> incoming "decaying" particle

        self.topology = (
            # s-channels first:
            base_objects.VertexList([
                base_objects.Vertex({
                    'id': DUMMY, # Irrelevant
                    'legs': base_objects.LegList([
                        base_objects.Leg({
                            'id': 21,
                            'number': 5,
                            'state': FINAL,
                        }),
                        base_objects.Leg({
                            'id': 21,
                            'number': 6,
                            'state': FINAL,
                        }),
                        base_objects.Leg({
                            'id': 21,
                            'number': -1,
                            'state': FINAL,
                        })
                    ])
                }),
                base_objects.Vertex({
                    'id': DUMMY, # Irrelevant
                    'legs': base_objects.LegList([
                        base_objects.Leg({
                            'id': 21,
                            'number': -1,
                            'state': FINAL,
                        }),
                        base_objects.Leg({
                            'id': 1,
                            'number': 4,
                            'state': FINAL,
                        }),
                        base_objects.Leg({
                            'id': 1,
                            'number': -2,
                            'state': FINAL,
                        })
                    ])
                }),
                base_objects.Vertex({
                    'id': DUMMY, # Irrelevant
                    'legs': base_objects.LegList([
                        base_objects.Leg({
                            'id': 1,
                            'number': -2,
                            'state': FINAL,
                        }),
                        base_objects.Leg({
                            'id': 1,
                            'number': 3,
                            'state': FINAL,
                        }),
                        base_objects.Leg({
                            'id': 22,
                            'number': -3,
                            'state': FINAL,
                        })
                    ])
                }),
            ]),
            # t-channels then:
            base_objects.VertexList([
                # The dummy vertex below is just to close the diagram and connect
                # with the number -3 which is to be understood here as the initial state #2.
                base_objects.Vertex({
                    'id': DUMMY, # Irrelevant
                    'legs': base_objects.LegList([
                        base_objects.Leg({
                            'id': 11,
                            'number': 1,
                            'state': INITIAL,
                        }),
                        base_objects.Leg({
                            'id': 22,
                            'number': -3,
                            'state': FINAL,
                        }),
                        base_objects.Leg({
                            'id': 11,
                            'number': -4,
                            'state': INITIAL,
                        })
                    ])
                }),
            ])
        )

        E_cm = math.sqrt(max(sum(LorentzVector(vec) for vec in hyperparameters['CrossSection']['incoming_momenta']).square(),0.0))
        self.E_cm = E_cm
        if debug: logger.debug("E_cm detected=%.16e"%self.E_cm)

        self.generator = PS.SingleChannelPhasespace([0.]*2, [0.]*(SG_info['topo']['n_loops']+1),
                beam_Es =(self.E_cm/2.,self.E_cm/2.), beam_types=(0,0), 
                model=model, topology=self.topology, #path = [[0,1], []],
                dimensions = self.dimensions
        )

        logger.debug("Considering the following topology:")
        logger.debug("-"*10)
        logger.debug(self.generator.get_topology_string(self.generator.topology, path_to_print=self.generator.path))
        logger.debug("-"*10)
        logger.debug("With the following path = %s"%str(self.generator.path))

#       For testing we can substitute here a flat PS generator.
#        self.generator = PS.FlatInvertiblePhasespace(
#            [0.]*2, [0.]*(SG_info['topo']['n_loops']+1), 
#            beam_Es =(E_cm/2.,E_cm/2.), beam_types=(0,0),
#            dimensions = self.dimensions
#        )

        # Variables for testing that the inverse jacobian works well.
        # Sett test_inverse_jacobian to True in order to debug the inverse jacobian.
        self.test_inverse_jacobian = False
        if self.test_inverse_jacobian:
            logger.debug("WARNING: Running in debug mode 'test_inverse_jacobian'.")
            time.sleep(1.0)
        self.rv = None
        self.i_count = 0

    def __call__(self, random_variables):
        """ 
        Generates a point using a sampling reflecting the topolgy of SG_QG36.
        It will return the point directly in x-space as well as the final weight to combine
        the result with.
        """

        #random_variables=[4.27617846e-01, 5.55958133e-01, 1.57593910e-01, 1.13340867e-01, 2.74746432e-01, 2.16284116e-01, 1.99368173e-04, 8.08726361e-02, 1.40842072e-01]


        if self.test_inverse_jacobian:
            if self.i_count==0:
                self.rv = copy.deepcopy(random_variables)
            else:
                random_variables = self.rv
            self.i_count += 1

        x_t = random_variables[self.generator.dim_name_to_position['t']]

        if self.debug: logger.debug('x_t=%s'%str(x_t))

        # We must now promote the point to full 3^L using the causal flow.
        # First transform t in [0,1] to t in [0,\inf]
        rescaling_t, wgt_t = self.h_function(x_t)

        if self.test_inverse_jacobian:
            rescaling_t /= float(self.i_count)

        if self.debug: logger.debug('t, wgt_t=%s, %s'%(rescaling_t, wgt_t))

        PS_point, wgt_param, x_1, x_2 = self.generator.get_PS_point(random_variables)

        if self.debug: logger.debug("PS point from param:\n%s"%str(LorentzVectorList(PS_point)))
        
        if self.debug: logger.debug("orig wgt_param: %s"%str(wgt_param))
        # Correct for the 1/(2E) of each cut propagator that alphaLoop includes but which is already included in the normal PS parameterisation
        for v in PS_point[2:]:
            wgt_param *= 2*v[0]
        if self.debug: logger.debug("after wgt_param: %s"%str(wgt_param))

        # Now rotate this point to the lMB of the topology.
        # In this case, since we generate flat, we don't care and simply take the first :-1 vectors
        # (and remove the initial external states)
        # and strip off the energies

        # The process is e+(1) e-(2) > a > d(3) d~(4) g(5) g(6)
        # The SG_QG36 topoology assigns externals as follows:
        # LMB vec #1 = pq3 = -d(3)
        # LMB vec #2 = pq6 = +d~(4)
        # LMB vec #3 = pq10 = +g(5)
        # dependnet momentum = pq9 = -(-#1 +#2 +#3 +Q) 
        PS_point_LMB = vectors.LorentzVectorList([])
        PS_point_LMB.append(-PS_point[2])
        PS_point_LMB.append(PS_point[3])
        PS_point_LMB.append(PS_point[4])
        PS_point_LMB.append(-PS_point[2]-PS_point[3]-PS_point[4]-vectors.LorentzVector([self.E_cm,0.,0.,0.]))

        PS_point_LMB_vec = [ list(v[1:]) for v in PS_point_LMB ]
        if self.debug: logger.debug("PS point in LMB:\n%s"%str(pformat(PS_point_LMB_vec)))

        # Apply the rescaling
        rescaled_PS_point_LMB_vec = [ [ rescaling_t*ki for ki in k ] for k in PS_point_LMB_vec]

        k_norms = [ math.sqrt(sum([ ki**2 for ki in k ])) for k in PS_point_LMB_vec ]

        # Exclude normalising func when testing the inverse jacobian.
        if not self.test_inverse_jacobian:
            normalising_func = self.h_function.PDF(rescaling_t)
        else:
            normalising_func = 1.0
        
        if self.debug: logger.debug('normalising_func=%s'%str(normalising_func))

        # Now we inverse parameterise this point using alphaLoop (non-multichanneled internal param)
        aL_xs = []

        if self.debug: logger.debug('h(t)=%s'%str(self.h_function.PDF(rescaling_t)))
        if self.debug: logger.debug('1/t=%s'%str(1./rescaling_t))
        if self.debug: logger.debug('h(1/t)=%s'%str(self.h_function.PDF(1./rescaling_t)))
        if self.debug: logger.debug('rescaling_t**(3*(len(PS_point_LMB)-1))=%s'%str(rescaling_t**(3*(len(PS_point_LMB)-1))))
        if self.debug: logger.debug('sum(k_norms)=%s'%str(sum(k_norms)))
        if self.debug: logger.debug('1./(rescaling_t*sum(k_norms))=%s'%str(1./(rescaling_t*sum(k_norms))))

        inv_aL_jac = ( 1.0 / self.h_function.PDF(1./rescaling_t) ) * (rescaling_t**((len(PS_point_LMB)-1)*3)) * ( rescaling_t * sum(k_norms) )

        for i_v, v in enumerate(rescaled_PS_point_LMB_vec[:-1]):
            kx, ky, kz, inv_jac = self.rust_worker.inv_parameterize(v, i_v, self.E_cm**2)
            inv_aL_jac *= inv_jac
            aL_xs.extend([kx, ky, kz])

        if self.debug: logger.debug('inv_aL_jac=%s'%str(inv_aL_jac))

        # The final jacobian must then be our param. jac together with that of t divided by the one from alphaloop.
        final_jacobian = wgt_param * wgt_t * normalising_func * inv_aL_jac

        if self.debug: logger.debug('xs=%s'%str(aL_xs))
        if self.debug: logger.debug('jacobian=%s'%str(final_jacobian))

        return aL_xs, final_jacobian