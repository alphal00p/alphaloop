#!/usr/bin/env python3

import sys
import os
from multiprocessing import Value 
import time
from pprint import pprint, pformat
import time
import random

import progressbar
import yaml
import shutil
import alpha_loop.utils as utils
import random
from scipy import optimize
import warnings
import traceback

import alpha_loop.integrator.phase_space_generators as PS
import madgraph.core.base_objects as base_objects
import madgraph.various.cluster as cluster
import madgraph.various.misc as misc
from madgraph import InvalidCmd, MadGraph5Error, MG5DIR, ReadWrite
import re
import math

import LTD.vectors as vectors
import LTD.ltd_utils as ltd_utils
from LTD.vectors import LorentzVector, LorentzVectorList

import alpha_loop.integrator.integrands as integrands
import alpha_loop.integrator.vegas3_integrator as vegas3

try:
    import scipy.optimize as optimize
    scipy_optimize_loaded = True
    import warnings
    warnings.filterwarnings('error')
except:
    scipy_optimize_loaded = False

#import matplotlib.pyplot as plt
import numpy as np

import logging

sampler_path = os.path.dirname(os.path.realpath( __file__ ))

pjoin = os.path.join

import copy
import math
import random

FINAL=True
INITIAL=False

DUMMY=99

logger = logging.getLogger("alpha_loop.sampler")

def transform_spherical_to_cartesian(rtp):
    cartesian = np.zeros(rtp.shape)
    jacobians = np.zeros(len(rtp))
    cartesian[:,0] = rtp[:,0]*np.sin(rtp[:,1])*np.cos(rtp[:,2])
    cartesian[:,1] = rtp[:,0]*np.sin(rtp[:,1])*np.sin(rtp[:,2])
    cartesian[:,2] = rtp[:,0]*np.cos(rtp[:,1])
    jacobians[:] = rtp[:,0]**2*np.sin(rtp[:,1])
    return cartesian, jacobians 

def transform_cartesian_to_spherical(xyz):
    spherical = np.zeros(xyz.shape)
    inv_jacobians = np.zeros(len(xyz))
    xy = xyz[:,0]**2 + xyz[:,1]**2
    spherical[:,0] = np.sqrt(xy + xyz[:,2]**2)
    spherical[:,1] = np.arctan2(np.sqrt(xy), xyz[:,2])
    spherical[:,2] = np.arctan2(xyz[:,1], xyz[:,0])
    inv_jacobians[:] = 1./(spherical[:,0]**2*np.sin(spherical[:,1]))
    return spherical, inv_jacobians 

def lin_conformal(xs, scale):
    ks = np.zeros(xs.shape)
    jacobians = np.zeros(xs.shape)
    ks[:] = scale*(xs[:]/(1.-xs[:]))
    jacobians[:] = scale*(1./(1.-xs[:])**2)
    return ks, jacobians

def lin_inverse_conformal(ks, scale):
    xs = np.zeros(ks.shape)
    inv_jacobians = np.zeros(xs.shape)
    xs[:] = ks[:]/(ks[:]+scale)
    inv_jacobians[:] = 1./(scale*(1./(1.-xs[:])**2))
    return xs, inv_jacobians

def log_conformal(xs,scale):
    ks = np.zeros(xs.shape)
    jacobians = np.zeros(xs.shape)
    ks[:] = scale*np.log(1./(1.-xs[:]))
    jacobians[:] = scale/(1.-xs[:])
    return ks, jacobians

def log_inverse_conformal(ks,scale):
    xs = np.zeros(ks.shape)
    inv_jacobians = np.zeros(xs.shape)
    exponentials = np.exp(ks[:]/scale)
    xs[:] = (exponentials[:]-1.)/exponentials[:]
    inv_jacobians[:] = (1.-xs[:])/scale
    return xs, inv_jacobians

class SamplerError(MadGraph5Error):
    pass

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

        self.exact_PDF = self.PDF

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

    def inverse_sampling(self, t):
        """ Inverse sampling, which is trivial in this case given that we do an exact inversion, but it may not always be so.
        Returns the corresponding sampling x in [0,1] and the inverse of the Jacobian.
        """

        return self.CDF(t, 0), self.PDF(t)

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

        PDF = self.PDF(root_result.root)
        if self.debug>2 or not root_result.converged or PDF==0.:
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
        
        return root_result.root, 1./PDF

class FlatHFunction(object):

    def __init__(self, *args, debug=0):

        self.debug = debug
        
        # Constant PDF in x gives 1/(1+t)^2 PDF in t
        self.normalisation = 1

        self.PDF = (
            lambda t : 1./(1.+t)**2
        )
        #self.exact_PDF = self.PDF

        normalisation = 3.5748532344428743120454466520221
        self.exact_PDF = (
            lambda t : math.exp( - (1. + t**2) / t ) * normalisation
        )

    def inverse_sampling(self, t):
        """ Inverse sampling, which is trivial in this case given that we do an exact inversion, but it may not always be so.
        Returns the corresponding sampling x in [0,1] and the inverse of the Jacobian.
        """

        x_t = (t / (1. + t) )

        return x_t, (1./(1.+t)**2)

    def __call__(self, x):

        t = x/(1.-x)

        return t, 1./(1./(1.+t)**2)

class DiscreteExponentialHFunction(object):
    r"""
    Implements sampling from the noramlised distibution:

            PDF(t) = exp(-(t^2+1)/t) * (1 / 2 BesselK[1, 2])
 
    The resulting CDF cannot be computed analytically, so we instead pick here a discretised version, with a discretisation 
    that can be specified upon instantiation
    """

    _DISCRETE_BINS_FILE = 'ExponentialHFunctionDiscretizationData_300bins.csv'
    #_DISCRETE_BINS_FILE = 'ExponentialHFunctionDiscretizationData_10bins.csv'

    def __init__(self, sigma, debug=0):

        self.debug = debug

        self.sigma = sigma

        if self.sigma != 1:
            raise BaseException("Currently the DiscreteExponentialHFunction is only implemented for sigma=1.")

        
        # Hard-coding 1./(2 BesselK[1, 2]) here.
        self.normalisation = 3.5748532344428743120454466520221

        self.exact_PDF = (
            lambda t : math.exp( - (1. + t**2) / t ) * self.normalisation
        )

        discretization_data_file = os.path.join(sampler_path,self._DISCRETE_BINS_FILE)
        self.discrete_bins = []
        try:
            headers = []
            for i_line, line in enumerate(open(discretization_data_file,'r').readlines()):
                if i_line == 0:
                    headers = [ header.replace('"','') for header in line.strip('\n').split(',') ]
                    continue
                self.discrete_bins.append({key : float(value) for key, value in zip(headers, line.strip('\n').split(',')) })
        except Exception as e:
            logger.info("Could not read exponential H-function discretization in file %s. Error: %s"%(discretization_data_file, str(e)))

        logger.info("Successfully loaded %d discretized bin data from file %s."%(
            len(self.discrete_bins), discretization_data_file
        ))

    def inverse_sampling(self, t):
        """ Inverse sampling, which is trivial in this case given that we do an exact inversion, but it may not always be so.
        Returns the corresponding sampling x in [0,1] and the inverse of the Jacobian.
        """

        x_t = (t / (1. + t) )

        selected_bin = None
        for a_bin in self.discrete_bins:
            if a_bin['left_edge'] <= x_t <= a_bin['right_edge']:
                selected_bin = a_bin
                break
        if selected_bin is None:
            raise Exception("Could not find a discretized bin for t=%.16f, x_t=%.16f"%(t, x_t))

        a = 0.5*selected_bin['a']
        b = selected_bin['b']
        c = -0.5*selected_bin['a']*selected_bin['left_edge']**2 - selected_bin['b']*selected_bin['left_edge'] + selected_bin['CDF_left'] 

        CDF_t = a*(x_t**2) + b*(x_t) + c 

        PDF_t = ( selected_bin['a'] * x_t + selected_bin['b'] ) / ((1+t)**2)

        return CDF_t, PDF_t

    def __call__(self, x):

        # First find the relevant bin
        selected_bin = None
        for a_bin in self.discrete_bins:
            if a_bin['CDF_left'] <= x <= a_bin['CDF_right']:
                selected_bin = a_bin
                break
        if selected_bin is None:
            raise Exception("Could not find a discretized bin for random variable %.16f."%x)

        a = 0.5*selected_bin['a']
        b = selected_bin['b']
        c = -0.5*selected_bin['a']*selected_bin['left_edge']**2 - selected_bin['b']*selected_bin['left_edge'] + selected_bin['CDF_left'] - x

        x_t = ( -b + math.sqrt(b**2-4*a*c) ) / ( 2* a )

        t = x_t / ( 1 - x_t )

        PDF_t = ( selected_bin['a'] * x_t + selected_bin['b'] ) / ((1+t)**2)

        return t, 1./PDF_t

class TestHFuncIntegrand(integrands.VirtualIntegrand):
    """An integrand for this phase-space volume test."""

    def __init__(self, h_function, debug=0):
    
        super(TestHFuncIntegrand, self).__init__( integrands.DimensionList([ integrands.ContinuousDimension('x_1',lower_bound=0.0, upper_bound=1.0), ]) )

        self.h_function = h_function

        # Integrate trivial one
        self.trial_normalised_integrand = lambda t: self.h_function.exact_PDF(t)

        # Integrate similarly shaped integrand
        #HF2 = HFunction(2)
        #HF3 = HFunction(3)
        #self.trial_normalised_integrand = lambda t: (HF2.PDF(t)+HF3.PDF(t))/2.

        self.debug = debug

        self.n_evals = Value('i', 0)

    def __call__(self, continuous_inputs, discrete_inputs, **opts):

        x = list(continuous_inputs)[0]
        t, wgt = self.h_function(x)

        #final_res = (1./(1+t**2))*(1./(math.pi/2.)) * wgt
        final_res = self.trial_normalised_integrand(t) * wgt

        if self.debug: logger.debug("Final wgt returned to integrator: %s"%str(final_res))
        if self.debug > 2: time.sleep(self.debug)

        self.update_evaluation_statistics([x,],final_res)

        return final_res

class DefaultALIntegrand(integrands.VirtualIntegrand):
    """An integrand for this phase-space volume test."""

    def __init__(self, rust_worker, generator, debug=0, phase='real', frozen_momenta=None):

        self.generator = generator
        self.dimensions = generator.dimensions
        self.debug = debug
        self.frozen_momenta = frozen_momenta

        super(DefaultALIntegrand, self).__init__( self.dimensions )

        self.rust_worker = rust_worker
        self.phase = phase
        self.n_evals = Value('i', 0)

        #if type(self.phase_space_generator).__name__ == 'SingleChannelPhasespace':
        #    self.my_random_path = self.phase_space_generator.generate_random_path()
    
        self.first_final_res = None

    def __call__(self, continuous_inputs, discrete_inputs, **opts):

        xs, wgt = self.generator(continuous_inputs)

        # ks = []
        # for i_v, v in enumerate([xs[0:3],xs[3:6],xs[6:9]]):
        #     kx, ky, kz, jac = self.rust_worker.parameterize(list(v), i_v, 125.0**2)
        #     ks.append([kx, ky, kz])

        # all_aL_xs = []
        # for signs in [
        #         (1.,1.,1.),
        #         (1.,1.,-1.),
        #         (1.,-1.,1.),
        #         (1.,-1.,-1.),
        #         (-1.,1.,1.),
        #         (-1.,1.,-1.),
        #         (-1.,-1.,1.),
        #         (-1.,-1.,-1.),
        #     ]:
        #     all_aL_xs.append( ( 1.0, []) )
        #     for i_v, v in enumerate([ 
        #             [ signs[0]*ks[0][0], signs[0]*ks[0][1], signs[0]*ks[0][2] ],
        #             [ signs[1]*ks[1][0], signs[1]*ks[1][1], signs[1]*ks[1][2] ],
        #             [ signs[2]*ks[2][0], signs[2]*ks[2][1], signs[2]*ks[2][2] ],
        #         ]):
        #         kx, ky, kz, inv_jac = self.rust_worker.inv_parameterize(list(v), i_v, 125.0**2)
        #         all_aL_xs[-1][-1].extend([kx, ky, kz])

        all_aL_xs = [(1.,xs),]
        res = complex(0., 0.)
        for overall_sign, aL_xs in all_aL_xs:
            if self.debug: logger.debug("aL input xs:: %s"%str(aL_xs))
            aL_res = complex( *self.rust_worker.evaluate_integrand( aL_xs ) )
            if self.debug: logger.debug("aL res :: %s"%str(aL_res))
            res += overall_sign*aL_res
        res /= float(len(all_aL_xs))

        re, im = res.real, res.imag

        if self.phase=='real':
            final_res = re*wgt
        else:
            final_res = im*wgt
        
        if self.debug: logger.debug("Final wgt returned to integrator: %.16e"%final_res)

        # Remove the padded xs for frozen momenta if applicable
        if self.frozen_momenta is not None:
            xs_without_frozen_momenta = xs[:-3*len(self.frozen_momenta['out'])]
        else:
            xs_without_frozen_momenta = xs

        self.update_evaluation_statistics(xs_without_frozen_momenta,final_res)

        return {'I':final_res}


class DummyFrozenPSGenerator(object):

    def __init__(self, beam_Es, frozen_momenta, rust_worker, n_loops, E_cm, *args, **opts):
        self.frozen_momenta = dict(frozen_momenta)

        self.dummy_initial_momenta = [ 
            vectors.LorentzVector([beam_Es[0],0.,0.,beam_Es[0]]),
            vectors.LorentzVector([beam_Es[0],0.,0.,-beam_Es[0]])
        ]
    
        # Pad xs with frozen momenta if necessary
        self.overall_inv_jac = 1.
        self.PS_xs = []
        for i_v, v in enumerate([v[1:] for v in self.frozen_momenta['out']]):
            xx, xy, xz, inv_jac = rust_worker.inv_parameterize(list(v), n_loops+i_v, E_cm**2)
            self.PS_xs.extend([xx, xy, xz])
            self.overall_inv_jac *= inv_jac*(2.0*math.pi)**4

        self.frozen_momenta['in'] = [ vectors.LorentzVector(v) for v in self.frozen_momenta['in'] ]
        self.frozen_momenta['out'] = [ vectors.LorentzVector(v) for v in self.frozen_momenta['out'] ]
        self.frozen_momenta['out'].append( sum(self.frozen_momenta['in'])-sum(self.frozen_momenta['out']) )

    def get_PS_point(self,*args, **opts):
        return vectors.LorentzVectorList(self.dummy_initial_momenta+self.frozen_momenta['out']), self.overall_inv_jac, 1., 1.

    def invertKinematics(self, *args, **opts):
        return self.PS_xs, self.overall_inv_jac

class DummyFrozenHFunction(object):

    def __init__(self, *args, **opts):

        normalisation = 3.5748532344428743120454466520221
        self.exact_PDF = (
            lambda t :  1.
        )

    def inverse_sampling(self, t):

        return 1., 1.

    def __call__(self, x):

        return 1., 1.

class AdvancedIntegrand(integrands.VirtualIntegrand):

    def __init__(self, rust_worker, SG, model, h_function, hyperparameters, channel_for_generation = None, external_phase_space_generation_type="flat", debug=0, phase='real', 
        selected_cut_and_side=None, selected_LMB=None, show_warnings=True, return_individual_channels=False, frozen_momenta=None,
        **opts):
        """ Set channel_for_generation to (selected_cut_and_side_index, selected_LMB_index) to disable multichaneling and choose one particular channel for the parameterisation. """

        self.debug = debug
        self.show_warnings = show_warnings
        self.rust_worker = rust_worker
        self.model = model
        self.SG = SG
        self.phase = phase
        self.hyperparameters = hyperparameters
        self.conformal_map = lin_conformal
        self.inverse_conformal_map = lin_inverse_conformal
        self.channel_for_generation = channel_for_generation
        self.return_individual_channels = return_individual_channels
        self.frozen_momenta = frozen_momenta

        if self.frozen_momenta is None:
            self.h_function = h_function
        else:
            self.h_function = DummyFrozenHFunction()

        self.selected_cut_and_side=selected_cut_and_side
        self.selected_LMB=selected_LMB

        self.n_aborted_evals = Value('i', 0)

        if external_phase_space_generation_type not in ['flat','advanced']:
            raise InvalidCmd("External phase-space generation mode not supported: %s"%external_phase_space_generation_type)
        self.external_phase_space_generation_type = external_phase_space_generation_type

        # Built external momenta. Remember that they appear twice.
        self.n_initial = len(hyperparameters['CrossSection']['incoming_momenta'])
        self.external_momenta = [ vectors.Vector(v[1:]) for v in self.hyperparameters['CrossSection']['incoming_momenta'] ]
        self.external_momenta.extend(self.external_momenta)

        self.dimensions = integrands.DimensionList()
        
        if self.frozen_momenta is None:
            # The first dimension will always be for the upscaling 't' variable
            self.dimensions.append(integrands.ContinuousDimension('x_t',lower_bound=0.0, upper_bound=1.0))
            # Then the remaining ones will first be for the final-state mappings (the number depends on the particular Cutkosky cut)
            # and of the remaining loop variables
            self.dimensions.extend([ 
                integrands.ContinuousDimension('x_%d'%i_dim,lower_bound=0.0, upper_bound=1.0) 
                for i_dim in range(1,SG['topo']['n_loops']*3) 
            ])
        else:
            self.dimensions.extend([ 
                integrands.ContinuousDimension('x_%d'%i_dim,lower_bound=0.0, upper_bound=1.0) 
                for i_dim in range(1,(SG['topo']['n_loops']-len(self.frozen_momenta['out']))*3+1) 
            ])

        super(AdvancedIntegrand, self).__init__( self.dimensions )

        E_cm = self.SG.get_E_cm(self.hyperparameters)        
        self.E_cm = E_cm
        if debug: logger.debug("E_cm detected=%.16e"%self.E_cm)

        # We must now build the reduced topology of the supergraph to obtain channels
        self.SG.set_integration_channels()
        if self.debug>2: 
            selected_channels = [(i_channel, channel_info) for i_channel, channel_info in enumerate(self.SG['SG_multichannel_info']) if 
                                    self.external_phase_space_generation_type != 'flat' or channel_info['side']=='left']
            logger.debug("\n>>> Integration channel analysis (%d channels):"%len(selected_channels))
            for i_cut, cut_info in enumerate(self.SG['cutkosky_cuts']):
                logger.info("\n>> Multichannel info for cut #%d:\n%s"%(i_cut,pformat(cut_info['multichannel_info'])))
            logger.info("\n>>> Overall multichannel info (%d channels):"%len(selected_channels))
            for i_channel, channel_info in selected_channels:
                logger.info("\n>>Info of channel #%d:\n%s"%(
                    i_channel, pformat(channel_info)
                ))

        # First we must regenerate a TopologyGenerator instance for this supergraph
        # We actually won't need it!
        #edges_list = SG['topo_edges']
        #SG_topo_gen = ltd_utils.TopologyGenerator(
        #    [e[:-1] for e in edges_list],
        #    powers = { e[0] : e[-1] for e in edges_list }
        #)

        # Now assign a tree-level structure generator to each channel 
        self.multichannel_generators = {}

        self.edges_PDG = { e[0] : e [1] for e in self.SG['edge_PDGs'] }

        unknown_particles = set([ PDG for edge_name, PDG in self.edges_PDG.items() if self.model.get_particle(PDG) is None ])
        if len(unknown_particles)>0:
            raise SamplerError("The particle PDGs %s are unknown in the model '%s' loaded. Make sure you loaded the same model that was used for generating the output."%(
                str(unknown_particles), self.model.get('name')
            ))
        try:
            self.edge_masses = { edge_name : self.model['parameter_dict'][
                    self.model.get_particle(self.edges_PDG[edge_name]).get('mass')
            ].real for edge_name in self.edges_PDG }
        except Exception as e:
            raise SamplerError("Some of the particle masses used in that supergraph are not defined in the model. Make sure you loaded the proper model and with the proper restrict card.")

        for i_channel, SG_channel_info in enumerate(self.SG['SG_multichannel_info']):

            tree_topology_info = self.SG['cutkosky_cuts'][
                    SG_channel_info['cutkosky_cut_id']
                ]['multichannel_info']['tree_topologies'][
                    SG_channel_info['side']
                ]
            initial_state_masses = [
                self.edge_masses[edge_name] for edge_name, direction in tree_topology_info['incoming_edges']
            ]
            final_state_masses = [
                self.edge_masses[edge_name] for edge_name, direction in tree_topology_info['outgoing_edges']
            ]
            
            n_external_state_dofs = len(SG_channel_info['independent_final_states'])*3-1

            if self.n_initial==1:
                # For 1>N topologies we fake a 2>N massless one:
                generator_initial_state_masses = [0., 0.]
                generator_beam_Es = (initial_state_masses[0]/2.,initial_state_masses[0]/2.)
            else:
                generator_initial_state_masses = initial_state_masses
                generator_beam_Es = (E_cm/2.,E_cm/2.)
            
            if self.frozen_momenta is not None:
                
                self.multichannel_generators[i_channel] = DummyFrozenPSGenerator(
                    generator_beam_Es, self.frozen_momenta, self.rust_worker, 
                    (SG['topo']['n_loops']-len(self.frozen_momenta['out'])), self.E_cm)

            elif self.external_phase_space_generation_type == 'flat':
                
                self.multichannel_generators[i_channel] = PS.FlatInvertiblePhasespace(
                        generator_initial_state_masses, final_state_masses, 
                        beam_Es=generator_beam_Es, beam_types=(0,0),
                        dimensions = self.dimensions.get_dimensions(['x_%d'%i_dim for i_dim in range(1,n_external_state_dofs+1)])
                )

            elif self.external_phase_space_generation_type == 'advanced':
        
                # TODO in complicated cases of Breit-Wigner competition we want to also MC over the different paths of 
                # building the *same* topology (i.e. think g g > h > w w > 4-leptons).
                # But for now we can hard-code this to be the first possible path detected, which we get by specifying None below.
                selected_generator_path = None

                #generator_initial_state_masses = [125.0,]
                #generator_beam_Es = tuple([125.0,])
                if self.n_initial==1:
                    # For emulating a 2>1 in the SingleChannelPhasespace class, we must move the existing "closing t-channel" into an s-channel 
                    # and add a fake "closing" t-channel from 2 photons for instance going into the one particle decaying 
                    selected_topology = list(copy.deepcopy(tree_topology_info['s_and_t_propagators']))
                    initial_state_leg = None
                    for leg in tree_topology_info['s_and_t_propagators'][1][0].get('legs'):
                        if leg.get('state') == INITIAL:
                            initial_state_leg = leg
                            initial_state_leg.set('state',FINAL)
                            break
                    selected_topology[0].append(tree_topology_info['s_and_t_propagators'][1][0])
                    selected_topology[1] = base_objects.VertexList([base_objects.Vertex({
                        'id': DUMMY, # Irrelevant
                        'legs': base_objects.LegList([
                            base_objects.Leg({
                                'id': 22,
                                'number': 1,
                                'state': INITIAL,
                            }),
                            base_objects.Leg({
                                'id': 22,
                                'number': initial_state_leg.get('number'),
                                'state': FINAL,
                            }),
                            base_objects.Leg({
                                'id': 22,
                                'number': initial_state_leg.get('number')-1,
                                'state': INITIAL,
                            })
                        ])
                    })])
                    selected_topology = tuple(selected_topology)
                else:
                    selected_topology = tree_topology_info['s_and_t_propagators']

                self.multichannel_generators[i_channel] =  PS.SingleChannelPhasespace(
                        generator_initial_state_masses, final_state_masses, 
                        beam_Es=generator_beam_Es, beam_types=(0,0),
                        model=model, topology=selected_topology, path = selected_generator_path,
                        dimensions = self.dimensions.get_dimensions(['x_%d'%i_dim for i_dim in range(1,n_external_state_dofs+1)])
                    )

    def loop_parameterisation(self, xs):

        loop_jac = 1.0
        n_loops = len(xs)//3
        if n_loops == 0:
            return [], loop_jac

        # First the radii
        radii, radii_jacobians = self.conformal_map(np.array([
            xs[i] for i in range(0, len(xs), 3)
        ]), self.E_cm)
        loop_jac *= np.prod(radii_jacobians)
        # The theta angles then
        thetas = [ math.pi*xs[i+1] for i in range(0, len(xs), 3) ]
        loop_jac *= (math.pi)**n_loops
        # The phis angles then
        phis = [ 2.0*math.pi*xs[i+2] for i in range(0, len(xs), 3) ]
        loop_jac *= (2.0*math.pi)**n_loops
        
        ks, jacobians = transform_spherical_to_cartesian(np.array(
            list(zip( radii, thetas, phis ))
        ))
        loop_jac *= np.prod(jacobians)

        ks = [e for k in ks for e in k]
        return ks, loop_jac

    def loop_inv_parameterisation(self, ks):

        loop_inv_jac = 1.0
        n_loops = len(ks)//3
        if n_loops == 0:
            return [], loop_inv_jac

        cartesian_inputs = [ks[i:i+3] for i in range(0, len(ks),3)]
        spherical_coordinates, inv_jacobians = transform_cartesian_to_spherical(
            np.array(cartesian_inputs)
        )
        loop_inv_jac *= np.prod(inv_jacobians)

        # Now remap the spherical coordinates

        # First the radii
        xs_radii , radii_inv_jacobians = self.inverse_conformal_map(np.array([
            sph_coords[0] for sph_coords in spherical_coordinates]), self.E_cm)
        loop_inv_jac *= np.prod(radii_inv_jacobians)
        # The theta angles then
        xs_thetas = [ sph_coords[1]/math.pi for sph_coords in spherical_coordinates ]
        loop_inv_jac /= (math.pi)**n_loops
        # The phis angles then
        xs_phis = [ sph_coords[2]/(2.0*math.pi) for sph_coords in spherical_coordinates ]
        loop_inv_jac /= (2.0*math.pi)**n_loops

        #xs = sum((list(e) for e in list(zip(xs_radii,xs_thetas,xs_phis)),[])
        xs = [e for x in zip(xs_radii,xs_thetas,xs_phis) for e in x]
        return xs, loop_inv_jac

    def numerical_solve_soper_flow(self, final_state_configurations, Q0):

        if self.debug:
            logger.info("Numerical Soper Flow solving with Q0=%.3e and inputs:\n%s"%(Q0,pformat(final_state_configurations)))

        if not scipy_optimize_loaded:
            raise SamplerError("The numerical solution of Soper flow requires scipy. Install it.") 
           
        square_roots = [ {
                'k.k' : cfg['k'].square(),
                'k.p' : cfg['k'].dot(cfg['shift']),
                'm2'  : cfg['mass']**2
            } for cfg in  final_state_configurations ]

        def energy_sum(t):
            return sum(math.sqrt((t**2)*sqr['k.k'] + t*sqr['k.p'] + sqr['m2']) for sqr in square_roots)/Q0 - 1.

        def energy_sum_squared(t):
            return energy_sum(t)**2

        if energy_sum(0.) > 0.:
            raise SamplerError("The Soper flow could not be solved for these square roots because it will "+
                "have no solutions or several positive ones which is not supported at the moment:\n%s "%pformat(final_state_configurations))

        upper_bound = Q0/max(sqr['k.k'] for sqr in square_roots)
        n_try = 1
        while energy_sum(upper_bound) < 0. and n_try < 100:
            upper_bound *= 10.
            n_try += 1

        if energy_sum(upper_bound) < 0.:
            raise SamplerError("The Soper flow could not be solved for these square roots because it could not find a point "+
                "for which the Cutkosky surface evaluates to a positive quantity :\n%s "%pformat(final_state_configurations))

        if self.debug:
            logger.info("Detected bounds for Numerical Soper Flow solving: (%.16e, %.16e)"%(0., upper_bound))

        scipy_res = None
        try:
            scipy_res = optimize.minimize_scalar(
                energy_sum_squared,
                bounds=(0.,upper_bound),
                #tol=1.0e-10,
                method='bounded',
                options={
                    'xatol': 1.0e-16,
                    'maxiter' : 1000
                }
            )
            fsolve_success = True
        except Warning as e:
            logger.info("Numerical Soper flow failed with Warning:\n%s"%str(e))
            fsolve_success = False
        
        if scipy_res is None or not fsolve_success:
             raise SamplerError("Could not numerically solve Soper Flow for:\n%s"%pformat(final_state_configurations))

        t_solution = scipy_res.x
        if self.debug:
            logger.info("Numerical Soper Flow solving returned:\n%s"%str(scipy_res))
            logger.info("Target function evaluation for t=%.16e found: %.16e"%(t_solution,energy_sum(t_solution)))

        _numerical_soper_flow_threshold = 1.0e-7
        if t_solution<=0. or abs(energy_sum(t_solution))>_numerical_soper_flow_threshold:
             raise SamplerError("The numerical solver for Soper Flow found a solution which is however unsatisfactory t=%.16e, abs(E(t))=%.3e>%.3e."%(
                 t_solution, abs(energy_sum(t_solution)), _numerical_soper_flow_threshold
             ))

        return t_solution

    def solve_soper_flow(self, final_state_configurations, external_momenta):
        """ Numerically solve for the specified final_state_configurations which is a list with entries:
            'k', 'shift', 'direction' and 'mass'.
            When all masses and shifts are zero, an anlytic solution can be found, otherwise use optimize.root_scalar
            to numerically find the solution.
        """
        
        # Make sure we are in the rest frame of the collision
        if  ( len(external_momenta)==1 and math.sqrt(vectors.Vector(external_momenta[0][1:]).square())/external_momenta[0][0] > 1.0e-10 ) or \
            ( len(external_momenta)>1 and math.sqrt(sum(vectors.Vector(q[1:]) for q in external_momenta).square()/vectors.Vector(external_momenta[0][1:]).square()) > 1.0e-10 ):
            raise SamplerError("The Soper flow can currently only be solved for kinematic configuration in the collision rest frame.")

        Q0 = sum(q[0] for q in external_momenta)

        if all( (cfg['mass']==0. and cfg['shift'].square()==0.) for cfg in final_state_configurations):
            return Q0/sum( math.sqrt(cfg['k'].square()) for cfg in final_state_configurations )
        else:
            return self.numerical_solve_soper_flow(final_state_configurations, Q0)

    def __call__(self, continuous_inputs, discrete_inputs, selected_cut_and_side=None, selected_LMB=None, **opts):

        start_time = time.time()

        self.n_evals.value += 1

        if selected_cut_and_side is None:
            selected_cut_and_side = self.selected_cut_and_side 
        if selected_LMB is None:
            selected_LMB = self.selected_LMB

        if self.debug: logger.debug('\n\n%s {{{ New integrand evaluation with xs = [ %s ] }}}%s\n\n'%(
            utils.bcolors.GREEN, ' '.join('%.16f'%x for x in continuous_inputs) ,utils.bcolors.ENDC
        ))

        final_res = 0.
        if self.debug or self.return_individual_channels:
            all_channel_weights = {}

            # Set a dummy result with all entries set to zero when we encounter a numerical crash
            all_channel_weights_set_to_none = {}
            for i_channel, channel_info in enumerate(self.SG['SG_multichannel_info']):
                if selected_cut_and_side is not None and i_channel not in selected_cut_and_side:
                    continue
                if self.external_phase_space_generation_type == 'flat' and channel_info['side']!='left':
                    continue
                for i_LMB, LMB_info in enumerate(channel_info['loop_LMBs']):
                    if selected_LMB is not None and i_LMB not in selected_LMB:
                        continue
                    all_channel_weights_set_to_none[(i_channel,i_LMB)] = None

        if self.channel_for_generation is None:

            # Now we loop over channels 
            for i_channel, channel_info in enumerate(self.SG['SG_multichannel_info']):
                if selected_cut_and_side is not None and i_channel not in selected_cut_and_side:
                    continue
                
                # When generating the phase-space flat, we always ever only have a single channel per cut.
                if self.external_phase_space_generation_type == 'flat' and channel_info['side']!='left':
                    continue

                for i_LMB, LMB_info in enumerate(channel_info['loop_LMBs']):
                    if selected_LMB is not None and i_LMB not in selected_LMB:
                        continue
                    with warnings.catch_warnings():
                        warnings.filterwarnings('error')

                        try:
                            this_channel_wgt = self.evaluate_channel(continuous_inputs, discrete_inputs, i_channel, i_LMB, multi_channeling=True, **opts)
                        except RuntimeWarning as w:
                            if self.show_warnings: logger.warning("The following numerical issue occurred when evaluating channel (%d,%d) of the integrand: %s\n%s"%(
                                i_channel, i_LMB, str(w), str(traceback.format_exc()) ))
                            if self.show_warnings: logger.warning("Aborting point now.")
                            this_channel_wgt = None
                        if this_channel_wgt is None:
                            self.n_aborted_evals.value += 1
                            if self.show_warnings: logger.warning("Exceptional sample point encountered at xs = [ %s ]"%( ' '.join('%.16f'%x for x in continuous_inputs) ))
                            if self.show_warnings: logger.warning("%sNumber of exceptional sample points %d/%d (%.3f%%)%s"%( 
                                utils.bcolors.RED,
                                self.n_aborted_evals.value, self.n_aborted_evals.value+self.n_evals.value, (self.n_aborted_evals.value/float(self.n_aborted_evals.value+self.n_evals.value))*100.,
                                utils.bcolors.ENDC
                             ))
                            final_res = None
                            break

                    if self.debug or self.return_individual_channels:
                        all_channel_weights[(i_channel,i_LMB)] = this_channel_wgt
                    final_res += this_channel_wgt

                if final_res is None:
                    break
            if final_res is not None and self.debug: logger.debug('%sAll channel weights:\n%s%s'%(utils.bcolors.GREEN,pformat(all_channel_weights),utils.bcolors.ENDC))

        else:
            try:
                this_channel_wgt = self.evaluate_channel(
                    continuous_inputs, discrete_inputs, self.channel_for_generation[0], self.channel_for_generation[1], multi_channeling=False, **opts)
            except RuntimeWarning as w:
                if self.show_warnings: logger.warning("The following numerical issue occurred when evaluating channel (%d,%d) of the integrand: %s\n%s"%(
                    self.channel_for_generation[0], self.channel_for_generation[1], str(w), str(traceback.format_exc()) ))
                if self.show_warnings: logger.warning("Aborting point now.")
                this_channel_wgt = None
            if this_channel_wgt is None:
                self.n_aborted_evals.value += 1
                if self.show_warnings: logger.warning("Exceptional sample point encountered at xs = [ %s ]"%( ' '.join('%.16f'%x for x in continuous_inputs) ))
                if self.show_warnings: logger.warning("%sNumber of exceptional sample points %d/%d (%.3f%%)%s"%( 
                    utils.bcolors.RED,
                    self.n_aborted_evals.value, self.n_aborted_evals.value+self.n_evals.value, (self.n_aborted_evals.value/float(self.n_aborted_evals.value+self.n_evals.value))*100.,
                    utils.bcolors.ENDC
                    ))
                final_res = None

            if final_res is not None:
                if self.debug or self.return_individual_channels:
                    all_channel_weights[(self.channel_for_generation[0],self.channel_for_generation[1])] = this_channel_wgt
                final_res += this_channel_wgt

        if self.debug: logger.debug('Python integrand evaluation time = %.2f ms.'%((time.time()-start_time)*1000.0))

        self.update_evaluation_statistics(list(continuous_inputs),final_res)

        if self.return_individual_channels:
            if final_res is None:
                all_channel_weights = all_channel_weights_set_to_none
            res_dic = { 'I' : final_res if final_res is not None else 0. }
            for (i_channel,i_LMB), wgt in all_channel_weights.items():
                res_dic['channel #(%d,%d) [ CC #%d, side=%s, LMB #%d (%s) ]'%(
                    i_channel, i_LMB, 
                    self.SG['SG_multichannel_info'][i_channel]['cutkosky_cut_id'], 
                    self.SG['SG_multichannel_info'][i_channel]['side'],
                    i_LMB,
                    ','.join(self.SG['SG_multichannel_info'][i_channel]['loop_LMBs'][i_LMB]['loop_edges'])
                )] = (wgt if wgt is not None else 0.)
            return res_dic
        else:
            return (final_res if final_res is not None else 0.)

    def evaluate_channel(self, continuous_inputs, discrete_inputs, selected_cut_and_side, selected_LMB, multi_channeling=True, **opts):
        """ The 'selected_cut_and_side' and 'selected_LMB' give the option of specifying a particle (subset of) all integration channels."""

        random_variables = list(continuous_inputs)
        
        if self.debug: logger.debug('%s>>> Starting new integrand evaluation using CC cut and side #%d and LMB index #%d for sampling and with random variables: %s%s'%(
            utils.bcolors.GREEN,selected_cut_and_side, selected_LMB, str(random_variables),utils.bcolors.ENDC
        ))

        if self.debug: logger.debug('Edge masses:\n%s'%pformat(self.edge_masses))

        if self.frozen_momenta is None:
            x_t = random_variables[self.dimensions.get_position('x_t')]
        else:
            x_t = 1.
        if self.debug: logger.debug('x_t=%s'%x_t)

        # We must now promote the point to full 3^L using the causal flow, and thus determine the rescaling_t with the right density.
        rescaling_t, wgt_t = self.h_function(x_t)
        if self.debug: logger.debug('t, 1/t, wgt_t=%s, %s, %s'%(rescaling_t, 1./rescaling_t, wgt_t))

        normalising_func = self.h_function.exact_PDF(rescaling_t)
        if self.debug: logger.debug('normalising_func=%s'%str(normalising_func))
        if self.debug: logger.debug('normalising_func * wgt-t =%s'%str(normalising_func*wgt_t))

        # Note that normalising_func*wgt_t = 1 when the sampling PDF matches exactly the h function, but it won't if approximated.

        i_channel = selected_cut_and_side

        channel_info = self.SG['SG_multichannel_info'][i_channel]

        generator = self.multichannel_generators[i_channel]
        if self.debug and self.frozen_momenta is None and self.external_phase_space_generation_type == 'advanced':
            logger.debug("Considering the following topology:")
            logger.debug("-"*10)
            logger.debug(generator.get_topology_string(generator.topology, path_to_print=generator.path))
            logger.debug("-"*10)

        topology_info = self.SG['cutkosky_cuts'][
                channel_info['cutkosky_cut_id']
            ]['multichannel_info']['tree_topologies'][
                channel_info['side']                    
            ]

        # Get the external phase-space random variables
        if self.frozen_momenta is None:
            external_phase_space_xs = [ 
                random_variables[self.dimensions.get_position('x_%d'%i_dim)] 
                for i_dim in range(1,len(channel_info['independent_final_states'])*3) 
            ]
        else:
            external_phase_space_xs = generator.PS_xs

        # build the external kinematics 
        try:
            PS_point, PS_jac, _, _ = generator.get_PS_point(external_phase_space_xs)
        except Exception as e:
            if self.show_warnings: logger.warning("Exception during phase-space generation: %s\n%s"%(str(e),str(traceback.format_exc())))
            PS_point = None

        if PS_point is None:
            if self.show_warnings: logger.warning(("\nPhase-space generator failed to generate a kinematic configuration for the following input xs:\n%s"+
                            "\ncorresponding to the following input random variables of the complete integrand:\n%s\n"+
                            "Immediately returning zero weight for this point.")%(
                                external_phase_space_xs, random_variables
                            ))
            return None

        if self.debug: logger.debug('PS xs=%s'%str(external_phase_space_xs))
        if self.debug:
            incoming_edge_names = [e for (e,d) in topology_info['incoming_edges']]
            if len(incoming_edge_names)==1:
                incoming_edge_names *= 2
            outgoing_edge_names = [e for (e,d) in topology_info['outgoing_edges']]
            logger.debug('PS_point=\n%s'%LorentzVectorList(PS_point).__str__(n_initial=2, leg_names=incoming_edge_names+outgoing_edge_names))

        inv_aL_jacobian = 1.
        if self.frozen_momenta is None:
            # Correct for the 1/(2E) of each cut propagator that alphaLoop includes but which is already included in the normal PS parameterisation
            for v in PS_point[2:]:
                inv_aL_jacobian *= 2*v[0]
        if self.debug: logger.debug('(-2*pi*I)/Es=%s'%str( (complex(0., -2.*math.pi)**len(PS_point[2:]) / inv_aL_jacobian)) )

        # We should undo the mock-up 2>N if n_initial is 1
        PS_initial_state_four_momenta = PS_point[:2] 
        PS_point = { edge_name : PS_point[2+i_edge].space()*edge_direction for i_edge, (edge_name, edge_direction) in enumerate(topology_info['outgoing_edges']) } 

        i_LMB = selected_LMB
        LMB_info = channel_info['loop_LMBs'][i_LMB]

        CMB_edges = list(LMB_info['loop_edges'])
        CMB_edges.extend([edge_name for edge_name, edge_direction in channel_info['independent_final_states']])

        if self.debug: logger.debug('independent final states=%s'%str(channel_info['independent_final_states']))
        if self.debug: logger.debug('dependent final states=%s'%str( [e for (e,d) in topology_info['outgoing_edges'] if (e,d) not in channel_info['independent_final_states']][0] ))
        if self.debug: logger.debug('loop edges=%s'%str(LMB_info['loop_edges']))
        if self.debug: logger.debug('all CMB_edges=%s'%str(CMB_edges))

        if self.frozen_momenta is None:
            loop_phase_space_xs = [ 
                random_variables[self.dimensions.get_position('x_%d'%i_dim)] 
                for i_dim in range(
                    len(channel_info['independent_final_states'])*3,
                    self.SG['topo']['n_loops']*3
                ) 
            ]
        else:
            loop_phase_space_xs = [ 
                random_variables[self.dimensions.get_position('x_%d'%i_dim)] 
                for i_dim in range(1,(self.SG['topo']['n_loops']-len(self.frozen_momenta['out']))*3+1) 
            ]

        ks, loop_jac = self.loop_parameterisation(loop_phase_space_xs)
        
        PS_point.update( { edge_name : vectors.Vector(ks[3*i_LMB_entry:3*(i_LMB_entry+1)]) for i_LMB_entry, edge_name in enumerate(LMB_info['loop_edges']) } )

        if self.debug: logger.debug('loop_phase_space_momenta=%s'%str(ks))
        if self.debug: logger.debug('loop_xs=%s'%str(loop_phase_space_xs))

        downscaled_input_momenta_in_CMB = [ PS_point[edge_name] for edge_name in CMB_edges ]
        if self.debug: logger.debug('downscaled_input_momenta_in_CMB=%s'%str(downscaled_input_momenta_in_CMB))

        transformation_matrix, parametric_shifts = LMB_info['transformation_to_defining_LMB']
        shifts = [ sum([self.external_momenta[i_shift]*shift 
                    for i_shift, shift in enumerate(parametric_shift)]) for parametric_shift in parametric_shifts ]
        downscaled_input_momenta_in_LMB = [ vectors.Vector(v) for v in transformation_matrix.dot(
            [list(rm+shift) for rm, shift in zip(downscaled_input_momenta_in_CMB, shifts)] ) ]

        if self.debug:
            all_downscaled_momenta = { edge_name : 
                    ( sum([ downscaled_input_momenta_in_LMB[i_k]*wgt for i_k, wgt in enumerate(edge_info[0]) ]) + 
                      sum([ self.external_momenta[i_shift]*wgt for i_shift, wgt in enumerate(edge_info[1]) ]) )
                for edge_name, edge_info in self.SG['edge_signatures'].items() }
            logger.debug('Downscaled momenta of all supergraph edges:\n%s'%pformat(all_downscaled_momenta))

        upscaled_input_momenta_in_LMB = [ k*rescaling_t for k in downscaled_input_momenta_in_LMB ]
        if self.debug: logger.debug('upscaled_input_momenta_in_LMB=%s'%str(upscaled_input_momenta_in_LMB))

        # Applying the rescaling to embed the momenta in the full R^(3L) space
        # and thus effectively compute the inverse of the jacobian that will be included by the full integrand in rust

        # First the delta function jacobian
        if self.frozen_momenta is None:
            delta_jacobian = 0.
            inv_rescaling_t = 1./rescaling_t
            for edge_name, edge_direction in topology_info['outgoing_edges']:
                
                k = sum([ upscaled_input_momenta_in_LMB[i_k]*wgt for i_k, wgt in enumerate(self.SG['edge_signatures'][edge_name][0]) ])
                shift = sum([ self.external_momenta[i_shift]*wgt for i_shift, wgt in enumerate(self.SG['edge_signatures'][edge_name][1]) ])

                delta_jacobian += (
                    ( 2. * inv_rescaling_t * k.square() + k.dot(shift) ) / 
                    ( 2. * math.sqrt( (inv_rescaling_t*k + shift).square() + self.edge_masses[edge_name]**2 ) )
                )

            delta_jacobian *= (inv_rescaling_t)**2
        else:
            delta_jacobian = 1.

        if self.debug: logger.debug('delta_jacobian=%s'%str(delta_jacobian))
        if self.debug: logger.debug('1./delta_jacobian=%s'%str(1./delta_jacobian))
        inv_aL_jacobian *= delta_jacobian

        if self.frozen_momenta is None:
            # Then the inverse of the H-function and of the Jacobian of the causal flow change of variables
            # WARNING: since h(t) = h (1/t) I am not 100% sure if 1/h(t) or 1/h(1/t) is logically he right thing to do below, but it would be the same anyway.
            try:
                inv_aL_jacobian *=  ( 1. / self.h_function.exact_PDF(1./rescaling_t) ) * (rescaling_t**(len(CMB_edges)*3)) 
            except ZeroDivisionError:
                if self.show_warnings: logger.warning("H-function evaluated to 0 for t=%.16f. Aborting point now."%(1./rescaling_t))
                return None

        # The final jacobian must then be our param. jac together with that of t divided by the one from alphaloop.
        final_jacobian = PS_jac * loop_jac * wgt_t * inv_aL_jacobian * normalising_func
        if math.isinf(final_jacobian):
            if self.show_warnings: logger.warning("final_jacobian evaluated to infinity. Aborting point now."+
                "\n(PS_jac=%.16e * loop_jac=%.16e * wgt_t=%.16e * inv_aL_jacobian=%.16e * normalising_func=%.16e)"%(
                    PS_jac , loop_jac , wgt_t , inv_aL_jacobian , normalising_func
                ))
            return None
        if self.debug: logger.debug('PS_jac=%s'%PS_jac)
        if self.debug: logger.debug('loop_jac=%s'%loop_jac)
        if self.debug: logger.debug('wgt_t=%s'%wgt_t)
        if self.debug: logger.debug('inv_aL_jacobian=%s'%inv_aL_jacobian)
        if self.debug: logger.debug('normalising_func=%s'%normalising_func)
        if self.debug: logger.debug('final jacobian=%s'%final_jacobian)


        # We must undo aL parameterisation as we must provide the inputs to the aL rust integrand in x-space.
        aL_xs = []
        undo_aL_parameterisation = 1.
        for i_v, v in enumerate(upscaled_input_momenta_in_LMB):
            if self.frozen_momenta is not None and i_v >= (len(upscaled_input_momenta_in_LMB)-len(self.frozen_momenta['out'])):
                continue
            kx, ky, kz, inv_jac = self.rust_worker.inv_parameterize(list(v), i_v, self.E_cm**2)
            undo_aL_parameterisation *= inv_jac
            aL_xs.extend([kx, ky, kz])
        # When there are frozen momenta we must pad the input xs of alphaloop with the frozen xs
        if self.frozen_momenta is not None:
            aL_xs.extend(external_phase_space_xs)

        if self.debug: logger.debug('aL xs=%s'%str(aL_xs))
        if self.debug: logger.debug('undo_aL_parameterisation=%s'%str(undo_aL_parameterisation))
        # Finally actually call alphaLoop full integrand
        rust_start_time = time.time()
        re, im = self.rust_worker.evaluate_integrand( aL_xs )
        if self.debug: logger.debug('Rust integrand evaluation time = %.2f ms.'%((time.time()-rust_start_time)*1000.0))

        aL_wgt = complex(re, im)
        if self.debug: logger.debug('aL res=%s'%str(aL_wgt))
        aL_wgt *= undo_aL_parameterisation
        if self.debug: logger.debug('aL*undo_aL_parameterisation res=%s'%str(aL_wgt))

        if self.debug and self.frozen_momenta is None:
            reconstituted_res = complex(0., 0.)
            for cut_ID, cut_info in enumerate(self.SG['cutkosky_cuts']):
                logger.debug("    Result for cut_ID #%d with external edges %s"%(cut_ID, ', '.join(c['name'] for c in cut_info['cuts']) ))
                scaling_solutions = list(self.rust_worker.get_scaling(upscaled_input_momenta_in_LMB, cut_ID))
                scaling, scaling_jacobian = scaling_solutions.pop(0)
                while scaling < 0.0:
                    if len(scaling_solutions)==0:
                        break
                    scaling, scaling_jacobian = scaling_solutions.pop(0)
                logger.debug("    > t-scaling, t-scaling jacobian = %.16f, %.16f"%(scaling, scaling_jacobian))
                re, im = self.rust_worker.evaluate_cut_f128(upscaled_input_momenta_in_LMB,cut_ID,scaling,scaling_jacobian)
                this_cut_res = complex(re, im)
                if self.hyperparameters['CrossSection']['picobarns'] or False:
                    this_cut_res /= 0.389379304e9
                logger.debug("    > aL res for this cut = %s"%str(this_cut_res))
                reconstituted_res += this_cut_res
            logger.debug("    Reconstituted aL res = %s (rel. diff = %s )"%(
                str(reconstituted_res), complex(
                    reconstituted_res.real/aL_wgt.real-1. if abs(aL_wgt.real) > 0. else reconstituted_res.real,
                    reconstituted_res.imag/aL_wgt.imag-1. if abs(aL_wgt.imag) > 0. else reconstituted_res.imag,
                )
            ))

        if self.debug: logger.debug('normalising_func=%s'%normalising_func)

        # And accumulate it to what will be the final result
        this_channel_jacobian = final_jacobian
        #this_channel_jacobian *= normalising_func
        final_res = this_channel_jacobian * aL_wgt


        #HACK
        #final_jacobian /= normalising_func * wgt_t * PS_jac * loop_jac

        if multi_channeling:

            # And accumulate this jacobian into the denominator of the multichanneling factor
            multichannel_factor_denominator = 1./final_jacobian

            multichanneling_cache = {'cutkosky_cut_and_side_cache': {} }
            # Uncomment the line below to disable this dynamic cache
            # multichanneling_cache = None

            if self.debug:
                all_multi_channel_denominator_weights = {
                    (i_channel, i_LMB) : final_jacobian
                } 
            for MC_i_channel, MC_channel_info in enumerate(self.SG['SG_multichannel_info']):
                # of the LMB chosen already. I leave this to future refinements.
                # When generating the phase-space flat, we always ever only have a single channel per cut.
                if self.external_phase_space_generation_type == 'flat' and MC_channel_info['side']!='left':
                    continue
                for MC_i_LMB, MC_LMB_info in enumerate(MC_channel_info['loop_LMBs']):
                    # Note that for debugging it is useful to uncomment the three lines below and verify explicitly that 
                    #         MC_final_jacobian = final_jacobian
                    # for (MC_i_channel, MC_i_LMB) == (i_channel, i_LMB)
                    # because this is indeed a non-trivial check of all the inversion of the parameterisations used here.
                    if (self.debug<=0) and (MC_i_channel, MC_i_LMB) == (i_channel, i_LMB):
                        # This contributions was already accounted for as part of final_jacobian
                        continue
                    
                    # Simply set the option multichanneling_cache to None below if you want to disable it.
                    #multichanneling_cache = None
                    MC_final_jacobian = self.get_jacobian_for_channel(
                        MC_i_channel, MC_i_LMB, upscaled_input_momenta_in_LMB, PS_initial_state_four_momenta, multichanneling_cache = multichanneling_cache)

                    # Check for a numerical issue in the computation of the MC_final_jacobian
                    if MC_final_jacobian is None:
                        return None

                    if (MC_i_channel, MC_i_LMB) == (i_channel, i_LMB):
                        logger.debug("%sFinal jacobian vs MC Final jacobian test for channel (%d, %d): %.16e, %.16e (rel. dif: %.16e)%s"%(
                            utils.bcolors.GREEN, i_channel, i_LMB, final_jacobian, MC_final_jacobian, (final_jacobian-MC_final_jacobian)/final_jacobian, utils.bcolors.ENDC
                        ))
                        continue

                    if self.debug:
                        all_multi_channel_denominator_weights[(MC_i_channel, MC_i_LMB)] = 1./MC_final_jacobian
                    # Aggregate that final jacobian for this channel to the multichanneling denominator
                    multichannel_factor_denominator += 1./MC_final_jacobian
        
            # Include the multichanneling factor
            if self.debug: logger.debug("%sAll multi-channeling denominator components:\n%s%s"%(
                utils.bcolors.BLUE,pformat(all_multi_channel_denominator_weights),utils.bcolors.ENDC
            ))
            if self.debug: logger.debug('normalising_func =%s'%( normalising_func ))
            if self.debug: logger.debug('Jacobian of this channel=%s'%( this_channel_jacobian ))
            if self.debug: logger.debug('alphaLoop integrand weight=%s'%( aL_wgt ))
            if self.debug: logger.debug('Multi-channeling factor numerator=%s'%( 1./final_jacobian ))
            if self.debug: logger.debug('Multi-channeling factor denominator=%s'%( multichannel_factor_denominator ))
            if self.debug: logger.debug('Integrand without jacobian / Multi-channeling denominator=%s'%( aL_wgt / multichannel_factor_denominator ))
            # The quantity below is interesting because for some reason it is very often of order O(1) for the max weights found, hinting at how to possibly cure these large weight perhaps.
            if self.debug: logger.debug('(Integrand without jacobian / Multi-channeling denominator)*normalising_func=%s'%( (aL_wgt / multichannel_factor_denominator) * normalising_func ))
            if self.debug: logger.debug('Jacobian of this channel * Multi-channeling numerator=%s'%( this_channel_jacobian * (1./final_jacobian) ))

            final_res *= 1./final_jacobian / multichannel_factor_denominator

        if self.debug: logger.debug( 'Final weight returned = %s (phase = %s)'%( final_res, self.phase) )

        if self.phase=='real':
            return final_res.real
        else:
            return final_res.imag

    def get_jacobian_for_channel(self, MC_i_channel, MC_i_LMB, upscaled_input_momenta_in_LMB, PS_initial_state_four_momenta, multichanneling_cache=None):

        if self.debug: logger.debug('%s> Computing jacobian for CC cut and side #%d and LMB index #%d from the follwing upscaled LMB momenta:\n%s%s'%(
                    utils.bcolors.BLUE, MC_i_channel, MC_i_LMB, str(upscaled_input_momenta_in_LMB), utils.bcolors.ENDC
                ))

        if multichanneling_cache is None:
            multichanneling_cache = {'cutkosky_cut_and_side_cache': {} }

        MC_channel_info = self.SG['SG_multichannel_info'][MC_i_channel]

        MC_LMB_info = MC_channel_info['loop_LMBs'][MC_i_LMB]

        MC_topology_info = self.SG['cutkosky_cuts'][
            MC_channel_info['cutkosky_cut_id']
        ]['multichannel_info']['tree_topologies'][
            MC_channel_info['side']                    
        ]

        MC_CMB_edges = list(MC_LMB_info['loop_edges'])
        MC_CMB_edges.extend([edge_name for edge_name, edge_direction in MC_channel_info['independent_final_states']])

        MC_generator = self.multichannel_generators[MC_i_channel]
        if self.debug and self.frozen_momenta is None and self.external_phase_space_generation_type == 'advanced':
            logger.debug("MC Considering the following topology:")
            logger.debug("-"*10)
            logger.debug(MC_generator.get_topology_string(MC_generator.topology, path_to_print=MC_generator.path))
            logger.debug("-"*10)

        if self.debug: logger.debug('MC independent final states=%s'%str(MC_channel_info['independent_final_states']))
        if self.debug: logger.debug('MC dependent final states=%s'%str( [e for (e,d) in MC_topology_info['outgoing_edges'] if (e,d) not in MC_channel_info['independent_final_states']][0] ))
        if self.debug: logger.debug('MC loop edges=%s'%str(MC_LMB_info['loop_edges']))
        if self.debug: logger.debug('MC all CMB_edges=%s'%str(MC_CMB_edges))

        MC_upscaled_input_momenta_in_CMB = {}
#        MC_aL_xs = []
#        MC_inv_aL_jacobian = 1.
        for i_edge, edge_name in enumerate(MC_CMB_edges):
            k = sum([ upscaled_input_momenta_in_LMB[i_k]*wgt for i_k, wgt in enumerate(self.SG['edge_signatures'][edge_name][0]) ])
            shift = sum([ self.external_momenta[i_shift]*wgt for i_shift, wgt in enumerate(self.SG['edge_signatures'][edge_name][1]) ])
            MC_upscaled_input_momenta_in_CMB[edge_name] = k+shift

#            kx, ky, kz, inv_jac = self.rust_worker.inv_parameterize(
#                list(MC_upscaled_input_momenta_in_CMB[edge_name]), i_edge, self.E_cm**2
#            )
#            MC_aL_xs.extend([kx, ky, kz])
#            MC_inv_aL_jacobian *= inv_jac

        if self.debug: logger.debug('MC_upscaled_input_momenta_in_CMB=%s'%str(MC_upscaled_input_momenta_in_CMB))

        # The quantities below depend only on the selected CC cut and side but *not* on the particular LMB chosen, so we cache them.
        if MC_i_channel in multichanneling_cache['cutkosky_cut_and_side_cache']:
            MC_inv_rescaling_t = multichanneling_cache['cutkosky_cut_and_side_cache'][MC_i_channel]['MC_inv_rescaling_t']
            MC_PS_jac = multichanneling_cache['cutkosky_cut_and_side_cache'][MC_i_channel]['MC_PS_jac']
            MC_wgt_t = multichanneling_cache['cutkosky_cut_and_side_cache'][MC_i_channel]['MC_wgt_t']
            MC_inv_aL_jacobian = multichanneling_cache['cutkosky_cut_and_side_cache'][MC_i_channel]['MC_inv_aL_jacobian']
            MC_normalising_func = multichanneling_cache['cutkosky_cut_and_side_cache'][MC_i_channel]['MC_normalising_func']
        else:
            MC_upscaled_final_state_configurations = [ ]
            MC_upscaled_final_state_momenta = [ ]
            for edge_name, edge_direction in MC_topology_info['outgoing_edges']:
                k = sum([ upscaled_input_momenta_in_LMB[i_k]*wgt for i_k, wgt in enumerate(self.SG['edge_signatures'][edge_name][0]) ])
                shift = sum([ self.external_momenta[i_shift]*wgt for i_shift, wgt in enumerate(self.SG['edge_signatures'][edge_name][1]) ])
                MC_upscaled_final_state_momenta.append( (k+shift)*edge_direction )
                MC_upscaled_final_state_configurations.append({
                    'k' : k,
                    'shift' : shift,
                    'direction' : edge_direction,
                    'mass' : self.edge_masses[edge_name]
                })
            if self.frozen_momenta is None:
                MC_inv_rescaling_t = self.solve_soper_flow(
                    MC_upscaled_final_state_configurations, self.hyperparameters['CrossSection']['incoming_momenta'] )
            else:
                MC_inv_rescaling_t = 1.

            MC_rescaling_t = 1./MC_inv_rescaling_t

            if self.debug: logger.debug('MC_rescaling_t=%s'%str(MC_rescaling_t))

            MC_x_t, MC_inv_wgt_t = self.h_function.inverse_sampling(MC_rescaling_t)
            MC_wgt_t = 1./MC_inv_wgt_t

            if self.debug: logger.debug('MC_downscaled_input_momenta_in_CMB=%s'%str({e: v*MC_inv_rescaling_t for e,v in MC_upscaled_input_momenta_in_CMB.items()}))

            # Note that MC_normalising_func*MC_wgt_t = 1 when the sampling PDF matches exactly the h function, but it won't if approximated.
            MC_normalising_func = self.h_function.exact_PDF(MC_rescaling_t)
            if self.debug: logger.debug('MC_normalising_func=%s'%str(MC_normalising_func))

            MC_downscaled_final_state_momenta = []
            MC_downscaled_final_state_four_momenta = []
            for cfg in MC_upscaled_final_state_configurations:
                MC_downscaled_final_state_momenta.append(
                    (cfg['k']*MC_inv_rescaling_t + cfg['shift'])*cfg['direction']
                )
                MC_downscaled_final_state_four_momenta.append(vectors.LorentzVector(
                    [ math.sqrt(MC_downscaled_final_state_momenta[-1].square()+cfg['mass']**2), ]+list(MC_downscaled_final_state_momenta[-1])
                ))

            # Warning! It is irrelevant for our usecase here, but the inversion below is only valid for the jacobian and the xs that control invariant masses 
            # when using the SingleChannelPhaseSpaceGenerator (for the FlatInvertiblePhaseSpaceGenerator all inverted xs would be correct)
            # We need to pad with two dummy initial states
            if self.debug:
                MC_incoming_edge_names = [e for (e,d) in MC_topology_info['incoming_edges']]
                if len(MC_incoming_edge_names)==1:
                    MC_incoming_edge_names *= 2
                MC_outgoing_edge_names = [e for (e,d) in MC_topology_info['outgoing_edges']]
                logger.debug('MC PS_point=\n%s'%LorentzVectorList(
                    list(PS_initial_state_four_momenta)+MC_downscaled_final_state_four_momenta
                ).__str__(n_initial=2, leg_names=MC_incoming_edge_names+MC_outgoing_edge_names))

            try:
                MC_xs, MC_PS_jac= MC_generator.invertKinematics( self.E_cm, vectors.LorentzVectorList(list(PS_initial_state_four_momenta)+MC_downscaled_final_state_four_momenta) )
            except Exception as e:
                if self.show_warnings: logger.warning("Exception during phase-space generation inversion: %s\n%s"%(str(e),str(traceback.format_exc())))
                MC_xs = None
            if MC_xs is None:
                return None
    
            if self.debug: logger.debug('MC PS xs=%s'%str(MC_xs))

            # Correct for the 1/(2E) of each cut propagator that alphaLoop includes but which is already included in the normal PS parameterisation
            MC_inv_aL_jacobian = 1.
            if self.frozen_momenta is None:
                for v in MC_downscaled_final_state_four_momenta:
                    MC_inv_aL_jacobian *= 2*v[0]

            if self.frozen_momenta is None:
                MC_delta_jacobian = 0.
                for edge_name, edge_direction in MC_topology_info['outgoing_edges']:

                    k = sum([ upscaled_input_momenta_in_LMB[i_k]*wgt for i_k, wgt in enumerate(self.SG['edge_signatures'][edge_name][0]) ])
                    shift = sum([ self.external_momenta[i_shift]*wgt for i_shift, wgt in enumerate(self.SG['edge_signatures'][edge_name][1]) ])

                    MC_delta_jacobian += (
                        ( 2. * MC_inv_rescaling_t * k.square() + k.dot(shift) ) / 
                        ( 2. * math.sqrt( (MC_inv_rescaling_t*k + shift).square() + self.edge_masses[edge_name]**2 ) )
                    )
                MC_delta_jacobian *= (MC_inv_rescaling_t)**2
            else:
                MC_delta_jacobian = 1.
            if self.debug: logger.debug('MC_delta_jacobian=%s'%str(MC_delta_jacobian))
            MC_inv_aL_jacobian *= MC_delta_jacobian

            # Then the inverse of the H-function and of the Jacobian of the causal flow change of variables
            HFunction_wgt = self.h_function.exact_PDF(1./MC_rescaling_t)
            if HFunction_wgt == 0.:
                if self.show_warnings: logger.warning("H-function evaluated to 0 for t=%.16f. Aborting point now."%(1./MC_rescaling_t))
                return None

            # WARNING: since h(t) = h (1/t) I am not 100% sure if 1/h(t) or 1/h(1/t) is logically he right thing to do below, but it would be the same anyway.
            MC_inv_aL_jacobian *=  ( 1. / self.h_function.exact_PDF(1./MC_rescaling_t) ) * (MC_rescaling_t**(len(MC_CMB_edges)*3)) 

            # Store all these quantities being always identical independently of the particular LMB chosen for future use.
            multichanneling_cache['cutkosky_cut_and_side_cache'][MC_i_channel] = {}
            multichanneling_cache['cutkosky_cut_and_side_cache'][MC_i_channel]['MC_inv_rescaling_t'] = MC_inv_rescaling_t
            multichanneling_cache['cutkosky_cut_and_side_cache'][MC_i_channel]['MC_PS_jac'] = MC_PS_jac
            multichanneling_cache['cutkosky_cut_and_side_cache'][MC_i_channel]['MC_wgt_t'] = MC_wgt_t
            multichanneling_cache['cutkosky_cut_and_side_cache'][MC_i_channel]['MC_inv_aL_jacobian'] = MC_inv_aL_jacobian
            multichanneling_cache['cutkosky_cut_and_side_cache'][MC_i_channel]['MC_normalising_func'] = MC_normalising_func

        loop_phase_space_momenta = []
        for edge_name in MC_LMB_info['loop_edges']:
            loop_phase_space_momenta.extend(list(MC_upscaled_input_momenta_in_CMB[edge_name]*MC_inv_rescaling_t))

        MC_loop_xs, MC_loop_inv_jac = self.loop_inv_parameterisation(loop_phase_space_momenta)
        MC_loop_jac = 1./MC_loop_inv_jac
        if self.debug: logger.debug('MC_loop_phase_space_momenta=%s'%str(loop_phase_space_momenta))
        if self.debug: logger.debug('MC_loop_xs=%s'%str(MC_loop_xs))

        # Return the final jacobian of the parameterisation
        try:
            MC_final_jacobian = MC_PS_jac * MC_loop_jac * MC_wgt_t * MC_inv_aL_jacobian * MC_normalising_func
        except RuntimeWarning as w:
            if self.show_warnings: logger.warning("The following numerical issue occurred when computing a piece of the Multi-channeling denominator:"+
                " %s\nMultiplicative pieces are: %s\n                           %s"%(
                str(w),
                ' * '.join('%-25s'%factor for factor in ["MC_PS_jac", "MC_loop_jac", "MC_wgt_t", "MC_inv_aL_jacobian","MC_normalising_func"]),
                ' * '.join('%-25.16e'%factor for factor in [MC_PS_jac, MC_loop_jac, MC_wgt_t, MC_inv_aL_jacobian,MC_normalising_func])))
            if self.show_warnings: logger.warning("Aborting point now.")
            return None

        if self.debug: logger.debug('MC_PS_jac=%s'%MC_PS_jac)
        if self.debug: logger.debug('MC_loop_jac=%s'%MC_loop_jac)
        if self.debug: logger.debug('MC_wgt_t=%s'%MC_wgt_t)
        if self.debug: logger.debug('MC_inv_aL_jacobian=%s'%MC_inv_aL_jacobian)
        if self.debug: logger.debug('MC_normalising_func=%s'%MC_normalising_func)
        if self.debug: logger.debug('MC_final_jacobian=%s'%MC_final_jacobian)

        #HACK
        #MC_final_jacobian /= MC_normalising_func * MC_wgt_t * MC_PS_jac * MC_loop_jac

        return MC_final_jacobian

class CustomGenerator(object):
    pass

class generator_aL(CustomGenerator):

    def __init__(self, dimensions, rust_worker, SG_info, model, h_function, hyperparameters, debug=0, frozen_momenta=None, **opts):

        self.rust_worker = rust_worker
        self.model = model
        self.SG_info = SG_info
        self.hyperparameters = hyperparameters
        self.dimensions = dimensions
        self.debug = debug
        self.frozen_momenta = frozen_momenta

        E_cm = self.SG_info.get_E_cm(self.hyperparameters)        
        self.E_cm = E_cm
        if debug: logger.debug("E_cm detected=%.16e"%self.E_cm)

    def __call__(self, random_variables, **opts):
        """ 
        Generates a point using a sampling reflecting the topolgy of SG_QG0.
        It will return the point directly in x-space as well as the final weight to combine
        the result with.
        """
        xs = list(random_variables)

        # Pad xs with frozen momenta if necessary
        overall_inv_jac = 1.
        if self.frozen_momenta is not None:
            n_loop_vs = len(random_variables)//3
            for i_v, v in enumerate([v[1:] for v in self.frozen_momenta['out']]):
                xx, xy, xz, inv_jac = self.rust_worker.inv_parameterize(list(v), n_loop_vs+i_v, self.E_cm**2)
                xs.extend([xx, xy, xz])
                overall_inv_jac *= inv_jac*(2.0*math.pi)**4

        return xs, overall_inv_jac