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
import random
from scipy import optimize

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


        self.update_evaluation_statistics(xs,final_res)

        return final_res


class AdvancedIntegrand(integrands.VirtualIntegrand):

    def __init__(self, rust_worker, SG, model, h_function, hyperparameters, channel_for_generation = None, external_phase_space_generation_type="flat", debug=0, phase='real', 
        selected_cut_and_side=None, selected_LMB=None,
        **opts):
        """ Set channel_for_generation to (selected_cut_and_side_index, selected_LMB_index) to disable multichaneling and choose one particular channel for the parameterisation. """

        self.debug = debug
        self.rust_worker = rust_worker
        self.model = model
        self.SG = SG
        self.phase = phase
        self.hyperparameters = hyperparameters
        self.h_function = h_function
        self.conformal_map = lin_conformal
        self.inverse_conformal_map = lin_inverse_conformal
        self.channel_for_generation = channel_for_generation
        
        self.selected_cut_and_side=selected_cut_and_side
        self.selected_LMB=selected_LMB

        if external_phase_space_generation_type not in ['flat','advanced']:
            raise InvalidCmd("External phase-space generation mode not supported: %s"%external_phase_space_generation_type)
        self.external_phase_space_generation_type = external_phase_space_generation_type

        # Built external momenta. Remember that they appear twice.
        self.n_initial = len(hyperparameters['CrossSection']['incoming_momenta'])
        self.external_momenta = [ vectors.Vector(v[1:]) for v in self.hyperparameters['CrossSection']['incoming_momenta'] ]
        self.external_momenta.extend(self.external_momenta)

        self.dimensions = integrands.DimensionList()
        # The first dimension will always be for the upscaling 't' variable
        self.dimensions.append(integrands.ContinuousDimension('x_t',lower_bound=0.0, upper_bound=1.0))
        # Then the remaining ones will first be for the final-state mappings (the number depends on the particular Cutkosky cut)
        # and of the remaining loop variables
        self.dimensions.extend([ 
            integrands.ContinuousDimension('x_%d'%i_dim,lower_bound=0.0, upper_bound=1.0) 
            for i_dim in range(1,SG['topo']['n_loops']*3) 
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
        self.edge_masses = { edge_name : self.model['parameter_dict'][
                    self.model.get_particle(self.edges_PDG[edge_name]).get('mass')
            ].real for edge_name in self.edges_PDG }

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
            
            if self.external_phase_space_generation_type == 'flat':
                
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
            raise SamplerError("The numerical solution of Soper flow with optimize.root_scalar is not implemented yet.")

    def __call__(self, continuous_inputs, discrete_inputs, selected_cut_and_side=None, selected_LMB=None, **opts):

        start_time = time.time()

        self.n_evals.value += 1

        if selected_cut_and_side is None:
            selected_cut_and_side = self.selected_cut_and_side 
        if selected_LMB is None:
            selected_LMB = self.selected_LMB

        final_res = 0.
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

                    final_res += self.evaluate_channel(continuous_inputs, discrete_inputs, i_channel, i_LMB, multi_channeling=True, **opts)
        else:
            final_res += self.evaluate_channel(
                continuous_inputs, discrete_inputs, self.channel_for_generation[0], self.channel_for_generation[1], multi_channeling=False, **opts)

        if self.debug: logger.debug('Python integrand evaluation time = %.2f ms.'%((time.time()-start_time)*1000.0))

        self.update_evaluation_statistics(list(continuous_inputs),final_res)

        return final_res

    def evaluate_channel(self, continuous_inputs, discrete_inputs, selected_cut_and_side, selected_LMB, multi_channeling=True, **opts):
        """ The 'selected_cut_and_side' and 'selected_LMB' give the option of specifying a particle (subset of) all integration channels."""

        random_variables = list(continuous_inputs)
        
        if self.debug: logger.debug('>>> Starting new integrand evaluation using CC cut and side #%d and LMB index #%d for sampling and with random variables: %s'%(
            selected_cut_and_side, selected_LMB, str(random_variables)
        ))

        if self.debug: logger.debug('Edge masses:\n%s'%pformat(self.edge_masses))

        x_t = random_variables[self.dimensions.get_position('x_t')]
        if self.debug: logger.debug('x_t=%s'%x_t)

        # We must now promote the point to full 3^L using the causal flow, and thus determine the rescaling_t with the right density.
        rescaling_t, wgt_t = self.h_function(x_t)
        if self.debug: logger.debug('t, 1/t, wgt_t=%s, %s, %s'%(rescaling_t, 1./rescaling_t, wgt_t))

        normalising_func = self.h_function.PDF(rescaling_t)
        if self.debug: logger.debug('normalising_func=%s'%str(normalising_func))

        # Note that normalising_func*wgt_t = 1 when the sampling PDF matches exactly the h function, but it won't if approximated.

        i_channel = selected_cut_and_side

        channel_info = self.SG['SG_multichannel_info'][i_channel]

        generator = self.multichannel_generators[i_channel]
        if self.debug and self.external_phase_space_generation_type == 'advanced':
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
        external_phase_space_xs = [ 
            random_variables[self.dimensions.get_position('x_%d'%i_dim)] 
            for i_dim in range(1,len(channel_info['independent_final_states'])*3) 
        ]

        # build the external kinematics 
        PS_point, PS_jac, _, _ = generator.get_PS_point(external_phase_space_xs)
        if PS_point is None:
            logger.warning(("\nPhase-space generator failed to generate a kinematic configuration for the following input xs:\n%s"+
                            "\ncorresponding to the following input random variables of the complete integrand:\n%s\n"+
                            "Immediately returning zero weight for this point.")%(
                                external_phase_space_xs, random_variables
                            ))
            return 0.0

        if self.debug: logger.debug('PS xs=%s'%str(external_phase_space_xs))
        if self.debug: 
            incoming_edge_names = [e for (e,d) in topology_info['incoming_edges']]
            if len(incoming_edge_names)==1:
                incoming_edge_names *= 2
            outgoing_edge_names = [e for (e,d) in topology_info['outgoing_edges']]
            logger.debug('PS_point=\n%s'%LorentzVectorList(PS_point).__str__(n_initial=2, leg_names=incoming_edge_names+outgoing_edge_names))

        inv_aL_jacobian = 1.
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

        loop_phase_space_xs = [ 
            random_variables[self.dimensions.get_position('x_%d'%i_dim)] 
            for i_dim in range(
                len(channel_info['independent_final_states'])*3,
                self.SG['topo']['n_loops']*3
            ) 
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

        upscaled_input_momenta_in_LMB = [ k*rescaling_t for k in downscaled_input_momenta_in_LMB ]
        if self.debug: logger.debug('upscaled_input_momenta_in_LMB=%s'%str(upscaled_input_momenta_in_LMB))

        # Applying the rescaling to embed the momenta in the full R^(3L) space
        # and thus effectively compute the inverse of the jacobian that will be included by the full integrand in rust

        # First the delta function jacobian
        delta_jacobian = 0.
        inv_rescaling_t = 1./rescaling_t
        for edge_name, edge_direction in topology_info['outgoing_edges']:
            
            k = sum([ upscaled_input_momenta_in_LMB[i_k]*wgt for i_k, wgt in enumerate(self.SG['edge_signatures'][edge_name][0]) ])
            shift = sum([ self.external_momenta[i_shift]*wgt for i_shift, wgt in enumerate(self.SG['edge_signatures'][edge_name][1]) ])

            delta_jacobian += (
                ( 2. * inv_rescaling_t * k.square() + k.dot(shift) ) / 
                ( 2. * math.sqrt( (inv_rescaling_t*k + shift).square() - self.edge_masses[edge_name]**2 ) )
            )

        if self.debug: logger.debug('delta_jacobian=%s'%str(delta_jacobian))
        if self.debug: logger.debug('1./delta_jacobian=%s'%str(1./delta_jacobian))
        inv_aL_jacobian *= delta_jacobian

        # Then the inverse of the H-function and of the Jacobian of the causal flow change of variables
        inv_aL_jacobian *=  ( 1. / self.h_function.PDF(1./rescaling_t) ) * (rescaling_t**(len(CMB_edges)*3)) 

        # The final jacobian must then be our param. jac together with that of t divided by the one from alphaloop.
        final_jacobian = PS_jac * loop_jac * wgt_t * inv_aL_jacobian

        if self.debug: logger.debug('PS_jac=%s'%PS_jac)
        if self.debug: logger.debug('loop_jac=%s'%loop_jac)
        if self.debug: logger.debug('wgt_t=%s'%wgt_t)
        if self.debug: logger.debug('inv_aL_jacobian=%s'%inv_aL_jacobian)
        if self.debug: logger.debug('final jacobian=%s'%final_jacobian)


        # We must undo aL parameterisation as we must provide the inputs to the aL rust integrand in x-space.
        aL_xs = []
        undo_aL_parameterisation = 1.
        for i_v, v in enumerate(upscaled_input_momenta_in_LMB):
            kx, ky, kz, inv_jac = self.rust_worker.inv_parameterize(list(v), i_v, self.E_cm**2)
            undo_aL_parameterisation *= inv_jac
            aL_xs.extend([kx, ky, kz])
        if self.debug: logger.debug('aL xs=%s'%str(aL_xs))
        if self.debug: logger.debug('undo_aL_parameterisation=%s'%str(undo_aL_parameterisation))
        # Finally actually call alphaLoop full integrand
        rust_start_time = time.time()
        re, im = self.rust_worker.evaluate_integrand( aL_xs )
        if self.debug: logger.debug('Rust integrand evaluation time = %.2f ms.'%((time.time()-rust_start_time)*1000.0))

        aL_wgt = complex(re, im) * undo_aL_parameterisation
        if self.debug: logger.debug('aL res=%s'%str(aL_wgt))

        if self.debug:
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
                    reconstituted_res.real/aL_wgt.real if abs(aL_wgt.real) > 0. else reconstituted_res.real,
                    reconstituted_res.imag/aL_wgt.imag if abs(aL_wgt.imag) > 0. else reconstituted_res.imag,
                )
            ))

        if self.debug: logger.debug('normalising_func=%s'%normalising_func)
        # And accumulate it to what will be the final result
        final_res = final_jacobian * normalising_func * aL_wgt

        if multi_channeling:

            # And accumulate this jacobian into the denominator of the multichanneling factor
            multichannel_factor_denominator = final_jacobian

            multichanneling_cache = {'cutkosky_cut_and_side_cache': {} }
            # Uncomment the line below to disable this dynamic cache
            # multichanneling_cache = None

            for MC_i_channel, MC_channel_info in enumerate(self.SG['SG_multichannel_info']):
                # WARNING TODO doing the double for-loop here is not so optimal because a lot of information can be computed independently
                # of the LMB chosen already. I leave this to future refinements.
                # When generating the phase-space flat, we always ever only have a single channel per cut.
                if self.external_phase_space_generation_type == 'flat' and MC_channel_info['side']!='left':
                    continue
                for MC_i_LMB, MC_LMB_info in enumerate(MC_channel_info['loop_LMBs']):
                    # Note that for debugging it is useful to uncomment the three lines below and verify explicitly that 
                    #         MC_final_jacobian = final_jacobian
                    # for (MC_i_channel, MC_i_LMB) == (i_channel, i_LMB)
                    # because this is indeed a non-trivial check of all the inversion of the parameterisations used here.
                    #if (MC_i_channel, MC_i_LMB) == (i_channel, i_LMB):
                    #    # This contributions was already accounted for as part of final_jacobian
                    #    continue
                    
                    MC_final_jacobian = self.get_jacobian_for_channel(
                        MC_i_channel, MC_i_LMB, upscaled_input_momenta_in_LMB, PS_initial_state_four_momenta, multichanneling_cache = multichanneling_cache)

                    # Aggregate that final jacobian for this channel to the multichanneling denominator
                    multichannel_factor_denominator += MC_final_jacobian
        
            # Include the multichanneling factor
            if self.debug: logger.debug('Multi-channeling factor=%s'%( final_jacobian / multichannel_factor_denominator ))
            final_res *= final_jacobian / multichannel_factor_denominator

        if self.debug: logger.debug( 'Final weight returned = %s (phase = %s)'%( final_res, self.phase) )

        if self.phase=='real':
            return final_res.real
        else:
            return final_res.imag

    def get_jacobian_for_channel(self, MC_i_channel, MC_i_LMB, upscaled_input_momenta_in_LMB, PS_initial_state_four_momenta, multichanneling_cache=None):

        if self.debug: logger.debug('> Computing jacobian for CC cut and side #%d and LMB index #%d from the follwing upscaled LMB momenta:\n%s'%(
                    MC_i_channel, MC_i_LMB, str(upscaled_input_momenta_in_LMB)
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
        if self.debug and self.external_phase_space_generation_type == 'advanced':
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
            MC_inv_rescaling_t = self.solve_soper_flow(
                MC_upscaled_final_state_configurations, self.hyperparameters['CrossSection']['incoming_momenta'] )
            MC_rescaling_t = 1./MC_inv_rescaling_t

            if self.debug: logger.debug('MC_rescaling_t=%s'%str(MC_rescaling_t))

            MC_x_t, MC_inv_wgt_t = self.h_function.inverse_sampling(MC_rescaling_t)
            MC_wgt_t = 1./MC_inv_wgt_t

            if self.debug: logger.debug('MC_downscaled_input_momenta_in_CMB=%s'%str({e: v*MC_inv_rescaling_t for e,v in MC_upscaled_input_momenta_in_CMB.items()}))

            # Note that MC_normalising_func*MC_wgt_t = 1 when the sampling PDF matches exactly the h function, but it won't if approximated.
            #MC_normalising_func = self.h_function.PDF(MC_rescaling_t)

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

            MC_xs, MC_PS_jac= MC_generator.invertKinematics( self.E_cm, vectors.LorentzVectorList(list(PS_initial_state_four_momenta)+MC_downscaled_final_state_four_momenta) )
            if self.debug: logger.debug('MC PS xs=%s'%str(MC_xs))

            # Correct for the 1/(2E) of each cut propagator that alphaLoop includes but which is already included in the normal PS parameterisation
            MC_inv_aL_jacobian = 1.
            for v in MC_downscaled_final_state_four_momenta:
                MC_inv_aL_jacobian *= 2*v[0]

            MC_delta_jacobian = 0.
            for edge_name, edge_direction in MC_topology_info['outgoing_edges']:

                k = sum([ upscaled_input_momenta_in_LMB[i_k]*wgt for i_k, wgt in enumerate(self.SG['edge_signatures'][edge_name][0]) ])
                shift = sum([ self.external_momenta[i_shift]*wgt for i_shift, wgt in enumerate(self.SG['edge_signatures'][edge_name][1]) ])

                MC_delta_jacobian += (
                    ( 2. * MC_inv_rescaling_t * k.square() + k.dot(shift) ) / 
                    ( 2. * math.sqrt( (MC_inv_rescaling_t*k + shift).square() - self.edge_masses[edge_name]**2 ) )
                )
            if self.debug: logger.debug('MC_delta_jacobian=%s'%str(MC_delta_jacobian))
            MC_inv_aL_jacobian *= MC_delta_jacobian

            # Then the inverse of the H-function and of the Jacobian of the causal flow change of variables
            MC_inv_aL_jacobian *=  ( 1. / self.h_function.PDF(1./MC_rescaling_t) ) * (MC_rescaling_t**(len(MC_CMB_edges)*3)) 

            # Store all these quantities being always identical independently of the particular LMB chosen for future use.
            multichanneling_cache['cutkosky_cut_and_side_cache'][MC_i_channel] = {}
            multichanneling_cache['cutkosky_cut_and_side_cache'][MC_i_channel]['MC_inv_rescaling_t'] = MC_inv_rescaling_t
            multichanneling_cache['cutkosky_cut_and_side_cache'][MC_i_channel]['MC_PS_jac'] = MC_PS_jac
            multichanneling_cache['cutkosky_cut_and_side_cache'][MC_i_channel]['MC_wgt_t'] = MC_wgt_t
            multichanneling_cache['cutkosky_cut_and_side_cache'][MC_i_channel]['MC_inv_aL_jacobian'] = MC_inv_aL_jacobian

        loop_phase_space_momenta = []
        for edge_name in MC_LMB_info['loop_edges']:
            loop_phase_space_momenta.extend(list(MC_upscaled_input_momenta_in_CMB[edge_name]*MC_inv_rescaling_t))

        MC_loop_xs, MC_loop_inv_jac = self.loop_inv_parameterisation(loop_phase_space_momenta)
        MC_loop_jac = 1./MC_loop_inv_jac
        if self.debug: logger.debug('MC_loop_phase_space_momenta=%s'%str(loop_phase_space_momenta))
        if self.debug: logger.debug('MC_loop_xs=%s'%str(MC_loop_xs))

        # Return the final jacobian of the parameterisation
        MC_final_jacobian = MC_PS_jac * MC_loop_jac * MC_wgt_t * MC_inv_aL_jacobian # * MC_normalising_func

        if self.debug: logger.debug('MC_PS_jac=%s'%MC_PS_jac)
        if self.debug: logger.debug('MC_loop_jac=%s'%MC_loop_jac)
        if self.debug: logger.debug('MC_wgt_t=%s'%MC_wgt_t)
        if self.debug: logger.debug('MC_inv_aL_jacobian=%s'%MC_inv_aL_jacobian)
        if self.debug: logger.debug('MC_final_jacobian=%s'%MC_final_jacobian)
        return MC_final_jacobian

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