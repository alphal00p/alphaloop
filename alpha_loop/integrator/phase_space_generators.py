##########################################################################################
#
# Copyright (c) 2017 The MadGraph5_aMC@NLO Development team and Contributors
#
# This file is a part of the MadGraph5_aMC@NLO project, an application which
# automatically generates Feynman diagrams and matrix elements for arbitrary
# high-energy processes in the Standard Model and beyond.
#
# It is subject to the MadGraph5_aMC@NLO license which should accompany this
# distribution.
#
# For more information, visit madgraph.phys.ucl.ac.be and amcatnlo.web.cern.ch
#
##########################################################################################

import sys
import os
from pprint import pprint, pformat

if __name__ == '__main__':
    sys.path.append(os.path.join(
        os.path.dirname(os.path.realpath(__file__)), os.path.pardir, os.path.pardir ))

import logging
import math
import random
from LTD.vectors import Vector, LorentzVector, LorentzVectorList

MADEVENT= False
import madgraph.various.misc as misc
import alpha_loop.integrator.integrands as integrands
import madgraph.core.base_objects as base_objects
from madgraph import InvalidCmd, MadGraph5Error

logger = logging.getLogger('madgraph.PhaseSpaceGenerator')

class PhaseSpaceGeneratorError(Exception):
    """Exception raised if an exception is triggered in integrators.""" 
    pass

#=========================================================================================
# Phase space generation
#=========================================================================================

class VirtualPhaseSpaceGenerator(object):

    def __init__(self, initial_masses, final_masses,
                 beam_Es, 
                 beam_types=(1,1),
                 is_beam_factorization_active=(False, False),
                 correlated_beam_convolution = False,
                 dimensions=None,
                 **opts
                ):
        
        self.initial_masses  = initial_masses
        self.masses          = final_masses
        self.n_initial       = len(initial_masses)
        self.n_final         = len(final_masses)
        self.beam_Es         = beam_Es
        self.collider_energy = sum(beam_Es)
        self.beam_types      = beam_types
        self.is_beam_factorization_active = is_beam_factorization_active
        self.correlated_beam_convolution = correlated_beam_convolution
        # Sanity check
        if self.correlated_beam_convolution and self.is_beam_factorization_active != (True, True):
            raise PhaseSpaceGeneratorError(
                'The beam convolution cannot be set to be correlated if it is one-sided only')
        if dimensions is None:
            self.dimensions      = self.get_dimensions()
        else:
            self.dimensions = dimensions
        self.dim_ordered_names = [d.name for d in self.dimensions]

        self.dim_name_to_position = dict((d.name,i) for i, d in enumerate(self.dimensions))
        self.position_to_dim_name = dict((v,k) for (k,v) in list(self.dim_name_to_position.items()))
        
    def generateKinematics(self, E_cm, random_variables):
        """Generate a phase-space point with fixed center of mass energy."""

        raise NotImplementedError
    
    def get_PS_point(self, random_variables, **opts):
        """Generate a complete PS point, including Bjorken x's,
        dictating a specific choice of incoming particle's momenta."""

        raise NotImplementedError

    def boost_to_lab_frame(self, PS_point, xb_1, xb_2):
        """Boost a phase-space point from the COM-frame to the lab frame, given Bjorken x's."""
        
        if self.n_initial == 2 and (xb_1!=1. or xb_2!=1.):
            ref_lab = (PS_point[0]*xb_1 + PS_point[1]*xb_2)
            if ref_lab.rho2() != 0.:
                lab_boost = ref_lab.boostVector()
                for p in PS_point:
                    p.boost(-lab_boost)

    def boost_to_COM_frame(self, PS_point):
        """Boost a phase-space point from the lab frame to the COM frame"""
        
        if self.n_initial == 2:
            ref_com = (PS_point[0] + PS_point[1])
            if ref_com.rho2() != 0.:
                com_boost = ref_com.boostVector()
                for p in PS_point:
                    p.boost(-com_boost)

    def nDimPhaseSpace(self):
        """Return the number of random numbers required to produce
        a given multiplicity final state."""

        if self.n_final == 1:
            return 0
        return 3*self.n_final - 4

    def get_dimensions(self):
        """Generate a list of dimensions for this integrand."""
        
        dims = integrands.DimensionList()

        # Add the PDF dimensions if necessary
        if self.beam_types[0]==self.beam_types[1]==1:
            dims.append(integrands.ContinuousDimension('ycms',lower_bound=0.0, upper_bound=1.0))
            # The 2>1 topology requires a special treatment
            if not (self.n_initial==2 and self.n_final==1):
                dims.append(integrands.ContinuousDimension('tau',lower_bound=0.0, upper_bound=1.0)) 

        # Add xi beam factorization convolution factors if necessary
        if self.correlated_beam_convolution:
            # A single convolution factor xi that applies to both beams is needed in this case
            dims.append(integrands.ContinuousDimension('xi',lower_bound=0.0, upper_bound=1.0))
        else:
            if self.is_beam_factorization_active[0]:
                dims.append(integrands.ContinuousDimension('xi1',lower_bound=0.0, upper_bound=1.0))             
            if self.is_beam_factorization_active[1]:
                dims.append(integrands.ContinuousDimension('xi2',lower_bound=0.0, upper_bound=1.0))

        # Add the phase-space dimensions
        dims.extend([ integrands.ContinuousDimension('x_%d'%i,lower_bound=0.0, upper_bound=1.0) 
                                     for i in range(1, self.nDimPhaseSpace()+1) ])
        
        return dims

class MultiChannelPhasespace(VirtualPhaseSpaceGenerator):
    """ A phase space generator for a channel in multi-channel integration."""
        
    def __init__(self,*args,**opts):
        
        if 'model' not in opts:
            raise PhaseSpaceGeneratorError("A model must be specified with the option 'model' when"+
                " instantiating the class %s."%self.__class__.__name__)
        self.model = opts.pop('model')
        
        if 'topologies' not in opts:
            raise PhaseSpaceGeneratorError("A list of topologies must be specified with the "+
                "option 'topologies' when instantiating the class %s."%self.__class__.__name__)
        self.topologies = opts.pop('topologies')
        
        super(MultiChannelPhasespace, self).__init__(*args, **opts)
        
        self.channels = []
        for topology in self.topologies:
            opts = {'initial_masses': self.initial_masses, 'final_masses': self.masses,
                    'beam_Es': self.beam_Es, 'beam_types': self.beam_types, 'model': self.model}
            self.channels.append(SingleChannelPhasespace(topology=topology,**opts))
    
    def get_PS_point(self, random_variables, adaptive_wgts=None,channel_nr = None):
        """Provides a momentum configuration according to the right phase space parameterization
        and the multi-channel weight (the Jacobians from all channels with their channel weights (alpha)). """
        """ adaptive_wgts = channel wgts (alpha) """
    
        #TODO: flattening_strategy = {channel_id: {'flattening_technique': 'diagram_fuction to call'/'jacobians', 'use_alphas': True/False}}
        #TODO: get_flattener_for_channel(id)
                                
        if random_variables is None:
            random_variables = self.dimensions.random_sample()
        if adaptive_wgts is None:
            raise PhaseSpaceGeneratorError('Specify the channel weights.')
            #adaptive_wgts = [1./len(self.channels)]*len(self.channels)
        if channel_nr is None:
            raise PhaseSpaceGeneratorError('Specify the channel number.')
            #randomly pick a channel
            #channel_nr = random.randint(0,len(self.channels)-1)
        PS_random_variables = random_variables
        
        channel_wgts = [None]*len(self.channels)

        PS_point, channel_wgts[channel_nr], xb_1, xb_2 = self.channels[channel_nr].get_PS_point(PS_random_variables)#, path=None)
        if PS_point == None:
            return None, 0., 1., 1.
        
        for i, channel in enumerate(self.channels):
            if i != channel_nr:
                variables, channel_wgts[i] = channel.get_PS_point(PS_point)
                if channel_wgts[i] == 0.:
                    return None, 0, 1., 1.
            
        multi_channel_wgt = 1./sum( alpha / wgt_jac for alpha, wgt_jac in zip(adaptive_wgts,channel_wgts))
        
        #if channel_nr == 1:
        #    multi_channel_wgt*= 1e-9
    
        return PS_point, multi_channel_wgt, xb_1, xb_2
        
        #matrix element multi-channeling idea
        """
            full_path = '/Users/Dario/Desktop/Thesis/HighEnergySoftwares/MadGraph5/MadEvent6/mytest_MCPS_full/SubProcesses/P1_epem_mupmuma'
            sys.path.append(full_path)
            import matrix6py as full_matrix_element
            full_matrix_element.initialisemodel(os.path.abspath(os.path.join(full_path,os.pardir,os.pardir,'Cards','param_card.dat')))

            def invert_momenta(p):
                #fortran/C-python do not order table in the same order
                new_p = []
                for i in range(len(p[0])):  new_p.append([0]*len(p))
                for i, onep in enumerate(p):
                    for j, x in enumerate(onep):
                        new_p[j][i] = x
                return new_p
            
            P = invert_momenta(PS_point)
            matrix_elements = [None]*len(self.channels)
            for i in xrange(len(self.channels)):
                # diagram i might not correspond to function i+1, CHECK
                matrix_elements[i] = full_matrix_element.smatrix(P,i+1)
            
            multi_channel_wgt = channel_wgts[channel_nr]*matrix_elements[channel_nr]
            multi_channel_wgt *= 1./sum( alpha*me for alpha, me in zip(alphas,matrix_elements))
            return PS_point,multi_channel_wgt, xb_1, xb_2
        """
            

class SingleChannelPhasespace(VirtualPhaseSpaceGenerator):
    """Implementation of a phase-space generator that lines up integration variables
    with s- and t-channels specifying a particular provided topology.
    """
    
    # The lowest value that the center of mass energy can take.
    # Below 1 GeV non-perturbative effects dominate and factorization does not apply
    absolute_Ecm_min = 1.
    
    def __init__(self, *args, **opts):
        
        if 'model' not in opts:
            raise PhaseSpaceGeneratorError("A model must be specified with the option 'model' when"+
                " instantiating the class %s."%self.__class__.__name__)
        self.model = opts.pop('model')
        
        if 'topology' not in opts:
            raise PhaseSpaceGeneratorError("A list of s- and t-channels must be specified with the "+
                "option 'topology' when instantiating the class %s."%self.__class__.__name__)
        self.topology = opts.pop('topology')

        super(SingleChannelPhasespace, self).__init__(*args, **opts)

        if 'path' not in opts or opts['path'] is None:
            # In order to make resutls deterministic, force using the first available path as opposed 
            # as a random one.
            path = self.get_random_path(select_first=True)
        else:
            path = opts.pop('path')

        self.path = path
        
        ##topology_string = self.get_topology_string(self.topology,path_to_print=self.path)
        ##misc.sprint(topology_string)
        
        """
        N_t = len(self.topology[1]) - 1 # number of t-channels
        N_s = len(self.topology[0]) # number of s-channels
        
        min_index = 0
        max_index = N_s
        for i in range(min_index,max_index):
            self.dimensions[self.dim_name_to_position['x_%d' % (i+1)]].name = 'distr_s_%d' % (i+1-min_index)
        min_index = max_index
        max_index += N_t -1 
        for i in range(min_index,max_index):
            self.dimensions[self.dim_name_to_position['x_%d' % (i+1)]].name = 'unif_s_%d' % (i+1-min_index)
        min_index = max_index
        max_index += N_t
        for i in range(min_index,max_index):
            self.dimensions[self.dim_name_to_position['x_%d' % (i+1)]].name = 't_%d' % (i+1-min_index)
        min_index = max_index
        max_index += N_s
        for i in range(min_index,max_index):
            self.dimensions[self.dim_name_to_position['x_%d' % (i+1)]].name = 'cos_theta_%d' % (i+1-min_index)
        min_index = max_index
        max_index += N_t + N_s
        for i in range(min_index,max_index):
            self.dimensions[self.dim_name_to_position['x_%d' % (i+1)]].name = 'phi_%d' % (i+1-min_index)
        """
        
        # One can do additional business here upon instantiating this PS generator, like renaming
        # the random variables to names describing the "propagators generated" with them. 
         
    def get_topology_string(self,topology_to_print,path_to_print=None):
        """Example of a nice way to printout what these topologies are:"""
        if  topology_to_print[0] == None:
            topology_string = '\n no s-channels \nand t-channels: %s'%\
                (', '.join('%s > %d(%d)'%(
                    ' '.join('%d(%d)'%(leg['number'],leg['id']) for leg in vertex['legs'][:-1]),
                    vertex['legs'][-1]['number'],vertex['legs'][-1]['id']) for vertex in topology_to_print[1]))
        else:
            topology_string = '\ns-channels:     %s\nand t-channels: %s'%\
                (', '.join('%s > %d(%d)'%(
                    ' '.join('%d(%d)'%(leg['number'],leg['id']) for leg in vertex['legs'][:-1]),
                    vertex['legs'][-1]['number'],vertex['legs'][-1]['id']) for vertex in topology_to_print[0]), 
                  ', '.join('%s > %d(%d)'%(
                    ' '.join('%d(%d)'%(leg['number'],leg['id']) for leg in vertex['legs'][:-1]),
                    vertex['legs'][-1]['number'],vertex['legs'][-1]['id']) for vertex in topology_to_print[1]))
        if path_to_print != None:
            topology_string += '\nselected path:  %s'%path_to_print
        return topology_string
    
    def uniform_distr(self,r,min,max):
        """distributes r uniformly within (min, max), with jacobian dvariable"""
        dvariable = (max-min)
        variable = min + dvariable*r
        return variable, dvariable
    
    def inv_uniform_distr(self,variable,min,max):
        """inverse of uniform_distr, obtain r and inverse jacobian dvariable from variable"""
        assert(min<=variable<= max)
        dvariable = (max-min)
        r = (variable-min)/dvariable
        return r, dvariable
    
    def massless_distr(self,r,min,max,nu=1.1,m2 = 0):
        """distributes r within (min, max), with jacobian dvariable \propto variable^(nu)
        for cross-section \propto 1/(s-m)^2
        m2 is a small parameter < 0, that fixes numerical problems when min=0,
        instead of setting min=small number, this method still allows to map to min=0"""
        if min == 0 and m2==0 or min==m2:
            m2 = -self.absolute_Ecm_min**2
        if nu != 1:            
            variable = ((max-m2)**(1.-nu)*r+(min-m2)**(1.-nu)*(1.-r))**((1.-nu)**(-1.))+m2
            dvariable = ((variable-m2)**nu)*((max-m2)**(1.-nu)-(min-m2)**(1.-nu))/(1.-nu)
        else:
            variable = math.exp(r*math.log(max-m2)+(1.-r)*math.log(min-m2))+m2
            dvariable = (math.log(max-m2)-math.log(min-m2))*(variable-m2)
        assert(dvariable>=0)
        return variable, dvariable
    
    def inv_massless_distr(self,variable,min,max,nu=1.1,m2=0):
        """inverse of massless_distr, obtain r and inverse jacobian dvariable from variable"""
        assert(min<=variable<=max)
        if (min == 0 and m2==0) or (min==m2):
            m2 = -self.absolute_Ecm_min**2
        if nu != 1:
            r = ((variable-m2)**(1.-nu)-(min-m2)**(1.-nu))/((max-m2)**(1.-nu)-(min-m2)**(1.-nu))
            dvariable = ((variable-m2)**nu)*((max-m2)**(1.-nu)-(min-m2)**(1.-nu))/(1.-nu)
        else:
            r = (math.log(variable-m2)-math.log(min-m2))/(math.log(max-m2)-math.log(min-m2))
            dvariable = (math.log(max-m2)-math.log(min-m2))*(variable-m2)
        assert(dvariable>=0)
        return r, dvariable
    
    def massive_distr(self,r,mass,width,min,max):
        """distributes r within (min, max), with jacobian dvariable \propto 1/Breit-Winger
        for cross-section \propto Breit-Wigner"""
        mass = mass.real
        width = width.real
        y_1 = math.atan((min-mass**2)/(mass*width))
        y_2 = math.atan((max-mass**2)/(mass*width))               
        variable = mass**2 + mass*width*math.tan(y_1+(y_2-y_1)*r)
        dvariable = (y_2-y_1)*((variable-mass**2)**2 + (mass*width)**2)/(mass*width)
        # equivalently: dvariable = mass*width*(y_2-y_1)*math.cos(y_1+(y_2-y_1)*r)**(-2)
        assert(dvariable>=0)
        return variable, dvariable
    
    def inv_massive_distr(self,variable,mass,width,min,max):
        """inverse of massive_distr, obtain r and inverse jacobian dvariable from variable"""
        assert(min<=variable<=max)
        mass = mass.real
        width = width.real
        y_1 = math.atan((min-mass**2)/(mass*width))
        y_2 = math.atan((max-mass**2)/(mass*width))
        r = (math.atan((variable-mass**2)/(mass*width))-y_1)/(y_2-y_1)
        dvariable = (y_2-y_1)*((variable-mass**2)**2 + (mass*width)**2)/(mass*width)
        assert(dvariable>=0)
        return r, dvariable
    
    def inv_t(self,p_2,p1_2,p2_2,p3_2,p4_2,cos_theta):
        """Mandelstam invariant t=(p1-p3)^2 formula C21 in https://arxiv.org/pdf/hep-ph/0008033.pdf
        p=p1+p2 is at rest;
        p1, p2 are opposite along z-axis
        p3, p4 are opposite along the direction defined by theta
        theta is the angle in the center of mass frame between p1 & p3"""
        nom = (p_2+p3_2-p4_2)*(p_2+p1_2-p2_2) - math.sqrt(self.Lambda(p_2,p3_2,p4_2))*math.sqrt(self.Lambda(p_2,p1_2,p2_2))*cos_theta
        t = p3_2+p1_2 - nom/(2*p_2)
        if t>0:
            t = 0
        assert(t<=0)
        return t
    
    def cos_theta_from_inv_t(self,p_2,p1_2,p2_2,p3_2,p4_2,t):
        """https://arxiv.org/pdf/hep-ph/0008033.pdf forula C21
        invert t=(p1-p3)^2 to cos_theta = ..."""
        nom = (t-p3_2-p1_2)*2*p_2 + (p_2+p3_2-p4_2)*(p_2+p1_2-p2_2)
        denom = math.sqrt(self.Lambda(p_2,p3_2,p4_2))*math.sqrt(self.Lambda(p_2,p1_2,p2_2))
        cos_theta = nom/denom
        assert(-1<=cos_theta<=1)
        return cos_theta
    
    def angles_to_rotate_along_z(self,p):
        """gives the angles phi and theta,
        so that \vec{p} can be parametrized as |\vec{p}|*(cos phi*sin theta, sin phi*sin theta, cos theta)"""
        theta = math.acos(p[3]/p.rho())
        if p[1]> 0:
            phi = math.atan(p[2]/p[1])
        elif p[1]< 0:
            phi = math.atan(p[2]/p[1]) + math.pi
        else:
            phi = 0
        return phi,theta
    
    def rotate_along_z_inv(self,p,phi,theta):
        """rotates a 4 vector p in space
        the z-axis (0,0,1) rotates to (cos phi*sin theta, sin phi*sin theta, cos theta)"""
        l = LorentzVector([p[0],0,0,0])
        l[1] = math.cos(theta)*math.cos(phi)*p[1] -math.sin(phi)*p[2] +math.sin(theta)*math.cos(phi)*p[3]
        l[2] = math.cos(theta)*math.sin(phi)*p[1] + math.cos(phi)*p[2] + math.sin(theta)*math.sin(phi)*p[3]
        l[3] = -math.sin(theta)*p[1] +math.cos(theta)*p[3]
        return l
    
    def Lambda(self,x,y,z):
        return x**2+y**2+z**2-2*x*y-2*x*z-2*y*z
    
    def get_cm_momenta(self,s,p3_2,p4_2,cos_theta,phi):
        """generates p3, p4 in the center of mass frame, with E_cm^2 = s
        the direction is defined by the angles cos_theta and phi"""
        assert(-1 <= cos_theta <= 1.)
        assert(s>0)
        q = 1./(2*math.sqrt(s))
        sin_theta = math.sqrt(1-cos_theta**2)
        p3 = LorentzVector([1.,sin_theta*math.cos(phi),sin_theta*math.sin(phi), cos_theta])
        p4 = LorentzVector([1.,-sin_theta*math.cos(phi),-sin_theta*math.sin(phi), -cos_theta])
        p3[0] *= (s+p3_2-p4_2)*q
        p4[0] *= (s+p4_2-p3_2)*q
        assert(self.Lambda(s,p3_2,p4_2)>0)
        rho = math.sqrt(self.Lambda(s,p3_2,p4_2))
        for i in range(1,4):
            p3[i] *= rho*q
            p4[i] *= rho*q
        wgt_PS = rho/(8*s)
        return p3,p4,wgt_PS
    
    def get_two_body_PS_wgt(self,s,p3_2,p4_2):
        """returns the two body PS_wgt that is also generated in generate_cm_momenta"""
        assert(s>0)
        assert(self.Lambda(s,p3_2,p4_2)>0)
        rho = math.sqrt(self.Lambda(s,p3_2,p4_2))
        wgt_PS = rho/(8*s)
        return wgt_PS
    
    def get_random_path(self, select_first=False):
        """Generates a random path in which the (distributed and uniform) invariants are sampled.
        For now, only paths are generated that go from two known outer legs to the inner leg.
        Breit-Wigner competition (in distributed invariants) is not taken into account yet.
        When setting the flag 'select_first' to True, this is made deterministic and the first possible
        path is selected instead."""

        # TODO generalise for self.n_initial=1

        max_leg_nr = self.n_initial+ self.n_final
        min_leg_nr = self.topology[1][-1].get('legs')[-1].get('number')
        numbers = list(range(min_leg_nr,0))+list(range(self.n_initial+1,max_leg_nr+1))+[1]

        """distr_inv stores vertices for distributed (=s-channel) invariants that
        can be used to generate a next invariant (available)
        were already used to generate an invariant (finished)
        third category: vertices that haven't been used already but are not available yet to generate invariants
        distr_inv['finished'] will eventually be ordered from outer (final) legs to inner legs"""
        
        distr_inv = {'available': [],'finished': []}
        # kinematics stores if a leg is available (True/False)
        kinematics = dict((nr,{'is_available': False}) for nr in numbers)
        # final legs are available
        for nr in kinematics:
            if nr > self.n_initial:
                kinematics[nr]['is_available'] = True

        # when no t-channels, first s-channel invariant is fixed
        if len(self.topology[1]) == 1:
            last_vertex = self.topology[0][-1]
            leg1_nr = last_vertex.get('legs')[0].get('number')
            leg2_nr = last_vertex.get('legs')[1].get('number')
            kinematics[leg1_nr]['is_available'] = False
            kinematics[leg2_nr]['is_available'] = False
            
        # find out which s-channel vertices are available
        for i, vertex in enumerate(self.topology[0]):
            leg1_nr = vertex.get('legs')[0].get('number')
            leg2_nr = vertex.get('legs')[1].get('number') 
            if kinematics[leg1_nr]['is_available'] and kinematics[leg2_nr]['is_available']:
                distr_inv['available'].append(i)
        
        # generate a path in which distributed (=s-channel) invariants are generated
        # in general the ordering is relevant since it changes the boundaries
        # TODO: Breit-Wigner competition, for now, (final) outer to inner vertices
        if len(self.topology[0]) > 0:
            nr_available_vertices = len(distr_inv['available'])
            while (nr_available_vertices > 0):
                if not select_first:
                    i = random.choice(distr_inv['available'])                    
                else:
                    i = distr_inv['available'][0]                    
                vertex = self.topology[0][i]
                leg1_nr = vertex.get('legs')[0].get('number')
                leg2_nr = vertex.get('legs')[1].get('number')
                leg3_nr = vertex.get('legs')[-1].get('number')
                kinematics[leg1_nr]['is_available'] = False
                kinematics[leg2_nr]['is_available'] = False
                kinematics[leg3_nr]['is_available'] = True
                # distr_inv['finished'] is ordered from outer (final) legs to inner legs
                distr_inv['finished'].append(i)
                distr_inv['available'].remove(i)   
                    
                # check if another vertex is available now
                # if no t-channels (i.e. first s-channel invariant fixed), skip the fixed (=last) vertex
                if len(self.topology[1])==1:
                    s_channel_range= len(self.topology[0][:-1])
                else:
                    s_channel_range = len(self.topology[0])
                for m in range(i,s_channel_range):
                    if m in distr_inv['finished']:
                        continue
                    vertex = self.topology[0][m]
                    leg1_nr = vertex.get('legs')[0].get('number')
                    leg2_nr = vertex.get('legs')[1].get('number')
                    if leg3_nr == leg1_nr or leg3_nr == leg2_nr:
                        if kinematics[leg2_nr+leg1_nr-leg3_nr]['is_available']:
                            distr_inv['available'].append(m)
                            break
                nr_available_vertices= len(distr_inv['available'])
        
        t_channel_path = []
        
        """find a path to generate the uniformly sampled invariants
        idea: save pairs of vertices that correspond to a uniform invariant
        randomly pick a vertex i<last, save in t_channel_path [(0,i),(i+1,last)]
        for pair in t_channel_path pick random vertex j<last, do the same thing (-> i.e. (0,j) (j+1,i))
        save in t_channel_path except either j+1 = i, or 0=j
        to find possible pairs of vertices, one needs to go from outermost to innermost t-channels."""
        
        if len(self.topology[1]) > 2:
            start = 0
            end = len(self.topology[1])-1
            if not select_first:
                i = random.choice(list(range(start,end)))
            else:
                i = list(range(start,end))[0]
            if i != 0:
                t_channel_path.append((0,i))
            if i+1 != len(self.topology[1])-1:
                t_channel_path.append((i+1,len(self.topology[1])-1))
            for tuple in t_channel_path:
                start = tuple[0]
                end = tuple[1]
                if not select_first:
                    i = random.choice(list(range(start,end)))
                else:
                    i = list(range(start,end))[0]
                if i != start:
                    t_channel_path.append((start,i))
                if i+1 != end: 
                    t_channel_path.append((i+1,end))
        
        """paths are ordered lists, their order will determine in what order the invariants will be generated
        elements of s_channel_path are the indices of the s-channel vertices
        e.g. s_channel_path = [0,4,1,2,...]
        elements of t_channel_path are tuples of indices of the t-channel vertices
        e.g. t_channel_path = [(0,4),(5,6),(0,2),(3,4),...])"""
        s_channel_path = distr_inv['finished']
        t_channel_path.reverse()
        
        return [s_channel_path,t_channel_path]
        
    def invertKinematics(self, E_cm, momenta, **opts):

        return self.get_PS_point(
            LorentzVectorList([LorentzVector(v) for v in momenta]), **opts
        )

    def get_PS_point(self, input_variables, path=None, **opts):
        """ Generates a complete PS point, including Bjorken x's, dictating a specific choice
        of incoming particle's momenta,"""
        
        # PROBLEMS:
        #    0) line 644, tolerance=...
        #    1) line 1067, s-channel momentum. pt_cut, boost to c
        #    2) line 1044,1092,1120, momentum conservation, t-channel, s-channel, tot
        #    3) line 795, 961, no_phase_space
        #    4) line 632, randomly pick one path every time, or stick with one
   
        if path is None:
            path = self.path
            #path = self.get_random_path()
                
        # kinematics is the object with all important info it it
        max_leg_nr = self.n_initial+ self.n_final
        min_leg_nr = self.topology[1][-1].get('legs')[-1].get('number')
        numbers = list(range(min_leg_nr,0))+list(range(3,max_leg_nr+1))+[1]
        kinematics = dict((nr,{'inv_mass': None, 'momentum': None, 'is_available': False, 'inv_mass_limits': None}) for nr in numbers)

        wgt = 1. #store the weight, coming from the importance sampling and reparameterization (reconstruct mode)
        tolerance = 1e-3
        
        if isinstance(input_variables, LorentzVectorList):
            reconstruct = True
            PS_point = input_variables
            output_variables = self.get_dimensions()
            kinematics[min_leg_nr]['momentum'] = PS_point[1]
            kinematics[1]['momentum'] = PS_point[0]
            for nr in range(3,max_leg_nr+1):
                kinematics[nr]['momentum'] = PS_point[nr-1]
        else:
            reconstruct = False
            random_variables = input_variables
            if random_variables is None:
                random_variables = self.dimensions.random_sample()  
            output_momenta = [LorentzVector()]*max_leg_nr
        variable_index = 0 #index to pick random variable, there might be a nicer way to do it
        
        # fill final and initial mass info into kinematics
        kinematics[min_leg_nr]['inv_mass'] = self.initial_masses[1]**2
        kinematics[1]['inv_mass'] = self.initial_masses[0]**2
        for nr in range(3,max_leg_nr+1):
            kinematics[nr]['inv_mass'] = self.masses[nr-3]**2
            kinematics[nr]['is_available'] = True #flag is necessary to calculate limits
        
        # TODO: Breit-Wigner competition: tau is always generated first for now: can lead to problems in BW comp.
        # get the bjorken x's in case of pp-collison
        if self.beam_types[0]==self.beam_types[1]==1:
            # definition: ycm = 1/2*log(xb_1/xb_2)
            # definition: tau = xb_1*xb_2
            if reconstruct:
                E_cm = math.sqrt((kinematics[1]['momentum']+kinematics[min_leg_nr]['momentum']).square())
                tau = E_cm**2/self.collider_energy**2
                ycm = None # need Bjorken x's as inputs for that
            else:
                x_ycm = random_variables[self.dim_name_to_position['ycms']]
                x_tau = random_variables[self.dim_name_to_position['tau']]
            tot_final_state_masses = sum(self.masses)
            tau_min = (max(tot_final_state_masses, self.absolute_Ecm_min)/self.collider_energy)**2
            tau_max = 1.
            if len(self.topology[1]) == 1: #if there are no t-channels, tau is distributed
                last_vertex = self.topology[0][-1]
                id = last_vertex.get('legs')[-1].get('id')
                particle = self.model.get_particle(id)
                mass_param = particle.get('mass')
                if mass_param.lower() == 'zero':
                    if reconstruct:
                        x_tau, wgt_jac = self.inv_massless_distr(tau, tau_min, tau_max)
                    else:
                        tau, wgt_jac = self.massless_distr(x_tau,tau_min,tau_max)
                else:
                    mass = self.model.get('parameter_dict')[mass_param]/self.collider_energy
                    width = self.model.get('parameter_dict')[particle.get('width')]/self.collider_energy
                    if reconstruct:
                        x_tau, wgt_jac = self.inv_massive_distr(tau, mass, width, tau_min, tau_max)
                    else:
                        tau, wgt_jac = self.massive_distr(x_tau,mass,width, tau_min, tau_max)
            else: #if there are t-channels, tau is uniform
                if reconstruct:
                    x_tau, wgt_jac = self.inv_uniform_distr(tau,tau_min,tau_max)
                else:
                    tau, wgt_jac = self.uniform_distr(x_tau,tau_min,tau_max)
            wgt *= wgt_jac
            # ycm always sampled uniformly
            ycm_min = 0.5 * math.log(tau)
            ycm_max = -ycm_min
            if reconstruct:
                x_ycm, wgt_jac = None, (ycm_max-ycm_min)
                output_variables[self.dim_name_to_position['tau']] = x_tau
                output_variables[self.dim_name_to_position['ycms']] = x_ycm
            else:
                ycm, wgt_jac = self.uniform_distr(x_ycm,ycm_min,ycm_max)
                xb_1 = math.sqrt(tau)*math.exp(ycm)
                xb_2 = math.sqrt(tau)*math.exp(-ycm)
                E_cm = math.sqrt(tau)*self.collider_energy
            wgt *= wgt_jac
        elif self.beam_types[0]==self.beam_types[1]==0: # set default in case of ll-collision
            xb_1 = 1.
            xb_2 = 1.
            E_cm = self.collider_energy
        else:
            raise InvalidCmd("This basic PS generator does not yet support collider mode (%d,%d)."%self.beam_types)  

        # Also generate the ISR collinear factorization convolutoin variables xi<i> if
        # necessary. In order for the + distributions of the PDF counterterms and integrated
        # collinear ISR counterterms to hit the PDF only (and not the matrix elements or
        # observables functions), a change of variable is necessary: xb_1' = xb_1 * xi1
        if self.correlated_beam_convolution:
            # Both xi1 and xi2 must be set equal then
            xi1 = random_variables[self.dim_name_to_position['xi']]
            xi2 = random_variables[self.dim_name_to_position['xi']]
        else:
            if self.is_beam_factorization_active[0]:
                xi1 = random_variables[self.dim_name_to_position['xi1']]
            else:
                xi1 = None
            if self.is_beam_factorization_active[1]:
                xi2 = random_variables[self.dim_name_to_position['xi2']]
            else:
                xi2 = None

        if not reconstruct:
            """generate initial momenta in CENTER OF MASS frame, along z-axis
            -> ALL output_momenta are going to be in this frame
            BUT: cuts are dependent on lab frame, therefore return Bjorken x's too"""
            leg1_nr = 1
            leg2_nr = min_leg_nr
            kinematics[leg1_nr]['momentum'] = E_cm/2.*LorentzVector([1.,0,0, 1.])
            kinematics[leg2_nr]['momentum'] = E_cm/2.*LorentzVector([1.,0,0, -1.])
        
        # if there are no t-channels do the first s-channel, its invariant mass is fixed
        if len(self.topology[1])==1:
            # set the momentum of the very first s-channel propagator (p_initial1+p_initial2)
            last_vertex = self.topology[0][-1]
            leg1_nr = 1
            leg2_nr = min_leg_nr
            leg3_nr = last_vertex.get('legs')[-1].get('number')
            kinematics[leg3_nr]['inv_mass'] = E_cm**2
            kinematics[leg3_nr]['momentum'] = kinematics[leg1_nr]['momentum'] + kinematics[leg2_nr]['momentum']
          
        # set the invariant masses for the s-channels, all distributed
        # remember: if there is only one s-channel, there is no invariant to generate
        if len(self.topology[0]) > 0:
            for i in path[0]:
                vertex = self.topology[0][i]
                leg1_nr = vertex.get('legs')[0].get('number')
                leg2_nr = vertex.get('legs')[1].get('number')
                leg3_nr = vertex.get('legs')[-1].get('number')

                kinematics[leg1_nr]['is_available'] = False #flag is necessary to calculate limits
                kinematics[leg2_nr]['is_available'] = False
                sum_masses = sum([math.sqrt(kinematics[nr]['inv_mass']) for nr in list(kinematics.keys()) if kinematics[nr]['is_available']])
                s_max = (E_cm - sum_masses)**2
                m1 = math.sqrt(kinematics[leg1_nr]['inv_mass'])
                m2 = math.sqrt(kinematics[leg2_nr]['inv_mass'])
                s_min = (m1 + m2)**2
                kinematics[leg3_nr]['inv_mass_limits'] = (s_min,s_max)
                #use importance sampling
                id = vertex.get('legs')[-1].get('id')
                particle = self.model.get_particle(id)
                mass_param = particle.get('mass')
                variable_index += 1
                if reconstruct:
                    assert(kinematics[leg3_nr]['momentum']==None)
                    kinematics[leg3_nr]['momentum'] = kinematics[leg1_nr]['momentum']+kinematics[leg2_nr]['momentum']
                    s = kinematics[leg3_nr]['momentum'].square()
                    if 0. < s_min-s < tolerance**2*self.absolute_Ecm_min**2:
                        s = s_min
                    if 0. < s-s_max < tolerance**2*s_max:
                        s = s_max
                    kinematics[leg3_nr]['inv_mass'] = s
                    if mass_param.lower() == 'zero':
                        x_s, wgt_jac = self.inv_massless_distr(s,s_min,s_max)
                    else:
                        mass = self.model.get('parameter_dict')[mass_param]
                        width = self.model.get('parameter_dict')[particle.get('width')]
                        if width == 0:
                            x_s, wgt_jac = self.inv_massless_distr(s, s_min, s_max, m2=mass**2)
                        else:
                            x_s, wgt_jac = self.inv_massive_distr(s,mass,width, s_min, s_max)
                    output_variables[self.dim_name_to_position['x_%d' % variable_index]] = x_s
                else:
                    assert(kinematics[leg3_nr]['inv_mass'] == None)
                    random_variable = random_variables[self.dim_name_to_position['x_%d' % variable_index]]
                    if mass_param.lower() == 'zero':
                            s, wgt_jac = self.massless_distr(random_variable,s_min,s_max)
                    else:
                        mass = self.model.get('parameter_dict')[mass_param]
                        width = self.model.get('parameter_dict')[particle.get('width')]
                        if width == 0:
                            s, wgt_jac = self.massless_distr(random_variable, s_min, s_max, m2=mass**2)
                        else:
                            s, wgt_jac = self.massive_distr(random_variable,mass,width, s_min, s_max)
                    # think about a better solution to handle if particle is produced at rest
                    no_phase_space = 0.*max([s,m1**2,m2**2])**2 >= self.Lambda(s,m1**2,m2**2) >= -tolerance**2*s**2
                    if no_phase_space:
                        misc.sprint('This kinematic configuration leaves no phase space ', self.Lambda(s,m1**2,m2**2))
                        return None,0.,(xb_1, xi1) , (xb_2, xi2)
                    kinematics[leg3_nr]['inv_mass'] = s

                wgt *= wgt_jac
                kinematics[leg3_nr]['is_available'] = True
        
        # if there is more than one t-channel, there are uniform invariants to be generated
        uniform_inv = []
        if len(self.topology[1]) > 2:
            # go from innermost to outermost vertex pair in t-channel path
            for item in path[1]:
                # extract the path-info into uniform_inv dictionary
                uniform_inv.append({'vertices': (item[0],item[1]), 'inv_mass': None, 'is_available': False, 'inv_mass_limits': None, 'momentum': None})
            for k,item in enumerate(uniform_inv):
                start = item['vertices'][0]
                finish = item['vertices'][1]
                if (start+1) == finish:
                    vertex1 = self.topology[1][start]
                    vertex2 = self.topology[1][finish]
                    leg1_nr = vertex1.get('legs')[1].get('number')
                    leg2_nr = vertex2.get('legs')[1].get('number')
                    m1 = math.sqrt(kinematics[leg1_nr]['inv_mass'])
                    m2 =  math.sqrt(kinematics[leg2_nr]['inv_mass'])
                    kinematics[leg1_nr]['is_available'] = False #flag is necessary to calculate limits
                    kinematics[leg2_nr]['is_available'] = False
                    if reconstruct:
                        p1 = kinematics[leg1_nr]['momentum']
                        p2 = kinematics[leg2_nr]['momentum']
                else:
                    for l in reversed(range(0,k)):
                        if uniform_inv[l]['vertices'][0] == start:
                            m1 = math.sqrt(uniform_inv[l]['inv_mass'])
                            if reconstruct:
                                p1 = uniform_inv[l]['momentum']
                            uniform_inv[l]['is_available'] = False
                            break
                    else:
                        vertex1 = self.topology[1][start]
                        leg1_nr = vertex1.get('legs')[1].get('number')
                        m1 = math.sqrt(kinematics[leg1_nr]['inv_mass'])
                        if reconstruct:
                                p1 = kinematics[leg1_nr]['momentum']
                        kinematics[leg1_nr]['is_available'] = False
                    for l in reversed(range(0,k)):
                        if uniform_inv[l]['vertices'][1] == finish:
                            m2 = math.sqrt(uniform_inv[l]['inv_mass'])
                            if reconstruct:
                               p2 = uniform_inv[l]['momentum']
                            uniform_inv[l]['is_available'] = False
                            break
                    else:
                        vertex2 = self.topology[1][finish]
                        leg2_nr = vertex1.get('legs')[1].get('number')
                        m2 = math.sqrt(kinematics[leg1_nr]['inv_mass'])
                        if reconstruct:
                            p2 = kinematics[leg2_nr]['momentum']
                        kinematics[leg2_nr]['is_available'] = False
                sum_masses = sum([math.sqrt(kinematics[nr]['inv_mass']) for nr in list(kinematics.keys()) if kinematics[nr]['is_available']])
                sum_masses += sum([math.sqrt(inv['inv_mass']) for inv in uniform_inv if inv['is_available']])
                s_max = (E_cm - sum_masses)**2
                s_min = (m1 + m2)**2
                item['inv_mass_limits']= (s_min,s_max)
                variable_index += 1
                if reconstruct:
                    item['momentum']=p1+p2
                    s = item['momentum'].square()
                    if 0. < s_min-s < tolerance**2*self.absolute_Ecm_min**2:
                        s = s_min
                    if 0. < s-s_max < tolerance**2*s_max:
                        s = s_max
                    item['inv_mass'] = s
                    x_s, wgt_jac = self.inv_uniform_distr(s,s_min,s_max)
                    output_variables[self.dim_name_to_position['x_%d' % variable_index]] = x_s
                else:
                    random_variable = random_variables[self.dim_name_to_position['x_%d' % variable_index]]
                    s, wgt_jac = self.uniform_distr(random_variable,s_min,s_max)
                    item['inv_mass'] = s
                wgt *= wgt_jac
                item['is_available'] = True
        
        # all invariant masses are generated, so for consistency set all available-flags to False
        # REMARK: there are some that say True but are not available i.e.
        #     1) if no t-channels and more than one s-channel -> last vertex final leg 1,2 say True
        #     2) if there's t-channels:
        #            one single: final leg 1,2 say True;
        #            more than one: either final leg 1 or 2 say True
        for item in uniform_inv:
            item['is_available'] = False
        for nr in list(kinematics.keys()):
            kinematics[nr]['is_available'] = False
        
        # sample the t-variables (uniform if massive, 1/x^nu if massless)
        # go from outermost to innermost t-channel
        # structure: leg1 + leg2 -> leg3 + leg4
        # if there are t-channels
        if len(self.topology[1]) > 1:
            uniform_inv.append({'vertices': (0,len(self.topology[1])-1), 'inv_mass': E_cm**2, 'is_available': False, 'momentum': None})
            uniform_inv.reverse()
            for k,item in enumerate(uniform_inv):
                start = item['vertices'][0]
                finish = item['vertices'][1]
                leg1_nr = self.topology[1][start].get('legs')[0].get('number')
                leg2_nr = self.topology[1][finish].get('legs')[-1].get('number')
                p1 = kinematics[leg1_nr]['momentum'].copy()
                p2 = kinematics[leg2_nr]['momentum'].copy()
                if leg1_nr < 0:
                    # t-channel propagator momenta always point upwards
                    # this code procedure requires p1 to point downwards
                    p1 = -p1.copy()
                # p is the sum of incoming momenta
                p = p1+p2
                if (start + 1) ==finish:
                    leg3_nr = self.topology[1][start].get('legs')[1].get('number')
                    leg4_nr = self.topology[1][finish].get('legs')[1].get('number')
                    if reconstruct:
                        p3 = kinematics[leg3_nr]['momentum']
                        p4 = kinematics[leg4_nr]['momentum']
                    p3_2 = kinematics[leg3_nr]['inv_mass']
                    p4_2 = kinematics[leg4_nr]['inv_mass']
                    prop_leg = self.topology[1][start].get('legs')[-1]
                else:
                    for inv in uniform_inv[k+1:]:
                        if inv['vertices'][0] == start:
                            leg3_nr = None
                            if reconstruct:
                                p3 = inv['momentum']
                            p3_2 = inv['inv_mass']
                            prop_leg = self.topology[1][inv['vertices'][1]].get('legs')[-1]
                            break
                        elif inv['vertices'][0] == start + 1:
                            leg3_nr = self.topology[1][start].get('legs')[1].get('number')
                            if reconstruct:
                                p3 = kinematics[leg3_nr]['momentum']
                            p3_2 = kinematics[leg3_nr]['inv_mass']
                            prop_leg = self.topology[1][inv['vertices'][0]].get('legs')[0]
                            break
                    for inv in uniform_inv[k+1:]:
                        if inv['vertices'][1] == finish:
                            leg4_nr = None
                            if reconstruct:
                                p4 = inv['momentum']
                            p4_2 = inv['inv_mass']
                            #not needed because already determined above
                            #prop_leg = self.topology[1][inv['vertices'][0]].get('legs')[0]
                            break
                        elif inv['vertices'][1] == finish - 1:
                            leg4_nr = self.topology[1][finish].get('legs')[1].get('number')
                            if reconstruct:
                                p4 = kinematics[leg4_nr]['momentum']
                            p4_2 = kinematics[leg4_nr]['inv_mass']
                            #not needed because already determined above
                            #prop_leg = self.topology[1][inv['vertices'][1]].get('legs')[-1]
                            break
                prop_id = prop_leg.get('id')
                s = p.square()
                assert(s >= 0)
                if not reconstruct:
                    # think about a better solution to handle if particle is produced at rest
                    no_phase_space = 0.*max([s,p3_2,p4_2])**2 >= self.Lambda(s,p3_2,p4_2) >= -tolerance**4*max([s,p3_2,p4_2])**2
                    if no_phase_space:
                        misc.sprint('This kinematic configuration leaves no phase space ', self.Lambda(s,p3_2,p4_2))
                        return None,0.,(xb_1, xi1) , (xb_2, xi2)
                t_max = self.inv_t(s, p1.square(), p2.square(), p3_2, p4_2, 1.)
                if abs(t_max) < tolerance**2*self.absolute_Ecm_min**2:
                    t_max = 0. 
                t_min = self.inv_t(s, p1.square(), p2.square(), p3_2, p4_2, -1.)
                assert(t_min<=t_max<=0.)
                particle = self.model.get_particle(prop_id)
                mass_param = particle.get('mass')
                prop_nr = prop_leg.get('number')
                assert(prop_nr<0)
                variable_index += 1
                if reconstruct:
                    kinematics[prop_nr]['momentum'] = p3-p1
                    t = kinematics[prop_nr]['momentum'].square()
                    if abs(t) < tolerance**2*self.absolute_Ecm_min**2:
                        t = t_max
                    if 0. < t_min-t < tolerance**2*t_min:
                        t = t_min
                    kinematics[prop_nr]['inv_mass'] = t
                    if mass_param.lower() == 'zero':
                        # t is always negative, distribute abs(t)
                        x_abs_t,wgt_jac = self.inv_massless_distr(abs(t), abs(t_max), abs(t_min))
                    else:
                        x_abs_t, wgt_jac = self.inv_uniform_distr(abs(t), abs(t_max), abs(t_min))
                    output_variables[self.dim_name_to_position['x_%d' % variable_index]] = x_abs_t
                else:
                    random_variable = random_variables[self.dim_name_to_position['x_%d' % variable_index]]
                    if mass_param.lower() == 'zero':
                        # t is always negative, distribute abs(t)
                        abs_t,wgt_jac = self.massless_distr(random_variable, abs(t_max), abs(t_min))
                        t = -abs_t
                    else:
                        abs_t, wgt_jac = self.uniform_distr(random_variable, abs(t_max), abs(t_min))
                        t = -abs_t
                wgt *= wgt_jac
                cos_theta = self.cos_theta_from_inv_t(s, p1.square(), p2.square(), p3_2, p4_2, t)
                variable_index += 1
                if reconstruct:
                    x_phi, wgt_jac = None, 2*math.pi
                    output_variables[self.dim_name_to_position['x_%d' % variable_index]] = x_phi
                else:
                    random_variable = random_variables[self.dim_name_to_position['x_%d' % variable_index]]
                    phi, wgt_jac = self.uniform_distr(random_variable, 0, 2*math.pi)
                wgt *= wgt_jac
                if reconstruct:
                    wgt_PS = self.get_two_body_PS_wgt(s, p3_2, p4_2)
                else:
                    p3, p4, wgt_PS = self.get_cm_momenta(s,p3_2,p4_2,cos_theta,phi)
                #t-channel only factor:
                wgt_PS *= 2.*s/(math.sqrt(self.Lambda(s,p3_2,p4_2))*math.sqrt(self.Lambda(s,p1.square(),p2.square())))
                wgt *= wgt_PS
                
                if not reconstruct:
                    p1_com = p1.copy()
                    p1_com.boost(-p.boostVector())
                    
                    # test if boost is correct
                    """
                    p2_com = p2.copy()
                    p2_com.boost(-p.boostVector())
                    assert(all(abs(x) < tolerance for x in (p1_com+p2_com)[1:]))
                    """
                    
                    phi, theta = self.angles_to_rotate_along_z(p1_com)
                    
                    # this is to test if rotation works as expected
                    """
                    phi1,theta1 = (math.pi/2,math.pi/6)
                    ptest = LorentzVector([10,math.cos(phi1)*math.sin(theta1),math.sin(phi1)*math.sin(theta1),math.cos(theta1)])
                    ptest_z = LorentzVector([10,0,0,1])
                    misc.sprint('\n',ptest)
                    phi2,theta2 = self.angles_to_rotate_along_z(ptest)
                    misc.sprint(phi1-phi2,theta1-theta2) #should be 0,0
                    misc.sprint('\n',self.rotate_along_z_inv(ptest_z,phi2, theta2)-ptest) #should be [0,0,0,0]
                    """
                    
                    p3 = self.rotate_along_z_inv(p3,phi, theta)
                    p4 = self.rotate_along_z_inv(p4,phi, theta)
                    assert(abs((p3-p1_com).square()-t) < tolerance**2*s) 
                    
                    p3.boost(p.boostVector())
                    p4.boost(p.boostVector())
                    
                    if leg3_nr != None:
                        assert(kinematics[leg3_nr]['momentum'] == None)
                        kinematics[leg3_nr]['momentum'] = p3
                    if leg4_nr != None:
                        assert(kinematics[leg4_nr]['momentum'] == None)
                        kinematics[leg4_nr]['momentum'] = p4
                        
                    assert(kinematics[prop_nr]['momentum'] == None)
                    # t-propagator-momentum always points upwards
                    kinematics[prop_nr]['momentum'] = p3-p1
                    assert(abs(kinematics[prop_nr]['momentum'].square()-t) < tolerance**2*s)
                    assert(all(abs(x) < tolerance*math.sqrt(s) for x in (p-p3-p4)))
        
        # s-channel momentum generation
        # start at the very end of topology[0] and go upwards, ordering doesn't matter
        # structure leg0 -> leg1 + leg2
        for vertex in reversed(self.topology[0]):
            leg0_nr = vertex.get('legs')[-1].get('number')
            leg1_nr = vertex.get('legs')[0].get('number')
            leg2_nr = vertex.get('legs')[1].get('number')
            s = kinematics[leg0_nr]['inv_mass']
            p = kinematics[leg0_nr]['momentum']
            p1_2 = kinematics[leg1_nr]['inv_mass']
            p2_2 = kinematics[leg2_nr]['inv_mass']

            variable_index += 1
            if reconstruct:
                x_cos_theta, wgt_jac = None, 2
                output_variables[self.dim_name_to_position['x_%d' % variable_index]] = x_cos_theta
            else:
                # this is for the case where a massless propagator is produced almost on shell, this leads
                # to a boost to velocity c which is numerically unstable
                # it's not an actual problem since it can only happen when massless splits into 2 massless,
                # e.g. e- > e- a, photon or gluon radiation. This would eventually be treated by pt_cuts anyways
                if s < tolerance**2*self.absolute_Ecm_min**2:
                    #misc.sprint(s)
                    return None,0.,(xb_1, xi1) , (xb_2, xi2)
                random_variable = random_variables[self.dim_name_to_position['x_%d' % variable_index]]
                cos_theta, wgt_jac = self.uniform_distr(random_variable, -1, 1)
            wgt *= wgt_jac
            variable_index += 1
            if reconstruct:
                x_phi, wgt_jac = None, 2*math.pi
                output_variables[self.dim_name_to_position['x_%d' % variable_index]] = x_phi
            else:
                random_variable = random_variables[self.dim_name_to_position['x_%d' % variable_index]]
                phi, wgt_jac = self.uniform_distr(random_variable, 0, 2*math.pi)
            wgt *= wgt_jac
            if reconstruct:
                wgt_PS = self.get_two_body_PS_wgt(s, p1_2, p2_2)
            else:
                p1,p2,wgt_PS = self.get_cm_momenta(s, p1_2, p2_2, cos_theta, phi)
            wgt *= wgt_PS
            if not reconstruct:
                p1.boost(p.boostVector())
                p2.boost(p.boostVector())
                assert(kinematics[leg1_nr]['momentum'] == None and kinematics[leg2_nr]['momentum'] == None)
                kinematics[leg1_nr]['momentum'] = p1
                kinematics[leg2_nr]['momentum'] = p2
                
                if not (all(abs(x) < tolerance*math.sqrt(s) for x in (p-p1-p2))):
                    logger.debug('Possible loss of precision with tolerance %.2e: %s'%(tolerance*math.sqrt(s), str(p-p1-p2)))
                #assert(all(abs(x) < tolerance*math.sqrt(s) for x in (p-p1-p2)))
        
        if not reconstruct:
            # sanity check if all random variables were actually used
            # Do no count the t scaling random variable
            offset = 1 if 't' in self.dim_name_to_position else 0
            if self.beam_types == (1,1):
                assert(variable_index +2 == len(random_variables)-offset)
            if self.beam_types == (0,0):
                assert(variable_index== len(random_variables)-offset)
        
        if not reconstruct:
            # save final momenta into output_momenta -> PS_point
            for nr in kinematics:
                if nr > 0:
                    output_momenta[nr-1] = kinematics[nr]['momentum']
                if nr == min_leg_nr:
                    output_momenta[1] = kinematics[nr]['momentum']
            PS_point = LorentzVectorList(output_momenta)
            
            # test if total momentum is conserved
            p_start = sum(PS_point[0:2])
            p_finish = sum(PS_point[2:])
            assert(all(abs(x) < tolerance*self.collider_energy for x in (p_start-p_finish)))
        
        assert(wgt>0)
        
        if reconstruct:
            return output_variables, wgt
        else:
            return LorentzVectorList(PS_point), wgt, (xb_1, xi1) , (xb_2, xi2)

class FlatInvertiblePhasespace(VirtualPhaseSpaceGenerator):
    """Implementation following S. Platzer, arxiv:1308.2922"""

    # This parameter defines a thin layer around the boundary of the unit hypercube
    # of the random variables generating the phase-space,
    # so as to avoid extrema which are an issue in most PS generators.
    epsilon_border = 1e-10

    # The lowest value that the center of mass energy can take.
    # We take here 1 GeV, as anyway below this non-perturbative effects dominate
    # and factorization does not make sense anymore
    absolute_Ecm_min = 1.

    # For reference here we put the flat weights that Simon uses in his
    # Herwig implementation. I will remove them once I will have understood
    # why they don't match the physical PS volume.
    # So these are not used for now, and get_flatWeights() is used instead.
    flatWeights =  { 2 :  0.039788735772973833942,
                     3 :  0.00012598255637968550463,
                     4 :  1.3296564302788840628e-7,
                     5 :  7.0167897579949011130e-11,
                     6 :  2.2217170114046130768e-14 
                   }

    def __init__(self, *args, **opts):
        super(FlatInvertiblePhasespace, self).__init__(*args, **opts)
        if self.n_initial == 1:
            raise InvalidCmd("This basic generator does not support decay topologies.")

    def get_dimensions(self):
        """ Make sure the collider setup is supported."""

        # Check if the beam configuration is supported
        if (not abs(self.beam_types[0])==abs(self.beam_types[1])==1) and \
           (not self.beam_types[0]==self.beam_types[1]==0):
            raise InvalidCmd(
                "This basic generator does not support the collider configuration: (lpp1=%d, lpp2=%d)"%
                             (self.run_card['lpp1'], self.run_card['lpp2']))
        
        if self.beam_Es[0]!=self.beam_Es[1]:
            raise InvalidCmd(
                "This basic generator only supports colliders with incoming beams equally energetic.")

        return super(FlatInvertiblePhasespace,self).get_dimensions()

    @staticmethod
    def get_flatWeights(E_cm, n, mass=None):
        """ Return the phase-space volume for a n massless final states.
        Vol(E_cm, n) = (pi/2)^(n-1) *  (E_cm^2)^(n-2) / ((n-1)!*(n-2)!)
        """
        if n==1: 
            # The jacobian from \delta(s_hat - m_final**2) present in 2->1 convolution
            # must typically be accounted for in the MC integration framework since we
            # don't have access to that here, so we just return 1.
            return 1.

        vol = math.pow((math.pi/2.0),n-1)*\
            (math.pow((E_cm**2),n-2)/(math.factorial(n-1)*math.factorial(n-2)))

        return vol

    @staticmethod
    def bisect(v, n, target=1.e-16, maxLevel=80):
        """Solve v = (n+2) * u^(n+1) - (n+1) * u^(n+2) for u."""
        
        if (v == 0. or v == 1.): return v

        level = 0
        left  = 0.
        right = 1.
            
        checkV = -1.
        u = -1.

        while (level < maxLevel):
            u = (left + right) * (0.5**(level + 1))
            checkV = (u**(n+1)) * (n+2.-(n+1.)*u)
            error = abs(1. - checkV / v)
            if (error == 0. or error <= target):
                break
            left *= 2.
            right *= 2.
            if (v <= checkV ): right -= 1.
            else: left += 1.
            level += 1

        return u
    
    @staticmethod
    def rho(M, N, m):
        """Returns sqrt((sqr(M)-sqr(N+m))*(sqr(M)-sqr(N-m)))/(8.*sqr(M))"""

        Msqr = M**2
        return ((Msqr-(N+m)**2) * (Msqr-(N-m)**2) )**0.5 / (8.*Msqr)

    def setInitialStateMomenta(self, output_momenta, E_cm):
        """Generate the initial state momenta."""

        if self.n_initial not in [1,2]:
            raise InvalidCmd(
               "This PS generator only supports 1 or 2 initial states")

        if self.n_initial == 1:
            if self.initial_masses[0]==0.:
                raise PhaseSpaceGeneratorError(
                    "Cannot generate the decay phase-space of a massless particle.")
            if self.E_cm != self.initial_masses[0]:
                raise PhaseSpaceGeneratorError(
                    "Can only generate the decay phase-space of a particle at rest.")

        if self.n_initial == 1:
            output_momenta[0] = LorentzVector([self.initial_masses[0] , 0., 0., 0.])
            return

        elif self.n_initial == 2:
            if self.initial_masses[0] == 0. or self.initial_masses[1] == 0.:
                output_momenta[0] = LorentzVector([E_cm/2.0 , 0., 0., +E_cm/2.0])
                output_momenta[1] = LorentzVector([E_cm/2.0 , 0., 0., -E_cm/2.0])
            else:
                M1sq = self.initial_masses[0]**2
                M2sq = self.initial_masses[1]**2
                E1 = (E_cm**2+M1sq-M2sq)/ E_cm
                E2 = (E_cm**2-M1sq+M2sq)/ E_cm
                Z = math.sqrt(E_cm**4 - 2*E_cm**2*M1sq - 2*E_cm**2*M2sq + M1sq**2 - 2*M1sq*M2sq + M2sq**2) / E_cm
                output_momenta[0] = LorentzVector([E1/2.0 , 0., 0., +Z/2.0])
                output_momenta[1] = LorentzVector([E2/2.0 , 0., 0., -Z/2.0])
        return

    def get_PS_point(self, random_variables, **opts):
        """Generate a complete PS point, including Bjorken x's,
        dictating a specific choice of incoming particle's momenta.
        """

        # if random_variables are not defined, than just throw a completely random point
        if random_variables is None:
            random_variables = self.dimensions.random_sample()
        
        # Check the sensitivity of te inputs from the integrator
        if any(math.isnan(r) for r in random_variables):
            logger.warning('Some input variables from the integrator are malformed: %s'%
                ( ', '.join( '%s=%s'%( name, random_variables[pos]) for name, pos in 
                                                     list(self.dim_name_to_position.items()) ) ))
            logger.warning('The PS generator will yield None, triggering the point to be skipped.')
            return None, 0.0, (0., 0.), (0., 0.)
        
        # Phase-space point weight to return
        wgt = 1.0
        
        #if any(math.isnan(r) for r in random_variables):
        #    misc.sprint(random_variables)
        
        # Avoid extrema since the phase-space generation algorithm doesn't like it
        random_variables = [min(max(rv,self.epsilon_border),1.-self.epsilon_border) for rv in random_variables]

        # Assign variables to their meaning.
        if 'ycms' in self.dim_name_to_position:
            PDF_ycm = random_variables[self.dim_name_to_position['ycms']]
        else:
            PDF_ycm = None
        if 'tau' in self.dim_name_to_position:
            PDF_tau = random_variables[self.dim_name_to_position['tau']]
        else:
            PDF_tau = None

        PS_random_variables  = [rv for i, rv in enumerate(random_variables) if self.position_to_dim_name[i].startswith('x_') ]

        # Also generate the ISR collinear factorization convolutoin variables xi<i> if
        # necessary. In order for the + distributions of the PDF counterterms and integrated
        # collinear ISR counterterms to hit the PDF only (and not the matrix elements or
        # observables functions), a change of variable is necessary: xb_1' = xb_1 * xi1
        if self.correlated_beam_convolution:
            # Both xi1 and xi2 must be set equal then
            xi1 = random_variables[self.dim_name_to_position['xi']]
            xi2 = random_variables[self.dim_name_to_position['xi']]
        else:
            if self.is_beam_factorization_active[0]:
                xi1 = random_variables[self.dim_name_to_position['xi1']]
            else:
                xi1 = None
            if self.is_beam_factorization_active[1]:
                xi2 = random_variables[self.dim_name_to_position['xi2']]
            else:
                xi2 = None

        # Now take care of the Phase-space generation:
        # Set some defaults for the variables to be set further
        xb_1 = 1.
        xb_2 = 1.
        E_cm = self.collider_energy
        
        # We generate the PDF from two variables \tau = x1*x2 and ycm = 1/2 * log(x1/x2), so that:
        #  x_1 = sqrt(tau) * exp(+ycm)
        #  x_2 = sqrt(tau) * exp(-ycm)
        # The jacobian of this transformation is 1.
        if abs(self.beam_types[0])==abs(self.beam_types[1])==1:
            
            tot_final_state_masses = sum(self.masses)
            if tot_final_state_masses > self.collider_energy:
                raise PhaseSpaceGeneratorError("Collider energy is not large enough, there is no phase-space left.")
            
            # Keep a hard cut at 1 GeV, which is the default for absolute_Ecm_min
            tau_min = (max(tot_final_state_masses, self.absolute_Ecm_min)/self.collider_energy)**2
            tau_max = 1.0

            if self.n_initial == 2 and self.n_final == 1:
                # Here tau is fixed by the \delta(xb_1*xb_2*s - m_h**2) which sets tau to 
                PDF_tau = tau_min
                # Account for the \delta(xb_1*xb_2*s - m_h**2) and corresponding y_cm matching to unit volume
                wgt *= (1./self.collider_energy**2)
            else:
                # Rescale tau appropriately
                PDF_tau = tau_min+(tau_max-tau_min)*PDF_tau
                # Including the corresponding Jacobian
                wgt *= (tau_max-tau_min)

            # And we can now rescale ycm appropriately
            ycm_min = 0.5 * math.log(PDF_tau)
            ycm_max = -ycm_min
            PDF_ycm = ycm_min + (ycm_max - ycm_min)*PDF_ycm            
            # and account for the corresponding Jacobian
            wgt *= (ycm_max - ycm_min)

            xb_1 = math.sqrt(PDF_tau) * math.exp(PDF_ycm)
            xb_2 = math.sqrt(PDF_tau) * math.exp(-PDF_ycm)
            # /!\ The mass of initial state momenta is neglected here.
            E_cm = math.sqrt(xb_1*xb_2)*self.collider_energy

        elif self.beam_types[0]==self.beam_types[1]==0:
            xb_1 = 1.
            xb_2 = 1.
            E_cm = self.collider_energy
        else:
            raise InvalidCmd("This basic PS generator does not yet support collider mode (%d,%d)."%self.beam_types)

        # Now generate a PS point
        PS_point, PS_weight = self.generateKinematics(E_cm, PS_random_variables)
        
        # Apply the phase-space weight
        wgt *= PS_weight
        
        return LorentzVectorList(PS_point), wgt, (xb_1, xi1) , (xb_2, xi2)

    def generateKinematics(self, E_cm, random_variables):
        """Generate a self.n_initial -> self.n_final phase-space point
        using the random variables passed in argument.
        """

        # Make sure the right number of random variables are passed
        assert (len(random_variables)==self.nDimPhaseSpace())

        # Make sure that none of the random_variables is NaN.
        if any(math.isnan(rv) for rv in random_variables):
            raise PhaseSpaceGeneratorError("Some of the random variables passed "+
              "to the phase-space generator are NaN: %s"%str(random_variables))

        # The distribution weight of the generate PS point
        weight = 1.
        
        output_momenta = []

        mass = self.masses[0]
        if self.n_final == 1:
            if self.n_initial == 1:
                raise InvalidCmd("1 > 1 phase-space generation not supported.")
            if mass/E_cm < 1.e-7 or ((E_cm-mass)/mass) > 1.e-7:
                raise PhaseSpaceGeneratorError("1 > 2 phase-space generation needs a final state mass equal to E_c.o.m.")
            output_momenta.append(LorentzVector([mass/2., 0., 0., +mass/2.]))
            output_momenta.append(LorentzVector([mass/2., 0., 0., -mass/2.]))
            output_momenta.append(LorentzVector([mass   , 0., 0.,       0.]))
            weight = self.get_flatWeights(E_cm, 1)
            return output_momenta, weight
  
        M    = [ 0. ]*(self.n_final-1)
        M[0] = E_cm

        weight *= self.generateIntermediatesMassive(M, E_cm, random_variables)
        M.append(self.masses[-1])

        Q     = LorentzVector([M[0], 0., 0., 0.])
        nextQ = LorentzVector()

        for i in range(self.n_initial+self.n_final-1):
            
            if i < self.n_initial:
                output_momenta.append(LorentzVector())
                continue

            q = 4.*M[i-self.n_initial]*self.rho(
                M[i-self.n_initial],M[i-self.n_initial+1],self.masses[i-self.n_initial] )
            cos_theta = 2.*random_variables[self.n_final-2+2*(i-self.n_initial)]-1.
            sin_theta = math.sqrt(1.-cos_theta**2)
            phi = 2.*math.pi*random_variables[self.n_final-1+2*(i-self.n_initial)]
            cos_phi = math.cos(phi)
            sin_phi = math.sqrt(1.-cos_phi**2)

            if (phi > math.pi):
                sin_phi = -sin_phi
            
            p = LorentzVector([0., q*sin_theta*cos_phi, q*sin_theta*sin_phi, q*cos_theta])
            p.set_square(self.masses[i-self.n_initial]**2)
            p.boost(Q.boostVector())
            p.set_square(self.masses[i-self.n_initial]**2)
            output_momenta.append(p)

            nextQ = Q - p
            nextQ.set_square(M[i-self.n_initial+1]**2)

            Q = nextQ
       
        output_momenta.append(Q)

        self.setInitialStateMomenta(output_momenta, E_cm)

        return LorentzVectorList(output_momenta), weight

    def generateIntermediatesMassless(self, M, E_cm, random_variables):
        """Generate intermediate masses for a massless final state."""
        
        for i in range(2, self.n_final):
            u = self.bisect(random_variables[i-2], self.n_final-1-i)
            M[i-1] = math.sqrt(u*(M[i-2]**2))

        return self.get_flatWeights(E_cm,self.n_final)
   

    def generateIntermediatesMassive(self, M, E_cm, random_variables):
        """Generate intermediate masses for a massive final state."""

        K = list(M)
        K[0] -= sum(self.masses)

        weight = self.generateIntermediatesMassless(K, E_cm, random_variables)
        del M[:]
        M.extend(K)
        
        for i in range(1,self.n_final):
            for k in range(i,self.n_final+1):
                M[i-1] += self.masses[k-1]
        
        weight *= 8.*self.rho(
            M[self.n_final-2],
            self.masses[self.n_final-1],
            self.masses[self.n_final-2] )

        for i in range(2,self.n_final):
            weight *= (self.rho(M[i-2],M[i-1],self.masses[i-2]) / self.rho(K[i-2],K[i-1],0.)) * (M[i-1]/K[i-1])

        weight *= math.pow(K[0]/M[0],2*self.n_final-4)

        return weight

    def invertKinematics(self, E_cm, momenta):
        """ Returns the random variables that yields the specified momenta configuration."""

        # Make sure the right number of momenta are passed
        assert (len(momenta) == (self.n_initial + self.n_final) )
        moms = momenta.get_copy()

        # The weight of the corresponding PS point
        weight = 1.

        if self.n_final == 1:
            if self.n_initial == 1:
                raise PhaseSpaceGeneratorError("1 > 1 phase-space generation not supported.")
            return [], self.get_flatWeights(E_cm,1) 

        # The random variables that would yield this PS point.
        random_variables = [-1.0]*self.nDimPhaseSpace()
        
        M    = [ 0. ]*(self.n_final-1)
        M[0] = E_cm

        Q     = [LorentzVector(4*[0., ])]*(self.n_final-1)
        Q[0]  = LorentzVector([M[0],0.,0.,0.])

        for i in range(2,self.n_final):
            for k in range(i, self.n_final+1):
                Q[i-1] = Q[i-1] + moms[k+self.n_initial-1]
            M[i-1] = abs(Q[i-1])

        weight = self.invertIntermediatesMassive(M, E_cm, random_variables)

        for i in range(self.n_initial,self.n_final+1):
            # BALDY another copy? moms not used afterwards
            p = LorentzVector(moms[i])
            # Take the opposite boost vector
            boost_vec = -Q[i-self.n_initial].boostVector()
            p.boost(boost_vec)
            random_variables[self.n_final-2+2*(i-self.n_initial)] = (p.cosTheta()+1.)/2.
            phi = p.phi()
            if (phi < 0.):
                phi += 2.*math.pi
            random_variables[self.n_final-1+2*(i-self.n_initial)] = phi / (2.*math.pi)
        
        return random_variables, weight

    def invertIntermediatesMassive(self, M, E_cm, random_variables):
        """ Invert intermediate masses for a massive final state."""

        K = list(M)
        for i in range(1, self.n_final):
            K[i-1] -= sum(self.masses[i-1:])
        
        weight = self.invertIntermediatesMassless(K, E_cm, random_variables)
        weight *= 8.*self.rho(M[self.n_final-2],
                              self.masses[self.n_final-1],
                              self.masses[self.n_final-2])
        for i in range(2, self.n_final):
            weight *= (self.rho(M[i-2],M[i-1],self.masses[i-2])/self.rho(K[i-2],K[i-1],0.)) \
                                                                      * (M[i-1]/K[i-1])
        
        weight *= math.pow(K[0]/M[0],2*self.n_final-4)

        return weight

    def invertIntermediatesMassless(self, K, E_cm, random_variables):
        """ Invert intermediate masses for a massless final state."""

        for i in range(2, self.n_final):
            u = (K[i-1]/K[i-2])**2
            random_variables[i-2] = \
                (self.n_final+1-i)*math.pow(u,self.n_final-i) - \
                (self.n_final-i)*math.pow(u, self.n_final+1-i)
        
        return self.get_flatWeights(E_cm, self.n_final)

#=========================================================================================
# Standalone main for debugging / standalone trials
#=========================================================================================
if __name__ == '__main__':

    import random

    E_cm  = 5000.0

    # Try to run the above for a 2->8.
    my_PS_generator = FlatInvertiblePhasespace([0.]*2, [100. + 10.*i for i in range(8)],
                                            beam_Es =(E_cm/2.,E_cm/2.), beam_types=(0,0))
    # Try to run the above for a 2->1.    
    #    my_PS_generator = FlatInvertiblePhasespace([0.]*2, [5000.0])
    
    random_variables = [random.random() for _ in range(my_PS_generator.nDimPhaseSpace())]

    momenta, wgt = my_PS_generator.generateKinematics(E_cm, random_variables)
   
    print("\n =========================")
    print(" ||    PS generation    ||")
    print(" =========================")

    print("\nRandom variables :\n",random_variables)
    print("\n%s\n"%momenta.__str__(n_initial=my_PS_generator.n_initial))
    print("Phase-space weight : %.16e\n"%wgt, end=' ')

    variables_reconstructed, wgt_reconstructed = \
                                         my_PS_generator.invertKinematics(E_cm, momenta)

    print("\n =========================")
    print(" || Kinematic inversion ||")
    print(" =========================")
    print("\nReconstructed random variables :\n",variables_reconstructed)
    differences = [
        abs(variables_reconstructed[i]-random_variables[i])
        for i in range(len(variables_reconstructed))
    ]
    print("Reconstructed weight = %.16e"%wgt_reconstructed)
    if differences:
        print("\nMax. relative diff. in reconstructed variables = %.3e"%\
            max(differences[i]/random_variables[i] for i in range(len(differences))))
    print("Rel. diff. in PS weight = %.3e\n"%((wgt_reconstructed-wgt)/wgt))


    print(('-'*100))
    print('SIDE EXPERIMENT, TO REMOVE LATER')
    print(('-'*100))
    import copy    

    my_PS_generator = FlatInvertiblePhasespace([0.]*2, [100. + 10.*i for i in range(2)],
                                            beam_Es =(E_cm/2.,E_cm/2.), beam_types=(0,0))
    # Try to run the above for a 2->1.    
    #    my_PS_generator = FlatInvertiblePhasespace([0.]*2, [5000.0])
    
    random_variables = [random.random() for _ in range(my_PS_generator.nDimPhaseSpace())]

    momenta, wgt = my_PS_generator.generateKinematics(E_cm, random_variables)


    print("original momenta")
    print((str(momenta)))
    print("computing boost vector to lab frame with x1=0.25 and x2=0.6")
    boost_vector_to_lab_frame = None
    xb_1, xb_2 = 0.25, 0.6
    ref_lab = (momenta[0]*xb_1 + momenta[1]*xb_2)
    boost_vector_to_lab_frame = -ref_lab.boostVector()

    print("then further rescaling with xi1=0.12 and xi2=0.71")
    xi1, xi2 = 0.12, 0.71
    momenta[0] = xi1*momenta[0]
    momenta[1] = xi2*momenta[1]

    print("Boosting back to c.o.m:")
    ref_com = (momenta[0] + momenta[1])
    xi_boost = ref_com.boostVector()
    boost_vector_to_lab_frame += xi_boost
    momenta_boosted_to_com = copy.deepcopy(momenta)
    for p in momenta_boosted_to_com:
        p.boost(-xi_boost)
    test_momenta = copy.deepcopy(momenta_boosted_to_com)

    print("Final PS point should be in c.o.m:")
    print((str(momenta_boosted_to_com)))
    print("\nFinal check\n")
    print("Initial state momenta obtained from simple rescaling:")
    print(("1: %s"%str(momenta_boosted_to_com[0]*xb_1*xi1)))
    print(("2: %s"%str(momenta_boosted_to_com[1]*xb_2*xi2)))

    print("\n Initial state momenta obtained with xi boost:")
    print(("xi_boost=",str(xi_boost)))
    test_momenta[0].boost(xi_boost)
    test_momenta[1].boost(xi_boost)
    print(("1: %s"%str(test_momenta[0])))
    print(("2: %s"%str(test_momenta[1])))
    print(".vs.")
    print(("1: %s"%str(momenta[0])))
    print(("2: %s"%str(momenta[1])))

    print("\n Initial state momenta obtained with subsequent bjorken boost:")
    test_momenta[0].boost(-ref_lab.boostVector())
    test_momenta[1].boost(-ref_lab.boostVector())
    print(("1: %s"%str(test_momenta[0])))
    print(("2: %s"%str(test_momenta[1])))
    print(".vs.")
    print(("1: %s"%str(momenta[0]*xb_1)))
    print(("2: %s"%str(momenta[1]*xb_2)))

    print("\n Initial state momenta obtained with overall boost:")
    momenta[0].boost(boost_vector_to_lab_frame)
    momenta[1].boost(boost_vector_to_lab_frame)
    print(("1: %s"%str(momenta[0])))
    print(("2: %s"%str(momenta[1])))


