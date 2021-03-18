#####################################################
#                                                   #
#  Source file of the alphaLoop MG5aMC plugin.      #
#                                                   #
#####################################################

import os
import logging
import sys
import re
import random
import sympy
import math
import timeit
import functools
import copy
import resource
import traceback
import shutil
import math
import time
import numpy as np
import itertools
from argparse import ArgumentParser
from pprint import pprint, pformat
import progressbar
import glob
import threading
import networkx as nx

#import matplotlib.pyplot as plt
#from matplotlib.font_manager import FontProperties
mpl_logger = logging.getLogger("matplotlib")
mpl_logger.setLevel(logging.WARNING)

from distutils.version import LooseVersion, StrictVersion

plugin_path = os.path.dirname(os.path.realpath( __file__ ))
import madgraph

import multiprocessing

from madgraph import InvalidCmd, MadGraph5Error, MG5DIR, ReadWrite
import madgraph.interface.extended_cmd as cmd
import madgraph.interface.madgraph_interface as madgraph_interface
import madgraph.various.misc as misc
import madgraph.various.cluster as cluster
import madgraph.core.color_algebra as color
import madgraph.core.base_objects as base_objects
import models.model_reader as model_reader

import alpha_loop.utils as utils
from alpha_loop.ltd_commons import HyperParameters
from alpha_loop.ltd_commons import hyperparameters as default_hyperparameters
from LTD.vectors import Vector, LorentzVector, LorentzVectorList
import LTD.ltd_utils as ltd_utils
import alpha_loop.integrator.sampler as sampler
import alpha_loop.integrator.integrands as integrands
import alpha_loop.integrator.integrators as integrators
import alpha_loop.integrator.vegas3_integrator as vegas3_integrator
import alpha_loop.integrator.pyCubaIntegrator as pyCubaIntegrator

Colours = utils.bcolors

from madgraph.iolibs.files import cp, ln, mv

logger = logging.getLogger('alphaLoop.Interface')

pjoin = os.path.join
template_dir = pjoin(plugin_path, 'Templates')

FINAL=True
INITIAL=False

DUMMY=99

# # The class below is a trick to side-step the difficulty of multiprocessing pools in python 3.9+ to work properly
# # with attributes that are nested functions or a rust worker instance for example.
# # The wrapper below will store them in the global dictionary below thus avoiding the issue.
# _CALLABLE_INSTANCES_POOL = {}
# class CallableInstanceWrapper(object):
#     def __init__(self, callable_instance):
#         if len(_CALLABLE_INSTANCES_POOL)==0:
#             self.instance_ID = 1
#         else:
#             self.instance_ID = max(_CALLABLE_INSTANCES_POOL.keys())+1
#         _CALLABLE_INSTANCES_POOL[self.instance_ID] = callable_instance

#     def __call__(self, *args, **opts):
#         return _CALLABLE_INSTANCES_POOL[self.instance_ID](*args, **opts)

#     def __getattr__(self, name):
#         try:
#             return getattr(_CALLABLE_INSTANCES_POOL[self.instance_ID],name)
#         except Exception as e:
#             return getattr(self,name) 
#             logger.critical("Faced exception %s when attempting to access attribute '%s' from instance of type '%s'."%(
#                 str(e), name, type(_CALLABLE_INSTANCES_POOL[self.instance_ID])
#             ))
#     def __exit__(self, exc_type, exc_val, exc_tb):
#         _CALLABLE_INSTANCES_POOL.pop(self.instance_ID)

# Simply trade the version above for the function below in order to remove the above "hack"
def CallableInstanceWrapper(instance):
    return instance

class alphaLoopRunInterfaceError(MadGraph5Error):
    """ Error for the alphaLoop plugin """
    pass

class alphaLoopInvalidRunCmd(InvalidCmd):
    """ Invalid command issued to the alphaLoop interface. """
    pass

try:
    import yaml
    from yaml import Loader, Dumper
except ImportError:
    raise alphaLoopRunInterfaceError("Install yaml python module in order to import topologies from yaml.")

class RunHyperparameters(HyperParameters):
    """ Subclassing of default parameters so as to adjust them for the context of this run interface."""

    def __init__(self, dict_update=None, dir_path=None):
        self.update(default_hyperparameters, allow_undefined=True)

        # Add here default modifications tailored to the alphaLoop run interface
        alphaloop_run_interface_specific_params = {
            'Integrator.integrated_phase'                   : 'real',
            'General.deformation_strategy'                  : 'none',
            'General.debug'                                 : 0,
            'Integrator.state_filename_prefix'              : '',
            'Integrator.keep_state_file'                    : False,
            'Integrator.load_from_state_file'               : False,
            'Integrator.n_max'                              : int(1.0e9),
            'Integrator.n_new'                              : int(1.0e5),
            'Integrator.n_start'                            : int(1.0e6),
            'Integrator.n_increase'                         : int(1.0e5),
            'Selectors.active_selectors'                    : [],
            'Selectors.jet.min_jets'                        : 0,
            'Selectors.jet.dR'                              : 0.4,
            'Selectors.jet.min_jpt'                         : 10.0,
            'General.multi_channeling'                      : True,
            'General.multi_channeling_including_massive_propagators' : True,
            'CrossSection.incoming_momenta'                 : [[125.0, 0.0, 0.0 ,0.0],],
            'CrossSection.m_uv_sq'                          : 155.0**2,
            'CrossSection.mu_r_sq'                          : 155.0**2,
            'CrossSection.gs'                               : 1.2177157847767195,
            'CrossSection.NormalisingFunction.name'         : 'left_right_exponential',
            'CrossSection.NormalisingFunction.center'       : 1.0,
            'CrossSection.NormalisingFunction.spread'       : 1.0,
            'General.stability_checks'                      : [
                {
                    # number of samples to take for the numerical stability check
                    'n_samples': 3,
                    'use_f128': False,
                    'use_pf': True,
                    # number of digits that should be the same between rotated versions
                    'relative_precision': 4.0,
                    # force an upgrade when a new weight is this threshold times the current maximum weight
                    'escalate_for_large_weight_threshold': 0.8,
                    'minimal_precision_to_skip_further_checks': 99.0
                },
                {
                    'n_samples': 3,
                    'use_f128': True,
                    'use_pf': True,
                    'relative_precision': 8.0,
                    'escalate_for_large_weight_threshold': -1.,
                    'minimal_precision_to_skip_further_checks': 99.0
                }
            ],
            # Can be yaml, FORM, FORM_integrand
            'CrossSection.numerator_source'                      :   'FORM_integrand',
            # Can be LTD, PF
            'CrossSection.integrand_type'                        :   'PF',
            'CrossSection.picobarns'                             :   False,
            # evaluate the C expression for the sum of diagrams
            'CrossSection.sum_diagram_sets'                      :   False,
            # compare locally against the same topology written in another loop momentum basis
            'CrossSection.compare_with_additional_topologies'    :   False,
            'CrossSection.inherit_deformation_for_uv_counterterm':   True,
        }

        for param, value in alphaloop_run_interface_specific_params.items():
            self.set_parameter(param, value)

        if dict_update is not None:
            self.update(dict_update)

class CrossSectionSet(dict):
    """ Container of supergraphs. """

    def __init__(self, record):

        if isinstance(record, str):
            if '\n' in record:
                flat_record = yaml.load(record, Loader=Loader)
            else:
                flat_record = yaml.load(open(record,'r'), Loader=Loader)
        elif isinstance(record, dict):
            flat_record = record
        else:
            raise alphaLoopRunInterfaceError(
                "A CrossSectionSet instance can only be created from a path to a "
                "yaml record or a dictionary containing the data.")

        self.update(flat_record)
        if not self.check_validity():
            if isinstance(record, str) and '\n' not in record:
                raise alphaLoopInvalidRunCmd("Cross section set specified at '%s' appears invalid."%record)
            else:
                raise alphaLoopInvalidRunCmd("Cross section set specified from data directly appears invalid.")

    def check_validity(self):
        
        for attr in ['name','topologies']:
            if attr not in self:
                return False

        if not isinstance(self['topologies'], (list, tuple)):
            return False

        for SG in self['topologies']:
            for attr in ['multiplicity','name']:
                if attr not in SG:
                    return False

        return True

    def __str__(self):
        return pformat(self)

    def summary_str(self):
        res = []
        res.append(('Number of supergraphs','%d'%len(self['topologies'])))
        sorted_mult_SG = sorted(( (SG['multiplicity'],SG['name']) for SG in self['topologies']) ,key=lambda el:el[0])
        res.append(('Min multiplicity','%d (%s)'%(sorted_mult_SG[0][0],sorted_mult_SG[0][1])))
        res.append(('Max multiplicity','%d (%s)'%(sorted_mult_SG[-1][0],sorted_mult_SG[-1][1])))

        res_str = []
        for k, v in res:
            res_str.append('%s%-40s%s : %s'%(Colours.GREEN,k,Colours.END,v))
        return '\n'.join(res_str)

class SuperGraph(dict):
    """ Container for storing yaml dump of each SG."""
    
    def __init__(self, record):

        if isinstance(record, str):
            if '\n' in record:
                flat_record = yaml.load(record, Loader=Loader)
            else:
                flat_record = yaml.load(open(record,'r'), Loader=Loader)
        elif isinstance(record, dict):
            flat_record = record
        else:
            raise alphaLoopRunInterfaceError(
                "A SuperGraph instance can only be created from a path to a "
                "yaml record or a dictionary containing the data.")

        self.update(flat_record)
        if not self.check_validity():
            if isinstance(record, str) and '\n' not in record:
                raise alphaLoopInvalidRunCmd("Supergraph specified at '%s' appears invalid."%record)
            else:
                raise alphaLoopInvalidRunCmd("Supergraph specified from data directly appears invalid.")
    
    def get_E_cm(self, global_hyperparameters=None):
        """ Determines the relevant E_cm for this topology."""

        # First get the Q^2 collision energy
        if 'default_kinematics' in self:
            incoming_momenta_summed = sum(LorentzVector(k) for k in self['default_kinematics'][0])
        else:
            if global_hyperparameters is None:
                raise alphaLoopRunInterfaceError("When no default kinematics is specified in a supergraph, hyperparameters must be provided.")
            incoming_momenta_summed = sum(LorentzVector(k) for k in global_hyperparameters['CrossSection']['incoming_momenta'])
        return math.sqrt(abs(incoming_momenta_summed.square()))

    def get_random_x_input(self):
        """ Generates a random input with the right overall scale."""
        return [ random.random() for _ in range(3*self['topo']['n_loops']) ]

    def get_random_momenta_input(self, E_cm):
        """ Generates a random input with the right overall scale."""
        return [
                [ E_cm*(r/(1.-r)) for r in [random.random() for _ in range(3)] ]
            for _ in range(self['topo']['n_loops'])
        ]

    def check_validity(self):
        return True

    def __str__(self):
        return pformat(self)

    def get_cut_characteristics(self, cut):

        n_CC = len(cut['cuts'])
        n_diag_sets = len(cut['diagram_sets'])

        n_loops_left = 0
        n_loops_right = 0
        for diag_set in cut['diagram_sets']:
            sum_n_loops_left = 0
            sum_n_loops_right = 0
            for diag_info in diag_set['diagram_info']:
                if not diag_info['conjugate_deformation']:
                    sum_n_loops_left += diag_info['graph']['n_loops']
                else:
                    sum_n_loops_right += diag_info['graph']['n_loops']
            n_loops_left = max(n_loops_left,sum_n_loops_left)
            n_loops_right = max(n_loops_right,sum_n_loops_right)

        return {
            'n_loops_left'  : n_loops_left,
            'n_loops_right' : n_loops_right,
            'n_CC'          : n_CC,
            'n_diag_sets'   : n_diag_sets,
        }

    def summary_str(self):
        res = []
        res.append(('Total number of cuts','%d'%len(self['cutkosky_cuts'])))
        res.append(('Contributions distribution:',''))
        distribution = {}
        for i_cut, cut in enumerate(self['cutkosky_cuts']):
            cut_info = self.get_cut_characteristics(cut)
            signature = (   cut_info['n_loops_left']+cut_info['n_loops_right'], cut_info['n_loops_left'], 
                            cut_info['n_CC'],cut_info['n_loops_right']  )
            if signature not in distribution:
                distribution[signature] = {'n_diag_sets': cut_info['n_diag_sets'], 'n_cuts': 1, 'cuts': [i_cut]}
            else:
                distribution[signature]['n_diag_sets'] += cut_info['n_diag_sets']
                distribution[signature]['n_cuts'] += 1
                distribution[signature]['cuts'].append(i_cut)

        sorted_keys = sorted(distribution.keys())
        for n_loops_tot, n_loops_left, n_CC, n_loops_right in sorted_keys:
            entry = distribution[(n_loops_tot, n_loops_left, n_CC, n_loops_right)]
            res.append(('( n_loops_tot=%d, n_loops_left=%d, n_CC=%d, n_loops_right=%d )'%(
                n_loops_tot, n_loops_left, n_CC, n_loops_right
            ),'n_cuts=%d n_diag_sets=%d cuts=%s'%(entry['n_cuts'],entry['n_diag_sets'],str(entry['cuts']))))

        res_str = []
        for k, v in res:
            if v!='':
                res_str.append('%s%-80s%s : %s'%(Colours.GREEN,k,Colours.END,v))
            else:
                res_str.append('%s%s%s'%(Colours.BLUE,k,Colours.END))
        return '\n'.join(res_str)

    def show_timing_statistics(self):
        if ('DERIVED_timing_profile_f64' not in self) and ('DERIVED_timing_profile_f128' not in self):
            return "No timing profile information available."

        res = []
        for precision in ['f64','f128']:
            if 'DERIVED_timing_profile_%s'%precision not in self:
                continue
            res.append(('Statistics for %s precision:'%precision,''))
            res.append(('-'*len(res[-1][0]),''))
            res.append(('Timing for overall evaluate','%.3f ms'%(self['DERIVED_timing_profile_%s'%precision])))
            sorted_cut_evaluations = sorted(
                [(i,c['DERIVED_timing_profile_%s'%precision]) for i, c in enumerate(self['cutkosky_cuts'])],
                key = lambda el: el[1]
            )
            res.append(('Summed timing overall evaluate cuts','%.3f ms'%sum([c[1] for c in sorted_cut_evaluations])))
            res.append(("Slowest cut evaluation", '%.3f ms (cut #%d)'%( sorted_cut_evaluations[-1][1],sorted_cut_evaluations[-1][0] )))
            res.append(("Fastest cut evaluation", '%.3f ms (cut #%d)'%( sorted_cut_evaluations[0][1],sorted_cut_evaluations[0][0] )))
            for i_cut, t in sorted_cut_evaluations:
                cut = self['cutkosky_cuts'][i_cut]
                cut_info = {'t': t,'COLOUR':Colours.GREEN, 'END': Colours.END}
                cut_info.update(self.get_cut_characteristics(cut))
                res.append(("Timing for cut #%d"%i_cut, 
                    ('%(t).3f ms (n_loops_left=%(COLOUR)s%(n_loops_left)d%(END)s, n_CC=%(COLOUR)s%(n_CC)d%(END)s, n_loops_right=%(COLOUR)s%(n_loops_right)d%(END)s,'+
                    ' n_diag_sets=%(COLOUR)s%(n_diag_sets)d%(END)s)')%cut_info))

        res_str = []
        for k, v in res:
            if v!='':
                res_str.append('%s%-40s%s : %s'%(Colours.GREEN,k,Colours.END,v))
            else:
                res_str.append('%s%s%s'%(Colours.BLUE,k,Colours.END))
        return '\n'.join(res_str)

    def show_IR_statistics(self, show_momenta=True):

        # TODO improve rendering
        res_str = []

        res_str.append("All %d E-surfaces considered in tests:"%(len(self['E_surfaces_analysis'])))
        for i_surf, E_surface in enumerate(self['E_surfaces_analysis']):
            res_str.append("#%-3d: %s%d%s-loop E-surface %s : %s"%(E_surface['id'], Colours.GREEN,E_surface['n_loops'],Colours.END, '(pinched)' if  E_surface['pinched'] else ' '*9,
                ', '.join('%-21s'%('%s%s%-4s%s'%('%s%s%s'%(Colours.GREEN,'+',Colours.END) if term['energy_shift_sign']> 0 else '%s%s%s'%(Colours.RED,'-',Colours.END),  
                    Colours.BLUE, term['name'], Colours.END) ) 
                for term in  E_surface['onshell_propagators'])))
        
        E_surface_ID_to_E_surface = { E_surf['id'] : E_surf for E_surf in self['E_surfaces_analysis'] }

        for inter_length, all_results in self['E_surfaces_intersection_analysis'].items():
            res_str.append("All %s%d%s tests performed for the intersection of %s%d E-surfaces%s:"%(
                Colours.BLUE, len(all_results), Colours.END,
                Colours.GREEN, inter_length, Colours.END
            ))
            for E_surface_combination, results in all_results.items():
                res_str.append('Intersection of E-surfaces %s%s%s'%(
                    Colours.BLUE,
                    '^'.join('dE(%s)'%(','.join(osp['name'] for osp in E_surface_ID_to_E_surface[E_surf_ID]['onshell_propagators'])) for E_surf_ID in E_surface_combination),
                    Colours.END
                ))
                if show_momenta:
                    intersection_point = [ Vector(v) for v in results['intersection_point'] ]
                    res_str.append( '\n'.join('\n'.join('   %s%-5s%s : %-20s'%(Colours.BLUE, osp['name'], Colours.END, 
                            ', '.join( '%s%.10e'%('+' if vi>=0. else '', vi) for vi in list(sum( l*factor for l, factor in zip(intersection_point,osp['loop_sig']) )+Vector(osp['v_shift'])) )
                        ) for osp in E_surface_ID_to_E_surface[E_surf_ID]['onshell_propagators']
                    ) for E_surf_ID in E_surface_combination) )

                res_list = [('Complete integrand','%sPASS%s'%(Colours.GREEN, Colours.END) if results['dod_computed'][-1] else '%sFAIL%s'%(Colours.RED, Colours.END), {k:v for k,v in results.items() if k!='cut_results'}),]+\
                    [(' > Cut #%-2d (%s)'%(cut_ID, 
                    ','.join(c['name'] for c in self['cutkosky_cuts'][cut_ID]['cuts'])
                    ), '%sPASS%s'%(Colours.BLUE, Colours.END) if cut_res['dod_computed'][-1] else '%sFAIL%s'%(Colours.BLUE, Colours.END), cut_res) for cut_ID, cut_res in  sorted(list(results.get('cut_results',{}).items()),key = lambda k: k[0]) ]
                for lead, middle, info in res_list:
                    res_str.append('   %-35s: %s %s'%(
                        lead, middle, '%-7.4f +/- %-7.4f (max for lambda = %.2e: %s)'%(info['dod_computed'][0],info['dod_computed'][1],info['max_result'][0],str(complex(*info['max_result'][1])) ) 
                    ))
        return '\n'.join(res_str)

    def show_UV_statistics(self):
        if ('DERIVED_UV_dod' not in self):
            return "No UV profile information available."

        res = []
        if len(self['DERIVED_UV_dod'])>0 and 'DERIVED_dod_consistency' in self:
            res.append(("Consistency of dod between eval_integrand and eval_cuts", "%.4f (%.4f vs %.4f)"%self['DERIVED_dod_consistency']))
        for k, v in self['DERIVED_UV_dod'].items():
            res.append(('Sum over cuts dod for LMB=%s, UV_indices=%s'%(
                ','.join(k[0]),str(k[1])
            ),'%-7.4f +/- %-7.4f %s'%(
                v[0],v[1], '%sPASS%s'%(Colours.GREEN, Colours.END) if v[-1] else '%sFAIL%s'%(Colours.RED, Colours.END)
            )))
        sorted_cut_evaluations = sorted(
            [(i,k,v) for i, c in enumerate(self['cutkosky_cuts']) for k,v in c['DERIVED_UV_dod'].items()],
            key = lambda el: el[2][0], reverse=True
        )
        for cut_ID, k, v in sorted_cut_evaluations:
            res.append(('Cut dod for cut_ID=%d, LMB=%s, UV_indices=%s'%(
                cut_ID, ','.join(k[0]),str(k[1])
            ),'%-7.4f +/- %-7.4f %s'%(
                v[0],v[1], '%sPASS%s'%(Colours.GREEN, Colours.END) if v[-1] else '%sFAIL%s'%(Colours.RED, Colours.END)
            )))

        res_str = []
        max_len = (max(len(k) for k, v in res if v!='') if len(res)>0 else 45)+5
        for k, v in res:
            if v!='':
                res_str.append(('%s%-{:d}s%s : %s'.format(max_len))%(Colours.GREEN,k,Colours.END,v))
            else:
                res_str.append('%s%s%s'%(Colours.BLUE,k,Colours.END))
        return '\n'.join(res_str)

    def set_integration_channels(self):
        """ This function shrinks all loops within each left- and right- amplitude graph of each cutkosky cut and builds the corresponding
        tree topology associated with it."""

        edge_powers = { e[0] : e[-1] for e in self['topo_edges'] }
        
        edges_PDG = { e[0] : e [1] for e in self['edge_PDGs'] }

        for i_cut, cutkosky_cut in enumerate(self['cutkosky_cuts']):

            # Obtain the (non-UV) diagram sets composing the amplitudes to the left and right of this cutkosky cut.
            left_graph_loop_lines = []
            right_graph_loop_lines = []

            for diagram_set in cutkosky_cut['diagram_sets']:
                found_diag_set = False
                for diag_info in diagram_set['diagram_info']:
                    # Ignore UV diagram sets
                    if all(not prop['uv'] for ll in diag_info['graph']['loop_lines'] for prop in ll['propagators']):
                        found_diag_set = True
                        # Remove tree propagator loop lines
                        filtered_loop_lines = [ll for ll in diag_info['graph']['loop_lines'] if not all(lle==0 for lle in ll['signature']) ]
                        if not diag_info['conjugate_deformation']:
                            left_graph_loop_lines.extend( filtered_loop_lines )
                        else:
                            right_graph_loop_lines.extend( filtered_loop_lines )
                    else:
                        break
                if found_diag_set:
                    break
            
            # We must now build disjoint sets of these propagator to identify which ones need to be shrunk together
            # First assign an independent ID to each group and merge them when we find common loop lines 
            left_graph_loop_lines = [[ [loop_line,], i_group] for i_group, loop_line in enumerate(left_graph_loop_lines)]
            right_graph_loop_lines = [[ [loop_line,], i_group] for i_group, loop_line in enumerate(right_graph_loop_lines)]
            # Now iteratively group these propagators together whenever they have a loop in common
            # (using the advanced algo for disjoint set computations using linked lists would be overkill here)
            for all_loop_lines in [left_graph_loop_lines, right_graph_loop_lines]:
                while True:
                    for i_group_A, groupA in enumerate(all_loop_lines):
                        for groupB in all_loop_lines[i_group_A:]:
                            if any( any(llA['signature'][i_signature]!=0 for llA in groupA[0]) for 
                                    llB in groupB[0] for i_signature, signature_int in enumerate(llB['signature']) if signature_int!=0 ):
                                # These two groups belong to the same loop that needs to be shrunk at once
                                groupA[1] = min(groupA[1],groupB[1])
                                groupB[1] = min(groupA[1],groupB[1])
                    
                    # Now collect and fuse all identical groups
                    new_all_loop_lines =[]
                    for group_ID in sorted(list(set(group[1] for group in all_loop_lines))):
                        new_all_loop_lines.append(
                            [ sum([ll_group[0] for ll_group in all_loop_lines if ll_group[1]==group_ID],[]), group_ID]
                        )
                    if len(new_all_loop_lines)==len(all_loop_lines):
                        break
                    else:
                        all_loop_lines[:] = new_all_loop_lines                        

            # We know only care about the propagator names that have been grouped together on either side of the cutkosky cut
            # So let's retain this information only
            left_graph_loop_propagators = { i_group : 
                    sum([ [ prop['name'] for prop in loop_line['propagators'] ] for loop_line in loop_lines],[]) 
                for loop_lines, i_group in left_graph_loop_lines
            }
            right_graph_loop_propagators = { i_group : 
                    sum([ [ prop['name'] for prop in loop_line['propagators'] ] for loop_line in loop_lines],[]) 
                for loop_lines, i_group in right_graph_loop_lines
            }

            # We can now use these propagators to shrink the corresponding loops in the supergraph
            non_shrunk_edges_for_this_CC_cut = {e[0]:tuple(e[1:]) for e in self['topo_edges']}
            effective_vertices = {}
            for i_side, loop_propagators_groups in enumerate([left_graph_loop_propagators, right_graph_loop_propagators]):
                for group_ID, loop_propagators in loop_propagators_groups.items():
                    non_shrunk_edges_for_this_CC_cut, subgraph_info = self.shrink_edges(non_shrunk_edges_for_this_CC_cut, loop_propagators)
                    subgraph_info['side_of_cutkosky_cut'] = 'left' if i_side==0 else 'right'
                    effective_vertices[subgraph_info.pop('effective_node')] = subgraph_info
            
            # We must now list all LMBs for the shrunk effective nodes since we will use this information to multichannel over them
            for effective_node, effective_vertex_info in effective_vertices.items():
                all_internal_edges = list(effective_vertex_info['internal_edges'].items())
                SG_topo_gen = ltd_utils.TopologyGenerator(
                    [ [e_name,]+list(e_info[:-1]) for e_name, e_info in all_internal_edges ],
                    powers = { e_name : e_info[-1] for e_name, e_info in all_internal_edges }
                )
                all_LMBs_for_this_effective_vertex = SG_topo_gen.loop_momentum_bases()
                effective_vertex_info['loop_momenta_bases'] = [ 
                    [ all_internal_edges[e_index][0] for e_index in one_lmb] 
                        for one_lmb in all_LMBs_for_this_effective_vertex 
                ]

            cut_edge_names = [c['name'] for c in cutkosky_cut['cuts']]
            external_edges = {'left':[], 'right':[]}
            for e_name, e_sig in self['edge_signatures'].items():
                if not all(s==0 for s in e_sig[0]):
                    continue
                left_external_e_sig = e_sig[1][:len(e_sig[1])//2] 
                right_external_e_sig = e_sig[1][len(e_sig[1])//2:] 
                if all(s==0 for s in right_external_e_sig) and left_external_e_sig.count(1)==1 and left_external_e_sig.count(0)==len(left_external_e_sig)-1:
                    external_edges['left'].append((e_name, left_external_e_sig.index(1)))
                if all(s==0 for s in left_external_e_sig) and right_external_e_sig.count(1)==1 and right_external_e_sig.count(0)==len(right_external_e_sig)-1:
                    external_edges['right'].append((e_name, right_external_e_sig.index(1)))
            # Now reorder the vertices in the order of the signatures
            external_edges['left'] = [ edge[0] for edge in sorted(external_edges['left'], key=lambda e: e[1]) ]
            external_edges['right'] = [ edge[0] for edge in sorted(external_edges['right'], key=lambda e: e[1]) ]

            # Make sure there is at least one in each
            assert(len(external_edges['left'])>=1)
            assert(len(external_edges['right'])>=1)

            cut_tree = ltd_utils.TopologyGenerator(
                [ (e_name, edge[0],edge[1]) for e_name, edge in non_shrunk_edges_for_this_CC_cut.items() if e_name not in cut_edge_names],
                powers = { e_name: edge_powers[e_name] for e_name in non_shrunk_edges_for_this_CC_cut if e_name not in cut_edge_names }
            )
            tree_topologies = { 'left':{}, 'right':{} }
            for side in ['left', 'right']:
                sub_tree_indices = []
                cut_tree.generate_spanning_trees(sub_tree_indices, tree={non_shrunk_edges_for_this_CC_cut[external_edges[side][0]][0]})
                tree_topologies[side]['edges'] = {
                    cut_tree.edge_map_lin[i][0] : non_shrunk_edges_for_this_CC_cut[cut_tree.edge_map_lin[i][0]]
                    for i in sub_tree_indices[0] if cut_tree.edge_map_lin[i][0] not in external_edges[side] }

                #internal_edge_nodes = list(set(sum( [ [edge[0], edge[1]] for e_name, edge in tree_topologies[side]['edges'].items() ],[])))

                # Add incoming edges to the subtree
                for edge_name in external_edges[side]:
                    tree_topologies[side]['edges'][edge_name] = non_shrunk_edges_for_this_CC_cut[edge_name]
                # We want the nodes in subtree to *not* include the nodes from the final state cuts, so we compute them here.
                nodes_in_subtree = list(set(sum( [ [edge[0], edge[1]] for e_name, edge in tree_topologies[side]['edges'].items() ],[])))
                for edge_name in cut_edge_names:
                    # Lower the power of the cut propagators as they have been cut
                    edge_info = list(non_shrunk_edges_for_this_CC_cut[edge_name])
                    edge_info[-1] -= 1
                    tree_topologies[side]['edges'][edge_name] = tuple(edge_info)

                internal_edge_nodes = sum( [ [edge[0], edge[1]] for e_name, edge in tree_topologies[side]['edges'].items() ],[])
                # Remove all nodes appearing once only and outside of the original list of nodes appearing in the subtree of that side
                internal_edge_nodes = [e for e in internal_edge_nodes if internal_edge_nodes.count(e)>1 and e in nodes_in_subtree]
                # Make each node occurrence unique
                internal_edge_nodes = list(set(internal_edge_nodes))

                tree_topologies[side]['effective_vertices_contained'] = [node for node in internal_edge_nodes if node in effective_vertices]
                tree_topologies[side]['incoming_edges'] = [
                    (e_name, 1 if non_shrunk_edges_for_this_CC_cut[e_name][1] in internal_edge_nodes else -1) for e_name in external_edges[side]
                ]
                tree_topologies[side]['outgoing_edges'] = [
                    (e_name, 1 if non_shrunk_edges_for_this_CC_cut[e_name][0] in internal_edge_nodes else -1) for e_name in cut_edge_names
                ]

                # Now create the list of s- and t-channels for this tree structure
                external_outgoing_nodes = set([ non_shrunk_edges_for_this_CC_cut[e_name][1 if direction==1 else 0] for e_name, direction in tree_topologies[side]['outgoing_edges'] ])
                external_incoming_nodes = set([ non_shrunk_edges_for_this_CC_cut[e_name][0 if direction==1 else 0] for e_name, direction in tree_topologies[side]['incoming_edges'] ])
                remaining_internal_nodes = set(internal_edge_nodes)
                pos_leg_id_counter = 0
                edge_name_to_leg_number = {}
                sorted_incoming_edges = sorted(tree_topologies[side]['incoming_edges'],key=lambda e:e[0])
                for edge_name, edge_direction in sorted_incoming_edges:
                    pos_leg_id_counter += 1
                    edge_name_to_leg_number[edge_name] = pos_leg_id_counter
                # For a technical reason in the SingleChannelPhasespace generator, the final state legs must start at 3.
                pos_leg_id_counter = max(pos_leg_id_counter,2)
                for edge_name, edge_direction in sorted(tree_topologies[side]['outgoing_edges'],key=lambda e:e[0]):
                    pos_leg_id_counter += 1
                    edge_name_to_leg_number[edge_name] = pos_leg_id_counter
                neg_leg_id_counter = 0
                for edge_name in tree_topologies[side]['edges']:
                    if edge_name not in edge_name_to_leg_number:
                        neg_leg_id_counter -= 1
                        edge_name_to_leg_number[edge_name] = neg_leg_id_counter
                ancestor_legs_for_node = {}
                for e_name, direction in tree_topologies[side]['outgoing_edges']:
                    internal_node = non_shrunk_edges_for_this_CC_cut[e_name][1 if direction==1  else 0]
                    if internal_node in ancestor_legs_for_node:
                        ancestor_legs_for_node[internal_node].append((e_name,FINAL))
                    else:
                        ancestor_legs_for_node[internal_node] = [ (e_name,FINAL), ]
                for e_name, direction in tree_topologies[side]['incoming_edges']:
                    internal_node = non_shrunk_edges_for_this_CC_cut[e_name][1 if direction==1 else 0]
                    if internal_node in ancestor_legs_for_node:
                        ancestor_legs_for_node[internal_node].append((e_name,INITIAL))
                    else:
                        ancestor_legs_for_node[internal_node] = [ (e_name,INITIAL), ]
                tree_topologies[side]['edge_name_to_leg_number'] = edge_name_to_leg_number

                # Construct the edges necessary for building an nx graph reprepsenting this tree topology
                tree_topologies[side]['nx_graph_edges'] = []
                for edge_name, (u_node, v_node, power) in tree_topologies[side]['edges'].items():
                    nx_nodes = []
                    # Leg number for the nx graph will be -1 for all internal and a positive one for externals
                    nx_leg_number = -1
                    for incoming_edge_name, direction in tree_topologies[side]['incoming_edges']:
                        if incoming_edge_name==edge_name:
                            nx_nodes = (edge_name_to_leg_number[edge_name], -v_node) if direction==1 else (-u_node, edge_name_to_leg_number[edge_name])
                            nx_leg_number = edge_name_to_leg_number[edge_name]
                            break
                    else:
                        for outgoing_edge_name, direction in tree_topologies[side]['outgoing_edges']:
                            if outgoing_edge_name==edge_name:
                                nx_nodes = (-u_node, edge_name_to_leg_number[edge_name]) if direction==1 else (edge_name_to_leg_number[edge_name], -v_node)
                                nx_leg_number = edge_name_to_leg_number[edge_name]
                                break
                        else:
                            nx_nodes = (-u_node, -v_node)

                    tree_topologies[side]['nx_graph_edges'].append(
                        (nx_nodes[0], nx_nodes[1], {'power':power, 'pdg': edges_PDG[edge_name] , 'leg_number': nx_leg_number})
                    )

                s_channels = base_objects.VertexList([])
                t_channels = base_objects.VertexList([])
                sink_edge = sorted_incoming_edges[-1][0] 
                sink_internal_node = non_shrunk_edges_for_this_CC_cut[sink_edge][
                    1 if sorted_incoming_edges[-1][1]==1 else 0
                ]
                sink_external_node = non_shrunk_edges_for_this_CC_cut[sink_edge][
                    0 if sorted_incoming_edges[-1][1]==1 else 1
                ] 
                while len(remaining_internal_nodes)>0:
                    # Find one node if all connected nodes being external but one
                    for node in remaining_internal_nodes:
                        # We want to force to finish on the sink_internal_node, so we veto probing that one until it is last
                        if len(remaining_internal_nodes)>1 and node == sink_internal_node:
                            continue
                        connected_nodes = [
                            (e,nodes[0]) if nodes[1]==node else (e,nodes[1]) 
                            for e, nodes in tree_topologies[side]['edges'].items() if node in nodes
                        ]
                        connected_non_explored_edges = [ (e, a_node) for e, a_node in connected_nodes if a_node in remaining_internal_nodes ]
                        if len(connected_non_explored_edges)>1:
                            continue
                        if len(connected_non_explored_edges)==1:
                            connected_non_explored_edge = connected_non_explored_edges[0][0]
                            connected_non_explored_node = connected_non_explored_edges[0][1]
                            # Add an s- or t-channel

                            connected_external_nodes = set(a_node for e, a_node in connected_nodes if a_node not in remaining_internal_nodes)
                            ancestors_states = [ ancestor_legs_for_node[a_node] for a_node in connected_external_nodes ]
                            final_state_ancestors = []
                            initial_state_ancestors = []
                            for ancestors in ancestors_states:
                                for ancestor_edge_name, ancestor_state in ancestors:
                                    if ancestor_state == FINAL:
                                        final_state_ancestors.append( (ancestor_edge_name, ancestor_state) ) 
                                    else:
                                        initial_state_ancestors.append( (ancestor_edge_name, ancestor_state) )
                            
                            # Add the ancestors to this node:
                            #node_to_update_ancestors_for = connected_non_explored_node
                            node_to_update_ancestors_for = node
                            if node_to_update_ancestors_for in ancestor_legs_for_node:
                                ancestor_legs_for_node[node_to_update_ancestors_for].extend(final_state_ancestors+initial_state_ancestors)
                            else:
                                ancestor_legs_for_node[node_to_update_ancestors_for] = final_state_ancestors+initial_state_ancestors
                            new_vertex = base_objects.Vertex({
                                'id': DUMMY, # Irrelevant
                                'legs': base_objects.LegList(
                                    [
                                        base_objects.Leg({
                                            'id': edges_PDG[edge_name],
                                            'number': edge_name_to_leg_number[edge_name],
                                            'state': INITIAL if edge_name in tree_topologies[side]['incoming_edges'] else FINAL,
                                        }) for edge_name, a_node in connected_nodes if a_node!=connected_non_explored_node
                                    ] + [
                                        base_objects.Leg({
                                            'id': edges_PDG[connected_non_explored_edge],
                                            'number': edge_name_to_leg_number[connected_non_explored_edge],
                                            # The state of internal leg is irrelevant and we choose them to be final
                                            'state': FINAL,
                                        }),
                                    ]
                                )
                            })

                            if len(initial_state_ancestors)==0 or len(final_state_ancestors)==0:
                                # s-channel propagator
                                s_channels.append(new_vertex)
                            else:
                                # t-channel propagator
                                t_channels.append(new_vertex)

                            remaining_internal_nodes.remove(node)
                            break

                        elif len(connected_non_explored_edges)==0:
                            # We are done, simply add the connecting fake vertex now
                            # The algorithm should guarantee that this is the sink_internal_node
                            assert(node==sink_internal_node)

                            # By convention in these phase-space generators based on t- and s-channel decomposition
                            # the last leg of the "gluing" last t-channel is the last initial state with an ID given 
                            # to be the largest available negative number
                            edge_name_to_leg_number[sink_edge] = neg_leg_id_counter-1
                            t_channels.append(base_objects.Vertex({
                                'id': DUMMY, # Irrelevant
                                'legs': base_objects.LegList(
                                    [
                                        base_objects.Leg({
                                            'id': edges_PDG[edge_name],
                                            'number': edge_name_to_leg_number[edge_name],
                                            'state': INITIAL if edge_name in tree_topologies[side]['incoming_edges'] else FINAL,
                                        }) for edge_name, a_node in connected_nodes if a_node!=sink_external_node
                                    ] + [
                                        base_objects.Leg({
                                            'id': edges_PDG[sink_edge],
                                            'number': edge_name_to_leg_number[sink_edge], #edge_name_to_leg_number[sink_internal_node],
                                            # The state of internal leg is irrelevant and we choose them to be final
                                            'state': INITIAL,
                                        }),
                                    ]
                                )
                            }))

                            remaining_internal_nodes.remove(node)
                            break
                    else:
                        raise alphaLoopRunInterfaceError("Could not build kinematic topology for side %s of Cutkosky cut #%d of SG %s."%(
                            side, i_cut, self['name']
                        ))
                
                tree_topologies[side]['s_and_t_propagators'] = ( s_channels, t_channels )
            
            multichannel_info = {}
            multichannel_info['effective_vertices'] = effective_vertices
            multichannel_info['tree_topologies'] = tree_topologies

            # Finally assign the newly generated information to this cutkosky cut
            cutkosky_cut['multichannel_info'] = multichannel_info

        # Now analyse all multichaneling topologies and recognize identical ones and group all LMBs.
        SG_multichannel_info = []
        for i_cut, cutkosky_cut in enumerate(self['cutkosky_cuts']):
            
            nx_graphs = {} 
            all_effective_vertices = []
            for side in ['left','right']:
                all_effective_vertices.extend(cutkosky_cut['multichannel_info']['tree_topologies'][side]['effective_vertices_contained'])
                nx_graphs[side] = nx.MultiDiGraph()
                nx_graphs[side].add_edges_from(cutkosky_cut['multichannel_info']['tree_topologies'][side]['nx_graph_edges'])

            are_both_side_isomorphic = self.is_isomorphic_to(nx_graphs['left'],nx_graphs['right'])
            sides_to_consider_for_SG_multichanneling = ['left',] if are_both_side_isomorphic else ['left','right']

            for side in sides_to_consider_for_SG_multichanneling:

                independent_final_states = cutkosky_cut['multichannel_info']['tree_topologies'][side]['outgoing_edges'][:-1]

                SG_multichannel_info.append(
                    {
                        'cutkosky_cut_id' : i_cut,
                        'side': side,
                        'independent_final_states' : independent_final_states,
                        'loop_LMBs' : []
                    }
                )

                # Build the multichanneling with the cartesian product of all LMBs of all remaining effective vertices.
                for LMB_combination in itertools.product( *[ 
                        cutkosky_cut['multichannel_info']['effective_vertices'][effective_vertex]['loop_momenta_bases'] 
                        for effective_vertex in all_effective_vertices
                    ] ):
                    # Flatten the combined LMBs
                    if len(LMB_combination)>0:
                        LMB_edges_for_this_channel = sum(LMB_combination,[])
                    else:
                        LMB_edges_for_this_channel = []

                    # construct the channel basis to LTD loop momentum basis mapping
                    mat = [ self['edge_signatures'][edge_name][0] 
                                for edge_name in LMB_edges_for_this_channel + [e[0] for e in independent_final_states] ]
                    transformation_matrix = np.linalg.inv(np.array(mat))
                    parametric_shifts = [ [-p for p in self['edge_signatures'][edge_name][1]] 
                                for edge_name in LMB_edges_for_this_channel + [e[0] for e in independent_final_states] ]
                    transformation_to_defining_LMB = (transformation_matrix, parametric_shifts)
                    # The shifts and transformation matrix can be used to perform the transformation as follows:
                    # Built external momenta. Remember that they appear twice.
                    #external_momenta = [ Vector(v[1:]) for v in self.hyperparameters['CrossSection']['incoming_momenta'] ]
                    #external_momenta.extend(external_momenta)
                    #shifts = [ sum([external_momenta[i_shift]*shift 
                    #            for i_shift, shift in enumerate(parametric_shift)]) for parametric_shift in parametric_shifts ]
                    #transformed_momenta_in_LMB = transfo.dot(
                    #    [list(rm+shift) for rm, shift in zip(momenta_to_transform, shifts)] )
                    
                    SG_multichannel_info[-1]['loop_LMBs'].append(
                        {
                            'loop_edges' : LMB_edges_for_this_channel,
                            'transformation_to_defining_LMB' : transformation_to_defining_LMB
                        }
                    )

            self['SG_multichannel_info'] = SG_multichannel_info

    # Make a copy here of the version of this function used for renormalisation as we may want to save/organised different data for it.
    @classmethod
    def shrink_edges(cls, edges, edges_to_shrink):
        subgraph_info = {
            'internal_edges' : {},
            'internal_nodes' : [],
            'in_edges' : {},
            'out_edges' : {},
            'effective_node' : None,
        }
        # First collect all the subgraph nodes
        subgraph_info['internal_nodes'] = sorted(list(set(sum([list(edges[edge][:2]) for edge in edges_to_shrink],[]))))
        # Elect the effective node
        subgraph_info['effective_node'] = subgraph_info['internal_nodes'][0]
        subgraph_info['internal_edges'] = { edge_name: edges[edge_name] for edge_name in edges_to_shrink }
        new_graph_edges = {}
        for ek, ee in edges.items():
            if ek in edges_to_shrink:
                continue
            new_edge = tuple(ee)
            if ee[0] in subgraph_info['internal_nodes']:
                assert(ee[1] not in subgraph_info['internal_nodes'])
                subgraph_info['out_edges'][ek] = ee
                new_edge = (subgraph_info['effective_node'],ee[1],ee[2])
            elif ee[1] in subgraph_info['internal_nodes']:
                assert(ee[0] not in subgraph_info['internal_nodes'])
                subgraph_info['in_edges'][ek] = ee
                new_edge = (ee[0],subgraph_info['effective_node'],ee[2])
            new_graph_edges[ek] = new_edge

        return new_graph_edges, subgraph_info

    @classmethod
    def is_isomorphic_to(cls, directed_grah_A, directed_graph_B):
        """ Uses networkx to decide if the two graphs in argument are isomorphic."""

        def edge_match_function(e1,e2):
            # This function needs tu support multi-edges
            # For now all we do is making sure the set of PDGs of 
            # all edges connecting the two nodes match.
            # This is however fine since we only aim to compare trees here.
            # TODO: Fix ambiguity with fermion flow and part / antipart
            # For now consider two edges equal whenever the *abs* of PDGs matches.
            return set((abs(e['pdg']),e['power'], e['leg_number']) for e in e1.values()) == \
                set((abs(e['pdg']),e['power'], e['leg_number']) for e in e2.values())

        # We do not care about edge orientation for comparing tree kinematic topologies
        graphA, graphB = directed_grah_A.to_undirected(), directed_graph_B.to_undirected()

        # We do not need to consider vertex ID as they do not matter for kinematic topology comparison
        #def node_match_function(n1,n2):
        #    return n1['vertex_id']==n2['vertex_id']
        def node_match_function(n1,n2):
            return True

        return nx.is_isomorphic(graphA, graphB,
            edge_match=edge_match_function,
            node_match=lambda n1,n2: node_match_function
        )

    def export(self, SG_name, dir_path):
        with open(pjoin(dir_path,'PROCESSED_%s.yaml'%SG_name),'w') as f:
            f.write(yaml.dump(dict(self), Dumper=Dumper, default_flow_style=False))

class SuperGraphCollection(dict):
    
    def __init__(self, *args, **opts):
        super(SuperGraphCollection, self).__init__(*args, **opts)

    def export(self, dir_path):
        for SG_name, SG in self.items():
            SG.export(SG_name, dir_path)

    def summary_str(self):

        res_str = []
        res_str.append('')
        for SG_name in sorted(list(self.keys())):
            SG = self[SG_name]
            res_str.append("\nSummary of the supergraph %s%s%s:\n%s"%(Colours.GREEN, SG_name, Colours.END, SG.summary_str()))
        res_str.append('')

        res = []
        res.append(('Total number of supergraphs','%d'%len(self)))
        res.append(('Total number of cuts','%d'%(sum( len(sg['cutkosky_cuts']) for sg in self.values() ))))
        res.append(('Contributions distribution:',''))
        distribution = {}
        for SG_name in sorted(list(self.keys())):
            SG = self[SG_name]

            for i_cut, cut in enumerate(SG['cutkosky_cuts']):
                cut_info = SG.get_cut_characteristics(cut)
                signature = (   cut_info['n_loops_left']+cut_info['n_loops_right'], cut_info['n_loops_left'], 
                                cut_info['n_CC'],cut_info['n_loops_right']  )
                if signature not in distribution:
                    distribution[signature] = {'n_diag_sets': cut_info['n_diag_sets'], 'n_cuts': 1, 'cuts': [i_cut]}
                else:
                    distribution[signature]['n_diag_sets'] += cut_info['n_diag_sets']
                    distribution[signature]['n_cuts'] += 1
                    distribution[signature]['cuts'].append(i_cut)

        sorted_keys = sorted(distribution.keys())
        for n_loops_tot, n_loops_left, n_CC, n_loops_right in sorted_keys:
            entry = distribution[(n_loops_tot, n_loops_left, n_CC, n_loops_right)]
            res.append(('( n_loops_tot=%d, n_loops_left=%d, n_CC=%d, n_loops_right=%d )'%(
                n_loops_tot, n_loops_left, n_CC, n_loops_right
            ),'n_cuts=%d n_diag_sets=%d'%(entry['n_cuts'],entry['n_diag_sets'])))

        for k, v in res:
            if v!='':
                res_str.append('%s%-40s%s : %s'%(Colours.GREEN,k,Colours.END,v))
            else:
                res_str.append('%s%s%s'%(Colours.BLUE,k,Colours.END))
        res_str.append('')

        return '\n'.join(res_str)

    def show_IR_statistics(self, show_momenta=True):

        # TODO improve rendering
        res_str = []

        res_str.append('')
        res_str.append('%s%s%s'%(Colours.BLUE,'IR statistics for individual supergraphs',Colours.END))
        for SG_name in sorted(list(self.keys())):
            if 'E_surfaces_analysis' not in self[SG_name]:
                continue
            res_str.append("\nIR profile of %s%s%s:\n%s"%(Colours.GREEN,SG_name,Colours.END, self[SG_name].show_IR_statistics(show_momenta=show_momenta)))

        return '\n'.join(res_str)

    def show_UV_statistics(self):

        res = []

        all_SG_dods = sorted([
            (SG_name,k,v) for SG_name,SG in self.items() if 'DERIVED_UV_dod' in SG for k,v in SG['DERIVED_UV_dod'].items()
        ],key=lambda el: el[2][0], reverse=True)
        fail_SG_dods = [ (SG_name,k,v) for SG_name,k,v in all_SG_dods if not v[-1] ]
        all_SG_cut_dods = sorted([
            (SG_name, cut_ID, k, v) for SG_name,SG in self.items() for cut_ID, cut in enumerate(SG['cutkosky_cuts'])
            if 'DERIVED_UV_dod' in cut for k,v in cut['DERIVED_UV_dod'].items()
        ],key=lambda el: el[3][0], reverse=True)
        fail_SG_cut_dods = [ (SG_name, cut_ID, k, v) for SG_name,cut_ID,k,v in all_SG_cut_dods if not v[-1] ]
        if len(all_SG_dods)==0 and len(all_SG_cut_dods)==0:
            return "No UV profile information available."
        
        all_SG_dod_consistencies = sorted([
            (SG_name,SG['DERIVED_dod_consistency']) for SG_name,SG in self.items() if 'DERIVED_dod_consistency' in SG
        ],key=lambda el: abs(el[1][0]), reverse=True)

        res.append(('Overall UV profile',''))
        res.append(('-'*len(res[-1][0]),''))
        res.append(('Total number of UV tests','%d'%(len(all_SG_dods)+len(all_SG_cut_dods))))
        if len(all_SG_dod_consistencies)>0:
            res.append(("Maximum inconsistency between dod of eval_integrand and eval_cuts", "%.4f (%.4f vs %.4f) (%s)"%(
                all_SG_dod_consistencies[0][1][0],all_SG_dod_consistencies[0][1][1], 
                all_SG_dod_consistencies[0][1][2], all_SG_dod_consistencies[0][0])) )
        res.append(('Failed UV tests for summed cuts','%s%d%s%s'%(
                Colours.GREEN if len(fail_SG_dods)==0 else Colours.RED, len(fail_SG_dods), Colours.END,
                    (('\n'+'\n'.join('%-60s -> %s'%(
                        '%s%s @ LMB=%s @ UV_indices=%s'%(
                            Colours.RED, SG_name, ','.join(test_k[0]),str(test_k[1])
                        ),
                        'dod=%-6.4f +/- %-6.4f%s'%(
                            test_v[0], test_v[1], Colours.END
                        )
                    ) 
                    for SG_name, test_k, test_v in fail_SG_dods
                )) if len(fail_SG_dods)>0 else '')
            )
        ))
        res.append(('Failed UV tests for individual cuts','%s%d%s%s'%(
                Colours.GREEN if len(fail_SG_cut_dods)==0 else Colours.RED, len(fail_SG_cut_dods), Colours.END,
                    (('\n'+'\n'.join('%-60s -> %s'%(
                            '%s%s @ cut_ID=%d @ LMB=%s @ UV_indices=%s'%(
                                Colours.RED,SG_name, cut_ID, ','.join(test_k[0]),str(test_k[1])
                        ),
                            'dod=%-6.4f +/- %-64f%s'%(
                                test_v[0], test_v[1], Colours.END
                        )
                    )
                    for SG_name, cut_ID, test_k, test_v in fail_SG_cut_dods
                )) if len(fail_SG_cut_dods)>0 else '')
            )
        ))
        if len(all_SG_dods)>0:
            res.append(('Maximum dod for summed cuts over all SGs','%-6.4f +/- %-6.4f (%s @ LMB=%s @ UV_indices=%s)'%(
                all_SG_dods[0][2][0], all_SG_dods[0][2][1],all_SG_dods[0][0], 
                ','.join(all_SG_dods[0][1][0]), str(all_SG_dods[0][1][1]) )))
        if len(all_SG_cut_dods)>0:
            res.append(('Maximum dod for individual cuts','%-6.4f +/- %-6.4f (%s @ cut_ID=%d @ LMB=%s @ UV_indices=%s)'%(
                all_SG_cut_dods[0][3][0], all_SG_cut_dods[0][3][1],all_SG_cut_dods[0][0],all_SG_cut_dods[0][1],
                ','.join(all_SG_cut_dods[0][2][0]), str(all_SG_cut_dods[0][2][1]) )))
        
        allFailedSGNames = sorted(list(set([SG_name for SG_name,k,v in fail_SG_dods]+[SG_name for SG_name,cut_ID,k,v in fail_SG_cut_dods])))
        res.append(('Non-converging supergraphs','%s(%d/%d)%s'%(
            Colours.GREEN if len(allFailedSGNames)==0 else Colours.RED, len(allFailedSGNames),len(self),' : ' if len(allFailedSGNames)>0 else ''
        )+'%s%s%s'%(Colours.RED, ', '.join(allFailedSGNames), Colours.END) ))

        res_str = []

        res_str.append('')
        res_str.append('%s%s%s'%(Colours.BLUE,'UV statistics for individual supergraphs',Colours.END))
        for SG_name in sorted(list(self.keys())):
            res_str.append("\nUV profile of %s%s%s:\n%s"%(Colours.GREEN,SG_name,Colours.END, self[SG_name].show_UV_statistics()))

        res_str.append('')
        max_len = (max(len(k) for k, v in res if v!='') if len(res)>0 else 45)+5
        for k, v in res:
            if v!='':
                res_str.append(('%s%-{:d}s%s : %s'.format(max_len))%(Colours.GREEN,k,Colours.END,v))
            else:
                res_str.append('%s%s%s'%(Colours.BLUE,k,Colours.END))

        return '\n'.join(res_str)

    def show_timing_statistics(self):

        res = []
        for precision in ['f64','f128']:
            all_SG_times = sorted([
                (SG_name,SG['DERIVED_timing_profile_%s'%precision]) for SG_name,SG in self.items()
                if 'DERIVED_timing_profile_%s'%precision in SG
            ],key=lambda el: el[1])
            all_SG_summed_times = sorted([
                (SG_name,sum(c['DERIVED_timing_profile_%s'%precision] for c in SG['cutkosky_cuts'])) for SG_name,SG in self.items()
                if 'DERIVED_timing_profile_%s'%precision in SG
            ],key=lambda el: el[1])
            if len(all_SG_times)==0:
                continue
            res.append(('Overall statistics for %s precision:'%precision,''))
            res.append(('-'*len(res[-1][0]),''))
            res.append(('Summed timing for overall evaluate','%.3f ms'%sum(t[1] for t in all_SG_times)))
            res.append(('Summed timing for summed evaluate cuts','%.3f ms'%sum(t[1] for t in all_SG_summed_times)))
            res.append(('Average timing over all supergaphs','%.3f ms'%(sum(t[1] for t in all_SG_times)/float(len(all_SG_times)))))
            res.append(("Slowest supergraph evaluation", "%.3f ms (%s)"%( all_SG_times[-1][1],all_SG_times[-1][0] )))
            res.append(("Fastest supergraph evaluation", "%.3f ms (%s)"%( all_SG_times[0][1],all_SG_times[0][0] )))

        if len(res)==0:
            return "No timing profile information available."

        res_str = []

        res_str.append('')
        res_str.append('%s%s%s'%(Colours.BLUE,'Statistics for individual supergraphs',Colours.END))
        for SG_name in sorted(list(self.keys())):
            res_str.append("\nTiming profile for supergraph %s%s%s:\n%s"%(Colours.GREEN,SG_name,Colours.END, self[SG_name].show_timing_statistics()))

        res_str.append('')
        for k, v in res:
            if v!='':
                res_str.append('%s%-40s%s : %s'%(Colours.GREEN,k,Colours.END,v))
            else:
                res_str.append('%s%s%s'%(Colours.BLUE,k,Colours.END))

        return '\n'.join(res_str)

# To be used as a decorator
def with_tmp_hyperparameters(options_to_overwrite=None):
    """ Decorate a function and automatically temporarily overwrite hyperparameters.yaml"""
    def add_hyperparameters_in_function(f):
        
        def modified_function(*args, **opt):
            orig_hyperparameters = copy.deepcopy(args[0].hyperparameters)
            try:
                if options_to_overwrite is not None:
                    for param, value in options_to_overwrite.items():
                        args[0].hyperparameters.set_parameter(param, value)
                return f(*args, **opt)
            except:
                args[0].hyperparameters = orig_hyperparameters
                raise
            args[0].hyperparameters = orig_hyperparameters
        return modified_function
    
    return add_hyperparameters_in_function

# To be used as a decorator
def wrap_in_process():
    """ Decorate the function so as to automatically sandbox it in a separate Process."""
    def wrap_function_in_process(f):
        q=multiprocessing.Queue()
        def fowarad_function_output_to_queue(*args):
            q.put(f(*args[1:],**args[0]))
        def modified_function(*args, **opts):
            # Note: using the multiprocessing instead of the threading module breaks on MacOS with python 3.8+ because of
            # a backward incompatible change in the way processes are spawn by Python on MacOS. 
            # p = multiprocessing.Process(target=fowarad_function_output_to_queue, args=tuple([opts,]+list(args)))
            # p.start()
            # p.join()
            t = threading.Thread(target=fowarad_function_output_to_queue, args=tuple([opts,]+list(args)))
            t.start()
            return q.get()
        return modified_function

    return wrap_function_in_process

class alphaLoopRunInterface(madgraph_interface.MadGraphCmd, cmd.CmdShell):
    """ Interface for steering the running of an alphaLoop output.
    We make it inherit from CmdShell so that launch_ext_prog does not attempt to start in WebMode."""

    _rust_inputs_folder = 'Rust_inputs'
    _cross_section_set_yaml_name = 'all_QG_supergraphs.yaml'
    _run_workspace_folder = 'run_workspace'
    _FORM_folder = 'FORM'

    def __init__(self, dir_path, alphaLoop_interface, launch_options={}, *args, **opts):
        """ Define attributes of this class."""

        self.dir_path = dir_path
        self.alphaLoop_interface = alphaLoop_interface
        self.launch_options = launch_options

        if not os.path.isdir(pjoin(self.dir_path, self._run_workspace_folder)):
            os.makedirs(pjoin(self.dir_path, self._run_workspace_folder))
            os.makedirs(pjoin(self.dir_path, self._run_workspace_folder,'stats'))

        self.hyperparameters = RunHyperparameters(dir_path=self.dir_path)
        if self.launch_options['reuse']!='NO_REUSE':
            if self.launch_options['reuse']=='default':
                hyperparameters_path = pjoin(self.dir_path, self._run_workspace_folder, 'hyperparameters.yaml')
            else:
                hyperparameters_path = os.path.abspath(self.launch_options['reuse'])
            if not os.path.isfile(hyperparameters_path):
                raise alphaLoopInvalidRunCmd("Cannot reuse existing workspace hyperparameters as they are not found at:\n%s"%(
                    hyperparameters_path))
            else:
                self.hyperparameters.update(yaml.load(open(hyperparameters_path, 'r'), Loader=yaml.Loader))                

        self.hyperparameters.export_to(pjoin(self.dir_path, self._run_workspace_folder, 'hyperparameters.yaml'))

        if os.path.isfile(pjoin(self.dir_path, self._FORM_folder,'generation_statistics.txt')):
            self.generation_statistics = open(pjoin(self.dir_path, self._FORM_folder,'generation_statistics.txt'),'r').read()
        else:
            self.generation_statistics = {}

        # The cross-section set file may not be standard
        if os.path.isfile(pjoin(self.dir_path, self._rust_inputs_folder, self._cross_section_set_yaml_name)):
            cross_section_set_yaml_file_path = pjoin(self.dir_path, self._rust_inputs_folder, self._cross_section_set_yaml_name)
        else:
            candidates = [fp for fp in glob.glob(pjoin(self.dir_path, self._rust_inputs_folder,'*.yaml')) if 
                            not os.path.basename(fp).startswith('SG') and not os.path.basename(fp).startswith('PROCESSED_SG')]
            if len(candidates)!=1:
                raise alphaLoopInvalidRunCmd("Could not find cross-section set yaml file in path %s"%(pjoin(self.dir_path, self._rust_inputs_folder)))
            cross_section_set_yaml_file_path = candidates[0]

        self.cross_section_set = CrossSectionSet(cross_section_set_yaml_file_path)
        self.all_supergraphs = self.load_supergraphs()

        super(alphaLoopRunInterface, self).__init__(*args, **opts)

    def get_rust_worker(self, supergraph_name):
        """ Return a rust worker instance to evaluate the LU representation """

        try:
            # Import the rust bindings
            from ltd import CrossSection
        except ImportError:
            raise alphaLoopRunInterfaceError("ERROR: Could not import the rust back-end 'ltd' module. Compile it first with:\n"
                " ./make_lib\nfrom within the pyNLoop directory.")
    
        #os.environ['MG_NUMERATOR_PATH'] = proc_path if proc_path.endswith('/') else '%s/'%proc_path

        # Write hyperparameters to read in in a tmp file
        self.hyperparameters.export_to(pjoin(self.dir_path, self._run_workspace_folder, 'tmp_hyperparameters.yaml'))            

        try:
            rust_worker = CrossSection( 
                pjoin(self.dir_path, self._rust_inputs_folder, supergraph_name+'.yaml'), 
                pjoin(self.dir_path, self._run_workspace_folder, 'tmp_hyperparameters.yaml') 
            )
        except:
            os.remove(pjoin(self.dir_path, self._run_workspace_folder, 'tmp_hyperparameters.yaml'))
            raise

        return CallableInstanceWrapper(rust_worker)

    def load_supergraphs(self):
        
        SG_collection = SuperGraphCollection()

        for SG in self.cross_section_set['topologies']:
            yaml_path = pjoin(self.dir_path, self._rust_inputs_folder, 'PROCESSED_'+SG['name']+'.yaml')
            if not os.path.isfile(yaml_path):
                yaml_path = pjoin(self.dir_path, self._rust_inputs_folder, SG['name']+'.yaml')
                if not os.path.isfile(yaml_path):
                    raise alphaLoopInvalidRunCmd("Could not find yaml file at '%s' specifying the supergraph information."%yaml_path)
            SG_collection[SG['name']] = SuperGraph(yaml_path)
        
        return SG_collection

    #### TIMING PROFILE COMMAND
    timing_profile_parser = ArgumentParser(prog='timing_profile')
    timing_profile_parser.add_argument('SG_name', metavar='SG_name', type=str, nargs='?',
                    help='the name of a supergraph to display')
    timing_profile_parser.add_argument("-n","--n_points", dest='n_points', type=int, default=0,
                    help='force a certain number of points to be considered for the timing profile')
    timing_profile_parser.add_argument("-s","--seed", dest='seed', type=int, default=0,
                    help='specify random seed')
    timing_profile_parser.add_argument("-t","--time", dest='time', type=float, default=5.0,
                    help='target evaluation time per profile, in seconds.')
    timing_profile_parser.add_argument(
        "-f", "--f128", action="store_true", dest="f128", default=False,
        help="Enable timing profile of the f128 output too.")
    def help_timing_profile(self):
        self.timing_profile_parser.print_help()
        return
    # We must wrape this function in a process because of the border effects of the pyO3 rust Python bindings
    @wrap_in_process()
    @with_tmp_hyperparameters({
        'Integrator.dashboard': False,
        'General.multi_channeling' : False
    })
    def do_timing_profile(self, line):
        """ Automatically timing profile a process output."""

        if line=='help':
            self.timing_profile_parser.print_help()
            return 

        args = self.split_arg(line)
        args = self.timing_profile_parser.parse_args(args)

        if args.SG_name is None:
            selected_SGs = list(self.all_supergraphs.keys())
        else:
            selected_SGs = [args.SG_name,]
        
        if args.seed != 0:
            random.seed(args.seed)

        max_count = sum( len(self.all_supergraphs[SG_name]['cutkosky_cuts']) for SG_name in selected_SGs )*(
            2 if args.f128 else 1 )
        logger.info("Starting timing profile...")
        # WARNING it is important that the rust workers instantiated only go out of scope when this function terminates
        self.hyperparameters['General']['stability_checks']=[self.hyperparameters['General']['stability_checks'][0],]
        self.hyperparameters['General']['stability_checks'][0]['use_f128']=False
        self.hyperparameters['General']['stability_checks'][0]['n_samples']=1
        self.hyperparameters['General']['stability_checks'][0]['relative_precision']=1.0e-99
        self.hyperparameters['General']['stability_checks'][0]['escalate_for_large_weight_threshold']=-1.0
        self.hyperparameters['General']['stability_checks'][0]['minimal_precision_to_skip_further_checks']=1.0e-99
        rust_workers = {SG_name: self.get_rust_worker(SG_name) for SG_name in selected_SGs}
        if args.f128:
            hyperparameters_backup=copy.deepcopy(self.hyperparameters)
            for entry in self.hyperparameters['General']['stability_checks']:
                entry['use_f128'] = True
            rust_workers_f128 = {SG_name: self.get_rust_worker(SG_name) for SG_name in selected_SGs}
            self.hyperparameters = hyperparameters_backup

        t_start_profile = time.time()
        with progressbar.ProgressBar(
                prefix=("Timing profile: {variables.SG}/{variables.n_SG} ({variables.SG_name}), "+
                        "cut ID: {variables.cut}/{variables.n_cuts}, Avg t per cut: {variables.t} [ms]"), 
                max_value=max_count,variables={
                    'SG_name':'N/A', 'SG': '0', 'n_SG':len(selected_SGs),'cut':0, 
                    'n_cuts': len(self.all_supergraphs[selected_SGs[0]]['cutkosky_cuts']),
                    't' : 'N/A'
                }
            ) as bar:
            running_avg_time_per_cut = 0.0
            n_cuts = 0.0
            for i_SG, SG_name in enumerate(selected_SGs):

                bar.update(SG_name=SG_name)
                bar.update(SG=i_SG)

                #rust_worker = self.get_rust_worker(SG_name)
                rust_worker = rust_workers[SG_name]
                if args.f128:
                    rust_worker_f128 = rust_workers_f128[SG_name]
                SG = self.all_supergraphs[SG_name]
                E_cm = SG.get_E_cm(self.hyperparameters)
                funcs_to_test = [('f64', rust_worker.evaluate_integrand)]
                if args.f128:
                    funcs_to_test.append(('f128', rust_worker_f128.evaluate_integrand))
                for precision, rust_function in funcs_to_test:
                    if args.n_points == 0:
                        t_start = time.time()
                        _res = rust_function(SG.get_random_x_input())
                        delta_t = time.time()-t_start
                        n_points = int(args.time/delta_t)
                    else:
                        n_points = args.n_points
                    
                    t_start = time.time()
                    for _ in range(n_points):
                        _res = rust_function(SG.get_random_x_input())
                        #misc.sprint('A', _res)
                    delta_t = time.time()-t_start

                    SG['DERIVED_timing_profile_%s'%precision] = (delta_t/float(n_points))*1000.0
                
                for cut_ID, cut in enumerate(SG['cutkosky_cuts']):
                    bar.update(cut=(cut_ID+1))
                    funcs_to_test = [('f64', rust_worker.get_scaling, rust_worker.evaluate_cut)]
                    if args.f128:
                        funcs_to_test.append(('f128', rust_worker.get_scaling, rust_worker.evaluate_cut_f128))
                    for precision, get_scaling_function, cut_evaluate_function in funcs_to_test:
                        if args.n_points == 0:
                            t_start = time.time()
                            random_momenta = SG.get_random_momenta_input(E_cm)
                            scaling_solutions = list(get_scaling_function(random_momenta,cut_ID))
                            scaling, scaling_jacobian = scaling_solutions.pop(0)
                            while scaling < 0.0:
                                if len(scaling_solutions)==0:
                                    break
                                scaling, scaling_jacobian = scaling_solutions.pop(0)
                            _res = cut_evaluate_function(random_momenta,cut_ID,scaling,scaling_jacobian)
                            delta_t = time.time()-t_start
                            n_points = int(args.time/delta_t)
                        else:
                            n_points = args.n_points

                        t_start = time.time()
                        for _ in range(n_points):
                            random_momenta = SG.get_random_momenta_input(E_cm)
                            scaling_solutions = list(get_scaling_function(random_momenta,cut_ID))
                            scaling, scaling_jacobian = scaling_solutions.pop(0)
                            while scaling < 0.0:
                                if len(scaling_solutions)==0:
                                    break
                                scaling, scaling_jacobian = scaling_solutions.pop(0)
                            _res = cut_evaluate_function(random_momenta,cut_ID,scaling,scaling_jacobian)
                            #misc.sprint('B', _res)
                        delta_t = time.time()-t_start
                        t_cut = (delta_t/float(n_points))*1000.0
                        running_avg_time_per_cut += t_cut
                        n_cuts += 1
                        bar.update(t='%.3f'%(running_avg_time_per_cut/float(n_cuts)))
                        bar.update(bar.value+1)

                        cut['DERIVED_timing_profile_%s'%precision] = t_cut

        delta_t = time.time()-t_start_profile
        logger.info("Timing profile completed in %d [s]."%(int(delta_t)))
        # Write out the results into processed topologies
        logger.info("Writing out processed yaml supergaphs on disk...")
        self.all_supergraphs.export(pjoin(self.dir_path, self._rust_inputs_folder))
        if len(selected_SGs)==1:
            self.do_display('%s --timing'%selected_SGs[0])
        else:
            self.do_display('--timing')


    #### IR PROFILE COMMAND
    ir_profile_parser = ArgumentParser(prog='ir_profile')
    ir_profile_parser.add_argument('SG_name', metavar='SG_name', type=str, nargs='?',
                    help='the name of a supergraph to display')
    ir_profile_parser.add_argument("-n","--n_points", dest='n_points', type=int, default=20,
                    help='force a certain number of points to be considered for the ir profile')
    ir_profile_parser.add_argument("-max","--max_scaling", dest='max_scaling', type=float, default=1.0e-09,
                    help='maximum IR scaling to consider')
    ir_profile_parser.add_argument("-min","--min_scaling", dest='min_scaling', type=float, default=1.0e0,
                    help='minimum IR scaling to consider')
    ir_profile_parser.add_argument("-s","--seed", dest='seed', type=int, default=0,
                    help='specify random seed')
    ir_profile_parser.add_argument("-rp","--required_precision", dest='required_precision', type=float, default=None,
                    help='minimum required relative precision for returning a result.')
    ir_profile_parser.add_argument("-t","--target_scaling", dest='target_scaling', type=int, default=0,
                    help='set target IR scaling (default=0)')
    ir_profile_parser.add_argument(
        "-mm", "--use_mathematica", action="store_true", dest="mathematica", default=False,
        help="Use a mathematica analysis of the E-surfaces.")
    ir_profile_parser.add_argument(
        "-f", "--f128", action="store_true", dest="f128", default=False,
        help="Perfom the UV profile using f128 arithmetics.")
    ir_profile_parser.add_argument(
        "-nf", "--no_f128", action="store_true", dest="no_f128", default=False,
        help="Forbid automatic promotion to f128.")
    ir_profile_parser.add_argument(
        "-nw", "--no_warnings", action="store_false", dest="show_warnings", default=True,
        help="Do not show warnings about this profiling run.")
    ir_profile_parser.add_argument(
        "-srw", "--show_rust_warnings", action="store_true", dest="show_rust_warnings", default=True,
        help="Show rust warnings.")
    ir_profile_parser.add_argument(
        "-nsof", "--no_skip_once_failed", action="store_false", dest="skip_once_failed", default=True,
        help="Do not skip the probing of a supergraph once it failed.")
    ir_profile_parser.add_argument(
        "-nsf", "--no_show_fails", action="store_false", dest="show_fails", default=True,
        help="Show exhaustive information for each fail.")
    ir_profile_parser.add_argument(
        "-relevant_cuts", "--only_relevant_cuts", action="store_true", dest="only_relevant_cuts", default=False,
        help="Only explore the scaling of cuts relevant for a particular E-surface intersection configuration.")
    ir_profile_parser.add_argument("-n_max","--n_max", dest='n_max', type=int, default=-1,
                    help='Set the maximum number of IR tests to perform per SG (default: all)')
    ir_profile_parser.add_argument("-maxnE","--max_E_surfaces_in_intersections", dest='max_E_surfaces_in_intersections', type=int, default=3,
                    help='Set the maximum number of E-surfaces in an intersection (default: all)')
    ir_profile_parser.add_argument("-nshifts","--n_shifts_to_test_for_finding_intersection", dest='n_shifts_to_test_for_finding_intersection', type=int, default=10,
                    help='Set the maximum number of shifts to test for finding an intersection (default: %(default)s)')
    ir_profile_parser.add_argument("-e_surfaces","--e_surfaces", dest='selected_e_surfaces', type=str, nargs='*', default=None,
                    help='Set the particular E-surfaces to study by specifying their edge names. Example --e_surfaces ("pq1","pq3") ("pq3","pq5","pq8") (default: All)')
    ir_profile_parser.add_argument("-intersections","--intersections", dest='intersections', type=str, nargs='*', default=None,
                    help='Only when specifying E-surfaces, then this options allows to specify the intersection of interest. Example --intersections (0,1) (1,2) (default: All)')
    ir_profile_parser.add_argument("-intersection_point","--intersection_point", dest='intersection_point', type=str, default=None,
                    help='Specify the intersection point in the LMB (excluding frozen momenta). Example --intersection_point (0.1244323,2.432e+02,...") (default: Automatic)')
    ir_profile_parser.add_argument("-approach_direction","--approach_direction", dest='approach_direction', type=str, default=None,
                    help='Particular direction in LMB used for approaching the intersection point. Example --approach_direction "(0.1244323,2.432e+02,1.03,...") (default: random)')
    ir_profile_parser.add_argument("-reanalyze","--reanalyze_E_surfaces", action="store_true", dest="reanalyze_E_surfaces", default=False,
                    help='Force the re-analysis of E-surfaces even if result already found in cache (default: %(default)s)')
    ir_profile_parser.add_argument(
        "-sm","--show_momenta", action="store_true", dest="show_momenta", default=False,
        help="Show the momenta of the edges in the E-surfaces for the intersection point approached in the IR.")
    ir_profile_parser.add_argument(
        "-v", "--verbose", action="store_true", dest="verbose", default=False,
        help="Enable verbose output.")
    def help_ir_profile(self):
        self.ir_profile_parser.print_help()
        return
    # We must wrap this function in a process because of the border effects of the pyO3 rust Python bindings
    @wrap_in_process()
    @with_tmp_hyperparameters({
        'Integrator.dashboard': False,
        'General.minimal_precision_for_returning_result': 1.0,
        'CrossSection.NormalisingFunction.name'         : 'left_right_exponential',
        'CrossSection.NormalisingFunction.center'       : 1.0,
        'CrossSection.NormalisingFunction.spread'       : 1.0,
        'General.multi_channeling'                      : False
    })
    def do_ir_profile(self, line):
        """ Automatically probe all UV limits of a process output."""

        from alpha_loop.E_surface_intersection_finder import EsurfaceIntersectionFinder
        import cvxpy

        if line=='help':
            self.ir_profile_parser.print_help()
            return 

        args = self.split_arg(line)
        args = self.ir_profile_parser.parse_args(args)

        # We need to detect here if we are in the amplitude-mock-up situation with frozen external momenta.
        frozen_momenta = None
        if 'external_data' in self.cross_section_set:
            frozen_momenta = {
                'in' : self.cross_section_set['external_data']['in_momenta'],
                'out' : self.cross_section_set['external_data']['out_momenta'],
            }
            # Also force the specified incoming momenta specified in the hyperparameters to match the frozen specified ones.
            self.hyperparameters.set_parameter('CrossSection.incoming_momenta',frozen_momenta['in'])
            self.hyperparameters.set_parameter('CrossSection.do_rescaling',False)
            self.hyperparameters.set_parameter('CrossSection.fixed_cut_momenta',frozen_momenta['out'])

        if args.selected_e_surfaces is not None:
            args.selected_e_surfaces = [
                set(eval(e_surfs)) for e_surfs in args.selected_e_surfaces
            ]
        if args.intersections is not None:
            if args.selected_e_surfaces is None:
                raise alphaLoopInvalidRunCmd("The --intersections option can only be used together with the --e_surfaces one.")
            args.intersections = [eval(inter) for inter in args.intersections]

        if args.intersection_point is not None:
            args.intersection_point = eval(args.intersection_point)

        if args.approach_direction is not None:
            args.approach_direction = eval(args.approach_direction)

        if args.SG_name is None:
            selected_SGs = list(self.all_supergraphs.keys())
        else:
            selected_SGs = [args.SG_name,]

        if args.required_precision is None:
            self.hyperparameters['General']['stability_checks'][-1]['relative_precision']=1.0e-99
        else:
            for entry in self.hyperparameters['General']['stability_checks']:
                entry['relative_precision'] = args.required_precision

        logger.info("Starting IR profile...")

        # Prepare the run
        IR_info_per_SG_and_E_surfaces_set = {}

        n_intersections_found = 0
        n_intersections_rejected = 0

        with progressbar.ProgressBar(
                prefix=("IR preparation. E-surface intersections: {variables.intersection} {variables.i_comb}/{variables.n_comb} {variables.inter_found}\u2713, {variables.inter_failed}\u2717, SG: {variables.SG_name} "), 
                max_value=len(selected_SGs),variables={
                    'intersection': 'N/A', 'inter_found': 0, 'inter_failed': 0, 'SG_name': 'N/A', 'i_comb': 0, 'n_comb': 0
                }
            ) as bar:

            for i_SG, SG_name in enumerate(selected_SGs):
                
                bar.update(SG_name=SG_name)
                bar.update(i_SG)

                SG = self.all_supergraphs[SG_name]

                if not args.reanalyze_E_surfaces and all(entry in SG for entry in ['E_surfaces','E_surfaces_intersection']):
                    
                    user_E_surfaces = []
                    all_E_surfaces_found = True
                    if args.selected_e_surfaces is not None:
                        for e_surf in args.selected_e_surfaces:
                            for E_surf in SG['E_surfaces']:
                                if set([osp['name'] for osp in E_surf['onshell_propagators']]) == e_surf:
                                    user_E_surfaces.append(E_surf)
                                    break
                        all_E_surfaces_found = (len(user_E_surfaces)==len(args.selected_e_surfaces))

                    
                    user_intersections = {}
                    all_intersections_found = True
                    if args.intersections is not None:
                        n_inter_found = 0
                        for intersection in args.intersections:
                            if len(intersection) not in SG['E_surfaces_intersection']:
                                continue
                            for inter, inter_info in SG['E_surfaces_intersection'][len(intersection)].items():
                                if inter==tuple(sorted([user_E_surfaces[inter_E_id]['id'] for inter_E_id in intersection])):
                                    n_inter_found += 1
                                    if len(intersection) not in user_intersections:
                                        user_intersections[len(intersection)] = {}    
                                    user_intersections[len(intersection)][inter] = inter_info
                                    break
                        all_intersections_found = (n_inter_found == len(args.intersections))

                    elif args.selected_e_surfaces is not None:
                        for inter_length in SG['E_surfaces_intersection']:
                            for inter, inter_info in SG['E_surfaces_intersection'][inter_length].items():
                                if all( (E_surf_id in [user_E_surf['id'] for user_E_surf in user_E_surfaces]) for E_surf_id in inter):
                                    if inter_length not in user_intersections:
                                            user_intersections[inter_length] = {}    
                                    user_intersections[inter_length][inter] = inter_info

                    if args.selected_e_surfaces is None:
                        user_E_surfaces = SG['E_surfaces']
                        user_intersections = SG['E_surfaces_intersection']

                    if args.intersection_point is not None:
                        if len(user_intersections) != 1 or len(user_intersections[list(user_intersections.keys())[0]]) != 1:
                            raise alphaLoopInvalidRunCmd('An intersection point can only be provided if only one test is selected with the option --e_surfaces and --intersections.')
                        intersection_length = list(user_intersections.keys())[0]
                        intersection_key = list(user_intersections[intersection_length].keys())[0]
                        user_intersections[intersection_length][intersection_key]['intersection_point'] = [ args.intersection_point[i:i+3] for i in range(0,len(args.intersection_point),3) ]

                    if all_E_surfaces_found and all_intersections_found:
                        IR_info_per_SG_and_E_surfaces_set[SG_name] = {
                            'E_surfaces' : user_E_surfaces,
                            'E_surfaces_intersection' : user_intersections
                        }
                        bar.update(intersection='SKIPPED')
                        continue
                
                IR_info_per_SG_and_E_surfaces_set[SG_name] = {}

                if args.seed != 0:
                    random.seed(args.seed)

                E_cm = SG.get_E_cm(self.hyperparameters)

                # First we must regenerate a TopologyGenerator instance for this supergraph
                edges_list = SG['topo_edges']
                SG_topo_gen = ltd_utils.TopologyGenerator(
                    [e[:-1] for e in edges_list],
                    powers = { e[0] : e[-1] for e in edges_list }
                )
                loop_SG = ltd_utils.LoopTopology.from_flat_format(SG['topo'])

                # We must adjust entries in the loop_SG for the specific external momenta specified
                external_momenta = [ LorentzVector(v) for v in self.hyperparameters['CrossSection']['incoming_momenta'] ]
                external_momenta.extend(external_momenta)
                loop_SG.external_kinematics = LorentzVectorList([list(v) for v in external_momenta])
                for i_ll, ll in enumerate(loop_SG.loop_lines):
                    for i_p, p in enumerate(ll.propagators):
                        p.q = sum([ external_momenta[i_shift]*wgt for i_shift, wgt in enumerate(SG['edge_signatures'][p.name][1]) ])

                edge_signatures = SG['edge_signatures']

                cvxpy_source_coordinates = [cvxpy.Variable(3) for _ in range(loop_SG.n_loops)]

                pinched_E_surface_keys = []
                extra_info = {}
                consider_pinches = None
                #consider_pinches = pinched_E_surface_keys
                ellipsoids, ellipsoid_param, delta_param, expansion_threshold = loop_SG.build_existing_ellipsoids(
                    cvxpy_source_coordinates, pinched_E_surfaces=consider_pinches, extra_info=extra_info,allow_for_zero_shifts=True)
                prop_id_to_name = {
                    (ll_index, p_index) : p.name for ll_index, ll in enumerate(loop_SG.loop_lines) for p_index, p in enumerate(ll.propagators)
                }
                delta_param_per_prop = {}
                for i_ll, ll in enumerate(loop_SG.loop_lines):
                    for i_p, p in enumerate(ll.propagators):
                        delta_param_per_prop[prop_id_to_name[(i_ll,i_p)]] = delta_param[len(delta_param_per_prop)]

                E_surfaces = []
                for E_surface_key, expression in ellipsoids.items():
                    E_surfaces.append({
                        'E_surface_key' : E_surface_key,
                        'n_loops' : len(E_surface_key)-1,
                        'onshell_propagators' : sorted( [
                            {
                                'name': prop_id_to_name[os[0]],
                                'square_root_sign' : os[1],
                                'energy_shift_sign' : os[2],
                                'loop_sig' : loop_SG.loop_lines[os[0][0]].propagators[os[0][1]].signature,
                                'v_shift'  : loop_SG.loop_lines[os[0][0]].propagators[os[0][1]].q.space(),
                                'm_squared': loop_SG.loop_lines[os[0][0]].propagators[os[0][1]].m_squared
                            }
                            for os in E_surface_key], 
                                key=lambda el: el['name'] if re.match('^pq\d+$','pq1743') is None else int(el['name'][2:]) ),
                        'cxpy_expression' : expression,
                        'ellipsoid_param' : ellipsoid_param[E_surface_key],
                        'pinched' : (E_surface_key in pinched_E_surface_keys),
                        'E_shift' : sum( os[2]*loop_SG.loop_lines[os[0][0]].propagators[os[0][1]].q[0] for os in E_surface_key )
                    })
                E_surfaces.sort(key=lambda e: (
                    e['n_loops'],
                    (not e['pinched']),
                    tuple([os['name'] for os in e['onshell_propagators']]),
                    tuple([os['square_root_sign'] for os in e['onshell_propagators']]), 
                    tuple([os['energy_shift_sign'] for os in e['onshell_propagators']])
                ))
                for i_surf, E_surf in enumerate(E_surfaces):
                    E_surf['id'] = i_surf

                if args.selected_e_surfaces is not None:
                    E_surfaces = [ E_surf for E_surf in E_surfaces if set([os['name'] for os in E_surf["onshell_propagators"]]) in args.selected_e_surfaces ]
                    if len(E_surfaces)==0:
                        raise alphaLoopInvalidRunCmd("No E-surface found within the selected list: %s"%str(args.selected_e_surfaces))

                # We want to work in the convention were all E-surface square root signs are positive
                assert(all(all(os['square_root_sign']==1 for os in E_surf['onshell_propagators']) for E_surf in E_surfaces))

                # Assert that if some E surfaces are pinched then all internal masses are zero and the E surface shift is zero too.
                for E_surf in E_surfaces:
                    if not E_surf['pinched']:
                        continue
                    assert(E_surf['E_shift']==0.) 
                    assert(all(osp['m_squared']==0 for osp in E_surf['onshell_propagators']))
                    assert(all(list(osp['v_shift'])==[0.,0.,0.] for osp in E_surf['onshell_propagators']))

                IR_info_per_SG_and_E_surfaces_set[SG_name]['E_surfaces']=E_surfaces
                if args.verbose:
                    logger.info("All %d existing E-surfaces from supergraph %s:"%(len(E_surfaces),SG_name))
                    for i_surf, E_surface in enumerate(E_surfaces):
                        logger.info("#%-3d: %s%d%s-loop E-surface %s : %s"%(E_surface['id'], Colours.GREEN,E_surface['n_loops'],Colours.END, '(pinched)' if  E_surface['pinched'] else ' '*9,
                            ', '.join('%-21s'%('%s%s%-4s%s'%('%s%s%s'%(Colours.GREEN,'+',Colours.END) if term['energy_shift_sign']> 0 else '%s%s%s'%(Colours.RED,'-',Colours.END),  
                                Colours.BLUE, term['name'], Colours.END) ) 
                            for term in  E_surface['onshell_propagators'])))

                # Write all E-surface specifications to a mathematica output
                if args.mathematica:
                    with open(pjoin(self.dir_path, self._run_workspace_folder, "E_surfaces.m"),'w') as f:
                        f.write('<|\n')
                        def MMform(fl): return ('%.16e'%fl).replace('e','*^')
                        for i_surf, E_surf in enumerate(E_surfaces):
                            f.write('%d->%s%s\n'%(
                                E_surf['id'],
                                '<|"E_shift" -> %s, "onshell_propagators"->{%s}|>'%(
                                    MMform(E_surf['E_shift']),
                                    ','.join(
                                        '<|"loop_sig"->%s,"v_shift"->%s,"m_squared"->%s|>'%(
                                            '{%s}'%(','.join( '%d'%sig for sig in osp['loop_sig'] )),
                                            '{%s}'%(','.join( MMform(vs) for vs in osp['v_shift'] )),
                                            MMform(osp['m_squared']) 
                                        ) for osp in E_surf['onshell_propagators']
                                    )
                                ),
                                ',' if i_surf!=(len(E_surfaces)-1) else ''
                            ))
                        f.write('|>')

                IR_info_per_SG_and_E_surfaces_set[SG_name]['E_surfaces_intersection']={}
                E_surfaces_combinations = []
                if not args.intersections:
                    for n_E_surf_in_intersection in range(2,args.max_E_surfaces_in_intersections+1):
                        for E_surfaces_combination in itertools.combinations(range(0,len(E_surfaces)),n_E_surf_in_intersection):
                            E_surfaces_combinations.append(E_surfaces_combination)
                else:
                    E_surfaces_combinations = args.intersections

                for n_E_surf_in_intersection in range(2,max(max(len(comb) for comb in E_surfaces_combinations),args.max_E_surfaces_in_intersections)+1):
                    IR_info_per_SG_and_E_surfaces_set[SG_name]['E_surfaces_intersection'][n_E_surf_in_intersection] = {}

                bar.update(n_comb=len(E_surfaces_combinations))

                if args.intersection_point is not None:
                    if len(E_surfaces_combinations)!=1:
                        raise alphaLoopInvalidRunCmd("An intersection point can only be specified if there is a single interesection to sample.")
                    E_surfaces_combination_with_id = tuple(sorted([E_surfaces[E_surf_index]['id'] for E_surf_index in E_surfaces_combinations[0]]))
                    IR_info_per_SG_and_E_surfaces_set[SG_name]['E_surfaces_intersection'][len(E_surfaces_combination_with_id)][E_surfaces_combination_with_id] = {
                        'intersection_point' : [ args.intersection_point[i:i+3] for i in range(0,len(args.intersection_point),3) ]
                    }
                else:
                    for i_comb, E_surfaces_combination in enumerate(E_surfaces_combinations):
                        bar.update(i_comb=i_comb)
                        bar.update(intersection=str(E_surfaces_combination))
                        #if set(E_surfaces_combination) in [{0,1},{0,2}]: continue
                        # The finder is really super verbose, so edit the line below if you really want its debug output
                        finder_verbosity = args.verbose and False
                        a_finder = EsurfaceIntersectionFinder([E_surfaces[E_surf_id] for E_surf_id in E_surfaces_combination],cvxpy_source_coordinates, 
                                                                E_cm, debug=finder_verbosity, seed_point_shifts=args.n_shifts_to_test_for_finding_intersection)
                        intersection_point = a_finder.find_intersection()
                        if intersection_point is not None:
                            n_intersections_found += 1
                            bar.update(inter_found=n_intersections_found)
                            E_surfaces_combination_with_id = tuple(sorted([E_surfaces[E_surf_index]['id'] for E_surf_index in E_surfaces_combination]))
                            IR_info_per_SG_and_E_surfaces_set[SG_name]['E_surfaces_intersection'][len(E_surfaces_combination_with_id)][E_surfaces_combination_with_id] = {
                                'intersection_point' : intersection_point
                            }
                        else:
                            n_intersections_rejected += 1
                            bar.update(inter_failed=n_intersections_rejected)

                logger.info("The IR profiler found a total %d(=%s) intersections for SG %s to inspect."%(
                    sum(len(v) for v in IR_info_per_SG_and_E_surfaces_set[SG_name]['E_surfaces_intersection'].values()),
                    '+'.join(
                        '%d'%len(IR_info_per_SG_and_E_surfaces_set[SG_name]['E_surfaces_intersection'][len_comb]) 
                        for len_comb in range(2,args.max_E_surfaces_in_intersections+1)
                    ),
                    SG_name
                ))
                for len_comb in range(2,args.max_E_surfaces_in_intersections+1):
                    logger.info("%d combinations of %d-E_surface intersecting: %s"%(
                        len(IR_info_per_SG_and_E_surfaces_set[SG_name]['E_surfaces_intersection'][len_comb]),
                        len_comb,', '.join('-'.join('%d'%E_id for E_id in comb) for comb in 
                        sorted(IR_info_per_SG_and_E_surfaces_set[SG_name]['E_surfaces_intersection'][len_comb])
                    )))

                if args.selected_e_surfaces is None or any(key not in SG for key in ['E_surfaces','E_surfaces_intersection']):
                    # Save the completed preprocessing
                    SG['E_surfaces'] = IR_info_per_SG_and_E_surfaces_set[SG_name]['E_surfaces']
                    SG['E_surfaces_intersection'] = IR_info_per_SG_and_E_surfaces_set[SG_name]['E_surfaces_intersection']

                    logger.info("Writing out processed yaml for supergaph '%s' on disk..."%SG_name)
                    SG.export(SG_name, pjoin(self.dir_path, self._rust_inputs_folder))


        # Compute the log-spaced sequence of rescaling
        scalings = [ 10.**((math.log10(args.min_scaling)+i*((math.log10(args.max_scaling)-math.log10(args.min_scaling))/(args.n_points-1))))
                        for i in range(args.n_points) ]

        # WARNING it is important that the rust workers instantiated only go out of scope when this function terminates
        rust_workers = {SG_name: self.get_rust_worker(SG_name) for SG_name in selected_SGs}
        if args.f128 or not args.no_f128:
            hyperparameters_backup=copy.deepcopy(self.hyperparameters)
            for entry in self.hyperparameters['General']['stability_checks']:
                entry['use_f128'] = True
            rust_workers_f128 = {SG_name: self.get_rust_worker(SG_name) for SG_name in selected_SGs}
            self.hyperparameters = hyperparameters_backup
        t_start_profile = time.time()
        n_passed = 0
        n_failed = 0
        n_fit_failed = 0
        SG_passed = 0
        SG_failed = 0

        with progressbar.ProgressBar(
                prefix=("IR profiling. E-surface intersection: {variables.intersection} {variables.i_comb}/{variables.n_comb} {variables.passed}\u2713, {variables.failed}\u2717, SGs: {variables.SG_passed}\u2713, {variables.SG_failed}\u2717, SG: {variables.SG_name} "), 
                max_value=len(selected_SGs),variables={
                    'intersection': 'N/A', 'passed': 0, 'failed': 0, 'SG_name': 'N/A', 'i_comb': 0, 'n_comb': 0, 'SG_passed': 0, 'SG_failed': 0
                }
            ) as bar:

            for i_SG, SG_name in enumerate(selected_SGs):
                
                E_surfaces = IR_info_per_SG_and_E_surfaces_set[SG_name]['E_surfaces']
                SG['E_surfaces_analysis'] = E_surfaces
                SG['E_surfaces_intersection_analysis'] = {}

                E_surfaces_intersection_analysis = IR_info_per_SG_and_E_surfaces_set[SG_name]['E_surfaces_intersection']
            
                if args.seed != 0:
                    random.seed(args.seed)

                skip_furhter_tests_in_this_SG = False
                this_SG_failed = False

                E_cm = SG.get_E_cm(self.hyperparameters)

                if args.approach_direction is None:
                    approach_direction = [ Vector([random.random()*E_cm for i_comp in range(0,3)]) for i_vec in range(0, SG['n_loops']) ]
                else:                    
                    approach_direction = [ Vector(list(args.approach_direction[i:i+3])) for i in range(0,len(args.approach_direction),3) ]
                    if len(approach_direction)!=SG['n_loops']:
                        raise alphaLoopInvalidRunCmd("The specified approach direction does not specify %d*3 components."%SG['n_loops'])

                bar.update(SG_name=SG_name)
                bar.update(i_SG)
                bar.update(n_comb=sum(len(E_surfaces_intersection_analysis[len_comb]) for len_comb in E_surfaces_intersection_analysis))
                
                E_surface_ID_to_E_surface = { E_surf['id'] : E_surf for E_surf in E_surfaces }

                rust_worker = rust_workers[SG_name]
                if args.f128 or not args.no_f128:
                    rust_worker_f128 = rust_workers_f128[SG_name]

                i_test = 0
                #for len_comb in range(2,args.max_E_surfaces_in_intersections+1):
                #    for E_surface_combination, comb_info in SG['E_surfaces_intersection_analysis'].get(len_comb,{}).items():
                all_tests = sum([list(E_surfaces_intersection_analysis[len_comb].items()) for len_comb in E_surfaces_intersection_analysis],[])
                for E_surface_combination, comb_info in all_tests:

                    if len(E_surface_combination) not in SG['E_surfaces_intersection_analysis']:
                        SG['E_surfaces_intersection_analysis'][len(E_surface_combination)] = {}
                    SG['E_surfaces_intersection_analysis'][len(E_surface_combination)][tuple(E_surface_combination)] = {}

                    if skip_furhter_tests_in_this_SG:
                        continue

                    i_test += 1
                    bar.update(i_comb=i_test)
                    bar.update(intersection=str(E_surface_combination))

                    SG['E_surfaces_intersection_analysis'][len(E_surface_combination)][tuple(E_surface_combination)]['intersection_point'] = comb_info['intersection_point']
                    intersection_point = [ Vector(v) for v in comb_info['intersection_point'] ]
                    momenta_str = '\n'.join('\n'.join('%-5s : %-20s'%(osp['name'], 
                            ', '.join( '%.10e'%vi for vi in list(sum( l*factor for l, factor in zip(intersection_point,osp['loop_sig']) )+Vector(osp['v_shift'])) )
                        ) for osp in E_surface_ID_to_E_surface[E_surf_ID]['onshell_propagators']
                    ) for E_surf_ID in E_surface_combination)

                    if args.verbose:
                        logger.info("Now studying interesection of E-surfaces %s == %s for SG %s with the following intersection point:\n%sand momenta:\n%s"%(
                            str(E_surface_combination),
                            '^'.join('dE(%s)'%(','.join(osp['name'] for osp in E_surface_ID_to_E_surface[E_surf_ID]['onshell_propagators'])) for E_surf_ID in E_surface_combination),
                            SG_name,
                            str(intersection_point),
                            momenta_str
                        ))

                    # Cut_ID None means running over the full supergraph
                    for cut_ID, cuts_info in [(None, None)]+list(enumerate(SG['cutkosky_cuts'])):
                        
                        if cut_ID is not None:
                            # Skip Cutkosky cuts not matching any of the specified thresholds
                            if args.only_relevant_cuts and not any( 
                                set(os['name'] for os in E_surface_ID_to_E_surface[E_surf_id]['onshell_propagators'])==set([cut['name'] for cut in cuts_info['cuts']])
                                for E_surf_id in E_surface_combination):
                                continue

                        use_f128 = args.f128
                        while True:
                            results = []
                            for scaling in scalings:
                                rescaled_momenta = [ v+approach_direction[i_v]*scaling for i_v, v in enumerate(intersection_point)]
                                #misc.sprint(rescaled_momenta)

                                # Now map these momenta in the defining LMB into x variables in the unit hypercube
                                xs_in_defining_LMB = []
                                overall_jac = 1.0
                                for i_k, k in enumerate(rescaled_momenta):
                                    # This is cheap, always do in f128
                                    if use_f128 or True:
                                        x1, x2, x3, jac = rust_worker.inv_parameterize_f128( list(k), i_k, E_cm**2)
                                    else:
                                        x1, x2, x3, jac = rust_worker.inv_parameterize( list(k), i_k, E_cm**2)
                                    xs_in_defining_LMB.extend([x1, x2, x3])
                                    overall_jac *= jac

                                if cut_ID is None:

        #                            for k in rescaled_momenta_in_defining_LMB:
        #                                misc.sprint(['%.16f'%k_i for k_i in k])
        #                            misc.sprint(xs_in_defining_LMB,E_cm)
                                    #with utils.Silence(active=(not args.show_rust_warnings)):
                                    with utils.suppress_output(active=(not args.show_rust_warnings)):
                                        if use_f128:
                                            res_re, res_im = rust_worker_f128.evaluate_integrand(xs_in_defining_LMB)
                                        else:
                                            res_re, res_im = rust_worker.evaluate_integrand(xs_in_defining_LMB)
                                    # We do *not* want to include the inverse jacobian in this case here, so do *not* multiply by overall_jac
                                    results.append( (scaling, complex(res_re, res_im) ) )

                                else:

                                    # Now obtain the rescaling for these momenta
                                    LU_scaling_solutions = rust_worker.get_scaling(rescaled_momenta,cut_ID)
                                    if LU_scaling_solutions is None or len(LU_scaling_solutions)==0 or all(LU_scaling[0]<0. for LU_scaling in LU_scaling_solutions):
                                        if args.show_warnings:
                                            logger.warning("Could not find rescaling for IR profiling of SG '%s' and cut #%d for the E_surface intersection %s : %s\nInput LMB momenta: %s"%(
                                                SG_name, cut_ID, str(E_surface_combination), str(LU_scaling_solutions), str(rescaled_momenta) ))
                                        continue
                                    LU_scaling_solutions = list(LU_scaling_solutions)
                                    LU_scaling, LU_scaling_jacobian = LU_scaling_solutions.pop(0)
                                    if LU_scaling>0.0 and args.show_warnings:
                                        logger.warning("Found unexpected rescaling solutions for IR profiling of SG '%s' and cut #%d for the E_surface intersection %s : %s\nInput LMB momenta: %s"%(
                                                SG_name, cut_ID, str(E_surface_combination), str(LU_scaling_solutions), str(rescaled_momenta) ))

                                    while LU_scaling < 0.0:
                                        if len(LU_scaling_solutions)==0:
                                            break
                                        LU_scaling, LU_scaling_jacobian = LU_scaling_solutions.pop(0)

                                    with utils.suppress_output(active=(not args.show_rust_warnings)):
                                        if use_f128:
                                            res_re, res_im = rust_worker.evaluate_cut_f128(rescaled_momenta,cut_ID,LU_scaling,LU_scaling_jacobian)                            
                                        else:
                                            res_re, res_im = rust_worker.evaluate_cut(rescaled_momenta,cut_ID,LU_scaling,LU_scaling_jacobian)
                                    results.append( (scaling, complex(res_re, res_im)/overall_jac ) )

    #                        misc.sprint(results)
                            # Here we are in x-space, so the dod read is already the one we want.
                            dod, standard_error, number_of_points_considered, successful_fit = utils.compute_dod(results, threshold=0.2)
                            # Flip sign since we approach the singularity with a decreasing scaling
                            dod *= -1.
                            max_result = max( (r for r in results) , key=lambda el: abs(el[1]) )

                            if cut_ID is not None:
                                if args.verbose:
                                    logger.info("IR scaling detected for cut #%d: %.3f +/- %.4f over %d points (max: %s)."%(
                                        cut_ID, dod, standard_error, number_of_points_considered, str(max_result)
                                    ))
                                if 'cut_results' not in SG['E_surfaces_intersection_analysis'][len(E_surface_combination)][tuple(E_surface_combination)]:
                                    SG['E_surfaces_intersection_analysis'][len(E_surface_combination)][tuple(E_surface_combination)]['cut_results'] = {}
                                if cut_ID not in SG['E_surfaces_intersection_analysis'][len(E_surface_combination)][tuple(E_surface_combination)]['cut_results']:
                                    SG['E_surfaces_intersection_analysis'][len(E_surface_combination)][tuple(E_surface_combination)]['cut_results'][cut_ID] = {}
                                container = SG['E_surfaces_intersection_analysis'][len(E_surface_combination)][tuple(E_surface_combination)]['cut_results'][cut_ID]
                                test_passed = (dod < float(args.target_scaling)+min(max(10.0*abs(standard_error),0.005),0.2) )
                                container['dod_computed'] = (dod, standard_error, test_passed)
                                container['max_result'] = (max_result[0], (max_result[1].real, max_result[1].imag))
                                break

                            # Support for IR-safety isolation cuts:
                            if all(r==complex(0.,0.) for r in results[-3:]):
                                # Flag this divergence as being cutaway with dod=None and set the test as passed:
                                dod = None
                                test_passed = True
                                do_debug = False
                                break

                            # We expect dod of at most five sigma above 0.0 for the integral to be convergent.
                            test_passed = (dod < float(args.target_scaling)+min(max(10.0*abs(standard_error),0.005),0.2) )

                            if (successful_fit and test_passed) or use_f128:
                                break
                            elif not args.no_f128:
                                use_f128 = True
                            else:
                                break
                        
                        if cut_ID is not None:
                            continue

                        do_debug = False
                        if not successful_fit and dod > float(args.target_scaling)-1.:
                            if args.show_warnings:
                                logger.critical("The fit for the IR scaling of SG '%s' for the E_surface intersection %s is unsuccessful (unstable). Found: %.3f +/- %.4f over %d points."%(
                                    SG_name, str(E_surface_combination), dod, standard_error, number_of_points_considered
                                ))
                                do_debug = True
                            n_fit_failed += 1
                            #bar.update(fit_failed=n_fit_failed)

                        #logger.info("For SG '%s' and UV edges %s and fixed edges %s, dod measured is: %.3f +/- %.4f over %d points"%(
                        #    SG_name, UV_edges_str, fixed_edges_str, dod, standard_error, number_of_points_considered
                        #))
                        if (args.verbose or do_debug) or (not test_passed and args.show_fails):
                            logger.info('%s : IR profile of SG %s with intersection %s. Intersection point:\n%s\nand momenta:\n%s\nThis can be rerun with option -e_surfaces %s -intersections %s -intersection_point %s -approach_direction %s'%(
                                '%sPASS%s'%(Colours.GREEN, Colours.END) if test_passed else '%sFAIL%s'%(Colours.RED, Colours.END), SG_name, 
                                str(E_surface_combination), str(intersection_point),momenta_str,
                                ' '.join('(%s)'%(','.join(
                                    '"%s"'%os['name'] for os in E_surface_ID_to_E_surface[E_surf_id]['onshell_propagators']
                                )) for E_surf_id in E_surface_combination),
                                '(%s)'%(','.join('%d'%i_surf for i_surf in range(0, len(E_surface_combination)))),
                                '(%s)'%(','.join(','.join('%.16e'%v_i for v_i in list(v) ) for v in intersection_point)),
                                '(%s)'%(','.join(','.join('%.16e'%v_i for v_i in list(v) ) for v in approach_direction))
                            ))
                            if args.verbose or do_debug:
                                logger.info('\n'+'\n'.join('%-13.5e -> %-13.5e'%(
                                    r[0], abs(r[1])) for i_r, r in enumerate(results)))
                            logger.info('%sdod= %.3f +/- %.3f%s (max: %s)'%(Colours.GREEN if test_passed else Colours.RED, dod, standard_error, Colours.END, str(max_result)))

                        if test_passed:
                            n_passed += 1
                        else:
                            n_failed += 1
                        bar.update(passed=n_passed)
                        bar.update(failed=n_failed)
                        
                        SG['E_surfaces_intersection_analysis'][len(E_surface_combination)][tuple(E_surface_combination)]['dod_computed'] = (dod, standard_error,test_passed)
                        SG['E_surfaces_intersection_analysis'][len(E_surface_combination)][tuple(E_surface_combination)]['max_result'] = (max_result[0], (max_result[1].real, max_result[1].imag))

                        if not test_passed and args.skip_once_failed:
                            skip_furhter_tests_in_this_SG = True
                        if not test_passed and not this_SG_failed:
                            SG_failed += 1
                            bar.update(SG_failed=SG_failed)
                            this_SG_failed = True

                if not this_SG_failed:
                    SG_passed += 1
                    bar.update(SG_passed=SG_passed)

        delta_t = time.time()-t_start_profile
        logger.info("IR profile completed in %d [s]: %s%d%s tests passed and %s failed (and %s fit failed)."%(
            int(delta_t),
            Colours.GREEN,
            n_passed,
            Colours.END,
            'none' if n_failed == 0 else '%s%d%s'%(Colours.RED, n_failed, Colours.END),
            'no' if n_fit_failed == 0 else '%s%d%s'%(Colours.RED, n_fit_failed, Colours.END)
        ))

        logger.info("Writing out processed yaml supergaphs on disk...")
        # Write out the results into processed topologies
        self.all_supergraphs.export(pjoin(self.dir_path, self._rust_inputs_folder))
        if len(selected_SGs)==1:
            self.do_display('%s --ir%s'%(selected_SGs[0],' --show_momenta' if args.show_momenta else ''))
        else:
            self.do_display('--ir%s'%(' --show_momenta' if args.show_momenta else ''))

    #### UV PROFILE COMMAND
    uv_profile_parser = ArgumentParser(prog='uv_profile')
    uv_profile_parser.add_argument('SG_name', metavar='SG_name', type=str, nargs='?',
                    help='the name of a supergraph to display')
    uv_profile_parser.add_argument("-n","--n_points", dest='n_points', type=int, default=20,
                    help='force a certain number of points to be considered for the uv profile')
    uv_profile_parser.add_argument("-max","--max_scaling", dest='max_scaling', type=float, default=1.0e5,
                    help='maximum UV scaling to consider')
    uv_profile_parser.add_argument("-min","--min_scaling", dest='min_scaling', type=float, default=1.0e3,
                    help='minimum UV scaling to consider')
    uv_profile_parser.add_argument("-rp","--required_precision", dest='required_precision', type=float, default=None,
                    help='minimum required relative precision for returning a result.')
    uv_profile_parser.add_argument("-s","--seed", dest='seed', type=int, default=0,
                    help='specify random seed')
    uv_profile_parser.add_argument("-t","--target_scaling", dest='target_scaling', type=int, default=-1,
                    help='set target UV scaling (default=0)')
    uv_profile_parser.add_argument("-hp","--h_power", dest='h_power', type=int, default=7,
                    help='h function dampening power')
    uv_profile_parser.add_argument("--LMB", dest='LMB', type=str, nargs='*', default=None,
                    help='set LMB to consider')
    uv_profile_parser.add_argument("--UV_indices", dest='UV_indices', type=int, nargs='*', default=None,
                    help='set UV indices to consider')
    uv_profile_parser.add_argument(
        "-f", "--f128", action="store_true", dest="f128", default=False,
        help="Perfom the UV profile using f128 arithmetics.")
    uv_profile_parser.add_argument(
        "-nf", "--no_f128", action="store_true", dest="no_f128", default=False,
        help="Forbid automatic promotion to f128.")
    uv_profile_parser.add_argument(
        "-nw", "--no_warnings", action="store_false", dest="show_warnings", default=True,
        help="Do not show warnings about this profiling run.")
    uv_profile_parser.add_argument(
        "-srw", "--show_rust_warnings", action="store_true", dest="show_rust_warnings", default=False,
        help="Show rust warnings.")
    uv_profile_parser.add_argument(
        "-nsof", "--no_skip_once_failed", action="store_false", dest="skip_once_failed", default=True,
        help="Do not skip the probing of a supergraph once it failed.")
    uv_profile_parser.add_argument(
        "-nsf", "--no_show_fails", action="store_false", dest="show_fails", default=True,
        help="Show exhaustive information for each fail.")
    uv_profile_parser.add_argument(
        "-sc", "--scale_cuts", action="store_true", dest="scale_cuts", default=False,
        help="Include UV scaling of edges being cut.")
    uv_profile_parser.add_argument("-n_max","--n_max", dest='n_max', type=int, default=-1,
                    help='Set the maximum number of UV test to perform per cut (default: all)')
    uv_profile_parser.add_argument(
        "-v", "--verbose", action="store_true", dest="verbose", default=False,
        help="Enable verbose output.")
    def help_uv_profile(self):
        self.uv_profile_parser.print_help()
        return
    # We must wrape this function in a process because of the border effects of the pyO3 rust Python bindings
    @wrap_in_process()
    @with_tmp_hyperparameters({
        'Integrator.dashboard': False,
        'General.minimal_precision_for_returning_result': 1.0,
        'CrossSection.NormalisingFunction.name'         : 'left_right_exponential',
        'CrossSection.NormalisingFunction.center'       : 1.0,
        'CrossSection.NormalisingFunction.spread'       : 1.0,
        'General.multi_channeling'                      : False
    })
    def do_uv_profile(self, line):
        """ Automatically probe all UV limits of a process output."""

        if line=='help':
            self.uv_profile_parser.print_help()
            return 

        args = self.split_arg(line)
        args = self.uv_profile_parser.parse_args(args)
        
        if args.SG_name is None:
            selected_SGs = list(self.all_supergraphs.keys())
        else:
            selected_SGs = [args.SG_name,]

        #self.hyperparameters['CrossSection']['NormalisingFunction']['spread'] = 1. # args.h_power
        if args.required_precision is None:
            self.hyperparameters['General']['stability_checks'][-1]['relative_precision']=1.0e-99
        else:
            for entry in self.hyperparameters['General']['stability_checks']:
                entry['relative_precision'] = args.required_precision

        logger.info("Starting UV profile...")

        # Compute the log-spaced sequence of rescaling
        scalings = [ 10.**(math.log10(args.min_scaling)+i*((math.log10(args.max_scaling)-math.log10(args.min_scaling))/(args.n_points-1))) 
                        for i in range(args.n_points) ]

        # We need to detect here if we are in the amplitude-mock-up situation with frozen external momenta.
        frozen_momenta = None
        if 'external_data' in self.cross_section_set:
            frozen_momenta = {
                'in' : self.cross_section_set['external_data']['in_momenta'],
                'out' : self.cross_section_set['external_data']['out_momenta'],
            }
            # Also force the specified incoming momenta specified in the hyperparameters to match the frozen specified ones.
            self.hyperparameters.set_parameter('CrossSection.incoming_momenta',frozen_momenta['in'])
            self.hyperparameters.set_parameter('CrossSection.do_rescaling',False)
            self.hyperparameters.set_parameter('CrossSection.fixed_cut_momenta',frozen_momenta['out'])

        # Built external momenta. Remember that they appear twice.
        external_momenta = [ Vector(v[1:]) for v in self.hyperparameters['CrossSection']['incoming_momenta'] ]
        external_momenta.extend(external_momenta)

        # Prepare the run
        UV_info_per_SG_and_cut = {}
        for i_SG, SG_name in enumerate(selected_SGs):

            if args.seed != 0:
                random.seed(args.seed)

            UV_info_per_SG_and_cut[SG_name] = {}
            SG = self.all_supergraphs[SG_name]
            # First we must regenerate a TopologyGenerator instance for this supergraph
            edges_list = SG['topo_edges']
            SG_topo_gen = ltd_utils.TopologyGenerator(
                [e[:-1] for e in edges_list],
                powers = { e[0] : e[-1] for e in edges_list }
            )
            edge_signatures = SG['edge_signatures']

            if args.LMB is not None and args.UV_indices is not None:
                all_LMBs = [tuple(edge_name for edge_name in args.LMB)]
            else:
                # Store the signature of each edge part of the LMBs w.r.t the defining LMB
                all_LMBs = [ tuple(edges_list[i_edge][0] for i_edge in lmb) for lmb in SG_topo_gen.loop_momentum_bases()]
                # Filter out all LMBs with the same loop momenta signatures
                new_all_LMBS = []
                signatures_encountered = []
                for LMB in all_LMBs:
                    this_LMB_sig = sorted(edge_signatures[edge_name][0] for edge_name in LMB)
                    if this_LMB_sig not in signatures_encountered:
                        signatures_encountered.append(this_LMB_sig)
                        new_all_LMBS.append(LMB)
                all_LMBs = new_all_LMBS

            all_LMBs_and_UV_edge_indices_and_signatures = []
            # Then only include combination of UV_edges that yield different signature combinations
            UV_info_per_SG_and_cut[SG_name]['LMBs_to_defining_LMB_transfo'] = {}
            # We must keep *all* signatures actually as it matters which other ones we keep fixed.
            #signature_combinations_considered = []
            for LMB in all_LMBs:
                # construct the cut basis to LTD loop momentum basis mapping
                mat = [edge_signatures[edge_name][0] for edge_name in LMB]
                transfo = np.linalg.inv(np.array(mat))
                shifts = [[-p for p in edge_signatures[edge_name][1]] for edge_name in LMB]
                UV_info_per_SG_and_cut[SG_name]['LMBs_to_defining_LMB_transfo'][LMB] = (transfo, shifts)
                
                if args.LMB is not None and args.UV_indices is not None:
                    indices_list = [args.UV_indices,]
                else:
                    indices_list = itertools.chain(*[
                        itertools.combinations(list(range(len(LMB))),n_UV) for n_UV in range(1,len(LMB)+1)
                    ])
                for UV_edge_indices in indices_list:
                    UV_edge_signatures = tuple([edge_signatures[LMB[UV_edge_index]][0] for UV_edge_index in UV_edge_indices])
                    # We must keep *all* signatures actually as it matters which other ones we keep fixed.
                    #sorted_UV_edge_signatures = sorted(list(UV_edge_signatures))
                    #if sorted_UV_edge_signatures in signature_combinations_considered:
                    #    continue
                    #signature_combinations_considered.append(sorted_UV_edge_signatures)
                    all_LMBs_and_UV_edge_indices_and_signatures.append((LMB, UV_edge_indices, UV_edge_signatures))

            UV_info_per_SG_and_cut[SG_name]['UV_edges_to_probe_for_cuts'] = []
            UV_info_per_SG_and_cut[SG_name]['UV_edges_to_probe'] = []
            uv_probes_to_add_for_this_SG = []
            uv_probes_to_add_for_this_cut = {}
            # Combination of loop momenta signatures considered
            for LMB, UV_edge_indices, UV_edge_signatures in all_LMBs_and_UV_edge_indices_and_signatures:
                added_UV_edges_to_probe=False
                for cut_ID, cuts_info in enumerate(SG['cutkosky_cuts']):
                    if (args.LMB is None or args.UV_indices is None) and (not args.scale_cuts):
                        # Rescaling in the UV a momenta part of the cuts will never yield a singularity
                        # First do a simple heuristics test
                        if ( 
                            any( cut['signature'][0] in UV_edge_signatures for cut in cuts_info['cuts'] ) or
                            any( [-c for c in cut['signature'][0]] in UV_edge_signatures for cut in cuts_info['cuts'] )
                        ):
                            continue
                        # And do an exact test otherwise
                        transfo, parametric_shifts = UV_info_per_SG_and_cut[SG_name]['LMBs_to_defining_LMB_transfo'][LMB]
                        prime_components = [1,3,5,7,11,13,17,19,23]
                        UV_momenta = [ prime_components[UV_index] if UV_index in UV_edge_indices else 0 for UV_index in range(len(LMB)) ]
                        UV_momenta_in_defining_LMB = transfo.dot(UV_momenta)
                        if any( UV_momenta_in_defining_LMB.dot(cut['signature'][0])!=0 for cut in cuts_info['cuts'] ):
                            continue

#                    UV_info_per_SG_and_cut[SG_name]['UV_edges_to_probe_for_cuts'].append((cut_ID, LMB, UV_edge_indices))
                    if cut_ID in uv_probes_to_add_for_this_cut:
                        uv_probes_to_add_for_this_cut[cut_ID].append((cut_ID, LMB, UV_edge_indices))
                    else:
                        uv_probes_to_add_for_this_cut[cut_ID]=[(cut_ID, LMB, UV_edge_indices),]
                    if not added_UV_edges_to_probe:
#                        UV_info_per_SG_and_cut[SG_name]['UV_edges_to_probe'].append((LMB, UV_edge_indices))
                        uv_probes_to_add_for_this_SG.append((LMB, UV_edge_indices))
                        added_UV_edges_to_probe = True

            if args.n_max <= 0:
                if frozen_momenta is None:
                    for cut_ID in sorted(list(uv_probes_to_add_for_this_cut.keys())):
                        UV_info_per_SG_and_cut[SG_name]['UV_edges_to_probe_for_cuts'].extend(uv_probes_to_add_for_this_cut[cut_ID])
                UV_info_per_SG_and_cut[SG_name]['UV_edges_to_probe'].extend(uv_probes_to_add_for_this_SG)
            else:
                if frozen_momenta is None:
                    for cut_ID in sorted(list(uv_probes_to_add_for_this_cut.keys())):
                        if len(uv_probes_to_add_for_this_cut[cut_ID])<=args.n_max:
                            UV_info_per_SG_and_cut[SG_name]['UV_edges_to_probe_for_cuts'].extend(uv_probes_to_add_for_this_cut[cut_ID])
                        else:
                            UV_info_per_SG_and_cut[SG_name]['UV_edges_to_probe_for_cuts'].extend(random.sample(uv_probes_to_add_for_this_cut[cut_ID],args.n_max))                        
                if len(uv_probes_to_add_for_this_SG) <= args.n_max:
                    UV_info_per_SG_and_cut[SG_name]['UV_edges_to_probe'].extend(uv_probes_to_add_for_this_SG)
                else:
                    UV_info_per_SG_and_cut[SG_name]['UV_edges_to_probe'].extend(random.sample(uv_probes_to_add_for_this_SG,args.n_max))

        max_count = sum( len(UV_info_per_SG_and_cut[SG_name]['UV_edges_to_probe']) + 
                         len(UV_info_per_SG_and_cut[SG_name]['UV_edges_to_probe_for_cuts']) for SG_name in selected_SGs )

        # WARNING it is important that the rust workers instantiated only go out of scope when this function terminates
        rust_workers = {SG_name: self.get_rust_worker(SG_name) for SG_name in selected_SGs}
        if args.f128 or not args.no_f128:
            hyperparameters_backup=copy.deepcopy(self.hyperparameters)
            for entry in self.hyperparameters['General']['stability_checks']:
                entry['use_f128'] = True
            rust_workers_f128 = {SG_name: self.get_rust_worker(SG_name) for SG_name in selected_SGs}
            self.hyperparameters = hyperparameters_backup
        t_start_profile = time.time()
        n_passed = 0
        n_failed = 0
        n_fit_failed = 0
        SG_passed = 0
        SG_failed = 0
        with progressbar.ProgressBar(
                prefix=("UV profile: {variables.passed}\u2713, {variables.failed}\u2717, SG: {variables.SG_passed}\u2713, {variables.SG_failed}\u2717, fit fails: {variables.fit_failed}, {variables.SG}/{variables.n_SG} ({variables.SG_name}), "+
                        "cut ID: {variables.cut}, test #: {variables.i_test}/{variables.n_tests} (UV: {variables.uv_edge_names} fixed: {variables.fixed_edge_names}) "), 
                max_value=max_count,variables={
                    'passed': n_passed, 'failed': n_failed, 'fit_failed': n_fit_failed, 'SG_passed': SG_passed, 'SG_failed': SG_failed,
                    'SG_name':selected_SGs[0], 'SG': 0, 'n_SG':len(selected_SGs),'cut':'sum', 
                    'i_test': 0, 'n_tests': len(UV_info_per_SG_and_cut[selected_SGs[0]]['UV_edges_to_probe']),
                    'uv_edge_names': 'N/A',
                    'fixed_edge_names': 'N/A'
                }
            ) as bar:

            for i_SG, SG_name in enumerate(selected_SGs):
                skip_furhter_tests_in_this_SG = False
                this_SG_failed = False

                SG = self.all_supergraphs[SG_name]
                E_cm = SG.get_E_cm(self.hyperparameters)

                bar.update(SG_name=SG_name)
                bar.update(SG=i_SG+1)
                bar.update(cut='sum')
                bar.update(n_tests=len(UV_info_per_SG_and_cut[SG_name]['UV_edges_to_probe']))
                bar.update(i_test=0)
                bar.update(uv_edge_names='N/A')

                #rust_worker = self.get_rust_worker(SG_name)
                rust_worker = rust_workers[SG_name]
                if args.f128 or not args.no_f128:
                    rust_worker_f128 = rust_workers_f128[SG_name]
                SG = self.all_supergraphs[SG_name]  

                if args.seed != 0:
                    random.seed(args.seed)
                random_momenta = SG.get_random_momenta_input(E_cm)
                random_momenta = [ Vector(v) for v in random_momenta]

                SG['DERIVED_UV_dod'] = {}
                SG['DERIVED_dod_consistency'] = (0.0, 0.0, 0.0)
                for i_test, (LMB, UV_edge_indices) in enumerate(UV_info_per_SG_and_cut[SG_name]['UV_edges_to_probe']):
                    bar.update(i_test=i_test+1)
                    UV_edges_str = ','.join(LMB[edge_index] for edge_index in UV_edge_indices)
                    bar.update(uv_edge_names=UV_edges_str)
                    fixed_edges_str = ','.join(LMB[edge_index] for edge_index in range(len(LMB)) if edge_index not in UV_edge_indices)
                    bar.update(fixed_edge_names=fixed_edges_str)
                    if skip_furhter_tests_in_this_SG:
                        bar.update(cut='SKIP')
                        bar.update(bar.value+1)
                        continue
                    if fixed_edges_str == '': fixed_edges_str = '[]'

                    use_f128 = args.f128
                    while True:
                        results = []
                        for scaling in scalings:
                            rescaled_momenta = [ v*scaling if i_v in UV_edge_indices else v for i_v, v in enumerate(random_momenta)]
                            # Now transform the momentum configuration in this LMB into the defining LMB
                            transfo, parametric_shifts = UV_info_per_SG_and_cut[SG_name]['LMBs_to_defining_LMB_transfo'][LMB]
                            shifts = [ sum([external_momenta[i_shift]*shift 
                                        for i_shift, shift in enumerate(parametric_shift)]) for parametric_shift in parametric_shifts ]
                            rescaled_momenta_in_defining_LMB = transfo.dot(
                                [list(rm+shift) for rm, shift in zip(rescaled_momenta,shifts)] )

                            # Now map these momenta in the defining LMB into x variables in the unit hypercube
                            xs_in_defining_LMB = []
                            overall_jac = 1.0
                            for i_k, k in enumerate(rescaled_momenta_in_defining_LMB):
                                # This is cheap, always do in f128
                                if use_f128 or True:
                                    x1, x2, x3, jac = rust_worker.inv_parameterize_f128( list(k), i_k, E_cm**2)
                                else:
                                    x1, x2, x3, jac = rust_worker.inv_parameterize( list(k), i_k, E_cm**2)
                                xs_in_defining_LMB.extend([x1, x2, x3])
                                overall_jac *= jac
#                            for k in rescaled_momenta_in_defining_LMB:
#                                misc.sprint(['%.16f'%k_i for k_i in k])
#                            misc.sprint(xs_in_defining_LMB,E_cm)
                            #with utils.Silence(active=(not args.show_rust_warnings)):
                            with utils.suppress_output(active=(not args.show_rust_warnings)):
                                if use_f128:
                                    res_re, res_im = rust_worker_f128.evaluate_integrand(xs_in_defining_LMB)
                                else:
                                    res_re, res_im = rust_worker.evaluate_integrand(xs_in_defining_LMB)
                            results.append( (scaling, complex(res_re, res_im)*overall_jac ) )

#                        misc.sprint(results)
                        # Here we are in x-space, so the dod read is already the one we want.
                        dod, standard_error, number_of_points_considered, successful_fit = utils.compute_dod(results)
                        dod += 3*float(len(UV_edge_indices))


                        # We expect dod of at most five sigma above 0.0 for the integral to be convergent.
                        test_passed = (dod < float(args.target_scaling)+min(max(5.0*abs(standard_error),0.005),0.1) )

                        if (successful_fit and test_passed) or use_f128:
                            break
                        elif not args.no_f128:
                            use_f128 = True
                        else:
                            break

                    do_debug = False
                    if not successful_fit and dod > float(args.target_scaling)-1.:
                        if args.show_warnings:
                            logger.critical("The fit for the UV scaling of SG '%s' with UV edges %s and fixed edges %s is unsuccessful (unstable). Found: %.3f +/- %.4f over %d points."%(
                                SG_name, UV_edges_str, fixed_edges_str, dod, standard_error, number_of_points_considered
                            ))
                            do_debug = True
                        n_fit_failed += 1
                        bar.update(fit_failed=n_fit_failed)

                    #logger.info("For SG '%s' and UV edges %s and fixed edges %s, dod measured is: %.3f +/- %.4f over %d points"%(
                    #    SG_name, UV_edges_str, fixed_edges_str, dod, standard_error, number_of_points_considered
                    #))
                    if (args.verbose or do_debug) or (not test_passed and args.show_fails):
                        logger.info("%s : UV scaling of SG '%s' with UV edges %s and fixed edges %s (--LMB %s --UV_indices %s):"%(
                            '%sPASS%s'%(Colours.GREEN, Colours.END) if test_passed else '%sFAIL%s'%(Colours.RED, Colours.END), SG_name, 
                            UV_edges_str, fixed_edges_str, ' '.join(LMB), ' '.join('%d'%uv_index for uv_index in UV_edge_indices)
                        ))
                        if args.verbose or do_debug:
                            logger.info('\n'+'\n'.join('%-13.5e -> %-13.5e'%(
                                r[0], abs(r[1])) for i_r, r in enumerate(results)))
                        logger.info('%sdod= %.3f +/- %.3f%s'%(Colours.GREEN if test_passed else Colours.RED, dod, standard_error, Colours.END))

                    if test_passed:
                        n_passed += 1
                    else:
                        n_failed += 1
                    bar.update(passed=n_passed)
                    bar.update(failed=n_failed)
                    bar.update(bar.value+1)

                    SG['DERIVED_UV_dod'][(tuple(LMB), tuple(UV_edge_indices))] = (dod, standard_error, test_passed)

                    if not test_passed and args.skip_once_failed:
                        skip_furhter_tests_in_this_SG = True
                    if not test_passed and not this_SG_failed:
                        SG_failed += 1
                        bar.update(SG_failed=SG_failed)
                        this_SG_failed = True

                bar.update(n_tests=len(UV_info_per_SG_and_cut[SG_name]['UV_edges_to_probe_for_cuts']))
                for cut in SG['cutkosky_cuts']:
                    cut['DERIVED_UV_dod'] = {}
                for i_test, (cut_ID, LMB, UV_edge_indices) in enumerate(UV_info_per_SG_and_cut[SG_name]['UV_edges_to_probe_for_cuts']):
                    bar.update(cut=cut_ID)
                    bar.update(i_test=i_test+1)
                    UV_edges_str = ','.join(LMB[edge_index] for edge_index in UV_edge_indices)
                    bar.update(uv_edge_names=UV_edges_str)
                    fixed_edges_str = ','.join(LMB[edge_index] for edge_index in range(len(LMB)) if edge_index not in UV_edge_indices)
                    bar.update(fixed_edge_names=fixed_edges_str)
                    if fixed_edges_str == '': fixed_edges_str = '[]'

                    if skip_furhter_tests_in_this_SG:
                        bar.update(cut='SKIP')
                        bar.update(bar.value+1)
                        continue

                    use_f128 = args.f128
                    while True:
                        results = []
                        LU_scalings = []
                        for scaling in scalings:
                            rescaled_momenta = [ v*scaling if i_v in UV_edge_indices else v for i_v, v in enumerate(random_momenta) ]
                            # Now transform the momentum configuration in this LMB into the defining LMB
                            transfo, parametric_shifts = UV_info_per_SG_and_cut[SG_name]['LMBs_to_defining_LMB_transfo'][LMB]
                            shifts = [ sum([external_momenta[i_shift]*shift 
                                        for i_shift, shift in enumerate(parametric_shift)]) for parametric_shift in parametric_shifts ]
                            rescaled_momenta_in_defining_LMB = transfo.dot(
                                [list(rm+shift) for rm, shift in zip(rescaled_momenta,shifts)] )

                            # Now obtain the rescaling for these momenta
                            LU_scaling_solutions = rust_worker.get_scaling(rescaled_momenta_in_defining_LMB,cut_ID)
                            if LU_scaling_solutions is None or len(LU_scaling_solutions)==0 or all(LU_scaling[0]<0. for LU_scaling in LU_scaling_solutions):
                                if args.show_warnings:
                                    logger.warning("Could not find rescaling for UV profiling of SG '%s' with cut ID #%d with UV edges %s and fixed edges %s: %s\nInput LMB momenta: %s"%(
                                        SG_name, cut_ID, UV_edges_str, fixed_edges_str, str(LU_scaling_solutions), str(rescaled_momenta_in_defining_LMB) ))
                                continue
                            LU_scaling_solutions = list(LU_scaling_solutions)
                            LU_scaling, LU_scaling_jacobian = LU_scaling_solutions.pop(0)
                            if LU_scaling>0.0 and args.show_warnings:
                                logger.warning("Found unexpected rescaling solutions for UV profiling of SG '%s' with cut ID #%d with UV edges %s and fixed edges %s: %s\nInput LMB momenta: %s"%(
                                    SG_name, cut_ID, UV_edges_str, fixed_edges_str, str(LU_scaling_solutions), str(rescaled_momenta_in_defining_LMB) ))

                            while LU_scaling < 0.0:
                                if len(LU_scaling_solutions)==0:
                                    break
                                LU_scaling, LU_scaling_jacobian = LU_scaling_solutions.pop(0)
                            LU_scalings.append((LU_scaling, LU_scaling_jacobian))
                            if use_f128:
                                res_re, res_im = rust_worker.evaluate_cut_f128(rescaled_momenta_in_defining_LMB,cut_ID,LU_scaling,LU_scaling_jacobian)                            
                            else:
                                res_re, res_im = rust_worker.evaluate_cut(rescaled_momenta_in_defining_LMB,cut_ID,LU_scaling,LU_scaling_jacobian)
                            results.append( (scaling, complex(res_re, res_im) ) )

                        dod, standard_error, number_of_points_considered, successful_fit = utils.compute_dod(results)
                        dod += 3*float(len(UV_edge_indices))

                        # We expect dod of at most five sigma above 0.0 for the integral to be convergent.
                        test_passed = (dod < float(args.target_scaling)+min(max(5.0*abs(standard_error),0.005),0.1) )

                        if (successful_fit and test_passed) or use_f128:
                            break
                        elif not args.no_f128:
                            use_f128 = True
                        else:
                            break

                    do_debug = False
                    if not args.scale_cuts and len(LU_scalings)>0:
                        t_sols = [LU_scaling[0] for LU_scaling in LU_scalings]
                        t_variance = (max(t_sols)-min(t_sols))/((max(t_sols)+min(t_sols))/2.0)
                        if t_variance > 1.0e-4:
                            if args.show_warnings:
                                logger.critical("UV profiling of SG '%s' with cut ID #%d with UV edges %s and fixed edges %s gives non-constant LU scaling. Found rel. vaiance of: %.5f over %d points."%(
                                    SG_name, cut_ID, UV_edges_str, fixed_edges_str, t_variance, len(LU_scalings)
                                ))
                                do_debug = True

                    if not successful_fit and (dod > float(args.target_scaling)-1.):
                        if args.show_warnings:
                            logger.critical("The fit for the UV scaling of SG '%s' with cut ID #%d with UV edges %s and fixed edges %s is unsuccessful (unstable). Found: %.3f +/- %.4f over %d points."%(
                                SG_name, cut_ID, UV_edges_str, fixed_edges_str, dod, standard_error, number_of_points_considered
                            ))
                            do_debug = True
                        n_fit_failed += 1
                        bar.update(fit_failed=n_fit_failed)

                    if args.verbose or do_debug or (not test_passed and args.show_fails):
                        logger.info("%s : UV scaling of SG '%s' with cut ID #%d with UV edges %s and fixed edges %s (--LMB %s --UV_indices %s):"%(
                            '%sPASS%s'%(Colours.GREEN, Colours.END) if test_passed else '%sFAIL%s'%(Colours.RED, Colours.END), SG_name, cut_ID, UV_edges_str, fixed_edges_str,
                            ' '.join(LMB), ' '.join('%d'%uv_index for uv_index in UV_edge_indices)
                        ))
                        if args.verbose or do_debug:
                            logger.info('\n'+'\n'.join('%-13.5e -> %-13.5e (LU scaling: %-13.5e , %-13.5e )'%(
                                r[0], abs(r[1]), LU_scalings[i_r][0], LU_scalings[i_r][1]) for i_r, r in enumerate(results)))
                        logger.info('%sdod= %.3f +/- %.3f%s'%(Colours.GREEN if test_passed else Colours.RED, dod, standard_error, Colours.END))

                    if test_passed:
                        n_passed += 1
                    else:
                        n_failed += 1
                    bar.update(passed=n_passed)
                    bar.update(failed=n_failed)
                    bar.update(bar.value+1)

                    SG['cutkosky_cuts'][cut_ID]['DERIVED_UV_dod'][(tuple(LMB), tuple(UV_edge_indices))] = (dod, standard_error, test_passed)

                    if not test_passed and args.skip_once_failed:
                        skip_furhter_tests_in_this_SG = True
                    if not test_passed and not this_SG_failed:
                        SG_failed += 1
                        bar.update(SG_failed=SG_failed)
                        this_SG_failed = True

                evaluate_integrand_dod = max(v[0] for v in SG['DERIVED_UV_dod'].values()) if len(SG['DERIVED_UV_dod'])>0 else None
                evaluate_cut_dod = None
                for cut in SG['cutkosky_cuts']:
                    cut_dod = max(v[0] for v in cut['DERIVED_UV_dod'].values()) if len(cut['DERIVED_UV_dod'])>0 else None
                    if cut_dod is not None:
                        if evaluate_cut_dod is None:
                            evaluate_cut_dod = cut_dod
                        else:
                            evaluate_cut_dod = max(evaluate_cut_dod, cut_dod)
                if evaluate_integrand_dod is not None and evaluate_cut_dod is not None:
                    SG['DERIVED_dod_consistency'] = (evaluate_integrand_dod-evaluate_cut_dod,evaluate_integrand_dod,evaluate_cut_dod)
                else:
                    SG['DERIVED_dod_consistency'] = (0.0,0.0,0.0)

                if not this_SG_failed:
                    SG_passed += 1
                    bar.update(SG_passed=SG_passed)          

        delta_t = time.time()-t_start_profile
        logger.info("UV profile completed in %d [s]: %s%d%s passed and %s failed (and %s fit failed)."%(
            int(delta_t),
            Colours.GREEN,
            n_passed,
            Colours.END,
            'none' if n_failed == 0 else '%s%d%s'%(Colours.RED, n_failed, Colours.END),
            'no' if n_fit_failed == 0 else '%s%d%s'%(Colours.RED, n_fit_failed, Colours.END)

        ))
        logger.info("Writing out processed yaml supergaphs on disk...")
        # Write out the results into processed topologies
        self.all_supergraphs.export(pjoin(self.dir_path, self._rust_inputs_folder))
        if len(selected_SGs)==1:
            self.do_display('%s --uv'%selected_SGs[0])
        else:
            self.do_display('--uv')

    def do_refresh_derived_data(self, line):
        """ Remove all processed data, on disk as well."""

        for file_path in glob.glob(pjoin(self.dir_path, self._rust_inputs_folder, 'PROCESSED_*.yaml')):
            shutil.move(file_path, pjoin(self.dir_path, self._rust_inputs_folder,'BACKUP_%s'%(
                os.path.basename(file_path)
            )))
        self.all_supergraphs = self.load_supergraphs()

    #### DISPLAY COMMAND
    display_parser = ArgumentParser(prog='display')
    display_parser.add_argument('SG_name', metavar='SG_name', type=str, nargs='?',
                    help='the name of a supergraph to display')
    display_parser.add_argument(
        '-t','--timing',action="store_true", dest="timing", default=False,
        help="Show timing profile information")
    display_parser.add_argument(
        '--uv',action="store_true", dest="uv", default=False,
        help="Show UV profile information")
    display_parser.add_argument(
        '--ir',action="store_true", dest="ir", default=False,
        help="Show IR profile information")
    display_parser.add_argument(
        "-f","--full", action="store_true", dest="full", default=False,
        help="exhaustively show information")
    display_parser.add_argument(
        "-sm","--show_momenta", action="store_true", dest="show_momenta", default=False,
        help="Show the momenta of the edges in the E-surfaces for the intersection point approached in the IR.")
    def help_display(self):
        self.display_parser.print_help()
        return
    def do_display(self, line):
        """ display command """

        if line=='help':
            self.display_parser.print_help()
            return 

        args = self.split_arg(line)
        args = self.display_parser.parse_args(args)

        # First code for particular information
        if args.timing:
            if args.SG_name:
                logger.info("Timing profile for supergraph '%s':\n%s"%(
                    args.SG_name, self.all_supergraphs[args.SG_name].show_timing_statistics()
                ))
            else:
                logger.info("Overall timing profile for all supergraphs:\n%s"%(
                    self.all_supergraphs.show_timing_statistics()
                ))

        if args.uv:
            if args.SG_name:
                sg_collection = SuperGraphCollection()
                sg_collection[args.SG_name] = self.all_supergraphs[args.SG_name]
                logger.info("UV profile for supergraph '%s':\n%s"%(
                    args.SG_name, sg_collection.show_UV_statistics()
                ))
            else:
                logger.info("Overall UV profile for all supergraphs:\n%s"%(
                    self.all_supergraphs.show_UV_statistics()
                ))

        if args.ir:
            if args.SG_name:
                sg_collection = SuperGraphCollection()
                sg_collection[args.SG_name] = self.all_supergraphs[args.SG_name]
                logger.info("UV profile for supergraph '%s':\n%s"%(
                    args.SG_name, sg_collection.show_IR_statistics(show_momenta=args.show_momenta)
                ))
            else:
                logger.info("Overall UV profile for all supergraphs:\n%s"%(
                    self.all_supergraphs.show_IR_statistics(show_momenta=args.show_momenta)
                ))

        # Only show general statistics when not showing anything else
        if args.timing or args.uv:
            return

        if args.SG_name is None:
            if args.full:
                logger.info("Exhaustive content of the cross-section set:\n%s"%str(self.cross_section_set))
            else:
                logger.info("Summary of the cross-section set:\n%s"%self.cross_section_set.summary_str())
        elif args.SG_name == 'ALL':
            logger.info("General information about all supergraphs:\n%s"%self.all_supergraphs.summary_str())
        elif args.SG_name in ['gen','generation_statistics']:
            if len(self.generation_statistics)==0:
                logger.info("Generation statistics not available.")
            else:
                logger.info("Generation statistics:\n%s"%pformat(self.generation_statistics))
        else:
            if args.SG_name not in self.all_supergraphs:
                raise alphaLoopInvalidRunCmd("Supergraph named '%s' not found in the supergraphs loaded."%args.SG_name)
            if args.full:
                logger.info("Exhaustive content of the supergraph '%s%s%s':\n%s"%(Colours.GREEN, args.SG_name, Colours.END, str(self.all_supergraphs[args.SG_name])))
            else:
                logger.info("Summary of the supergraph %s%s%s:\n%s"%(Colours.GREEN, args.SG_name, Colours.END, self.all_supergraphs[args.SG_name].summary_str()))

    #### SET_HYPERPARAMETER COMMAND
    set_hyperparameter_parser = ArgumentParser(prog='set_hyperparameter')
    set_hyperparameter_parser.add_argument('param_value', metavar='param_value', type=str, nargs=2,
                    help='parameter name and value to set')
    set_hyperparameter_parser.add_argument(
        "-w", "--write", action="store_true", dest="write", default=False,
        help="Write hyperparameter to disk.")
    def help_set_hyperparameter(self):
        self.set_hyperparameter_parser.print_help()
        return
    def do_set_hyperparameter(self, line):
        """ set_hyperparameter command """

        if line=='help':
            self.set_hyperparameter_parser.print_help()
            return 

        args = self.split_arg(line)
        args = self.set_hyperparameter_parser.parse_args(args)
        try:
            parsed_value = eval(args.param_value[1])
        except Exception as e:
            raise alphaLoopInvalidRunCmd("Could not evaluate specified value '%s' for parameter '%s':\n%s"%(
                args.param_value[0], args.param_value[1],str(e)
            ))
        if not self.hyperparameters.set_parameter(args.param_value[0],parsed_value):
            raise alphaLoopInvalidRunCmd("Failed to set the following hyperparameter '%s'."%args.param_value[0])
        if args.write:
            self.hyperparameters.export_to(pjoin(self.dir_path, 'hyperparameters.yaml'))            

    #### INTEGRATE COMMAND
    integrate_parser = ArgumentParser(prog='integrate')
    integrate_parser.add_argument('SG_name', metavar='SG_name', type=str, nargs='?',
                    help='the name of a supergraph to display')
    integrate_parser.add_argument('-s','--sampling', metavar='sampling', type=str, default='xs', 
                    choices=('xs','flat', 'advanced', 'test_h_function'), help='Specify the sampling method (default: %(default)s)')
    integrate_parser.add_argument('-i','--integrator', metavar='integrator', type=str, default='vegas3', 
                    choices=('naive','vegas', 'vegas3', 'inspect'), help='Specify the integrator (default: %(default)s)')
    integrate_parser.add_argument('-hf','--h_function', metavar='h_function', type=str, default='left_right_polynomial', 
                    choices=('left_right_polynomial','left_right_exponential', 'flat'), help='Specify the h-function to use in sampling (default: %(default)s)')
    integrate_parser.add_argument('-hfs','--h_function_sigma', metavar='h_function_sigma', type=int, default=3,
                    help='Spread of the h-function in sampling, higher=steeper (default: %(default)s).')
    integrate_parser.add_argument('-alhf','--al_h_function', metavar='al_h_function', type=str, default='left_right_exponential', 
                    choices=('left_right_polynomial','left_right_exponential'), help='Specify the h-function to use in alphaLoop (default: %(default)s)')
    integrate_parser.add_argument('-alhfs','--al_h_function_sigma', metavar='al_h_function_sigma', type=int, default=1, 
                    help='Spread of the h-function in alphaLoop, higher=steeper (default: %(default)s).')
    integrate_parser.add_argument('-v','--verbosity', metavar='verbosity', type=int, default=0,choices=(0,1,2,3),
                    help='verbosity level (default: %(default)s).')
    integrate_parser.add_argument('-c','--n_cores', metavar='n_cores', type=int, default=multiprocessing.cpu_count(),
                    help='Number of cores to parallelize on (default: %(default)s).')
    integrate_parser.add_argument('-nis','--n_iterations_survey', metavar='n_iterations_survey', type=int, default=10,
                    help='Number of iteration for the survey stage (default: %(default)s).')
    integrate_parser.add_argument('-nir','--n_iterations_refine', metavar='n_iterations_refine', type=int, default=5,
                    help='Number of iteration for the refine stage (default: %(default)s).')
    integrate_parser.add_argument('-nps','--n_points_survey', metavar='n_points_survey', type=int, default=int(1.0e5),
                    help='Number of sample points per iteration for the survey stage (default: %(default)s).')
    integrate_parser.add_argument('-npr','--n_points_refine', metavar='n_points_refine', type=int, default=int(1.0e5),
                    help='Number of sample points per iteration for the refine stage (default: %(default)s).')
    integrate_parser.add_argument('--n_max', metavar='n_max', type=int, default=int(1.0e7),
                    help='Maximum number of sample points in Vegas (default: as per hyperparameters).')
    integrate_parser.add_argument('--n_max_survey', metavar='n_max_survey', type=int, default=-1,
                    help='Maximum number of sample points in Vegas for survey (default: no survey).')
    integrate_parser.add_argument('--target_accuracy_survey', metavar='target_accuracy_survey', type=float, default=1.0e-5,
                    help='Target accuracy for Vegas survey stage (default: %(default)f).')
    integrate_parser.add_argument('--target_accuracy', metavar='target_accuracy', type=float, default=1.0e-5,
                    help='Target accuracy for Vegas refine stage (default: %(default)f).')
    integrate_parser.add_argument('--load_grids', metavar='load_grids', type=str, default=None,
                    help='Specify a Vegas grid file to load from. (default: None).')
    integrate_parser.add_argument('--n_start', metavar='n_start', type=int, default=-1,
                    help='Starting number of sample points in Vegas (default: as per hyperparameters).')
    integrate_parser.add_argument('--n_increase', metavar='n_increase', type=int, default=-1,
                    help='Starting number of sample points in Vegas (default: as per hyperparameters).')
    integrate_parser.add_argument('-bs','--batch_size', metavar='batch_size', type=int, default=-1,
                    help='Batch size for parallelisation (default: as per hyperparameters n_vec).')
    integrate_parser.add_argument('--seed', metavar='seed', type=int, default=0,
                    help='Specify the random seed for the integration (default: %(default)s).')
    integrate_parser.add_argument('-hp','--hyperparameters', metavar='hyperparameters', type=str, default=[], nargs='+',
                    help='Specify particular hyperparameters to overwrite in pairs of form <hp_name> <hp_str_expression_value>  (default: %(default)s).')
    integrate_parser.add_argument('-ccs','--selected_cutkosky_cuts_and_sides', metavar='selected_cutkosky_cuts_and_sides', type=int, nargs='+', default=[-1,],
                    help='Selected cutkosky cut and sides for the multichanneling. [-1,] means sum over all. (default: %(default)s).')
    integrate_parser.add_argument('-lmbs','--selected_lmbs', metavar='selected_lmbs', type=int, nargs='+', default=[-1,],
                    help='Selected lmb indices for the multichanneling. [-1,] means sum over all. (default: %(default)s).')
    integrate_parser.add_argument('-xs','--xs', metavar='xs', type=float, nargs='+', default=[-1.,],
                    help='Selected random variables to probe with inspect.')
    integrate_parser.add_argument(
        '-mc','--multichanneling',action="store_true", dest="multichanneling", default=False,
        help="Enable multichanneling (default: as per hyperparameters)")
    integrate_parser.add_argument(
        '-no_mc','--no_multichanneling',action="store_false", dest="multichanneling", default=False,
        help="Disable multichanneling (default: as per hyperparameters)")
    integrate_parser.add_argument(
        '-nw','--no_warnings',action="store_false", dest="show_warnings", default=True,
        help="Disable the printout of numerical warnings during the integration.")

    def help_integrate(self):
        self.integrate_parser.print_help()
        return
    # We must wrape this function in a process because of the border effects of the pyO3 rust Python bindings
    #@wrap_in_process()
    @with_tmp_hyperparameters({
        'Integrator.dashboard': False
    })
    def do_integrate(self, line):
        """ Integrate a given (set of) supergraphs using different sampling strategies."""

        if line=='help':
            self.integrate_parser.print_help()
            return 

        args = self.split_arg(line)
        args = self.integrate_parser.parse_args(args)

        # We need to detect here if we are in the amplitude-mock-up situation with frozen external momenta.
        frozen_momenta = None
        if 'external_data' in self.cross_section_set:
            frozen_momenta = {
                'in' : self.cross_section_set['external_data']['in_momenta'],
                'out' : self.cross_section_set['external_data']['out_momenta'],
            }
            # Also force the specified incoming momenta specified in the hyperparameters to match the frozen specified ones.
            self.hyperparameters.set_parameter('CrossSection.incoming_momenta',frozen_momenta['in'])
            self.hyperparameters.set_parameter('CrossSection.do_rescaling',False)
            self.hyperparameters.set_parameter('CrossSection.fixed_cut_momenta',frozen_momenta['out'])

        self.hyperparameters.set_parameter('General.multi_channeling',args.multichanneling)

        args.hyperparameters = [( args.hyperparameters[i],eval(args.hyperparameters[i+1]) ) for i in range(0,len(args.hyperparameters),2)]
        for hp_name, value in args.hyperparameters:
            self.hyperparameters.set_parameter(hp_name,value)

        if args.h_function == 'left_right_polynomial':
            selected_h_function = CallableInstanceWrapper(sampler.HFunction(args.h_function_sigma, debug=args.verbosity))
        elif args.h_function == 'left_right_exponential':
            selected_h_function = CallableInstanceWrapper(sampler.DiscreteExponentialHFunction(args.h_function_sigma, debug=args.verbosity))
        elif args.h_function == 'flat':
            selected_h_function = CallableInstanceWrapper(sampler.FlatHFunction(args.h_function_sigma, debug=args.verbosity))
        else:
            raise alphaLoopInvalidRunCmd("Unsupported h-function specification: %s'."%args.h_function)

        self.hyperparameters.set_parameter('CrossSection.NormalisingFunction.name',args.al_h_function)
        self.hyperparameters.set_parameter('CrossSection.NormalisingFunction.spread',float(args.al_h_function_sigma))

        if args.n_cores > 1:
            runner = cluster.MultiCore(args.n_cores)
        else:
            runner = cluster.onecore

        if args.n_start < 0:
            args.n_start = self.hyperparameters['Integrator']['n_start']
        if args.n_max < 0:
            args.n_max = self.hyperparameters['Integrator']['n_max']
        if args.n_increase < 0:
            args.n_increase = self.hyperparameters['Integrator']['n_increase']
        if args.batch_size < 0:
            args.batch_size = self.hyperparameters['Integrator']['n_vec']

        if args.seed > 0:
            random.seed(args.seed)

        if args.integrator == 'naive':
            selected_integrator = integrators.SimpleMonteCarloIntegrator
            integrator_options = {
                 'n_iterations' : args.n_iterations_refine,
                 'n_points_per_iterations' : args.n_points_refine,
                 'verbosity' : args.verbosity+1,
                 'seed' : args.seed
            }
        elif args.integrator == 'vegas':
            selected_integrator = pyCubaIntegrator.pyCubaIntegrator
            integrator_options = {
                 'cluster' : runner,
                 'max_eval' : args.n_max,
                 'n_start' : args.n_start,
                 'n_increase' : args.n_increase,
                 'n_start_survey' : args.n_start,
                 'n_increase_survey' : args.n_increase,
                 'n_batch' : args.batch_size,
                 'state_file_folder' : pjoin(self.dir_path, self._run_workspace_folder),
                 'max_eval_survey' : args.n_max_survey,
                 'target_accuracy_survey' :args.target_accuracy_survey,
                 'target_accuracy' :args.target_accuracy,
                 'load_grids' : args.load_grids,
                 'n_vec' : 1,
            }
        elif args.integrator == 'vegas3':
            selected_integrator = vegas3_integrator.Vegas3Integrator
            integrator_options = {
                 'cluster' : runner,
                 'verbosity' : args.verbosity+2,
                 'seed' : args.seed,
                 'n_iterations_survey' : args.n_iterations_survey,
                 'n_points_survey' : args.n_points_survey,
                 'n_iterations_refine' : args.n_iterations_refine,
                 'n_points_refine' : args.n_points_refine,
                 'batch_size' : args.batch_size,
            }
        elif args.integrator == 'inspect':
            # Inspection requires no integrator
            selected_integrator = None
        else:
            raise alphaLoopInvalidRunCmd("Unsupported integrator specification: %s'."%args.integrator)

        if args.sampling == 'test_h_function':

            logger.info("Dummy integrand for testing h-function with integrator '%s':"%args.integrator)

            my_integrand = sampler.TestHFuncIntegrand( selected_h_function, debug=args.verbosity )

            my_integrator = selected_integrator(my_integrand, **integrator_options)

            result = my_integrator.integrate()

            logger.info('')
            if my_integrator.n_evals_failed.value > 0:
                logger.info("Number of failed evaluations: %d (%.6g%%)"%(
                    my_integrator.n_evals_failed.value, (my_integrator.n_evals_failed.value/my_integrator.n_evals.value)*100.
                ))
            logger.info("Zero weights fraction: %.6g%%"%(
                float((my_integrator.n_zero_evals.value / my_integrator.n_evals.value))*100.
            ))
            if my_integrand.max_eval_positive.value != 0.0:
                logger.info("Maximum posivite weight found: %.6g (%.1e x central value) for xs=[%s]"%(
                    my_integrand.max_eval_positive.value, 
                    abs(my_integrand.max_eval_positive.value/result[0]),
                    ' '.join('%.16f'%x for x in my_integrand.max_eval_positive_xs) 
                ))
            if my_integrand.max_eval_negative.value != 0.0:
                logger.info("Maximum negative weight found: %.6g (%.1e x central value) for xs=[%s]"%(
                    my_integrand.max_eval_negative.value, 
                    abs(my_integrand.max_eval_negative.value/result[0]),
                    ' '.join('%.16f'%x for x in my_integrand.max_eval_negative_xs) 
                ))
            logger.info("Result of the test integration using h function %s with %d function calls (target is 1.0): %.4e +/- %.2e"%(
                args.h_function, my_integrator.tot_func_evals, result[0], result[1]))
            logger.info('')
            return
        else:
            SG_name = args.SG_name
            if SG_name not in self.all_supergraphs:
                raise alphaLoopInvalidRunCmd("Cannot find SG named in '%s' in the collection loaded."%SG_name)
            
            SG_info = self.all_supergraphs[SG_name]

            logger.info("Integrating SG '%s' with sampler '%s' and integrator '%s':"%(
                SG_name, args.sampling, args.integrator
            ))

            if args.sampling in ['xs',]:

                self.hyperparameters.set_parameter('General.multi_channeling',args.multichanneling)
                rust_worker = self.get_rust_worker(SG_name)
                
                n_integration_dimensions = SG_info['topo']['n_loops']*3
                if frozen_momenta is not None:
                    n_integration_dimensions -= 3*len(frozen_momenta['out'])

                if args.sampling == 'xs':
                    my_sampler = sampler.generator_aL( 
                        integrands.DimensionList([ 
                            integrands.ContinuousDimension('x_%d'%i_dim,lower_bound=0.0, upper_bound=1.0) 
                            for i_dim in range(1,n_integration_dimensions+1)
                        ]),
                        rust_worker,
                        SG_info,
                        self.alphaLoop_interface._curr_model,
                        selected_h_function,
                        self.hyperparameters,
                        debug=args.verbosity,
                        frozen_momenta=frozen_momenta
                    )

                my_integrand = sampler.DefaultALIntegrand( rust_worker, my_sampler, 
                    debug=args.verbosity, phase=self.hyperparameters['Integrator']['integrated_phase'], frozen_momenta=frozen_momenta)

            elif args.sampling in ['flat','advanced']:
                
                # The multichanneling is done in-house by our sampling for these args.sampling options
                self.hyperparameters.set_parameter('General.multi_channeling',False)

                self.hyperparameters.set_parameter('Parameterization.mapping','linear')
                self.hyperparameters.set_parameter('Parameterization.b',1.0)
                self.hyperparameters.set_parameter('Parameterization.mode','spherical')

                rust_worker = self.get_rust_worker(SG_name)

                # Specify particular channels options
                if args.multichanneling:
                    channel_for_generation = None
                    selected_cut_and_side = None if args.selected_cutkosky_cuts_and_sides[0]<0 else args.selected_cutkosky_cuts_and_sides
                    selected_LMB = None if args.selected_lmbs[0]<0 else args.selected_lmbs
                else:
                    channel_for_generation = ( 
                        args.selected_cutkosky_cuts_and_sides[0] if args.selected_cutkosky_cuts_and_sides[0]>=0 else 0,
                        args.selected_lmbs[0] if args.selected_lmbs[0]>=0 else 0 
                    )
                    selected_cut_and_side = None
                    selected_LMB = None                    

                computed_model = model_reader.ModelReader(self.alphaLoop_interface._curr_model)
                computed_model.set_parameters_and_couplings(
                                            pjoin(self.dir_path,'Source','MODEL','param_card.dat'))

                my_integrand = sampler.AdvancedIntegrand(
                    rust_worker,
                    SG_info,
                    computed_model,
                    selected_h_function,
                    self.hyperparameters,
                    debug=args.verbosity,
                    external_phase_space_generation_type = args.sampling,
                    channel_for_generation = channel_for_generation,
                    selected_cut_and_side=selected_cut_and_side, 
                    selected_LMB=selected_LMB,
                    phase=self.hyperparameters['Integrator']['integrated_phase'],
                    show_warnings=args.show_warnings,
                    return_individual_channels = (args.integrator in ['vegas3','inspect']),
                    frozen_momenta=frozen_momenta
                )

            # Inspect mode
            if selected_integrator is None:
                if len(args.xs)==0 or args.xs[0]<0.:
                    raise InvalidCmd("When running in inspect mode, the random variables must be supplied with -xs = x1 x2 x3 ...") 

                if len(args.xs) != (SG_info['topo']['n_loops'] - (0 if frozen_momenta is None else len(frozen_momenta['out'])))*3:
                    raise InvalidCmd("Expected %d random variables, but only %d were specified."%(SG_info['topo']['n_loops']*3,len(args.xs))) 

                res = my_integrand(args.xs, [])
                logger.info("Final weight for point xs=[%s] :\n %s"%(
                    ' '.join('%.16f'%x for x in args.xs), pformat(res) 
                ))
                return

                # import matplotlib
                # import matplotlib.pyplot as plt
                # mpl_logger = logging.getLogger("matplotlib")
                # mpl_logger.setLevel(logging.WARNING)
    
                # ts = []
                # wgts = []

                # min_t = 0.01
                # max_t = 0.99
                # t_incr = 0.01
                # all_ts = [min_t+i*t_incr for i in range(int((max_t-min_t)/t_incr))]

                # # min_t = 0.1
                # # max_t = 2.0
                # # t_incr = 0.05
                # # all_ts += [min_t+i*t_incr for i in range(int((max_t-min_t)/t_incr))]

                # # min_t = 2.0
                # # max_t = 10.0
                # # t_incr = 0.1
                # # all_ts += [min_t+i*t_incr for i in range(int((max_t-min_t)/t_incr))]

                # # min_t = 10.0
                # # max_t = 100.0
                # # t_incr = 1.0
                # # all_ts += [min_t+i*t_incr for i in range(int((max_t-min_t)/t_incr))]

                # for i_t, t_value in enumerate(all_ts):
                #     print("Currently at t #%d/%d"%(i_t, len(all_ts)),end='\r')
                #     #x_t, inv_wgt_t = my_integrand.h_function.inverse_sampling(t_value)
                #     #this_xs =[x_t,]+args.xs[1:] 
                #     this_xs = list(args.xs)
                #     tot_res = 0.
                #     this_xs[2] = t_value
                #     #this_xs[5] = t_value
                #     #this_xs[5] = args.xs[5]
                #     res = my_integrand(this_xs, [])
                #     tot_res += res['I']
                #     #this_xs[5] = (t_value-0.5) if t_value-0.5 > 0. else t_value+0.5
                #     #this_xs[8] = args.xs[8]
                #     #res = my_integrand(this_xs, [])
                #     #tot_res += res['I']
                #     # this_xs[5] = t_value
                #     # this_xs[8] = (args.xs[8]-0.5) if args.xs[8]-0.5 > 0. else args.xs[8]+0.5
                #     # res = my_integrand(this_xs, [])
                #     # tot_res += res['I']
                #     # this_xs[5] = (t_value-0.5) if t_value-0.5 > 0. else t_value+0.5
                #     # this_xs[8] = (args.xs[8]-0.5) if args.xs[8]-0.5 > 0. else args.xs[8]+0.5
                #     # res = my_integrand(this_xs, [])
                #     # tot_res += res['I']
                #     tot_res /= 1.
                #     ts.append(t_value)
                #     wgts.append(tot_res)
                #     #logger.info("Final weight for point xs=[%s] : %.16e"%(
                #     #    ' '.join('%.16f'%x for x in this_xs), res 
                #     #))
                #     #misc.sprint("Eval for t=%.16f = %.16e"%(t_value,res))

                # fig, ax = plt.subplots()
                # plt.yscale('linear')
                # plt.xscale('linear')
                # ax.plot(ts, wgts)

                # ax.set(xlabel='x_phi(pq2)', ylabel='I(x_phi(pq2))',
                #     title='Symmetrized integrand h > gg vs angle')
                # ax.grid()

                # #fig.savefig("test.png")
                # plt.show()

                # return

            my_integrator = selected_integrator(my_integrand, **integrator_options)

            result = my_integrator.integrate()

            if args.integrator == 'vegas3' and my_integrator.full_result is not None:
                logger.info('')
                logger.info('Summary of the integration:\n%s'%str(my_integrator.full_result.summary()))
                logger.info('')
                logger.info("Complete integration result for all integrand components:")
                max_key_length = max(len(k) for k in my_integrator.full_result.keys())
                for k in sorted(my_integrator.full_result.keys()):
                    logger.info(('%-{}s : %s%.6g +/- %.4g (%.2g%%)'.format(max_key_length))%(k, 
                        ' ' if my_integrator.full_result[k].mean>=0. else '',
                        my_integrator.full_result[k].mean, my_integrator.full_result[k].sdev,
                        abs(my_integrator.full_result[k].sdev/my_integrator.full_result[k].mean)*100. if my_integrator.full_result[k].mean!=0. else 0.
                    ))

            # RAvg.mean : mean 
            # RAvg.sdev : standard dev.
            # RAvg.chi2 : Chi^2 of the weighted average
            # RAvg.dof  : number of degreeas of freedom
            # RAvg.Q    : p-value of the weighted average 
            # RAvg.itn_results : list of the integral estimates for each iteration
            # RAvg.summary() : summary of the integration

            logger.info('')
            if my_integrand.n_evals_failed.value > 0:
                logger.info("Number of failed evaluations: %d (%.6g%%)"%(
                    my_integrand.n_evals_failed.value, (my_integrand.n_evals_failed.value/my_integrand.n_evals.value)*100.
                ))
            logger.info("Zero weights fraction: %.6g%%"%(
                float((my_integrand.n_zero_evals.value / my_integrand.n_evals.value))*100.
            ))
            if my_integrand.max_eval_positive.value != 0.0:
                logger.info("Maximum posivite weight found: %.6g (%.1e x central value) for xs=[%s]"%(
                    my_integrand.max_eval_positive.value, 
                    abs(my_integrand.max_eval_positive.value/result[0]),
                    ' '.join('%.16f'%x for x in my_integrand.max_eval_positive_xs) 
                ))
            if my_integrand.max_eval_negative.value != 0.0:
                logger.info("Maximum negative weight found: %.6g (%.1e x central value) for xs=[%s]"%(
                    my_integrand.max_eval_negative.value, 
                    abs(my_integrand.max_eval_negative.value/result[0]),
                    ' '.join('%.16f'%x for x in my_integrand.max_eval_negative_xs) 
                ))
            logger.info("Result of the cross-section for %s%s%s of %s%s%s with sampler %s%s%s and integrator %s%s%s, using %s%d%s function calls:\n%s %.6g +/- %.4g%s (%.2g%%)"%(
                utils.bcolors.GREEN, SG_name, utils.bcolors.ENDC,
                utils.bcolors.GREEN, os.path.basename(self.dir_path), utils.bcolors.ENDC,
                utils.bcolors.BLUE, args.sampling, utils.bcolors.ENDC,
                utils.bcolors.BLUE, args.integrator, utils.bcolors.ENDC,
                utils.bcolors.GREEN, my_integrator.tot_func_evals, utils.bcolors.ENDC,
                utils.bcolors.GREEN, result[0], result[1], utils.bcolors.ENDC,
                abs(result[1]/result[0])*100. if result[0]!=0. else 0.
            ))
            logger.info('')

            return

    #### EXPERIMENT COMMAND
    experiment_parser = ArgumentParser(prog='experiment')
    experiment_parser.add_argument('SG_name', metavar='SG_name', type=str, nargs='?',
                    help='The name of a supergraph to consider')
    experiment_parser.add_argument(
        '-e','--experiment', dest="experiment", type=str, default='default',
        help="Which experiment to run")
    def help_experiment(self):
        self.experiment_parser.print_help()
        return 
    # We must wrape this function in a process because of the border effects of the pyO3 rust Python bindings
    @wrap_in_process()
    @with_tmp_hyperparameters({
        'Integrator.dashboard': False
    })
    def do_experiment(self, line):
        """ Integrate a given (set of) supergraphs using different sampling strategies."""
        
        if line=='help':
            self.experiment_parser.print_help()
            return 

        args = self.split_arg(line)
        args = self.integrate_parser.parse_args(args)

        logger.critical("A PLACE TO PUT IN SOME EXPERIMENTS")

        SG = self.all_supergraphs[args.SG_name]
        E_cm = SG.get_E_cm(self.hyperparameters)
        rust_worker = self.get_rust_worker(args.SG_name)
        pass


    ######################################################################
    #
    # Example of the implementation of a trivial function 'hello_world'
    #
    ######################################################################

    def do_hello_world(self, line):
        """ Hello world command example."""

        logger.info('Hello World of alphaLoop run interface, and welcome you, %s%s%s!'%(utils.bcolors.GREEN, line,utils.bcolors.ENDC))

    def help_hello_world(self, line):
        """ Hello world command example."""

        logger.info('Contextual help for command hello world.')

    def complete_hello_world(self, text, line, begidx, endidx):
        """ Hello world command example."""

        return self.list_completion(text,['something', 'else'], line)

    def get_alpha_loop_banner(self):
        """ Returns a string of alpha loop banner."""

        res =[]
        res.append(( "%s"+"="*80+"=%s")%(utils.bcolors.GREEN, utils.bcolors.ENDC))
        res.append(( "%s||"+" "*36+u'\u03B1Loop'+" "*36+"||%s")%(utils.bcolors.GREEN, utils.bcolors.ENDC) )
        res.append(( "%s"+"="*80+"=%s")%(utils.bcolors.GREEN, utils.bcolors.ENDC) )
        return '\n'.join(res)

    #command to change the prompt 
    def preloop(self, *args, **opts):
        """only change the prompt after calling  the mother preloop command"""

        # The colored prompt screws up the terminal for some reason.
        #self.prompt = '\033[92mGGVV > \033[0m'
        self.prompt = "aLoop @ %s > "%os.path.basename(self.dir_path)

        # preloop mother
        madgraph_interface.CmdExtended.preloop(self)