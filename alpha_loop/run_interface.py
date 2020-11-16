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

#import matplotlib.pyplot as plt
#from matplotlib.font_manager import FontProperties

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
            'CrossSection.NormalisingFunction.name'         : 'left_right_polynomial',
            'CrossSection.NormalisingFunction.center'       : 1.0,
            'CrossSection.NormalisingFunction.spread'       : 3,
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
            # evaluate the C expression for the sum of diagrams
            'CrossSection.sum_diagram_sets'                      :   False,
            # compare locally against the same topology written in another loop momentum basis
            'CrossSection.compare_with_additional_topologies'    :   False,
            'CrossSection.inherit_deformation_for_uv_counterterm':   True
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

    def export(self, SG_name, dir_path):
        open(pjoin(dir_path,'PROCESSED_%s.yaml'%SG_name),'w').write(
            yaml.dump(dict(self), Dumper=Dumper, default_flow_style=False))

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

    def show_UV_statistics(self):

        res = []

        all_SG_dods = sorted([
            (SG_name,k,v) for SG_name,SG in self.items() for k,v in SG['DERIVED_UV_dod'].items() if 'DERIVED_UV_dod' in SG
        ],key=lambda el: el[2][0], reverse=True)
        fail_SG_dods = [ (SG_name,k,v) for SG_name,k,v in all_SG_dods if not v[-1] ]
        all_SG_cut_dods = sorted([
            (SG_name, cut_ID, k, v) for SG_name,SG in self.items() for cut_ID, cut in enumerate(SG['cutkosky_cuts'])
            for k,v in cut['DERIVED_UV_dod'].items() if 'DERIVED_UV_dod' in cut
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
        res_str.append('%s%s%s'%(Colours.BLUE,'Statistics for individual supergraphs',Colours.END))
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
            p = multiprocessing.Process(target=fowarad_function_output_to_queue, args=tuple([opts,]+list(args)))
            p.start()
            p.join()
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

        self.cross_section_set = CrossSectionSet(pjoin(self.dir_path,
                            self._rust_inputs_folder, self._cross_section_set_yaml_name))
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

        return rust_worker

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
    timing_profile_parser = ArgumentParser()
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
    # We must wrape this function in a process because of the border effects of the pyO3 rust Python bindings
    @wrap_in_process()
    @with_tmp_hyperparameters({
        'Integrator.dashboard': False,
        'General.multi_channeling' : False
    })
    def do_timing_profile(self, line):
        """ Automatically timing profile a process output."""
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

    #### UV PROFILE COMMAND
    uv_profile_parser = ArgumentParser()
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
        "-sof", "--skip_once_failed", action="store_true", dest="skip_once_failed", default=True,
        help="Skip the probing of a supergraph once it failed.")
    uv_profile_parser.add_argument(
        "-sf", "--show_fails", action="store_true", dest="show_fails", default=True,
        help="Show exhaustive information for each fail.")
    uv_profile_parser.add_argument(
        "-sc", "--scale_cuts", action="store_true", dest="scale_cuts", default=False,
        help="Include UV scaling of edges being cut.")
    uv_profile_parser.add_argument("-n_max","--n_max", dest='n_max', type=int, default=-1,
                    help='Set the maximum number of UV test to perform per cut (default: all)')
    uv_profile_parser.add_argument(
        "-v", "--verbose", action="store_true", dest="verbose", default=False,
        help="Enable verbose output.")
    # We must wrape this function in a process because of the border effects of the pyO3 rust Python bindings
    @wrap_in_process()
    @with_tmp_hyperparameters({
        'Integrator.dashboard': False,
        'General.minimal_precision_for_returning_result': 3,
        'CrossSection.NormalisingFunction.name'         : 'left_right_polynomial',
        'CrossSection.NormalisingFunction.center'       : 1.0,
        'General.multi_channeling'                      : False
    })
    def do_uv_profile(self, line):
        """ Automatically probe all UV limits of a process output."""
        args = self.split_arg(line)
        args = self.uv_profile_parser.parse_args(args)
        
        if args.SG_name is None:
            selected_SGs = list(self.all_supergraphs.keys())
        else:
            selected_SGs = [args.SG_name,]

        self.hyperparameters['CrossSection']['NormalisingFunction']['spread'] = args.h_power
        if args.required_precision is None:
            self.hyperparameters['General']['stability_checks'][-1]['relative_precision']=1.0e-99
        else:
            for entry in self.hyperparameters['General']['stability_checks']:
                entry['relative_precision'] = args.required_precision

        logger.info("Starting UV profile...")

        # Compute the log-spaced sequence of rescaling
        scalings = [ 10.**(math.log10(args.min_scaling)+i*((math.log10(args.max_scaling)-math.log10(args.min_scaling))/(args.n_points-1))) 
                        for i in range(args.n_points) ]

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
            UV_info_per_SG_and_cut['SG_topo_gen'] = SG_topo_gen
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
                for cut_ID in sorted(list(uv_probes_to_add_for_this_cut.keys())):
                    UV_info_per_SG_and_cut[SG_name]['UV_edges_to_probe_for_cuts'].extend(uv_probes_to_add_for_this_cut[cut_ID])
                UV_info_per_SG_and_cut[SG_name]['UV_edges_to_probe'].extend(uv_probes_to_add_for_this_SG)
            else:
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
                            # We must multiply by the overall jac to get a convergence UV scaling for dod = 0 as defined including line A).
                            results.append( (scaling, complex(res_re, res_im)*overall_jac ) )

#                        misc.sprint(results)
                        # Here we are in x-space, so the dod read is already the one we want.
                        dod, standard_error, number_of_points_considered, successful_fit = utils.compute_dod(results)

                        # A) The jacobian is not included in the evaluate_cut functions, so we must now account for it here for each UV edge
                        dod += 4*float(len(UV_edge_indices))
                        # Subtract one to the dod so that it's interpretation matches the power counting in the original Minkowski graph,
                        # i.e. -> dod=0 is a log divergence
                        dod -= 1

                        # We expect dod of at most five sigma above 0.0 for the integral to be convergent.
                        test_passed = (dod < float(args.target_scaling)+min(max(5.0*abs(standard_error),0.005),0.1) )

                        if (successful_fit and test_passed) or use_f128:
                            break
                        elif not args.no_f128:
                            use_f128 = True
                        else:
                            break

                    do_debug = False
                    if not successful_fit and dod > float(args.target_scaling)-2.0:
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
                            # Also divide the result by the UV scaling otherwise line A) below is incorrect
                            # when the scaling is not a constant.
                            #additional_t_scaling = LU_scaling**(-len(LMB))
                            additional_t_scaling = 1.0
                            results.append( (scaling, complex(res_re, res_im)*additional_t_scaling ) )
                        
                        dod, standard_error, number_of_points_considered, successful_fit = utils.compute_dod(results)

                        # A) The jacobian is not included in the evaluate_cut functions, so we must now account for it here for each UV edge
                        dod += 4*float(len(UV_edge_indices))
                        # Subtract one to the dod so that it's interpretation matches the power counting in the original Minkowski graph,
                        # i.e. -> dod=0 is a log divergence
                        dod -= 1

                        # We expect dod of at most five sigma above 0.0 for the integral to be convergent.
                        test_passed = (dod < float(args.target_scaling)+min(max(5.0*abs(standard_error),0.005),0.1) )

                        if (successful_fit and test_passed) or use_f128:
                            break
                        elif not args.no_f128:
                            use_f128 = True
                        else:
                            break

                    do_debug = False
                    if not args.scale_cuts:
                        t_sols = [LU_scaling[0] for LU_scaling in LU_scalings]
                        t_variance = (max(t_sols)-min(t_sols))/((max(t_sols)+min(t_sols))/2.0)
                        if t_variance > 1.0e-4:
                            if args.show_warnings:
                                logger.critical("UV profiling of SG '%s' with cut ID #%d with UV edges %s and fixed edges %s gives non-constant LU scaling. Found rel. vaiance of: %.5f over %d points."%(
                                    SG_name, cut_ID, UV_edges_str, fixed_edges_str, t_variance, len(LU_scalings)
                                ))
                                do_debug = True

                    if not successful_fit and (dod > float(args.target_scaling)-2.0):
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
    display_parser = ArgumentParser()
    display_parser.add_argument('SG_name', metavar='SG_name', type=str, nargs='?',
                    help='the name of a supergraph to display')
    display_parser.add_argument(
        '-t','--timing',action="store_true", dest="timing", default=False,
        help="Show timing profile information")
    display_parser.add_argument(
        '--uv',action="store_true", dest="uv", default=False,
        help="Show UV profile information")
    display_parser.add_argument(
        "-f","--full", action="store_true", dest="full", default=False,
        help="exhaustively show information")
    def do_display(self, line):
        """ display command """
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
    set_hyperparameter_parser = ArgumentParser()
    set_hyperparameter_parser.add_argument('param_value', metavar='param_value', type=str, nargs=2,
                    help='parameter name and value to set')
    set_hyperparameter_parser.add_argument(
        "-w", "--write", action="store_true", dest="write", default=False,
        help="Write hyperparameter to disk.")
    def do_set_hyperparameter(self, line):
        """ set_hyperparameter command """
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
    integrate_parser = ArgumentParser()
    integrate_parser.add_argument('SG_name', metavar='SG_name', type=str, nargs='?',
                    help='the name of a supergraph to display')
    integrate_parser.add_argument('-s','--sampling', metavar='sampling', type=str, default='flat', 
                    choices=('flat', 'advanced', 'test_h_function'), help='Specify the sampling method (default: %(default)s)')
    integrate_parser.add_argument('-i','--integrator', metavar='integrator', type=str, default='vegas3', 
                    choices=('naive','vegas', 'vegas3'), help='Specify the integrator (default: %(default)s)')
    integrate_parser.add_argument('-hf','--h_function', metavar='h_function', type=str, default='left_right_polynomial', 
                    choices=('left_right_polynomial',), help='Specify the h-function to use (default: %(default)s)')
    integrate_parser.add_argument('-hfs','--h_function_sigma', metavar='h_function_sigma', type=int, default=3, 
                    help='Spread of the h-function, higher=steeper (default: %(default)s).')
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
    integrate_parser.add_argument('--n_max', metavar='n_max', type=int, default=-1,
                    help='Maximum number of sample points in Vegas (default: as per hyperparameters).')
    integrate_parser.add_argument('--n_start', metavar='n_start', type=int, default=-1,
                    help='Starting number of sample points in Vegas (default: as per hyperparameters).')
    integrate_parser.add_argument('--n_increase', metavar='n_increase', type=int, default=-1,
                    help='Starting number of sample points in Vegas (default: as per hyperparameters).')
    integrate_parser.add_argument('-bs','--batch_size', metavar='batch_size', type=int, default=-1,
                    help='Batch size for parallelisation (default: as per hyperparameters n_vec).')
    integrate_parser.add_argument('--seed', metavar='seed', type=int, default=0,
                    help='Specify the random seed for the integration (default: %(default)s).')
    integrate_parser.add_argument(
        '-mc','--multichanneling',action="store_true", dest="multichanneling", default=False,
        help="Enable multichanneling (default: as per hyperparameters)")
    # We must wrape this function in a process because of the border effects of the pyO3 rust Python bindings
    @wrap_in_process()
    @with_tmp_hyperparameters({
        'Integrator.dashboard': False
    })
    def do_integrate(self, line):
        """ Integrate a given (set of) supergraphs using different sampling strategies."""
        
        args = self.split_arg(line)
        args = self.integrate_parser.parse_args(args)
        if args.h_function == 'left_right_polynomial':
            selected_h_function = sampler.HFunction(args.h_function_sigma, debug=args.verbosity)
        else:
            raise alphaLoopInvalidRunCmd("Unsupported h-function specification: %s'."%args.h_function)

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
                 'n_points_per_iterations' : args.n_points_survey,
                 'verbosity' : args.verbosity+1
            }
        elif args.integrator == 'vegas':
            selected_integrator = pyCubaIntegrator.pyCubaIntegrator
            integrator_options = {
                 'cluster' : runner,
                 'max_eval' : args.max_eval,
                 'n_start' : args.n_start,
                 'n_increase' : args.n_increase,
                 'n_batch' : args.batch_size
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

        else:
            raise alphaLoopInvalidRunCmd("Unsupported integrator specification: %s'."%args.integrator)

        if args.sampling == 'test_h_function':

            logger.info("Dummy integrand for testing h function:")

            my_integrand = sampler.TestHFuncIntegrand( selected_h_function, debug=args.verbosity )

            my_integrator = selected_integrator(my_integrand, **integrator_options)

            result = my_integrator.integrate()

            logger.info('')
            logger.info("Result of the test integration using h function '%s' with %d function calls (target is 1.0): %.4e +/- %.2e"%(
                args.h_function, my_integrator.tot_func_evals, result[0], result[1]))
            logger.info('')
            return

        dimensions = integrands.DimensionList(
                [ integrands.ContinuousDimension('x_%d'%i,lower_bound=0.0, upper_bound=1.0) for i in range(1, n_loops*3) ]+
                [ integrands.ContinuousDimension('t',lower_bound=0.0, upper_bound=1.0), ] 
                )
        my_PS_generator = None
        try:
            my_PS_generator = eval("generator_%s(dimensions, rust_worker, SG_info, model, selected_h_function, hyperparameters, debug=args.debug)"%args.run_mode)
        except Exception as e:
            logger.critical("Unreckognized run mode: '%s'"%args.run_mode)
            raise

        logger.info("Integrating '%s' with the parameterisation %s:"%(args.diag_name,args.run_mode))

        SG = self.all_supergraphs[args.SG_name]
        E_cm = SG.get_E_cm(self.hyperparameters)
        rust_worker = self.get_rust_worker(args.SG_name)

        my_integrand = integrator.DefaultALIntegrand(rust_worker, my_PS_generator, debug=args.verbosity)

        my_integrator = integrator.vegas3.Vegas3Integrator(my_integrand, 
                n_points_survey=args.n_points_survey, n_points_refine=args.n_points_refine, accuracy_target=None,
                verbosity=args.verbosity, cluster=runner
        )

        result = my_integrator.integrate()

        logger.info('')
        logger.info("Result of integration using the '%s' parameterisation with %d function calls: %.7e +/- %.2e"%(
            args.run_mode, my_integrator.tot_func_evals, result[0], result[1]))
        logger.info('')



    #### EXPERIMENT COMMAND
    experiment_parser = ArgumentParser()
    experiment_parser.add_argument('SG_name', metavar='SG_name', type=str, nargs='?',
                    help='The name of a supergraph to consider')
    experiment_parser.add_argument(
        '-e','--experiment', dest="experiment", type=str, default='default',
        help="Which experiment to run")
    # We must wrape this function in a process because of the border effects of the pyO3 rust Python bindings
    @wrap_in_process()
    @with_tmp_hyperparameters({
        'Integrator.dashboard': False
    })
    def do_experiment(self, line):
        """ Integrate a given (set of) supergraphs using different sampling strategies."""
        
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