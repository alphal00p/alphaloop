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
from argparse import ArgumentParser
from pprint import pprint, pformat
import progressbar

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
            'CrossSection.numerator_source'                      :   'yaml',
            # Can be LTD, PF
            'CrossSection.integrand_type'                        :   'PF',
            # evaluate the C expression for the sum of diagrams
            'CrossSection.sum_diagram_sets'                      :   True,
            # compare locally against the same topology written in another loop momentum basis
            'CrossSection.compare_with_additional_topologies'    :   False,
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
        return [
                [ random.random() for _ in range(3) ]
            for _ in range(self['topo']['n_loops'])
        ]

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
                res_str.append('%s%-40s%s : %s'%(Colours.GREEN,k,Colours.END,v))
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
            print(args[1:],args[0])
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
    timing_profile_parser.add_argument("-t","--time", dest='time', type=float, default=5.0,
                    help='target evaluation time per profile, in seconds.')
    timing_profile_parser.add_argument(
        "-f", "--f128", action="store_true", dest="f128", default=False,
        help="Enable timing profile of the f128 output too.")
    # We must wrape this function in a process because of the border effects of the pyO3 rust Python bindings
    @wrap_in_process()
    @with_tmp_hyperparameters()
    def do_timing_profile(self, line):
        """ Automatically timing profile a process output."""
        args = self.split_arg(line)
        args = self.timing_profile_parser.parse_args(args)
        
        if args.SG_name is None:
            selected_SGs = list(self.all_supergraphs.keys())
        else:
            selected_SGs = [args.SG_name,]
        
        max_count = sum( len(self.all_supergraphs[SG_name]['cutkosky_cuts']) for SG_name in selected_SGs )
        logger.info("Starting timing profile...")
        # WARNING it is important that the rust workers instantiated only go out of scope when this function terminates
        rust_workers = {SG_name: self.get_rust_worker(SG_name) for SG_name in selected_SGs}
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
                SG = self.all_supergraphs[SG_name]
                E_cm = SG.get_E_cm(self.hyperparameters)
                funcs_to_test = [('f64', rust_worker.evaluate)]
                if args.f128:
                    funcs_to_test.append(('f128', rust_worker.evaluate_f128))

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
                        rust_function(SG.get_random_x_input())
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
                            scaling_solutions = list(get_scaling_function(random_momenta,cut_ID,))
                            scaling, scaling_jacobian = scaling_solutions.pop(0)
                            while scaling < 0.0:
                                if len(scaling_solutions)==0:
                                    break
                                scaling, scaling_jacobian = scaling_solutions.pop(0)
                            cut_evaluate_function(random_momenta,cut_ID,scaling,scaling_jacobian)
                            delta_t = time.time()-t_start
                            n_points = int(args.time/delta_t)
                        else:
                            n_points = args.n_points

                        t_start = time.time()
                        for _ in range(n_points):
                            random_momenta = SG.get_random_momenta_input(E_cm)
                            scaling_solutions = list(get_scaling_function(random_momenta,cut_ID,))
                            scaling, scaling_jacobian = scaling_solutions.pop(0)
                            while scaling < 0.0:
                                if len(scaling_solutions)==0:
                                    break
                                scaling, scaling_jacobian = scaling_solutions.pop(0)
                            cut_evaluate_function(random_momenta,cut_ID,scaling,scaling_jacobian)
                        delta_t = time.time()-t_start
                        t_cut = (delta_t/float(n_points))*1000.0
                        running_avg_time_per_cut += t_cut
                        n_cuts += 1
                        bar.update(t='%.3f'%(running_avg_time_per_cut/float(n_cuts)))
                        bar.update(bar.value+1)

                        cut['DERIVED_timing_profile_%s'%precision] = t_cut

                misc.sprint("At iSG=%d"%i_SG)
                time.sleep(1.0)

        delta_t = time.time()-t_start_profile
        logger.info("Timing profile completed in %d [s]."%(int(delta_t)))
        # Write out the results into processed topologies
        self.all_supergraphs.export(pjoin(self.dir_path, self._rust_inputs_folder))
        if len(selected_SGs)==1:
            self.do_display('%s --timing'%selected_SGs[0])
        else:
            self.do_display('--timing')

    #### DISPLAY COMMAND
    display_parser = ArgumentParser()
    display_parser.add_argument('SG_name', metavar='SG_name', type=str, nargs='?',
                    help='the name of a supergraph to display')
    display_parser.add_argument(
        '-t','--timing',action="store_true", dest="timing", default=False,
        help="Show timing profile information")
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