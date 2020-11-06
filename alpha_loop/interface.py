#####################################################
#                                                   #
#  Source file of the alphaLoop MG5aMC plugin.      #
#                                                   #
#####################################################

import os
import logging
import sys
import shutil
import re
import random
import sympy
import math
import timeit
import functools
import copy
import resource
import traceback

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

_launch_usage = "launch [DIRPATH] [options]\n" + \
         "-- integrate or test the alphaLoop output at DIRPATH.\n" + \
         "   Example: launch MY_ALPHALOOP_OUTPUT \n"
_launch_parser = misc.OptionParser(usage=_launch_usage)
_launch_parser.add_option("-d","--debug", dest="debug", default=0, type=int,
        help="debug level")
_launch_parser.add_option("-r", "--reuse_hyperparametes", dest="reuse", type=str, default="NO_REUSE",
        help="Work from existing yaml hyperparameter (value: 'default') file in workspace or specify path to one.")

import alpha_loop.utils as utils
import alpha_loop.exporters as aL_exporters
import alpha_loop.helas_call_writers as aL_helas_call_writers
import alpha_loop.LTD_squared as LTD_squared
import alpha_loop.madgraph_patches as madgraph_patches
import alpha_loop.FORM_processing as FORM_processing
import alpha_loop.run_interface as run_interface

from madgraph.iolibs.files import cp, ln, mv

logger = logging.getLogger('alphaLoop.Interface')

pjoin = os.path.join
template_dir = pjoin(plugin_path, 'Templates')

class alphaLoopInterfaceError(MadGraph5Error):
    """ Error for the alphaLoop plugin """
    pass

class alphaLoopInvalidCmd(InvalidCmd):
    """ Invalid command issued to the alphaLoop interface. """
    pass

class alphaLoopInterface(madgraph_interface.MadGraphCmd, cmd.CmdShell):
    """ Interface for steering the generation/output of alphaLoop.
    We make it inherit from CmdShell so that launch_ext_prog does not attempt to start in WebMode."""

    _supported_FORM_output_formats = [None, 'c']

    def __init__(self, *args, **opts):
        """ Define attributes of this class."""
        
        self.alphaLoop_options = {
            # Set a soft limit for the vitual memory, default is no limit
            'virtual_memory' : None,
            'perturbative_orders' : { 'QCD': 0 },
            # Offers the possibility of overwriting the requirement of final_state_pdgs. Set None to have this automatically set.
            'final_state_pdgs' : None,
            # Offers the possibility of overwriting n_jets. Set to None to have this automatically set.
            'n_jets' : None,
            # By default we don't want denominators as those will be included in the rust backend
            'include_denominators' : False, 
            'use_physical_gluon_helicity_sum' : False,
            # Specify the number of rust inputs to generate, -1 considers all by default.
            # This is useful for debugging as this can be slow.
            'n_rust_inputs_to_generate' : -1,
            # We plan on include self-energies from two-point vertices so as to be able to apply our 
            # LTD^2 LSZ treatment. So 'include_self_energies_from_squared_amplitudes' should be kept
            # set to False, except for debugging.
            'include_self_energies_from_squared_amplitudes' : True,
            # We must not apply the graph ismorphism checks when attempting to compute LO cross-section from MG.
            'apply_graph_isomorphisms' : True,
            # Instead of None, the option below can be a list of supergraphs identified as a two-tuple of the 
            # of the two LTD2Diagram ids making them up on the left and right of the Cutkosky cut respectively.
            # For example [(2,2),(1,49)] will select only the supergraphs built from the 2x2 and 1x49 sewings of LT2Diagrams.
            'supergraphs_selection' : None,
            # One must differentiate particles from anti-particles in the isomorphism check.
            # However this is not properly working, at least not for self-energies, so we allow it
            # to be disabled with the option below.
            'differentiate_particle_from_antiparticle_in_graph_isomorphism' : False,
            # We do not want to group numerators when considering MG outputs
            'consider_edge_orientation_in_graph_isomorphism' : False,
            'consider_vertex_id_in_graph_isomorphism' : False,
            # Select which templete to use for qgraf generation:
            'qgraf_template_model': 'epem',
            # filter qgraf output based on cuts
            'qgraf_cut_filter': False,
            # Specify qgraf model
            'qgraf_model' : 'SM',
            # Loop induced processes will enforce the presence of at least two virtual corrections:
            'loop_induced': False,
            # Veto some field from the QGRAF generation
            'qgraf_vetoes': None,
            # Set the output processing format of Rust to `None` if you want to skip it.
            # Otherwise it can take values in ['rust',] for now.
            'FORM_processing_output_format' : None,
            # Select if the FORM integrand is the "PF" expression or "LTD" expression, both, or None
            'FORM_integrand_type' : "both",
            # Select if the sum of diagram sets should be generated with 'nosum', 'onlysum', 'both'
            'FORM_sum_diagram_sets': "nosum",
            'FORM_construct_numerator': False,
            # Limit nuber of cores used for compilation
            'FORM_compile_cores': None, 
            # Select what to compile of the FORM output. Default is ['all']
            'FORM_compile_arg' : ['all'],
            # Optimization level for the FORM processing and compilation. 0 to 3
            'FORM_compile_optimization' : 3
        }
        self.FORM_options=FORM_processing.FORM_processing_options
        self.plugin_output_format_selected = None

        self.model_backup_copy = None

        self.qgraf_exporter = None

        super(alphaLoopInterface, self).__init__(*args, **opts)

    def parse_set_option(self, args):
        """ Parsing arguments/options passed to the command set_alphaLoop option."""

        options = { }

        # First combine all value of the options (starting with '--') separated by a space
        opt_args = []
        new_args = []
        for arg in args:
            if arg.startswith('--'):
                opt_args.append(arg)
            elif len(opt_args) > 0:
                opt_args[-1] += ' %s' % arg
            else:
                new_args.append(arg)

        for arg in opt_args:
            try:
                key, value = arg.split('=')
            except:
                key, value = arg, None
            key = key[2:]

            # All options are declared valid in this contex
            options[key] = eval(str(value))

        return new_args, options

    def do_display_alphaLoop_option(self, line):
        """ Display alphaLoop options"""
        logger.info('%sGeneral alphaLoop options%s'%(utils.bcolors.GREEN, utils.bcolors.ENDC))
        logger.info('%s-----------------------%s'%(utils.bcolors.GREEN, utils.bcolors.ENDC))
        for opt in sorted(self.alphaLoop_options.keys()):
            logger.info('%-30s : %s'%(opt, str(self.alphaLoop_options[opt])))

    def do_display_FORM_option(self, line):
        """ Display FORM options"""
        logger.info('%sGeneral FORM processing options%s'%(utils.bcolors.GREEN, utils.bcolors.ENDC))
        logger.info('%s-----------------------%s'%(utils.bcolors.GREEN, utils.bcolors.ENDC))
        for opt in sorted(self.FORM_options.keys()):
            logger.info('%-30s : %s'%(opt, str(self.FORM_options[opt])))        

    def do_set_FORM_option(self, line):
        """ Logic for setting alphaLoop options."""
        args = self.split_arg(line)
        args, options = self.parse_set_option(args)
        key, value = args[:2]

        if key in self.FORM_options:
            if key=='renormalisation_finite_terms':
                processed_value=str(value)
            else:
                processed_value=eval(value)
            if isinstance(self.FORM_options[key],dict):
                self.FORM_options[key].update(processed_value)
            else:
                self.FORM_options[key] = processed_value
        else:
            raise alphaLoopInvalidCmd("Specified FORM option '%s' not reckognized."%key)

    def do_set_alphaLoop_option(self, line):
        """ Logic for setting alphaLoop options."""
        args = self.split_arg(line)
        args, options = self.parse_set_option(args)
        key, value = args[:2]
        
        if key == 'virtual_memory':
            vlimit_format = re.compile(r'^\d+(|B|K|M|G|T)$')
            vlimit_conv = {'B': 1, 'K': 1024, 'M': 1024**2, 'G': 1024**3, 'T': 1024**3}
            if vlimit_format.match(value) is None and value != '-1':
                raise alphaLoopInvalidCmd("Virtual memory limit had to be an integer in bytes (possible explicit units B, K, M, G, T)")
            if value[-1] in vlimit_conv.keys():
                vlimit = int(value[:-1]) * vlimit_conv[value[-1]]
            else:
                vlimit = int(value)
            # set soft limit, this applies to every subprocess as well
            # Not sefe enough for 'make -j NCORES'
            soft_limit,hard_limit=resource.getrlimit(resource.RLIMIT_AS)
            resource.setrlimit(resource.RLIMIT_AS, (vlimit,hard_limit))
            #print(resource.getrlimit(resource.RLIMIT_AS))
            #os.system('ulimit -v')
        elif key == 'perturbative_orders':
            try:
                orders_dict = eval(value)
            except:
                raise alphaLoopInvalidCmd("Specified perturbative orders should be a dictionary of orderrs, not '%s'."%value)
            for order,v in orders_dict.items():
                if order not in self._curr_model['order_hierarchy'].keys():
                    raise alphaLoopInvalidCmd("Specified order is not part of the model '%s'."%order)
                if v%2 != 0:
                    raise alphaLoopInvalidCmd("Specified perturbative orders shoudl be multiple of 2, not '%d'."%v)
            self.alphaLoop_options[key] = orders_dict
        elif key == 'n_jets':
            if value.upper()=='NONE':
                value = None
            else:
                try:
                    value = int(value)
                except ValueError:
                    raise alphaLoopInvalidCmd("alphaLoop option 'n_jets' should either be None or an integer, not '%s'."%value)
            self.alphaLoop_options[key] = value
        elif key == 'final_state_pdgs':
            if value.upper()=='NONE':
                fs_pdgs = None
            else:
                try:
                    fs_pdgs = eval(value)
                except:
                    raise alphaLoopInvalidCmd("alphaLoop option 'final_state_pdgs' should either be None or a list/tuple of integer, not '%s'."%value)
            self.alphaLoop_options[key] = tuple(fs_pdgs)
        elif key == 'include_denominators':
            if value.upper() not in ['TRUE','FALSE']:
                raise alphaLoopInvalidCmd("Specified value for 'include_denominators' should be 'True' or 'False', not '%s'."%value)
            
            bool_val = (value.upper()=='TRUE')
            self.alphaLoop_options['include_denominators'] = bool_val
            if bool_val:
                madgraph_patches.remove_propagator_denominators( restore=True )
            else:
                madgraph_patches.remove_propagator_denominators( restore=False )
        elif key == 'use_physical_gluon_helicity_sum':
            if value.upper() not in ['TRUE','FALSE']:
                raise alphaLoopInvalidCmd("Specified value for 'use_physical_gluon_helicity_sum' should be 'True' or 'False', not '%s'."%value)
            bool_val = (value.upper()=='TRUE')
            self.alphaLoop_options['use_physical_gluon_helicity_sum'] = bool_val
        elif key == 'loop_induced':
            if value.upper() not in ['TRUE','FALSE']:
                raise alphaLoopInvalidCmd("alphaLoop option 'loop_induced' should be 'True' or 'False', not %s"%value)
            self.alphaLoop_options['loop_induced'] = value.upper() == 'TRUE'
        elif key == 'qgraf_cut_filter':
            if value.upper() not in ['TRUE','FALSE']:
                raise alphaLoopInvalidCmd("Specified value for '%s' should be 'True' or 'False', not '%s'."%(key,value))
            bool_val = (value.upper()=='TRUE')
            self.alphaLoop_options[key] = bool_val
        elif key == 'qgraf_template_model':
            if value not in aL_exporters.HardCodedQGRAFExporter.qgraf_templates.keys():
                raise alphaLoopInvalidCmd("QGraf template model '{}' not supported.\nTry models (: example)\n{}"\
                    .format(value,"\n".join(["\t%s: %s"%(k,v['example']) for k,v in aL_exporters.HardCodedQGRAFExporter.qgraf_templates.items()]))
                )
            self.alphaLoop_options['qgraf_template_model'] = value
        elif key == 'qgraf_model':
            self.alphaLoop_options['qgraf_model'] = value
        elif key == 'FORM_compile_cores':
            if not value.isdigit():
                raise alphaLoopInvalidCmd("alphaLoop option 'FORM_compile_cores' takes values form 1 to %d. Not %s" % (multiprocessing.cpu_count(), value))
            cores = int(value)
            if cores > multiprocessing.cpu_count():
                cores = multiprocessing.cpu_count()
                logger.warning("Asking from more CORES then available. Overwriting the value to %s"%cores)
            self.FORM_options['cores'] = cores
        elif key == 'FORM_compile_arg':
            if not value in ['integrand', 'numerator', 'all']:
                raise alphaLoopInvalidCmd("alphaLoop option 'FORM_compile_output' should be one of 'numerator', 'integrand', 'all'")
            self.FORM_options['compilation-options'] += [value]
        elif key == 'FORM_compile_optimization':
            if not value in [str(opt_n) for opt_n in range(4)]:
                raise alphaLoopInvalidCmd("alphaLoop option 'FORM_compile_optimization' should be between 0 and 3, not %s"%value)
            self.FORM_options['compilation-options'] += ['-e', "OPTIMIZATION_LVL=%s"%value]
            self.FORM_options['extra-options']['OPTIMLVL'] = [1,2,4,4][int(value)]
        elif key == 'FORM_compile_optimization_f128':
            if not value in [str(opt_n) for opt_n in range(4)]:
                raise alphaLoopInvalidCmd("alphaLoop option 'FORM_compile_optimization_f128' should be between 0 and 3, not %s"%value)
            self.FORM_options['compilation-options'] += ['-e', "OPTIMIZATION_LVL_f128=%s"%value]
        elif key == 'FORM_integrand_type':
            if not value in ('PF', 'LTD', 'both', 'None'):
                raise alphaLoopInvalidCmd("alphaLoop option 'FORM_integrand_type' should be one of 'PF', 'LTD', 'None', not %s"%value)
            self.alphaLoop_options['FORM_integrand_type'] = value if value in ('PF', 'both', 'LTD') else None
        elif key == 'FORM_sum_diagram_sets':
            if not value in ('nosum', 'onlysum', 'both'):
                raise alphaLoopInvalidCmd("alphaLoop option 'FORM_sum_diagram_sets' should be one of 'nosum', 'onlysum', 'both', not %s"%value)
            self.alphaLoop_options['FORM_sum_diagram_sets'] = value if value in ('nosum', 'onlysum', 'both') else None
            self.FORM_options['extra-options']['SUMDIAGRAMSETS'] = value
        elif key == 'FORM_construct_numerator':
            if value.upper() not in ['TRUE','FALSE']:
                raise alphaLoopInvalidCmd("Specified value for '%s' should be 'True' or 'False', not '%s'."%(key,value))
            self.FORM_options['extra-options']['NUMERATOR'] = '1' if value.upper() == 'TRUE' else '0'
        elif key == 'include_self_energies_from_squared_amplitudes':
            if value.upper() not in ['TRUE','FALSE']:
                raise alphaLoopInvalidCmd("Specified value for 'include_self_energies_from_squared_amplitudes' should be 'True' or 'False', not '%s'."%value)
            bool_val = (value.upper()=='TRUE')
            if bool_val and len(self.alphaLoop_options['perturbative_orders'])>0 and sum(self.alphaLoop_options['perturbative_orders'].values())>0:
                logger.warning(
"""%sSelf-energies will be generated from two-point vertices so as to be able to apply the 
LTD^2 LSZ treatment. So 'include_self_energies_from_squared_amplitudes' should be kept
set to False, except for debugging, which seems to be what you are doing now, so we'll se it to True now.%s"""%(
utils.bcolors.RED,utils.bcolors.ENDC
))
            self.alphaLoop_options['include_self_energies_from_squared_amplitudes'] = bool_val   
        elif key in [
            'differentiate_particle_from_antiparticle_in_graph_isomorphism',
            'consider_edge_orientation_in_graph_isomorphism',
            'consider_vertex_id_in_graph_isomorphism',
            'apply_graph_isomorphisms'
            ]:
            if value.upper() not in ['TRUE','FALSE']:
                raise alphaLoopInvalidCmd("Specified value for '%s' should be 'True' or 'False', not '%s'."%(key,value))
            bool_val = (value.upper()=='TRUE')
            self.alphaLoop_options[key] = bool_val
        elif key == 'FORM_processing_output_format':
            if value.upper() in ['NONE']:
                value = None
            if value not in self._supported_FORM_output_formats:
                raise alphaLoopInvalidCmd("Specified value '%s' for 'FORM_processing_output_format' is not in %s."%(
                                                                        value,self._supported_FORM_output_formats))
            self.alphaLoop_options['FORM_processing_output_format'] = value
        elif key == 'supergraphs_selection':
            if value.upper() in ['NONE']:
                parsed_value = None
            else:
                try:
                    parsed_value = eval(value)
                except:
                    parsed_value = value
            if not parsed_value is None:
                if (not (isinstance(parsed_value,list) or isinstance(parsed_value,tuple))) or \
                    (not all(len(el)==2 for el in parsed_value)) or \
                    (not all(isinstance(el[0],int) and isinstance(el[1],int) for el in parsed_value)):
                    raise alphaLoopInvalidCmd("Invalid value for option supergraphs_selection: '%s'"%value)
            self.alphaLoop_options['supergraphs_selection'] = parsed_value

        elif key == 'n_rust_inputs_to_generate':
            try:
                self.alphaLoop_options['n_rust_inputs_to_generate'] = int(value)
            except ValueError:
                raise alphaLoopInvalidCmd("Specified value for 'n_rust_inputs_to_generate' should be an integer, not '%s'."%value)
        else:
            raise alphaLoopInvalidCmd("Unrecognized alphaLoop option: %s"%key)

    def do_prepare_model_for_self_energies(self, line):
        """Force the processing of the model to accomodate self-energy generation."""

        args = self.split_arg(line)
        if len(args)!=2:
            raise alphaLoopInvalidCmd("Supply exactly two arguments to the command 'do_prepare_model_for_self_energies'.")
        else: 
            if args[0]=='revert':
                if self.model_backup_copy is None:
                    logger.information("Nothing to revert to.")
                else:
                    self.prepare_model_for_self_energies(None,revert=True)
            else:
                try:
                    n_SE = int(args[0])
                except ValueError:
                    raise alphaLoopInvalidCmd("Supply a single integer argument to the command 'do_prepare_model_for_self_energies'"
                                             +" (number of self-energies to accomodate), and not '%s'."%args[0])
                try:
                    anchor_vertices = eval(args[1])
                except ValueError:
                    raise alphaLoopInvalidCmd("The second argument of 'do_prepare_model_for_self_energies' "+
                                                "must be a list of anchor vertices, not '%s'."%args[1])
                self.prepare_model_for_self_energies(anchor_vertices,revert=False, perturbative_order=n_SE)

    def prepare_model_for_self_energies(self, anchor_vertices, revert=False, perturbative_order=3):
        """Commands for adding the necessary TREE UV interactions to the current model.
        The anchor vertices are the selected location for attaching the bridge of each self-energy.
        """

        if perturbative_order<=0:
            return

        if revert:
            if self.model_backup_copy is None:
                raise MadGraph5Error("Cannot revert to a non-existing model.")
            logger.info('Self-energy patch for model %s successfully reverted.'%self.model_backup_copy.get('name'))
            self._curr_model = self.model_backup_copy
            return

        # First create a backup of the current model.
        self.model_backup_copy = copy.deepcopy(self._curr_model)        

        # Short-hand to the current model that will be modified
        model = self._curr_model

        lorentzStructClass = model['lorentz'][0].__class__
        new_lorentz_structures = {'SE_SSS':
            lorentzStructClass(name = 'SE_SSS',
               spins = [ 1, 1, 1 ],
               structure = '1')
        }

        pure_QCD_corrections_only = list(self.alphaLoop_options['perturbative_orders'].keys())==['QCD',]

        # First add the necessary self-energy coupling_orders to the available coupling_orders 
        # in the model
        for n_self_energy in range(1,perturbative_order+1):
            if (not model['order_hierarchy']) or any(order in [
                'SE_%d_ANCHOR'%n_self_energy,
                ] for order in model['order_hierarchy']):
                raise MadGraph5Error("Cannot prepare the currently active model for self-energies.")
            model['order_hierarchy']['SE_%d_ANCHOR'%n_self_energy]=0
            for order in self.alphaLoop_options['perturbative_orders']:
                model['order_hierarchy']['SE_%d_%s'%(n_self_energy,order)]=model['order_hierarchy'][order]
            for n_bridge in range(1,perturbative_order+1):
                if (not model['order_hierarchy']) or any(order in [
                    'SE_%d_BRIDGE_%d'%(n_self_energy,n_bridge),
                    'SE_%d_CONTACT_%d'%(n_self_energy,n_bridge)
                    ] for order in model['order_hierarchy']):
                    raise MadGraph5Error("Cannot prepare the currently active model for self-energies.")
                model['order_hierarchy']['SE_%d_BRIDGE_%d'%(n_self_energy,n_bridge)]=0
                model['order_hierarchy']['SE_%d_CONTACT_%d'%(n_self_energy,n_bridge)]=0

        new_order_hierarchy = dict(model['order_hierarchy'])
        new_perturbation_couplings = list(model['perturbation_couplings'])

        # Then add the necessary new particles
        new_particles = {}
        pdg_offset = LTD_squared.self_energy_global_id_offset
        bridge_offset = pdg_offset*100
        anchor_offset = pdg_offset*10
        for n_self_energy in range(1,perturbative_order+1):
            new_particles[n_self_energy] = {}

            # Let us create self-energy specific particles
            for particle in model['particles']:
                if pure_QCD_corrections_only:
                    if particle.get('color')==1:
                        continue
                new_particle=copy.deepcopy(particle)
                new_particle['pdg_code'] = particle['pdg_code'] + n_self_energy*pdg_offset
                new_particle['name'] = ('SE_%d_'%n_self_energy+particle['name']).lower()
                new_particle['antiname'] = ('SE_%d_'%n_self_energy+particle['antiname']).lower()
                new_particles[n_self_energy][particle.get_pdg_code()]=new_particle
                if not new_particle['self_antipart']:
                    new_anti_particle=copy.deepcopy(new_particle)
                    new_anti_particle['is_part'] = not new_particle['is_part']
                    new_particles[n_self_energy][particle.get_anti_pdg_code()]=new_anti_particle

            for n_anchor in range(0,perturbative_order+1):
                anchor = base_objects.Particle({
                    'pdg_code' : n_self_energy*pdg_offset+(n_anchor+1)*anchor_offset,
                    'name' : ('SE_%d_A_%d'%(n_self_energy,n_anchor)).lower(),
                    'antiname' : ('SE_%d_A_%d'%(n_self_energy,n_anchor)).lower(),
                    'spin': 1,
                    'color' : 1,
                    'is_part' : True,
                    'self_antipart' : True,
                    'charge' : 0.
                })
                new_particles[n_self_energy]['anchor_%d'%n_anchor] = anchor

            for n_bridge in range(1,perturbative_order+1):
                bridge = base_objects.Particle({
                    'pdg_code' : n_self_energy*pdg_offset+n_bridge*bridge_offset,
                    'name' : ('SE_%d_B_%d'%(n_self_energy,n_bridge)).lower(),
                    'antiname' : ('SE_%d_B_%d'%(n_self_energy,n_bridge)).lower(),
                    'spin': 1,
                    'color' : 1,
                    'is_part' : True,
                    'self_antipart' : True,
                    'charge' : 0.
                })
                new_particles[n_self_energy]['bridge_%d'%n_bridge] = bridge

        # Assign the new particles
        model.set('particles',base_objects.ParticleList(
            model.get('particles')+[p for se_parts in new_particles.values() 
                for p in se_parts.values() if p['is_part']] ))

        # Then create the new interactions
        interaction_offset = LTD_squared.self_energy_global_id_offset
        anchor_interaction_id_offset = interaction_offset*10
        dark_sector_id_offset = interaction_offset*100
        bridge_offset = interaction_offset*1000

        # First deterrmine that anchor vertices
        anchor_base_vertices = []
        all_anchor_vertices = []
        sorted_selected_anchor_vertices = [sorted(pdgs) for pdgs in anchor_vertices]
        for inter in model['interactions'].get_type('base'):
            particle_pdgs = sorted([p.get_pdg_code() for p in inter.get('particles')])
            if particle_pdgs not in sorted_selected_anchor_vertices:
                continue
            anchor_base_vertices.append(inter)

        new_interactions = []
        new_SE_sectors_interactions = {}
        for n_self_energy in range(1,perturbative_order+1):
            # Dubplicate base interactions so as to add the the anchor scalar #1 to it.
            new_SE_sectors_interactions[n_self_energy] = []
            for inter in model['interactions'].get_type('base'):
                if (set(inter['orders'].keys())-set(self.alphaLoop_options['perturbative_orders']))!=set([]):
                    continue
                new_interaction = copy.deepcopy(inter)
                new_interaction.set('color',inter.get('color'))
                new_orders = {}
                for order in new_interaction['orders']:
                    new_orders['SE_%d_%s'%(n_self_energy,order)]=new_interaction['orders'][order]
                new_orders['SE_%d_ANCHOR'%n_self_energy] = 1
                new_interaction.set('orders',new_orders)

                non_dark_particles_selected_already = set([])
                for i_leg in range(len(new_interaction['particles'])):
                    non_dark_pdg = new_interaction['particles'][i_leg].get_pdg_code()
                    if non_dark_pdg in non_dark_particles_selected_already:
                        continue
                    non_dark_particles_selected_already.add(non_dark_pdg)
                    a_new_interaction = copy.deepcopy(new_interaction)
                    a_new_interaction.set('color',inter.get('color'))
                    a_new_interaction.set('id',
                        n_self_energy*interaction_offset+
                        (i_leg+1)*anchor_interaction_id_offset+
                        inter.get('id')
                    )
                    a_new_interaction['particles'] = base_objects.ParticleList([
                       new_particles[n_self_energy][p.get_pdg_code()] if i_p!=i_leg else p
                            for i_p, p in enumerate(new_interaction['particles'])
                        ]+[new_particles[n_self_energy]['anchor_0'],]
                    )
                    a_new_interaction['lorentz'] = ['SE_%s'%l for l in new_interaction['lorentz']]
                    new_lorentz_names=[]
                    for l in new_interaction['lorentz']:
                        new_lorentz_name = 'SE_%s'%l
                        new_lorentz_names.append(new_lorentz_name)
                        if new_lorentz_name not in new_lorentz_structures:
                            new_lorentz_struct = copy.deepcopy(model.get_lorentz(l))
                            new_lorentz_struct.spins.append(1)
                            new_lorentz_struct.name = 'SE_%s'%new_lorentz_struct.name
                            new_lorentz_structures[new_lorentz_name] = new_lorentz_struct
                    a_new_interaction['lorentz'] = new_lorentz_names
                    new_interactions.append(a_new_interaction)

                # When using the "self-energy particles" above, then also copy down interactions
                # for the whole perturbed sector aming the self-energy particles for a given self-energy index.
                new_interaction = copy.deepcopy(inter)
                new_interaction.set('color',inter.get('color'))
                new_orders = {}
                for order in new_interaction['orders']:
                    new_orders['SE_%d_%s'%(n_self_energy,order)]=new_interaction['orders'][order]
                new_interaction.set('orders',new_orders)
                new_interaction.set('id',
                    n_self_energy*interaction_offset+
                    dark_sector_id_offset+
                    inter.get('id')
                )
                new_interaction['particles'] = base_objects.ParticleList([
                       new_particles[n_self_energy][p.get_pdg_code()] for p in inter['particles']
                        ]
                    )
                new_SE_sectors_interactions[n_self_energy].append(new_interaction)

            # Add an interaction between three scalars: Anchor_<i>_{n-1} | Bridge_<i>_{n} | Anchor_<i>_{n-1}
            for n_bridge in range(1,perturbative_order+1):
                new_interactions.append(base_objects.Interaction({
                    'id' : n_self_energy*interaction_offset+n_bridge*bridge_offset,
                    'particles' : base_objects.ParticleList([
                        new_particles[n_self_energy]['anchor_%d'%(n_bridge-1)],
                        new_particles[n_self_energy]['bridge_%d'%n_bridge],
                        new_particles[n_self_energy]['anchor_%d'%n_bridge],
                    ]),
                    'color': [color.ColorString(),],
                    'lorentz': [ 'SE_SSS',],
                    'couplings' : { (0, 0): 'CMPLX_ONE'},
                    'orders' : {'SE_%d_BRIDGE_%d'%(n_self_energy,n_bridge): 1}
                }))
                # This first bridge only connects to the core graph, only the subsequent ones connect to the inside
                # of the self-energies.
                if n_bridge == 1:
                    continue
                for inter in new_SE_sectors_interactions[n_self_energy]:
                    new_interaction = copy.deepcopy(inter)
                    new_interaction.set('color',inter.get('color'))
                    new_interaction.set('id',
                        n_self_energy*interaction_offset+
                        n_bridge*bridge_offset+
                        inter.get('id')
                    )
                    new_orders = new_interaction['orders']
                    new_orders['SE_%d_CONTACT_%d'%(n_self_energy,n_bridge)] = 1
                    new_interaction.set('orders',new_orders)
                    new_interaction['particles'] = base_objects.ParticleList([
                            p for p in new_interaction['particles']
                        ]+[new_particles[n_self_energy]['bridge_%d'%n_bridge],]
                    )
                    new_lorentz_names=[]
                    for l in new_interaction['lorentz']:
                        new_lorentz_name = 'SE_%s'%l
                        new_lorentz_names.append(new_lorentz_name)
                        if new_lorentz_name not in new_lorentz_structures:
                            new_lorentz_struct = copy.deepcopy(model.get_lorentz(l))
                            new_lorentz_struct.spins.append(1)
                            new_lorentz_struct.name = 'SE_%s'%new_lorentz_struct.name
                            new_lorentz_structures[new_lorentz_name] = new_lorentz_struct
                    new_interaction['lorentz'] = new_lorentz_names
                    new_interactions.append(new_interaction)

            # Now add an anchor for this self-energy
            new_anchor_vertices = []
            for inter in anchor_base_vertices:
                new_interaction = copy.deepcopy(inter)
                new_interaction.set('color',inter.get('color'))
                new_interaction.set('id',
                    n_self_energy*interaction_offset+
                    bridge_offset+
                    inter.get('id')
                )
                new_orders = new_interaction['orders']
                new_orders['SE_%d_CONTACT_%d'%(n_self_energy,1)] = 1
                new_interaction.set('orders',new_orders)
                new_interaction['particles'] = base_objects.ParticleList([
                        p for p in new_interaction['particles']
                    ]+[new_particles[n_self_energy]['bridge_1'],]
                )
                new_lorentz_names=[]
                for l in new_interaction['lorentz']:
                    new_lorentz_name = 'SE_%s'%l
                    new_lorentz_names.append(new_lorentz_name)
                    if new_lorentz_name not in new_lorentz_structures:
                        if l in new_lorentz_structures:
                            new_lorentz_struct = copy.deepcopy(new_lorentz_structures[l])
                        else:
                            new_lorentz_struct = copy.deepcopy(model.get_lorentz(l))
                        new_lorentz_struct.spins.append(1)
                        new_lorentz_struct.name = 'SE_%s'%new_lorentz_struct.name
                        new_lorentz_structures[new_lorentz_name] = new_lorentz_struct
                new_interaction['lorentz'] = new_lorentz_names 
                new_anchor_vertices.append(new_interaction)

            all_anchor_vertices.extend(new_anchor_vertices)
            anchor_base_vertices = new_anchor_vertices

        model.set('interactions',base_objects.InteractionList(
            model.get('interactions')+
            new_interactions+
            sum(new_SE_sectors_interactions.values(),[])+
            all_anchor_vertices
        ))

        model['lorentz'].extend(new_lorentz_structures.values())
        model.create_lorentz_dict()

        # Reset the model dictionaries to sync changes
        model.reset_dictionaries()
        # Refresh dictionaries manually (would be done automatically anyway)
        model.actualize_dictionaries()
    
        model.set('order_hierarchy', new_order_hierarchy)
        model.set('perturbation_couplings', new_perturbation_couplings)


        logger.info('Model %s successfully patched so as to accomodate self-energy generations.'%model.get('name'))
        return

    def check_launch(self, args, options):
        """ Specify here sanity check for launch options."""
        pass

    def do_launch(self, line):
        """ Start an alphaLoop run interface on an existing alphaLoop process output."""

        args = self.split_arg(line)
        # check argument validity and normalise argument
        (launch_options, args) = _launch_parser.parse_args(args)
        self.check_launch(args, launch_options)
        launch_options = launch_options.__dict__

        dir_path = os.path.abspath(pjoin(MG5DIR,args[0]))

        alphaLoop_run_interface = run_interface.alphaLoopRunInterface(
                            dir_path, self, launch_options=launch_options)
        stop = self.define_child_cmd_interface(alphaLoop_run_interface)
        return stop


    def complete_launch(self, text, line, begidx, endidx,formatting=True):
        """ complete the launch command"""
        args = self.split_arg(line[0:begidx])

        # Directory continuation
        if args[-1].endswith(os.path.sep):
            return self.path_completion(text,
                                        pjoin(*[a for a in args if a.endswith(os.path.sep)]),
                                        only_dirs = True)
        # Format
        if len(args) == 1:
            out = {'Path from ./': self.path_completion(text, '.', only_dirs = True)}
            if MG5DIR != os.path.realpath('.'):
                out['Path from %s' % MG5DIR] =  self.path_completion(text,
                                     MG5DIR, only_dirs = True, relative=False)

        #option
        if len(args) >= 2:
            out={}

        # Example of how to provide completion for options
        #if line[0:begidx].endswith('--laststep='):
        #    opt = ['parton', 'pythia', 'pgs','delphes','auto']
        #    out['Options'] = self.list_completion(text, opt, line)
        if False:
            pass
        else:
            # List here all options
            opt = []
            out['Options'] = self.list_completion(text, opt, line)


        return self.deal_multiple_categories(out,formatting)

    def do_qgraf_define(self, line):
        """ define specific multiparticles to be used at generation time. """
        self.do_define(line)
        self.alphaLoop_options['_multiparticles'] = self._multiparticles

    def do_qgraf_generate(self, line):
        """ qgraf-based output. """
 
        args = self.split_arg(line)
        process_definition = self.extract_process(" ".join(args), proc_number=0)
        
#        proc_output_path = pjoin(MG5DIR,args[1])

        # Check if qgraf has been compiled  
        QGRAF_path = pjoin(plugin_path,os.path.pardir,'libraries','QGRAF')
        if not os.path.isfile(pjoin(QGRAF_path,'qgraf')):
            logger.info("Compiling QGRAF (needs to be done once only)...")
            misc.compile(cwd=QGRAF_path)
            if not os.path.isfile(pjoin(QGRAF_path,'qgraf')):
                raise MadGraph5Error("Could not successfully compile QGRAF.")

        # TODO make this 82, -82 general
        self.alphaLoop_options['_jet_PDGs'] = tuple([
            pdg for pdg in self._multiparticles['j'] ]+
            [pdg for pdg in [82,-82] if self.alphaLoop_options['final_state_pdgs'] is None or 
                                        pdg not in self.alphaLoop_options['final_state_pdgs']])
        
        # Now generate the output directory structure:
        self.qgraf_exporter =  aL_exporters.HardCodedQGRAFExporter(
            process_definition, self._curr_model, 
            MG5aMC_options=self.options,alphaLoop_options=self.alphaLoop_options)

    def do_output(self, line):
        """ Wrapper to support the syntax output alphaLoop <args>.
        This just to add extra freedom in adding special action that may be needed at the output
        stage for these output formats.
        """
        args = self.split_arg(line)
        if len(args)>=1 and args[0]=='alphaLoop':
            self.plugin_output_format_selected = 'alphaLoop'
            self.do_output_alphaLoop(' '.join(args[1:]))    
        elif len(args)>=1 and args[0]=='qgraf':
            self.plugin_output_format_selected = 'qgraf'
            self.do_output_qgraf(' '.join(args[1:]))
        else:
            super(alphaLoopInterface,self).do_output(' '.join(args))

    def do_output_alphaLoop(self, line):
        args = self.split_arg(line)
        super(alphaLoopInterface,self).do_output(' '.join(['alphaLoop']+args))
    
    def do_output_qgraf(self, line):
        args = self.split_arg(line)
        print(' '.join(args))
        self.qgraf_exporter.output(' '.join(args))

    def export(self,*args,**opts):
        """Overwrite this so as to force a pythia8 type of output if the output mode is PY8MEs."""
        
        if self._export_format == 'plugin':
            # Also pass on the aloha model to the exporter (if it has been computed already)
            # so that it will be used when generating the model
            if self.plugin_output_format_selected == 'alphaLoop':
                # Set what are the jet pdgs in alphaLoop options
                if 'j' not in self._multiparticles:
                    raise alphaLoopInvalidCmd("alphaLoop requires the 'j' multiparticle label to be defined.")
                # TODO make this 82, -82 general
                self.alphaLoop_options['_jet_PDGs'] = tuple([
                    pdg for pdg in self._multiparticles['j']
                ]+[82,-82])
                self._curr_exporter = aL_exporters.alphaLoopExporter(self._export_dir,
                    alphaLoop_options=self.alphaLoop_options,
                    MG5aMC_options=self.options
                )
            else:
                raise alphaLoopInterfaceError("A plugin output format must have been specified at this stage.")

        super(alphaLoopInterface,self).export(*args, **opts)


    def do_output_LU_scalar(self, line):
        """ Generates and directly output a Local Unitarity computation of the specified scalar topology."""

        args = self.split_arg(line)

        if len(args)<=1:
            raise alphaLoopInvalidCmd("Missing mandatory first argument for output LU scalar which is the output directory name.")
        
        output_path = args.pop(0)

        processed_args = {
            'topology': None,
            'name': 'DefaultLUScalarName',
            'externals': None,
            'lmb' : None,
            'analytical_result': 0.0, 
        }
        mandatory_args = ['topology', 'externals']

        # Group arguments in between the '--' specifiers.
        # For example, it will group '--process=p' 'p' '>' 'd' 'd~' 'z'
        # into '--process=p p > d d~ z'.
        opt_args = []
        new_args = []
        for arg in args:
            if arg.startswith('--'):
                opt_args.append(arg)
            elif len(opt_args)>0:
                opt_args[-1] += ' %s'%arg
            else:
                new_args.append(arg)
        
        for arg in opt_args:
            try:
                key, value = arg.split('=',1)
                value = value.strip()
            except ValueError:
                key = arg
                value = None
        
            if key == '--topology':
                try:
                    topology=eval(value)
                except Exception as e:
                    traceback.print_exc()
                    raise alphaLoopInvalidCmd(
                        "Could not process specified topology: %s"%value+
                        "\n The following Python intepreter error occured then: %s"%str(e))
                if not isinstance(topology, (tuple, list)):
                    raise alphaLoopInvalidCmd("Specified topology must be a tuple or a list")
                if any(
                    ( (not isinstance(te,(tuple,list))) or (not isinstance(te[0],str)) or (not isinstance(te[1],int)) or (not isinstance(te[2],int)) )
                    for te in topology ):
                    raise alphaLoopInvalidCmd("Each element of the topology specified must be a tuple ('edge_name', node_int_id_left, node_int_id_right) ")                        

                processed_args['topology'] = topology

            if key == '--name':
                if "'" in value or '"' in value:
                    processed_args[key[2:]] = eval(value)
                else:
                    processed_args[key[2:]] = value
            
            if key == '--externals':
                try:
                    externals = eval(value)
                except Exception as e:
                    traceback.print_exc()
                    raise alphaLoopInvalidCmd(
                        "Could not process specified externals for this topology: %s"%value+
                        "\n The following Python intepreter error occured then: %s"%str(e))
                if not isinstance(externals, (tuple, list)) or any(not isinstance(el,str) for el in externals):
                    raise alphaLoopInvalidCmd("Each element of the list of externals must be a string identifying the edge name.")                        

                processed_args[key[2:]] = externals

            if key == '--lmb':
                try:
                    lmb = eval(value)
                except Exception as e:
                    traceback.print_exc()
                    raise alphaLoopInvalidCmd(
                        "Could not process specified lmb for this topology: %s"%value+
                        "\n The following Python intepreter error occured then: %s"%str(e))
                if (lmb is not None) and (not isinstance(lmb, (tuple, list)) or any(not isinstance(el,str) for el in lmb)):
                    raise alphaLoopInvalidCmd("Each element of the list of lmb must be a string identifying the edge name.")                        

                processed_args[key[2:]] = lmb

            if key == '--analytical_result':
                try:
                    analytical_result = complex(value)
                except Exception as e:
                    traceback.print_exc()
                    raise alphaLoopInvalidCmd(
                        "Specified analytical result cannot be parsed to a complex: %s"%value+
                        "\n The following Python intepreter error occured then: %s"%str(e))
                processed_args[key[2:]] = analytical_result

        for arg in mandatory_args:
            if processed_args[arg] is None:
                raise alphaLoopInvalidCmd("Missing mandatory option '%s' for command do_output_LU_scalar."%arg)


        lu_scalar_exporter = aL_exporters.LUScalarTopologyExporter(
            self,
            output_path, processed_args['topology'], 
            processed_args['externals'], processed_args['name'], 
            processed_args['lmb'], self._curr_model, 
            benchmark_result=processed_args['analytical_result'],
            alphaLoop_options=self.alphaLoop_options,
            MG5aMC_options=self.options
        )

        lu_scalar_exporter.output()

    def do_output_scalar_integral(self, line):
        """ Generates and directly output a scalar integal topology as a mock-up of a LU integrand (with or without numerator).
        """

        args = self.split_arg(line)

        if len(args)<=1:
            raise alphaLoopInvalidCmd("Missing mandatory first argument for output LU scalar which is the output directory name.")
        
        output_path = args.pop(0)

        processed_args = {
            'topology': None,
            'name': 'DefaultScalarIntegralName',
            'externals': None,
            'default_kinematics': None,
            'lmb' : None,
            'masses': None,
            'numerator': None,
            'analytical_result': 0.0, 
        }
        mandatory_args = ['topology', 'externals']

        # Group arguments in between the '--' specifiers.
        # For example, it will group '--process=p' 'p' '>' 'd' 'd~' 'z'
        # into '--process=p p > d d~ z'.
        opt_args = []
        new_args = []
        for arg in args:
            if arg.startswith('--'):
                opt_args.append(arg)
            elif len(opt_args)>0:
                opt_args[-1] += ' %s'%arg
            else:
                new_args.append(arg)
        
        for arg in opt_args:
            try:
                key, value = arg.split('=',1)
                value = value.strip()
            except ValueError:
                key = arg
                value = None
        
            if key == '--topology':
                try:
                    topology=eval(value)
                except Exception as e:
                    traceback.print_exc()
                    raise alphaLoopInvalidCmd(
                        "Could not process specified topology: %s"%value+
                        "\n The following Python intepreter error occured then: %s"%str(e))
                if not isinstance(topology, (tuple, list)):
                    raise alphaLoopInvalidCmd("Specified topology must be a tuple or a list")
                if any(
                    ( (not isinstance(te,(tuple,list))) or (not isinstance(te[0],str)) or (not isinstance(te[1],int)) or (not isinstance(te[2],int)) )
                    for te in topology ):
                    raise alphaLoopInvalidCmd("Each element of the topology specified must be a tuple ('edge_name', node_int_id_left, node_int_id_right) ")                        

                processed_args['topology'] = topology

            if key == '--name':
                if "'" in value or '"' in value:
                    processed_args[key[2:]] = eval(value)
                else:
                    processed_args[key[2:]] = value
            
            if key == '--numerator':
                if "'" in value or '"' in value:
                    processed_args[key[2:]] = eval(value)
                else:
                    processed_args[key[2:]] = value

            if key == '--externals':
                try:
                    externals = eval(value)
                except Exception as e:
                    traceback.print_exc()
                    raise alphaLoopInvalidCmd(
                        "Could not process specified externals for this topology: %s"%value+
                        "\n The following Python intepreter error occured then: %s"%str(e))
                if (not isinstance(externals, (tuple, list))) or (not len(externals)==2) or \
                    any(not isinstance(el,str) for el in externals[0]) or \
                    any(not isinstance(el,str) for el in externals[1]):
                    raise alphaLoopInvalidCmd("Each element of each of the two lists provided in the 2-tuple 'externals' must be a string identifying the edge name.")                        

                processed_args[key[2:]] = externals

            if key == '--default_kinematics':
                try:
                    default_kinematics = eval(value)
                except Exception as e:
                    traceback.print_exc()
                    raise alphaLoopInvalidCmd(
                        "Could not process specified default kinematics for this scalar integral: %s"%value+
                        "\n The following Python intepreter error occured then: %s"%str(e))

                if not isinstance(default_kinematics, dict) or any(not isinstance(k,str) for k in default_kinematics.keys()) or \
                    any( (not isinstance(v, (tuple, list)) or any( len(v)!=4 or not isinstance(vi,float) for vi in v) ) for v in default_kinematics.values()):
                    raise alphaLoopInvalidCmd("Keys of the default_kinematics provided must identify the edge name and values must be four-tuples specifying the kinematics.")                        

                processed_args[key[2:]] = default_kinematics

            if key == '--masses':
                try:
                    masses = eval(value)
                except Exception as e:
                    traceback.print_exc()
                    raise alphaLoopInvalidCmd(
                        "Could not process specified default kinematics for this scalar integral: %s"%value+
                        "\n The following Python intepreter error occured then: %s"%str(e))
                if not isinstance(masses, dict) or any(not isinstance(k,str) for k in masses.keys()) or \
                    any( not isinstance(mi,float) for mi in masses.values()):
                    raise alphaLoopInvalidCmd("Keys of the masses dictionary provided must identify the edge name and values must be a float specifying the mass.")                        

                processed_args[key[2:]] = masses

            if key == '--lmb':
                try:
                    lmb = eval(value)
                except Exception as e:
                    traceback.print_exc()
                    raise alphaLoopInvalidCmd(
                        "Could not process specified lmb for this topology: %s"%value+
                        "\n The following Python intepreter error occured then: %s"%str(e))
                if (lmb is not None) and (not isinstance(lmb, (tuple, list)) or any(not isinstance(el,str) for el in lmb)):
                    raise alphaLoopInvalidCmd("Each element of the list of lmb must be a string identifying the edge name.")                        

                processed_args[key[2:]] = lmb

            if key == '--analytical_result':
                try:
                    analytical_result = complex(value)
                except Exception as e:
                    traceback.print_exc()
                    raise alphaLoopInvalidCmd(
                        "Specified analytical result cannot be parsed to a complex: %s"%value+
                        "\n The following Python intepreter error occured then: %s"%str(e))
                processed_args[key[2:]] = analytical_result

        for arg in mandatory_args:
            if processed_args[arg] is None:
                raise alphaLoopInvalidCmd("Missing mandatory option '%s' for command do_output_LU_scalar."%arg)


        scalar_integral_exporter = aL_exporters.ScalarIntegralTopologyExporter(
            self,
            output_path, processed_args['topology'], 
            processed_args['externals'], 
            processed_args['default_kinematics'],
            processed_args['masses'],
            processed_args['name'], 
            processed_args['lmb'], self._curr_model, 
            numerator=processed_args['numerator'],
            benchmark_result=processed_args['analytical_result'],
            alphaLoop_options=self.alphaLoop_options,
            MG5aMC_options=self.options
        )

        scalar_integral_exporter.output()

    ######################################################################
    #
    # Example of the implementation of a trivial function 'hello_world'
    #
    ######################################################################

    def do_hello_world(self, line):
        """ Hello world command example."""

        logger.info('Hello World, and welcome you, %s%s%s!'%(utils.bcolors.GREEN, line,utils.bcolors.ENDC))

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
        self.prompt = "u'\u03B1Loop > "

        logger.info("\n\n%s\n"%self.get_alpha_loop_banner())
        logger.info("Loading default model for alphaLoop: sm")

        # By default, load the UFO Standard Model
        logger.info("Loading default model for alphaLoop: sm")
        self.exec_cmd('import model sm', printcmd=False, precmd=True)

        # preloop mother
        madgraph_interface.CmdExtended.preloop(self)
