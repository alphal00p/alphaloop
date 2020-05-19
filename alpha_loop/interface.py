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
import alpha_loop.exporters as aL_exporters
import alpha_loop.helas_call_writers as aL_helas_call_writers
import alpha_loop.LTD_squared as LTD_squared
import alpha_loop.madgraph_patches as madgraph_patches

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

    def __init__(self, *args, **opts):
        """ Define attributes of this class."""
        
        self.alphaLoop_options = {
            'perturbative_orders' : { 'QCD': 0 },
            # By default we don't want denominators as those will be included in the rust backend
            'include_denominators' : False, 
            'use_physical_gluon_helicity_sum' : False,
            # Specify the number of rust inputs to generate, -1 considers all by default.
            # This is useful for debugging as this can be slow.
            'n_rust_inputs_to_generate' : -1,
            # We plan on include self-energies from two-point vertices so as to be able to apply our 
            # LTD^2 LSZ treatment. So 'include_self_energies_from_squared_amplitudes' should be kept
            # set to False, except for debugging.
            'include_self_energies_from_squared_amplitudes' : False,
            # One must differentiate particles from anti-particles in the isomorphism check.
            # However this is not properly working, at least not for self-energies, so we allow it
            # to be disabled with the option below.
            'differentiate_particle_from_antiparticle_in_graph_isomorphism' : False
        }
        self.plugin_output_format_selected = None

        self.model_backup_copy = None

        super(alphaLoopInterface, self).__init__(*args, **opts)

    def parse_set_alphaLoop_option(self, args):
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

    def do_set_alphaLoop_option(self, line):
        """ Logic for setting alphaLoop options."""
        args = self.split_arg(line)
        args, options = self.parse_set_alphaLoop_option(args)
        key, value = args[:2]

        if key == 'perturbative_orders':
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
        elif key == 'include_self_energies_from_squared_amplitudes':
            if value.upper() not in ['TRUE','FALSE']:
                raise alphaLoopInvalidCmd("Specified value for 'include_self_energies_from_squared_amplitudes' should be 'True' or 'False', not '%s'."%value)
            bool_val = (value.upper()=='TRUE')
            if bool_val:
                logger.warning(
"""Self-energies will be generated from two-point vertices so as to be able to apply the 
LTD^2 LSZ treatment. So 'include_self_energies_from_squared_amplitudes' should be kept
set to False, except for debugging, which seems to be what you are doing now, so we'll se it to True now.
""")
            self.alphaLoop_options['include_self_energies_from_squared_amplitudes'] = bool_val   
        elif key == 'differentiate_particle_from_antiparticle_in_graph_isomorphism':
            if value.upper() not in ['TRUE','FALSE']:
                raise alphaLoopInvalidCmd("Specified value for 'differentiate_particle_from_antiparticle_in_graph_isomorphism' should be 'True' or 'False', not '%s'."%value)
            bool_val = (value.upper()=='TRUE')
            self.alphaLoop_options['differentiate_particle_from_antiparticle_in_graph_isomorphism'] = bool_val
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

    def OLD_prepare_model_for_self_energies(self, revert=False, perturbative_order=3):
        """Commands for adding the necessary TREE UV interactions to the current model."""
        
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

        pure_QCD_corrections_only = list(self.alphaLoop_options['perturbative_orders'].keys())==['QCD',]

        # First add the necessary Self-energt coupling_orders to the available coupling_orders 
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
        pdg_offset = 100000000
        for n_self_energy in range(1,perturbative_order+1):
            new_particles[n_self_energy] = {}

            # Let us not create self-energy specific particles just yet
            # TODO
#            for particle in model['particles']:
#                if pure_QCD_corrections_only:
#                    if particle.get('color')==1:
#                        continue
#                new_particle=copy.deepcopy(particle)
#                new_particle['pdg_code'] = particle['pdg_code'] + n_self_energy*pdg_offset
#                new_particle['name'] = ('SE_%d_'%n_self_energy+particle['name']).lower()
#                new_particle['antiname'] = ('SE_%d_'%n_self_energy+particle['antiname']).lower()
#                new_particles[n_self_energy][particle['pdg_code']]=new_particle

            for n_anchor in range(0,perturbative_order+1):
                anchor_offset = int(pdg_offset/10)
                anchor = base_objects.Particle({
                    'pdg_code' : n_self_energy*pdg_offset+n_anchor*anchor_offset,
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
                bridge_offset = int(pdg_offset/100)
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
            model.get('particles')+[p for se_parts in new_particles.values() for p in se_parts.values()] ))


        # Then create the new interactions
        interaction_offset = 1000000
        new_interactions = []
        for n_self_energy in range(1,perturbative_order+1):
            # Dubplicate base interactions so as to add the the anchor scalar #1 to it.
            id_offset =0
            for inter in model['interactions'].get_type('base'):
                if (set(inter['orders'].keys())-set(self.alphaLoop_options['perturbative_orders']))!=set([]):
                    continue
                id_offset += 1
                new_interaction = copy.deepcopy(inter)
                new_interaction.set('color',inter.get('color'))
                new_interaction.set('id',n_self_energy*interaction_offset+id_offset)
                new_orders = {}
                # Let's not use the "self-energy particles" just yet, but simply normal particles and therefore normal orders:
                # TODO
                for order in new_interaction['orders']:
                    new_orders[order]=new_interaction['orders'][order]
#                    new_orders['SE_%d_%s'%(n_self_energy,order)]=new_interaction['orders'][order]
                new_orders['SE_%d_ANCHOR'%n_self_energy] = 1
                new_interaction.set('orders',new_orders)
                # Let's not use the "self-energy particles" just yet, but simply normal particles:
                # TODO
                new_interaction['particles'] = base_objects.ParticleList([
#                       new_particles[n_self_energy][p['pdg_code']] for p in new_interaction['particles']
                        p for p in new_interaction['particles']
                    ]+[new_particles[n_self_energy]['anchor_0'],]
                )
                new_interaction['lorentz'] = ['SE_%s'%l for l in new_interaction['lorentz']]
                new_interactions.append(new_interaction)

                # TODO when using the "self-energy particles" above, then also copy down interactions
                # for the whole perturbed sector aming the self-energy particles for a given self-energy index.

            # Add an interaction between three scalars: Anchor_<i>_{n-1} | Bridge_<i>_{n} | Anchor_<i>_{n-1}
            for n_bridge in range(1,perturbative_order+1):
                bridge_offset = int(interaction_offset/10)
                new_interactions.append(base_objects.Interaction({
                    'id' : n_self_energy*interaction_offset+n_bridge*bridge_offset,
                    'particles' : base_objects.ParticleList([
                        new_particles[n_self_energy]['anchor_%d'%(n_bridge-1)],
                        new_particles[n_self_energy]['bridge_%d'%n_bridge],
                        new_particles[n_self_energy]['anchor_%d'%n_bridge],
                    ]),
                    'color': [color.ColorString(),],
                    'lorentz': [ 'SE_SSS',],
                    'couplings' : { (0, 0): '1'},
                    'orders' : {'SE_%d_BRIDGE_%d'%(n_self_energy,n_bridge): 1}
                }))
                id_offset = 0
                # TODO Make sure here to also roll over the interactions between "self-energy" particles.
                for inter in model['interactions'].get_type('base'):
                    id_offset += 1
                    new_interaction = copy.deepcopy(inter)
                    new_interaction.set('id',n_self_energy*interaction_offset+n_bridge*bridge_offset+id_offset)
                    new_orders = new_interaction['orders']
                    new_orders['SE_%d_CONTACT_%d'%(n_self_energy,n_bridge)] = 1
                    new_interaction.set('orders',new_orders)
                    new_interaction['particles'] = base_objects.ParticleList([
                            p for p in new_interaction['particles']
                        ]+[new_particles[n_self_energy]['bridge_%d'%n_bridge],]
                    )
                    new_interaction['lorentz'] = ['SE_%s'%l for l in new_interaction['lorentz']]
                    new_interactions.append(new_interaction)                    

        model.set('interactions',base_objects.InteractionList(model.get('interactions')+new_interactions))
        # Reset the model dictionaries to sync changes
        model.reset_dictionaries()
        # Refresh dictionaries manually (would be done automatically anyway)
        model.actualize_dictionaries()
    
        model.set('order_hierarchy', new_order_hierarchy)
        model.set('perturbation_couplings', new_perturbation_couplings)


        logger.info('Model %s successfully patched so as to accomodate self-energy generations.'%model.get('name'))
        return

    def do_output(self, line):
        """ Wrapper to support the syntax output alphaLoop <args>.
        This just to add extra freedom in adding special action that may be needed at the output
        stage for these output formats.
        """
        args = self.split_arg(line)
        if len(args)>=1 and args[0]=='alphaLoop':
            self.plugin_output_format_selected = 'alphaLoop'
            self.do_output_alphaLoop(' '.join(args[1:]))    
        else:
            super(alphaLoopInterface,self).do_output(' '.join(args))

    def do_output_alphaLoop(self, line):
        args = self.split_arg(line)
        super(alphaLoopInterface,self).do_output(' '.join(['alphaLoop']+args))

    def export(self,*args,**opts):
        """Overwrite this so as to force a pythia8 type of output if the output mode is PY8MEs."""
        
        if self._export_format == 'plugin':
            # Also pass on the aloha model to the exporter (if it has been computed already)
            # so that it will be used when generating the model
            if self.plugin_output_format_selected == 'alphaLoop':
                # Set what are the jet pdgs in alphaLoop options
                if 'j' not in self._multiparticles:
                    raise alphaLoopInvalidCmd("alphaLoop requires the 'j' multiparticle label to be defined.")
                self.alphaLoop_options['_jet_PDGs'] = tuple([
                    pdg for pdg in self._multiparticles['j']
                ])
                self._curr_exporter = aL_exporters.alphaLoopExporter(self._export_dir,
                    alphaLoop_options=self.alphaLoop_options,
                    MG5aMC_options=self.options
                )
            else:
                raise alphaLoopInterfaceError("A plugin output format must have been specified at this stage.")

        super(alphaLoopInterface,self).export(*args, **opts)

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
        self.prompt = 'alphaLoop > '

        logger.info("\n\n%s\n"%self.get_alpha_loop_banner())
        logger.info("Loading default model for alphaLoop: sm")

        # By default, load the UFO Standard Model
        logger.info("Loading default model for alphaLoop: sm")
        self.exec_cmd('import model sm', printcmd=False, precmd=True)

        # preloop mother
        madgraph_interface.CmdExtended.preloop(self)
