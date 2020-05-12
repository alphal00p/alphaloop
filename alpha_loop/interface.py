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

import alpha_loop.utils as utils
import alpha_loop.exporters as aL_exporters
import alpha_loop.helas_call_writers as aL_helas_call_writers
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
            'perturbative_order' : 'LO',
            # By default we don't want denominators as those will be included in the rust backend
            'include_denominators' : False, 
            'use_physical_gluon_helicity_sum' : False,
            # Specify the number of rust inputs to generate, -1 considers all by default.
            # This is useful for debugging as this can be slow.
            'n_rust_inputs_to_generate' : -1,
        }
        self.plugin_output_format_selected = None

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

        if key == 'perturbative_order':
            if not re.match(r'^(N*)LO$',value):
                raise alphaLoopInvalidCmd("Specified perturbative order should be of the form N*LO, not '%s'."%value)
            self.alphaLoop_options[key] = value
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
        elif key == 'n_rust_inputs_to_generate':
            try:
                self.alphaLoop_options['n_rust_inputs_to_generate'] = int(value)
            except ValueError:
                raise alphaLoopInvalidCmd("Specified value for 'n_rust_inputs_to_generate' should be an integer, not '%s'."%value)
        else:
            raise alphaLoopInvalidCmd("Unrecognized alphaLoop option: %s"%key)

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
