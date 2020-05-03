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
import alpha_loop.utils as utils
import madgraph.various.cluster as cluster

from madgraph.iolibs.files import cp, ln, mv

logger = logging.getLogger('alphaLoop.Interface')

pjoin = os.path.join
template_dir = pjoin(plugin_path, 'Templates')

class alphaLoopInterfaceError(MadGraph5Error):
    """ Error for the pyNLoop plugin """
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
        }

        super(alphaLoopInterface, self).__init__(*args, **opts)

    def parse_set_alphaLoop_option(self, args):
        """ Parsing arguments/options passed to the command set_pyNLoop option."""

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
        for opt in sorted(self.pyNLoop_options.keys()):
            logger.info('%-30s : %s'%(opt, str(self.pyNLoop_options[opt])))

    def do_set_alphaLoop_option(self, line):
        """ Logic for setting alphaLoop options."""
        args = self.split_arg(line)
        args, options = self.parse_set_pyNLoop_option(args)
        key, value = args[:2]

        if key == 'TO_IMPLEMENT':
            self.pyNLoop_options[key] = eval(value)
        else:
            raise alphaLoopInvalidCmd("Unrecognized alphaLoop option: %s"%key)

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
        logger.info("Loading default model for pyNLoop: sm")

        # By default, load the UFO Standard Model
        logger.info("Loading default model for alphaLoop: sm")
        self.exec_cmd('import model sm', printcmd=False, precmd=True)

        # preloop mother
        madgraph_interface.CmdExtended.preloop(self)
