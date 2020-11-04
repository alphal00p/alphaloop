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
from argparse import ArgumentParser
from pprint import pprint, pformat

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
        res = {}
        res['Number of supergraphs']='%d'%len(self)
        sorted_mult_SG = sorted(( (SG['multiplicity'],SG['name']) for SG in self['topologies']) ,key=lambda el:el[0])
        res['Min multiplicity']='%d (%s)'%(sorted_mult_SG[0][0],sorted_mult_SG[0][1])
        res['Max multiplicity']='%d (%s)'%(sorted_mult_SG[-1][0],sorted_mult_SG[-1][1])

        res_str = []
        for k, v in res.items():
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

    def check_validity(self):
        return True

    def __str__(self):
        return pformat(self)

    def summary_str(self):
        res = {}
        res['Number of cuts']='%d'%len(self['cutkosky_cuts'])

        res_str = []
        for k, v in res.items():
            res_str.append('%s%-40s%s : %s'%(Colours.GREEN,k,Colours.END,v))
        return '\n'.join(res_str)

class SuperGraphCollection(dict):
    
    def __init__(self, *args, **opts):
        super(SuperGraphCollection, self).__init__(*args, **opts)

class alphaLoopRunInterface(madgraph_interface.MadGraphCmd, cmd.CmdShell):
    """ Interface for steering the running of an alphaLoop output.
    We make it inherit from CmdShell so that launch_ext_prog does not attempt to start in WebMode."""

    _rust_inputs_folder = 'Rust_inputs'
    _cross_section_set_yaml_name = 'all_QG_supergraphs.yaml'

    def __init__(self, dir_path, alphaLoop_interface, launch_options={}, *args, **opts):
        """ Define attributes of this class."""

        self.dir_path = dir_path
        self.alphaLoop_interface = alphaLoop_interface
        self.launch_options = launch_options

        self.cross_section_set = CrossSectionSet(pjoin(self.dir_path,
                            self._rust_inputs_folder, self._cross_section_set_yaml_name))
        self.all_supergraphs = self.load_supergraphs()

        super(alphaLoopRunInterface, self).__init__(*args, **opts)

    def load_supergraphs(self):
        
        SG_collection = SuperGraphCollection()

        for SG in self.cross_section_set['topologies']:
            yaml_path = pjoin(self.dir_path, self._rust_inputs_folder, SG['name']+'.yaml')
            if not os.path.isfile(yaml_path):
                raise alphaLoopInvalidRunCmd("Could not find yaml file at '%s' specifying the supergraph information."%yaml_path)
            SG_collection[SG['name']] = SuperGraph(yaml_path)
        
        return SG_collection

    display_parser = ArgumentParser()
    display_parser.add_argument('SG_name', metavar='SG_name', type=str, nargs='?',
                    help='the name of a supergraph to display')
    display_parser.add_argument(
        "--full", action="store_true", dest="full", default=False,
        help="exhaustively show information")
    def do_display(self, line):

        args = self.split_arg(line)
        args = self.display_parser.parse_args(args)

        if args.SG_name is None:
            if args.full:
                logger.info("Exhaustive content of the cross-section set:\n%s"%str(self.cross_section_set))
            else:
                logger.info("Summary of the cross-section set:\n%s"%self.cross_section_set.summary_str())
        elif args.SG_name == 'ALL':
            for SG_name, SG in self.all_supergraphs.items():
                logger.info("Summary of the supergraph '%s%s%s':\n%s"%(Colours.GREEN, SG_name, Colours.END, SG.summary_str()))
        else:
            if args.SG_name not in self.all_supergraphs:
                raise alphaLoopInvalidRunCmd("Supergraph named '%s' not found in the supergraphs loaded."%args.SG_name)
            if args.full:
                logger.info("Exhaustive content of the supergraph '%s%s%s':\n%s"%(Colours.GREEN, args.SG_name, Colours.END, str(self.all_supergraphs[args.SG_name])))
            else:
                logger.info("Summary of the supergraph '%s%s%s':\n%s"%(Colours.GREEN, args.SG_name, Colours.END, self.all_supergraphs[args.SG_name].summary_str()))

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
        self.prompt = "u'\u03B1Loop @ %s > "%os.path.basename(self.dir_path)

        # preloop mother
        madgraph_interface.CmdExtended.preloop(self)