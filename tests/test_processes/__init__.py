import unittest
import os
import sys
import warnings
warnings._showwarnmsg = (lambda msg: 0)

import inspect
import shutil
from pprint import pformat
from copy import deepcopy

import logging
FORMAT = '%(name)s : %(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)

# Run the tests with `AL_TESTS_UPDATE=True python3.7 -m unittest -v` for having the tests return what are
# the new target values
UPDATE_TESTS = ('AL_TESTS_UPDATE' in os.environ)
    
DEBUG_MODE = int(os.environ.get('AL_DEBUG',0))
master_logger = logging.getLogger("aLTest")
master_logger.setLevel(logging.INFO)
# ===========================================================================
# Run the tests with `AL_DEBUG=3 python3.7 -m unittest -v` for debug mode
if DEBUG_MODE >= 1:
    master_logger.setLevel(logging.DEBUG)
# ===========================================================================

if UPDATE_TESTS:
    master_logger.info("AL_TESTS_UPDATE detected: will return new values for updating tests")

# Make sure madgraph remains silent
logging.getLogger('madevent').setLevel(logging.CRITICAL if DEBUG_MODE < 3 else logging.INFO)
logging.getLogger('madgraph').setLevel(logging.CRITICAL if DEBUG_MODE < 3 else logging.INFO)
logging.getLogger('cmdprint').setLevel(logging.CRITICAL if DEBUG_MODE < 3 else logging.INFO)

# Disable harmless warnings due to MadGraph
def ignore_warnings(test_func):
    def do_test(self, *args, **kwargs):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ResourceWarning)
            test_func(self, *args, **kwargs)
    return do_test

pjoin = os.path.join
_root_path = os.path.dirname(os.path.realpath( __file__ ))
_aL_dir = os.path.abspath(pjoin(_root_path,os.path.pardir,os.path.pardir))
_mg5_dir = os.path.abspath(pjoin(_aL_dir,os.path.pardir,os.path.pardir))

# If not in debug mode, disable progress bar too
if DEBUG_MODE < 3:
    from unittest.mock import MagicMock
    myprogressbar = MagicMock()
    def progressbar_mock(iterator, **kwargs):
        for result in iterator:
            yield result
    myprogressbar.progressbar = progressbar_mock
    sys.modules['progressbar'] = myprogressbar

# Important to make madgraph path have precedence for the models module loaded to be the madgaph one
if sys.path[:2] !=  [_mg5_dir,_aL_dir]:
    sys.path.insert(0, _aL_dir )
    sys.path.insert(0, _mg5_dir )

import alpha_loop
from alpha_loop.interface import alphaLoopInterface
from alpha_loop.run_interface import alphaLoopRunInterface, CrossSectionSet
from alpha_loop.ltd_commons import HyperParameters
from alpha_loop.ltd_commons import hyperparameters as default_hyperparameters
import alpha_loop.utils as utils

# Load Python API
from ltd import CrossSection

class ProcessTester(object):

    # To be specitied by daughter class
    _generation_card = None
    _global_hyperparameters = None

    @classmethod
    @ignore_warnings
    def setUpClass(cls):

        cls._logger = logging.getLogger('aLTest.%s'%cls.__name__)

        cls._logger.debug("Starting setUp")

        cls._workspace_dir = '%s_workspace'%cls.__name__

        recycled = False
        if os.path.isdir(pjoin(_root_path,cls._workspace_dir)):
            if DEBUG_MODE>=2:
                cls._logger.debug("Recycling existing process output at '%s'"%pjoin(_root_path,cls._workspace_dir))
                recycled = True
            else:
                cls._logger.debug("Cleaning up workspace '%s'"%pjoin(_root_path,cls._workspace_dir))
                shutil.rmtree(pjoin(_root_path,cls._workspace_dir))

        cls.proc_output_path = pjoin(_root_path,cls._workspace_dir,'LU_%s'%cls.__name__)
        
        # Now generate the process
        generation_card = cls._generation_card%({'output_path':cls.proc_output_path})

        if not recycled:
            cls._logger.debug("Generating process")
            with utils.suppress_output(active=(DEBUG_MODE>=3)):
                cls.cmd_line = alphaLoopInterface()
                for line in generation_card.split('\n'):
                    if line.strip()=='':
                        continue
                    cls.cmd_line.run_cmd(line.strip())

        cls.cross_section_set = CrossSectionSet(pjoin(cls.proc_output_path,'Rust_inputs','all_QG_supergraphs.yaml'))

        cls.hyperparams = HyperParameters()
        cls.hyperparams.update(default_hyperparameters, allow_undefined=True)
        for param, value in cls._global_hyperparameters.items():
            cls.hyperparams.set_parameter(param, value)

        cls._logger.debug("setUp finished")

    @classmethod
    def tearDownClass(cls):
        cls._logger.debug("Starting tearDown")
        if DEBUG_MODE<2:
            try:
                shutil.rmtree(pjoin(_root_path,cls._workspace_dir))
            except Exception as e:
                cls._logger.info("Failed to clean up directory '%s'. Exception: %s."%(cls.proc_output_path, str(e)))
        cls._logger.debug("tearDown finished")

    def run_local_evaluations(self, test_name, extra_hyperparameters, default_xs, targets ):
        
        hyperparams = deepcopy(self.hyperparams)
        for param, value in extra_hyperparameters.items():
            hyperparams.set_parameter(param, value)
        hyperparams.export_to(pjoin(self.proc_output_path,'hyperparameters_%s.yaml'%test_name))

        all_sg_names = sorted([sg['name'] for sg in self.cross_section_set['topologies']])
        if not UPDATE_TESTS:
            self.assertListEqual( sorted(list(targets.keys())), all_sg_names, "List of contributing supergraphs for test %s differ."%self.__class__.__name__)
        else:
            targets = { sg['name'] : [ targets.get(sg['name'],[None,])[0], sg['multiplicity'], None ] for sg in self.cross_section_set['topologies'] }

        multiplicities = {sg['name'] : sg['multiplicity'] for sg in self.cross_section_set['topologies']}
        for supergraph_name in all_sg_names:
            
            xs, multiplicity, result = targets[supergraph_name]
            if xs is None:
                xs = default_xs

            if not UPDATE_TESTS:
                self.assertAlmostEqual(multiplicity, multiplicities[supergraph_name], 15, "Multiplicity for %s differs"%supergraph_name)

            re, im = None, None
            with utils.suppress_output(active=(DEBUG_MODE>=3)):
                rust_worker = CrossSection(
                    pjoin(self.proc_output_path,'Rust_inputs',supergraph_name+'.yaml'), 
                    pjoin(self.proc_output_path,'hyperparameters_%s.yaml'%test_name)
                )
                re, im = rust_worker.evaluate_integrand(xs)

            if not UPDATE_TESTS:
                with self.subTest(SG=supergraph_name):
                    self.assertAlmostEqual(result.real, re, 13, "Real part for %s differs"%supergraph_name)
                    self.assertAlmostEqual(result.imag, im, 13, "Imaginary part for %s differs"%supergraph_name)
            else:
                targets[supergraph_name][2] = complex(re,im)

        if UPDATE_TESTS:
            self._logger.info("\nNew targets for test '%s%s%s' of class '%s%s%s':\n%s%s%s"%(
                utils.bcolors.BLUE, test_name, utils.bcolors.ENDC,
                utils.bcolors.BLUE, self.__class__.__name__, utils.bcolors.ENDC,
                utils.bcolors.GREEN, pformat(targets), utils.bcolors.ENDC
            ))