import unittest
import os
import sys
import warnings
warnings._showwarnmsg = (lambda msg: 0)

import inspect
import shutil
import math
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
from alpha_loop.run_interface import alphaLoopRunInterface, CrossSectionSet, SuperGraph
from alpha_loop.ltd_commons import HyperParameters
from alpha_loop.ltd_commons import hyperparameters as default_hyperparameters
import alpha_loop.utils as utils

Colours = utils.bcolors
#Colours = utils.NO_colors

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

        cls.cmd_line = alphaLoopInterface()
        cls._logger.debug("Generating process")
        with utils.suppress_output(active=(DEBUG_MODE>=3)):
            for line in generation_card.split('\n'):
                if line.strip()=='':
                    continue
                if recycled and 'output qgraf' in line:
                    continue
                cls.cmd_line.run_cmd(line.strip())

        cls.run_cmd_line = None

        cls.cross_section_set = CrossSectionSet(pjoin(cls.proc_output_path,'Rust_inputs','all_QG_supergraphs.yaml'))

        cls.hyperparams = HyperParameters()
        cls.hyperparams.update(default_hyperparameters, allow_undefined=True)
        for param, value in cls._global_hyperparameters.items():
            cls.hyperparams.set_parameter(param, value)
        cls.hyperparams.export_to(pjoin(cls.proc_output_path,'default_test_hyperparameters.yaml'))

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

    def run_UV_tests(self, warmup_cmds, uv_test_cmd):

        if self.run_cmd_line is None:
            self.run_cmd_line = alphaLoopRunInterface(
                self.proc_output_path, self.cmd_line, 
                launch_options={'reuse': pjoin(self.proc_output_path,'default_test_hyperparameters.yaml')}
            )
        
        res = None
        with utils.suppress_output(active=(DEBUG_MODE>=3)):
            # Always remove all previously generated test results
            self.run_cmd_line.run_cmd('refresh_derived_data')
            # Send warmup commands
            for line in warmup_cmds.split('\n'):
                if line.strip()=='':
                    continue
                self.run_cmd_line.run_cmd(line.strip())
            # Perform test
            res = self.run_cmd_line.run_cmd(uv_test_cmd.strip())
        

        for supergraph_name in sorted(list(res.keys())):
            sg_info = res[supergraph_name]

            if ('DERIVED_UV_dod' not in sg_info):
                continue
            for k, v in sg_info['DERIVED_UV_dod'].items():
                with self.subTest(SG=supergraph_name):
                    self.assertTrue(v[-1], ('%-s : %s')%('%s%s%s : Sum over cuts dod for LMB=%s, UV_indices=%s'%(
                        Colours.RED, supergraph_name, Colours.END,
                        ','.join(k[0]),str(k[1])
                    ),'%-7.4f +/- %-7.4f %s'%(
                        v[0],v[1], '%sPASS%s'%(Colours.GREEN, Colours.END) if v[-1] else '%sFAIL%s'%(Colours.RED, Colours.END)
                    )))

            sorted_cut_evaluations = sorted(
                [(i,k,v) for i, c in enumerate(sg_info['cutkosky_cuts']) for k,v in c['DERIVED_UV_dod'].items()],
                key = lambda el: el[2][0], reverse=True
            )
            for cut_ID, k, v in sorted_cut_evaluations:
                with self.subTest(SG=supergraph_name,cut_ID=cut_ID):
                    self.assertTrue(v[-1], ('%-s : %s')%('%s%s%s : Cut dod for cut_ID=%d, LMB=%s, UV_indices=%s'%(
                        Colours.RED, supergraph_name, Colours.END,
                        cut_ID, ','.join(k[0]),str(k[1])
                    ),'%-7.4f +/- %-7.4f %s'%(
                        v[0],v[1], '%sPASS%s'%(Colours.GREEN, Colours.END) if v[-1] else '%sFAIL%s'%(Colours.RED, Colours.END)
                    )))

    def run_IR_tests(self, warmup_cmds, ir_test_cmd):

        if self.run_cmd_line is None:
            self.run_cmd_line = alphaLoopRunInterface(
                self.proc_output_path, self.cmd_line, 
                launch_options={'reuse': pjoin(self.proc_output_path,'default_test_hyperparameters.yaml')}
            )
        
        res = None
        with utils.suppress_output(active=(DEBUG_MODE>=3)):
            # Always remove all previously generated test results
            self.run_cmd_line.run_cmd('refresh_derived_data')
            # Send warmup commands
            for line in warmup_cmds.split('\n'):
                if line.strip()=='':
                    continue
                self.run_cmd_line.run_cmd(line.strip())
            # Perform test
            res = self.run_cmd_line.run_cmd(ir_test_cmd.strip())


        for supergraph_name in sorted(list(res.keys())):
            sg_info = res[supergraph_name]

            if 'ir_limits_analysis' not in sg_info or len(sg_info['ir_limits_analysis'])==0:
                continue

            ir_limits_to_consider = sg_info['ir_limits_analysis']

            IR_limits_per_order = {}
            for ir_limit in ir_limits_to_consider:
                pert_order = SuperGraph.compute_ir_limit_perturbative_order(ir_limit)
                if pert_order in IR_limits_per_order:
                    IR_limits_per_order[pert_order].append(ir_limit)
                else:
                    IR_limits_per_order[pert_order] = [ir_limit,]

            for pert_order in sorted(list(IR_limits_per_order.keys())):
                max_IR_limit_str_len = max(len(SuperGraph.format_ir_limit_str(ir_limit,colored_output=False)) for ir_limit in ir_limits_to_consider )
                for ir_limit in sorted(IR_limits_per_order[pert_order]):
                    result = ir_limits_to_consider[ir_limit]
                    cuts_sum_dod_colour = Colours.GREEN if (result['cuts_sum']['dod']['central'] < -1 + max(min(max(10.0*abs(result['cuts_sum']['dod']['std_err']),0.05),0.2),sg_info['ir_limits_analysis_setup']['min_dod_tolerance'])) else Colours.RED
                    max_cut_dod = max(result['per_cut'].values(), key=lambda d: d['dod']['central'])
                    severity_central = (max_cut_dod['dod']['central']-result['complete_integrand']['dod']['central'])
                    severity_std_err = math.sqrt(result['complete_integrand']['dod']['std_err']**2+max_cut_dod['dod']['std_err']**2)
                    with self.subTest(SG=supergraph_name,ir_limit=SuperGraph.format_ir_limit_str(ir_limit,colored_output=False)):
                        self.assertTrue(result['status'][0],
                            ('> %-{}s : %s%s%s %-12s dod: itg = %-32s cuts_sum = %-32s severity = %-32s'.format(
                                    max_IR_limit_str_len+(len(SuperGraph.format_ir_limit_str(ir_limit,colored_output=True))-len(SuperGraph.format_ir_limit_str(ir_limit,colored_output=False)))
                                ))%(
                                SuperGraph.format_ir_limit_str(ir_limit),
                                Colours.GREEN if result['status'][0] else Colours.RED, 'PASSED' if result['status'][0] else 'FAILED', Colours.END,
                                '(%s)'%result['status'][1],
                                '%sN/A%s'%(Colours.BLUE, Colours.END) if result['complete_integrand']['dod']['status'] not in ['SUCESS','UNSTABLE'] else (
                                    '%s%s%.5f +/- %.1g%s'%(
                                        Colours.GREEN if result['status'][0] else Colours.RED,
                                        '+' if result['complete_integrand']['dod']['central']>=0. else '',
                                        result['complete_integrand']['dod']['central'], result['complete_integrand']['dod']['std_err'],
                                        Colours.END
                                )),
                                '%sN/A%s'%(Colours.BLUE, Colours.END) if result['complete_integrand']['dod']['status'] not in ['SUCESS','UNSTABLE'] else (
                                    '%s%s%.5f +/- %.1g%s'%(cuts_sum_dod_colour, '+' if result['cuts_sum']['dod']['central']>=0. else '',
                                    result['cuts_sum']['dod']['central'], result['cuts_sum']['dod']['std_err'],Colours.END)
                                ),
                                '%sN/A%s'%(Colours.BLUE, Colours.END) if result['complete_integrand']['dod']['status'] not in ['SUCESS','UNSTABLE'] else (
                                    '%s%s%.5f +/- %.1g%s'%(Colours.BLUE,'+' if severity_central>=0. else '', severity_central, severity_std_err, Colours.END)
                                )
                            )+'\n  %s'%result['command']
                        )