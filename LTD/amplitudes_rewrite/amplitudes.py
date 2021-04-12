import multiprocessing
import re
import copy
import logging
import os
from pathlib import Path
from pprint import pprint, pformat
import math
import igraph
import time
import numpy as np
from itertools import combinations_with_replacement
from collections import OrderedDict
import glob as glob_module

import progressbar
from itertools import chain, product
import sys
import subprocess
import argparse
import shutil
import py_compile
import sympy as sp
from sympy.simplify.simplify import simplify
from warnings import catch_warnings
import yaml
from yaml import Loader, Dumper
import glob

src_path = os.path.dirname(os.path.realpath(__file__))
pjoin = os.path.join


class AmpExporter():
    def __init__(self, out_dir=None, alphaloop_dir=None, process_name=None, amplitude_dict=None, amplitude_options=None, form_options=None):
        self.out_dir = out_dir
        self.alphaloop_dir = alphaloop_dir
        self.process_name = process_name
        self.amplitude_dict = amplitude_dict
        self.amplitude_options = amplitude_options
        self.form_options = form_options

    @classmethod
    def from_runcard(cls, amplitude_runcard_path):
        with open(amplitude_runcard_path) as f:
            runcard_options = yaml.load(f, Loader=yaml.FullLoader)
        # make proper path

        amplitude = runcard_options.get('amplitude')
        process_name = amplitude.get('name', 'my_amp')
        amplitude_options = runcard_options.get('amplitude_options', '')
        form_options = runcard_options.get('form_options', '')
        out_dir = os.path.abspath(pjoin(runcard_options.get('dirs').get('out_dir','./'),process_name))
        alphaloop_dir = os.path.abspath(runcard_options.get('dirs').get('alphaloop_dir'))
        amp_exporter = AmpExporter(out_dir=out_dir, amplitude=amplitude, process_name=process_name,
                                   form_options=form_options, amplitude_options=amplitude_options, alphaloop_dir=alphaloop_dir)
        return amp_exporter

    def generate_dir_structure(self, add_out_dir='', mode=''):
        # make proper path
        
        out_dir = pjoin(self.out_dir,add_out_dir)
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        if not mode == 'full':
            pass
            # create only top_dir
        else:            
            Path(pjoin(out_dir,'Rust_inputs')).mkdir(parents=True, exist_ok=True)
            Path(pjoin(out_dir,'lib')).mkdir(parents=True, exist_ok=True)
            Path(pjoin(out_dir,'FORM','workspace')).mkdir(parents=True, exist_ok=True)
        return out_dir

    def export(self):
        alphaloop_dir = self.alphaloop_dir
        out_dir = self.generate_dir_structure()
        amplitude = Amplitude.import_amplitude(self.amplitude_data)
        amplitude_list = amplitude.perform_color(
            out_dir=out_dir, alpha_dir=alphaloop_dir)

        for i, amp in amplitude_list:
            # TODO: do a proper export
            out_dir = self.generate_dir_structure(
                add_out_dir='color_struc_'+str(i), mode='full')
            topo_gen = DariosYamlFileGenerator(amp)
            topo_gen.export_yaml(out_dir)
            form_processor = FormProcessorAmp(amp)
            form_processor.generate_output(alphaloop_dir, out_dir)


class DariosYamlFileGenerator():
    def __init__(self):
        pass

    def export_yaml(self):
        pass


class AmplitueList():
    """ Holds a list of amplitudes, e.g. neccesary if we split by color
    """

    def __init__(self, amplitudes):
        self.amplitudes = amplitudes

    @classmethod
    def from_amp(amp):
        return AmplitueList([amp])

    def append_amp(self, amp):
        self.amplitudes += [amp]


class Amplitude():
    """This class holds amplitudes which consist of propagators and numerators
    """

    def __init__(self, name=None, diag_list=None, external_data=None, masses=None, additional_options={}):
        self.name = name
        self.diag_list = diag_list
        self.external_data = external_data
        self.additional_options = additional_options
        self.color_struc = None

    @classmethod
    def import_amplitude(self, amp, additional_options):
        if isinstance(amp.get('diag_list'), dict):
            return Amplitude(name=amp.get('name', 'myamp'),
                             diag_list=amp.get('diag_list'),
                             masses=amp.get('masses', {}),
                             external_data=amp.get('external_data'),
                             constants=amp.get('constants', {}),
                             additional_options=additional_options
                             )
        else:
            from pathlib import Path
            p = Path(amp.get('diag_list'))
            sys.path.insert(0, str(p.parent))
            m = __import__(p.stem)
            print("Imported {} diagrams.".format(len(m.diag_list)))

            try:
                amp["diag_list"] = m.diag_list
            except AttributeError:
                sys.exit("{} has no diag_list".format(amp.get('diag_list')))
            try:
                amp["masses"] = m.masses
            except AttributeError:
                pass
            try:
                amp["external_data"] = m.external_data
            except AttributeError:
                pass
            try:
                amp["constants"] = m.constants
            except AttributeError:
                pass
            try:
                amp["name"] = m.name
            except AttributeError:
                pass

            return Amplitude(
                name=amp.get('name', 'myamp'),
                diag_list=amp.get('diag_list'),
                masses=amp.get('masses', {}),
                external_data=amp.get('external_data'),
                constants=amp.get('constants', {}),
                additional_options=additional_options
            )

    def perform_color(self, out_dir='', alphaloop_dir=''):
        # TODO: implement. Is supposed to give amplitudeList
        self.color_struc = 'color(1)'
        return AmplitueList.from_amp(self)


class FormProcessorAmp():
    """This class allows for the generation of c-code for amplitude integrand evaluation
    """

    def __init__(self, amplitude, alphaloop_path, form_options={},form_wrk_space=None):

        plugin_path = alphaloop_path

        FORM_processing_options = {
            'FORM_path': str(Path(plugin_path).parent.joinpath('libraries', 'form', 'sources', 'form').resolve()),
            'tFORM_path': str(Path(plugin_path).parent.joinpath('libraries', 'form', 'sources', 'tform').resolve()),
            # Define the extra aguments for the compilation
            'compilation-options': [],
            'cores': 2,  # multiprocessing.cpu_count(),
            'extra-options': {'OPTIMITERATIONS': 1000},
            'optimisation_strategy': 'CSEgreedy',
            'FORM_setup': {
                #   'MaxTermSize':'100K',
                #   'Workspace':'1G'
            }
        }
        
        self.FORM_workspace = form_wrk_space
        self.alphaloop_dir = alphaloop_path

        if isinstance(amplitude, amplitude):
            self.additional_options = amplitude.additional_options
            # update FORM options
            self.form_options = copy(
                FORM_processing_options).update(form_options)
            # integrand list
            if self.additional_options.get('integrand_per_diag'):
                self.integrands = [diag.get('analytic_integrand')
                                   for diag in (amplitude.diag_list)]
            else:
                self.integrands = ''
                for diag in (amplitude.diag_list):
                    self.integrands += '+'+diag.get('analytic_integrand')
                self.integrands = [self.integrands]
            # constants
            self.additional_constants = copy.copy(amplitude.masses)
            self.additional_constants = self.additional_constants.update(
                amplitude.constants)
        else:
            sys.exit("not an amp")

    # TODO Remove once FORM will have fixed its C output bug

    def temporary_fix_FORM_output(self, FORM_output):

        new_output = []
        previous_line = None
        for line in FORM_output.split('\n'):
            if line.startswith('      _ +=  '):
                line = '      %s' % line[12:]
                if previous_line is not None:
                    new_output.append(previous_line[:-1])
                previous_line = line
            else:
                if previous_line is not None:
                    new_output.append(previous_line)
                previous_line = line
        if previous_line is not None:
            new_output.append(previous_line)

        return '\n'.join(new_output)

    def generate_c_files_from_proto_c(self, protoc_file, diagID):
        """Generates from a proto-c file obtained from FORM a valid c-routine for the evaluation of the integrand

        Args:
            protoc_file ([path to file]): [contains FORM optimized proto-c output]
            diagID ([integer]): [defines the diagram]
            out_dir (str, optional): [description]. Defaults to "../".
        """

        # SETUP UNIVERSAL REPLACEMENTS AND PARTS
        # regular expressions to be replaced in proto_c
        var_pattern = re.compile(r'Z\d*_')
        var_pattern_rest = re.compile(r'\b\w+Z\d*_')
        input_pattern = re.compile(r'lm(\d*)')
        norm_pattern = re.compile(r'norm')
        divide_pattern = re.compile(r'1/')
        wrapper_function_pattern = re.compile(r'\((\d*)\)')
        float_pattern = re.compile(r'((\d+\.\d*)|(\.\d+))')

        # header file inclusion
        header = """#include <tgmath.h>
# include <quadmath.h>
# include <signal.h>
# include \"numerator.h\""""

        # eval function (we will always have one config )
        eval_fct = """ void evaluate_PF_%(diagID)s(%(numbertype)s lm[], %(numbertype)s params[], int conf, %(numbertype)s* out) {
        switch(conf) {
            case 0: evaluate_PF_%(diagID)s_0%(mode)s(lm, params, out); return;
            default: *out = 0.;
        }
        }"""

        # GENERATION OF C-CODE in f64 and f128
        # read in
        with open(protoc_file, 'r') as f:
            form_output = f.read()
        form_output = self.temporary_fix_FORM_output(form_output)

        # determine temp vars from optimization
        temp_vars = list(set(var_pattern.findall(form_output)))
        temp_vars += list(set(var_pattern_rest.findall(form_output)))
        temp_vars = list(sorted(temp_vars))

        # prepend, header, function declarations and tempvar definition
        form_output = '\n%(header)s\n\nstatic inline void evaluate_PF_%(diagID)s_0%(mode)s\
        ( %(numbertype)s lm, %(numbertype)s params[], %(numbertype)s* out)' \
        + '{{\n\t{}\n\n{}\n}}'.format('%(numbertype)s {};'.format(','.join(temp_vars)) if len(temp_vars) > 0 else '', form_output
                                      )
        form_output = input_pattern.sub(r'lm[\1]', form_output)
        form_output = divide_pattern.sub(r'1./', form_output)
        form_output = norm_pattern.sub(r'1./', form_output)
        form_output = wrapper_function_pattern.sub(r'\1', form_output)
        # append eval:
        form_output += '\n\n'+eval_fct

        f64_code = form_output % {'mode': '',
                                  'numbertype': 'double complex', 'diagID': str(diagID), 'header': header}
        # f128 needs function replacements
        f128_code = form_output % {'mode': '_f128',
                                   'numbertype': '__complex128', 'diagID': diagID, 'header': header}
        f128_code = float_pattern.sub(r'\1q', f128_code)
        f128_code = f128_code.replace('pow', 'cpowq').replace(
            'sqrt', 'csqrtq').replace('log', 'clogq').replace('pi', 'M_PIq')

        # TODO: export the c-code

    def generate_c_header(self, additional_variables: dict, out_dir):
        """generates the c-header for amplitude integrands

        Args:
            additional_variables (dict): [contains the symbolic name and the numeric value of user-defined variables]
            out_dir (str, optional): [description]. Defaults to "../".
        """
        header_file = """# ifndef NUM_H
# define NUM_H

#define  pi M_PI
#define  mUV params[0]
#define  mu params[1]
#define  small_mass_sq params[3]
"""
        for key, value in additional_variables.items():
            header_file += '\n'+'#define '+key+' '+str(value)
        header_file += '\n'+'# endif'
        # export

    def generate_integrand_c_file(self):
        #TODO: implement
        pass

    def generate_form_integrand_files(self):
        """ generates input.h files used in integrand_amp.frm
        """

        # set additional symbols
        additional_constants = format('\nS {};'.format(','.join(list(
            self.additional_constants.keys())) if len(self.additional_constants) > 0 else ''))

        integrand = additional_constants + "L F ="
        for SGID, inte in self.integrands:
            dia_integrand = integrand + inte + ';'
            # TODO: export integrand with SGID

    def run_form(self):
        # TODO: derive form command and run
        pass

    def compile_integrand(self):
        # TODO: write routine for compilation
        pass

    def create_form_dir(self, alpha_loop_path):        
        Path(self.FORM_workspace).mkdir(parents=True, exist_ok=True)
        src_path=pjoin(self.alphaloop_dir,'LTD','amplitude_rewrite','form_codes_amplitudes')
        form_endings = [".h",".frm",".prc"]
        for ending in form_endings:
            flist = glob.glob(src_path + "/*"+ending)
            for ff in flist:
                shutil.copy(ff,self.FORM_workspace)

        

    def generate_output(self, alpha_loop_path, out_dir):
        # TODO: make this a sensible functions with some statistics similar to FORMPROCESSOR
        self.create_form_dir(out_dir)
        self.generate_form_integrand_files(out_dir)
        self.generate_c_header(out_dir)
        self.generate_c_files_from_proto_c(out_dir)
        self.generate_integrand_c_file(out_dir)
        self.compile_integrand(out_dir)
