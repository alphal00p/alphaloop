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
abspath = os.path.abspath
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
    def from_runcard(AmpExporter, amplitude_runcard_path):
        """Initializes and AmpExporter from an amplitude runcard

        Args:
            amplitude_runcard_path (path): Path to an amplitude runcard
        Returns:
            Class: AmpExporter
        """

        with open(amplitude_runcard_path) as f:
            runcard_options = yaml.load(f, Loader=yaml.FullLoader)
        amplitude = runcard_options.get('amplitude')
        process_name = amplitude.get('name', 'my_amp')
        amplitude_options = runcard_options.get('amplitude_options', '')
        form_options = runcard_options.get('form_options', '')
        out_dir = os.path.abspath(pjoin(runcard_options.get(
            'dirs').get('out_dir', './'), process_name))
        alphaloop_dir = os.path.abspath(
            runcard_options.get('dirs').get('alphaloop_dir'))
        amp_exporter = AmpExporter(out_dir=out_dir, amplitude_dict=amplitude, process_name=process_name,
                                   form_options=form_options, amplitude_options=amplitude_options, alphaloop_dir=alphaloop_dir)
        return amp_exporter

    def generate_dir_structure(self, add_out_dir='', mode=''):
        """Generates the output structure needed for and amplitude

        Args:
            add_out_dir (str, optional): Allows for a modification of the output directory. Defaults to ''.
            mode (str, optional): If 'full', then all neccessary sub-directories are build. Defaults to ''.

        Returns:
            out_path (str): path to output directory
        """        
        out_dir = pjoin(self.out_dir, add_out_dir)
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        if not mode == 'full':
            # create only top_dir
            pass
        else:
            Path(pjoin(out_dir, 'Rust_inputs')).mkdir(
                parents=True, exist_ok=True)
            Path(pjoin(out_dir, 'lib')).mkdir(parents=True, exist_ok=True)
            Path(pjoin(out_dir, 'FORM', 'workspace')).mkdir(
                parents=True, exist_ok=True)
        return out_dir

    def export(self):
        """ Produces the output for an amplitude
        """        
        alphaloop_dir = self.alphaloop_dir
        out_dir = self.generate_dir_structure()

        amplitude = Amplitude.import_amplitude(self.amplitude_dict,self.amplitude_options)
        amplitude_list = amplitude.perform_color(
            out_dir=out_dir, alpha_dir=alphaloop_dir)

        for i, amp in amplitude_list:
            out_dir = self.generate_dir_structure(
                add_out_dir='color_struc_'+str(i), mode='full')
            # topo_gen = DariosYamlFileGenerator(amp)
            # topo_gen.export_yaml(out_dir)
            form_processor = FormProcessorAmp(amp, alphaloop_dir, form_options=self.form_options, form_wrk_space=pjoin(out_dir,'FORM','workspace'))
            form_processor.generate_output()


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
        self.masses = masses

    @classmethod
    def import_amplitude(Amplitude, amp, additional_options):
        """Defines an amplitude by importing either directly form the runcard or alternatively from a specified .py file

        Args:
            amp (dict): an amplitude as defined in a amplitude runcard
            additional_options (dict): additional options for amplitude computations

        Returns:
            the imported amplitude
        """        
        #import directly from ymal
        if isinstance(amp.get('diag_list'), dict):
            return Amplitude(name=amp.get('name', 'myamp'),
                             diag_list=amp.get('diag_list'),
                             masses=amp.get('masses', {}),
                             external_data=amp.get('external_data'),
                             constants=amp.get('constants', {}),
                             additional_options=additional_options
                             )
        #import from .py
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
        # TODO: implement. Is supposed to give amplitudeList: needs implementing
        self.color_struc = 'color(1)'
        return AmplitueList.from_amp(self)


class FormProcessorAmp():
    """This class allows for the generation of c-code for amplitude integrand evaluation
    """

    def __init__(self, amplitude, alphaloop_path, form_options={}, form_wrk_space=None):

        plugin_path = alphaloop_path

        FORM_processing_options = {
            'FORM_path': str(Path(plugin_path).parent.joinpath('libraries', 'form', 'sources', 'form').resolve()),
            'tFORM_path': str(Path(plugin_path).parent.joinpath('libraries', 'form', 'sources', 'tform').resolve()),
            # Define the extra aguments for the compilation
            'cores': multiprocessing.cpu_count(),
            'extra-options': {'OPTIMITERATIONS': 1000},
            'OPTIMISATIONSTRATEGY': 'CSEgreedy',
            'FORM_setup': {
                #   'MaxTermSize':'100K',
                #   'Workspace':'1G'
            }
        }

        self.FORM_workspace = form_wrk_space
        self.alphaloop_dir = alphaloop_path

        if isinstance(amplitude, amplitude):
            # FORM variables
            ext_data = amplitude.external_data
            self.FORM_variables = {
                'NINITIALMOMENTA': ext_data.get("n_in"),
                'NFINALMOMENTA': ext_data.get("n_out")-1,
                'NLOOPMOMENTA':  ext_data.get("n_loop"),
                'NPOL': len(ext_data.get("pol", [])),
                'NCPOL': len(ext_data.get("cpol", [])),
                'NSPINV': len(ext_data.get("spinor_v", [])),
                'NSPINU': len(ext_data.get("spinor_u", [])),
                'NSPINVBAR': len(ext_data.get("spinor_vbar", [])),
                'NSPINUBAR': len(ext_data.get("spinor_ubar", [])),
                'INDSHIFT_LIST': []
            }
            self.additional_options = amplitude.additional_options
            # update FORM options
            self.form_options = copy(
                FORM_processing_options).update(form_options)
            # integrand list
            if self.additional_options.get('integrand_per_diag'):
                self.integrands = [diag.get('analytic_integrand')
                                   for diag in (amplitude.diag_list)]
                for diag in (amplitude.diag_list):
                    self.FORM_variables['INDSHIFT_LIST'] = self.FORM_variables['INDSHIFT_LIST'] + [
                        diag.get('index_shift', 0)]

            else:
                self.integrands = ''
                ind_shift = 0
                for diag in (amplitude.diag_list):
                    self.integrands += '+'+diag.get('analytic_integrand')
                    if self.diag.get('index_shift', 0) > ind_shift:
                        ind_shift = self.diag.get('index_shift', 0)
                self.integrands = [self.integrands]
                self.FORM_variables['INDSHIFT_LIST'] = [ind_shift]
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

    def generate_c_files_from_proto_c(self, diagID):
        """Generates from a proto-c file obtained from FORM a valid c-routine for the evaluation of the integrand

        Args:
            diagID ([integer]): [defines the diagram]
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
        protoc_file = pjoin(self.FORM_workspace,
                            'out_integrand_PF_'+str(diagID)+".proto_c")
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

        # export
        f64_file = pjoin(pjoin(self.FORM_workspace, '../'),
                         'integrand_PF_'+str(diagID)+'_f64.c')
        f128_file = pjoin(pjoin(self.FORM_workspace, '../'),
                          'integrand_PF_'+str(diagID)+'_f128.c')

        with open(f64_file, "w") as c_code:
            c_code.write(f64_code)
        with open(f128_file, "w") as c_code:
            c_code.write(f128_code)

    def generate_c_header(self):
        """generates the c-header for amplitude integrands
        """
        header_file = """# ifndef NUM_H
# define NUM_H

#define  pi M_PI
#define  mUV params[0]
#define  mu params[1]
#define  small_mass_sq params[3]
"""
        for key, value in self.additional_constants.items():
            header_file += '\n'+'#define '+key+' '+str(value)
        header_file += '\n'+'# endif'
        out_file = pjoin(pjoin(self.FORM_workspace, '../'), 'numerator.h')
        with open(out_file, "w") as numfile:
            numfile.write(header_file)

    def generate_integrand_c_files(self):
        out_file = pjoin("../", self.FORM_workspace)
        out_file = pjoin(out_file, 'numerator.c')

        header = "# include <tgmath.h>\n# include <quadmath.h>\# include <signal.h>\n\n // integrands \n"
        c_routines_f64 = '\n\\ f64 routines '
        c_routines_f128 = '\n\\ f64 routines '
        eval_f64 = "void evaluate_PF(double complex lm[], double complex params[], int diag, int conf, double complex* out) {\n\tswitch(diag) {"
        eval_f128 = "void evaluate_PF_f128(__complex128 lm[], __complex128 params[], int diag, int conf, __complex128* out) {\n\tswitch(diag) {"

        for i in range(len(self.integrands)):
            c_routines_f64 += "\nvoid evaluate_PF_{}(double complex[], double complex[], int conf, double complex* out);".format(
                i)
            c_routines_f128 += "\nvoid evaluate_PF_{}_f128(__complex128[], __complex128[], int conf, __complex128* out);".format(
                i)
            eval_f64 += "\n\t\tcase {}: evaluate_PF_{}(lm, params, conf, out); return;".format(
                i)
            eval_f128 += "\n\t\tcase {}: evaluate_PF_{}_f128(lm, params, conf, out); return;".format(
                i)
        eval_f64 += "\n\t\tdefault: raise(SIGABRT);\n\t}\n}"
        eval_f128 += "\n\t\tdefault: raise(SIGABRT);\n\t}\n}"

        full_integrand = header + c_routines_f64 + \
            c_routines_f128 + eval_f64 + eval_f128
        with open(out_file, "w") as c_code:
            c_code.write(full_integrand)

    def generate_form_integrand_files(self):
        """ generates input_`SGID'.h files used in integrand_amp.frm
        """

        # set additional symbols
        additional_constants = format('\nS {};'.format(','.join(list(
            self.additional_constants.keys())) if len(self.additional_constants) > 0 else ''))

        integrand = additional_constants + "\n\nL F ="
        for SGID, inte in self.integrands:
            dia_integrand = integrand + inte + ';'
            filename = pjoin(self.FORM_workspace, "input_"+str(SGID)+".h")
            with open(filename, "w") as form_file:
                form_file.write(dia_integrand)

    def run_form(self, diagID):
        """Runs form to create proto-c code for the integrand of diagram diagID

        Args:
            diagID (int): identifier of diagram
        """

        # verify file integrity
        FORM_source = pjoin(self.FORM_workspace, "integrand_amp.frm")
        FORM_target = pjoin(self.FORM_workspace, "input_"+str(diagID)+".h")
        if not os.path.isfile(FORM_source):
            self.create_form_dir()
        if not os.path.isfile(FORM_target):
            self.generate_form_integrand_files()

        form_vars = copy.deepcopy(self.FORM_variables)
        form_vars.update(self.form_options["extra-options"])
        workers = self.form_options['cores']
        if workers > 1:
            form_exec = self.form_options['tFORM_path']
            FORM_cmd = ' '.join([
                form_exec,
            ] +
                ['-D %s=%s' % (k, v) for k, v in form_vars.items()] +
                ['-w %s' % workers] +
                ['-M', '-l', '-C', 'numerator_%s.log' % diagID] +
                [FORM_source, ]
            )
        else:
            form_exec = self.form_options['FORM_path']
            FORM_cmd = ' '.join([
                form_exec,
            ] +
                ['-D %s=%s' % (k, v) for k, v in form_vars.items()] +
                ['-M', '-l', '-C', 'numerator_%s.log' % diagID] +
                [FORM_source, ]
            )

        with open(pjoin(self.FORM_workspace, 'FORM_run_cmd_%d.exe' % diagID), 'w') as f:
            f.write(FORM_cmd)

        r = subprocess.run(FORM_cmd,
                           shell=True,
                           cwd=self.FORM_workspace,
                           capture_output=True)
        if r.returncode != 0 or not os.path.isfile(pjoin(self.FORM_workspace, 'out_%d.proto_c' % diagID)):
            error_message = "FORM processing failed with error:\n%s\nFORM command to reproduce:\ncd %s; %s" % (
                r.stdout.decode('UTF-8'),
                self.FORM_workspace, FORM_cmd)
            sys.exit(error_message)

    def compile_integrand(self):
        """compiles the c-library for rust
        """

        root_output_path = pjoin(self.FORM_workspace,'../')

        if os.path.isfile(pjoin(root_output_path,'Makefile')):
            try:
                print("Now compiling FORM-generated numerators with options: %s ..."%(' '.join(self.form_options['compilation-options'])))

                t = time.time()
                cmd = 'make'
                arg = ' all'
                nb_core=self.form_options["cores"]
                cwd=root_output_path,mode='cpp'
                try:
                    if  nb_core > 1:
                        cmd.append('-j%s' % nb_core)
                    cmd += arg
                    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, 
                                        stderr=subprocess.STDOUT, cwd=cwd)
                    (out, err) = p.communicate()
                except OSError as error:
                    raise OSError('Directory %s doesn\'t exists. Impossible to run make' % cwd)
                print("Compilation time: {:.2}s".format(time.time() - t))
            except OSError as e:
                print("Compilation of FORM-generated numerator failed:\n%s"%(str(e)))        
        
    def create_form_dir(self):
        """copies the neccessary files for form
        """        
        Path(self.FORM_workspace).mkdir(parents=True, exist_ok=True)
        src_path = pjoin(self.alphaloop_dir, 'LTD',
                         'amplitude_rewrite', 'form_codes_amplitudes')
        form_endings = [".h", ".frm", ".prc"]
        for ending in form_endings:
            flist = glob.glob(src_path + "/*"+ending)
            for ff in flist:
                shutil.copy(ff, self.FORM_workspace)
        src_path = pjoin(src_path,'FORM_output_makefile')
        target = pjoin(pjoin(self.FORM_workspace,'../'),'Makefile')
        shutil.copy2(src_path,target)

    def generate_output(self):
        """generates a compiled c-integrand
        """                
        self.create_form_dir()
        self.generate_form_integrand_files()
        self.generate_c_header()
        # TODO: include with some statistics similar to FORMPROCESSOR
        for int_id in range(len(self.integrands)):
            print("form on integrand %s"%int_id)
            self.run_form(int_id)
            self.generate_c_files_from_proto_c(int_id)
        self.generate_integrand_c_files()
        self.compile_integrand()

