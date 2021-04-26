import multiprocessing
import re
import copy

import os
from pathlib import Path
import time

import progressbar

import sys
import subprocess
import shutil
import glob
import json

abspath = os.path.abspath
pjoin = os.path.join

    







class AmplitudeList():
    """ Holds a list of amplitudes, e.g. neccesary if we split by color
    """

    def __init__(self, amp_list):
        for amp in amp_list:
            if not isinstance(amp, Amplitude):
                sys.exit(
                    "Can not initialize AmplitudeList. Not all elements in list are of type Amplitude")
        self.amplitudes = amp_list

    @classmethod
    def from_amp(AmplitudeList, amp):
        return AmplitudeList([amp])

    def append_amp(self, amp):
        self.amplitudes += [amp]


class Amplitude():
    """This class holds amplitudes which consist of propagators and numerators
    """

    def __init__(self, name='my_amp', diagram_list=[], external_data={}, masses={}, additional_options={}, constants={}):

        _AMP_OPTS_DEFAULT = {
            # if true, every integrand will get its own c-routine, else one integrand for the complete amplitude
            "integrand_per_diag": True,
            # if true: the color structure per graph will be factorized: e.g. Tr[T^a T^b T^c]-Tr[T^a T^c T^b] will be treated as ``independent'' color-structure
            # if false: every structure is treated independently
            #   Example
            #  q q~ > g g:    True:     3 color-structures appear {T^a T^b T^c,T^a T^c T^b,T^a T^b T^c-T^a T^c T^b} where the last one is due to the 3g-vertex
            #                 False:    2 color-structures for all diagrams {T^a T^b T^c,T^a T^c T^b}
            "color_per_graph": False,
            "colorless": False
        }
        _AMP_OPTS_DEFAULT.update(additional_options)

        _DERIVABLE_EXT_DATA = {"n_in": 0, "n_out": 0, "n_loop": 0}
        _MANDATORY_EXT_DATA = ["in_momenta", "out_momenta"]
        _NON_MANDATORY_EXT_DATA = {
            "spinor_v": [],
            "spinor_vbar": [],
            "spinor_u": [],
            "spinor_ubar": [],
            "pol": [],
            "cpol": []}

        for wavefct, val in _NON_MANDATORY_EXT_DATA.items():
            if wavefct not in external_data:
                external_data.update({wavefct: val})
        for manda in _MANDATORY_EXT_DATA:
            assert manda in external_data, "missing mandatory entry %s in external_data" % manda

        if len(diagram_list) > 0:
            _DERIVABLE_EXT_DATA.update({
                "n_in": len(((diagram_list[0].get('propagators'))[0]).get('incoming_signature')),
                "n_out": len(((diagram_list[0].get('propagators'))[0]).get('outgoing_signature'))+1,
                "n_loop": len(((diagram_list[0].get('propagators'))[0]).get('loop_signature'))})

        external_data.update(_DERIVABLE_EXT_DATA)

        self.name = name
        self.diagram_list = diagram_list
        self.external_data = copy.copy(external_data)
        self.additional_options = copy.deepcopy(_AMP_OPTS_DEFAULT)
        self.color_strucs = []

        # treat color:
        if self.additional_options.get('colorless', 'False'):
            for diag in self.diagram_list:
                diag.update({'color_struc': 'color(1)'})
                self.color_strucs += ['color(1)']
        else:
            for diag in self.diagram_list:
                if not 'color_struc' in diag.keys():
                    diag.update({'color_struc': 'derive'})
                self.color_strucs += [diag.get('color_struc')]
            if len(set(self.color_strucs)) > 1 and 'derive' in set(self.color_strucs):
                sys.exit(
                    'The color-structures of the amplitude: \n %s \n contain a mix of determined and undetermined color structures' % self.color_strucs)

        self.constants = copy.copy(constants)
        self.masses = copy.copy(masses)

    @ classmethod
    def import_amplitude(Amplitude, amp, additional_options):
        """Defines an amplitude by importing either directly form the runcard or alternatively from a specified .py file

        Args:
            amp (dict): an amplitude as defined in a amplitude runcard
            additional_options (dict): additional options for amplitude computations

        Returns:
            the imported amplitude
        """
        # import directly from ymal
        if isinstance(amp.get('diagram_list'), list):
            return Amplitude(name=amp.get('name', 'myamp'),
                             diagram_list=amp.get('diagram_list'),
                             masses=amp.get('masses', {}),
                             external_data=amp.get('external_data'),
                             constants=amp.get('constants', {}),
                             additional_options=additional_options
                             )
        # import from .py:
        else:
            p = Path(amp.get('diagram_list'))
            filename, file_extension = os.path.splitext(p)
            if file_extension == '.json':
                with open(p) as f:
                    amp_json = json.load(f)
                    amp.update(amp_json)
            elif file_extension == '.py':

                sys.path.insert(0, str(p.parent))
                m = __import__(p.stem)
                print("Imported {} diagrams.".format(len(m.diagram_list)))

                try:
                    amp["diagram_list"] = m.diagram_list
                except AttributeError:
                    sys.exit("{} has no diagram_list".format(
                        amp.get('diagram_list')))
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
            else:
                sys.exit("Unknown extension for diagram list import in: %s" % p)

        return Amplitude(
            name=amp.get('name', 'myamp'),
            diagram_list=amp.get('diagram_list'),
            masses=amp.get('masses', {}),
            external_data=amp.get('external_data'),
            constants=amp.get('constants', {}),
            additional_options=additional_options
        )

    def append_diagram_list(self, diag_list):
        self.diagram_list = self.diagram_list + diag_list
        self.color_strucs = self.color_strucs + \
            [diag.get('color_struc', 'derive') for diag in diag_list]
        if len(diag_list) > 0:
            self.external_data.update({
                "n_in": len(((diag_list[0].get('propagators'))[0]).get('incoming_signature')),
                "n_out": len(((diag_list[0].get('propagators'))[0]).get('outgoing_signature'))+1,
                "n_loop": len(((diag_list[0].get('propagators'))[0]).get('loop_signature'))})

    def split_color(self, out_dir='', alphaloop_dir=''):

        # with color-computation
        if 'derive' in self.color_strucs:
            # we need that, since we are going to initialize a Form processor which needs to compute color per integrand
            amp_options = copy.deepcopy(self.additional_options)
            self.additional_options.update({"integrand_per_diag": True})

            # initialize empty amplitude: we will append diagrams and color-strucs to it
            color_ordered_amp = Amplitude(name=self.name, external_data=self.external_data,
                                          masses=self.masses, additional_options=amp_options, constants=self.constants)
            form_processor_color = FormProcessorAmp(
                self, alphaloop_path=alphaloop_dir, form_wrk_space=pjoin(out_dir, 'FORM', 'workspace'))
            color_decomp = form_processor_color.compute_color(
                {'color_per_graph': self.additional_options.get('color_per_graph')})
            diag_list = self.diagram_list

            for diagID, colorsplit_integrand in enumerate(color_decomp):
                for color, numerator in colorsplit_integrand.items():
                    diag = copy.deepcopy(diag_list[diagID])
                    diag.update(
                        {'color_struc': color, 'analytic_num': numerator})
                    color_ordered_amp.append_diagram_list(
                        diag_list=[diag])

            # create amplitude list:
            color_ordered_amps = {cc: []
                                  for cc in set(color_ordered_amp.color_strucs)}
            for i, cc in enumerate(color_ordered_amp.color_strucs):
                color_ordered_amps[cc] = color_ordered_amps[cc] + \
                    [color_ordered_amp.diagram_list[i]]
            amp_list = [Amplitude(name=color_ordered_amp.name + '_color_'+str(i),
                                  diagram_list=dia_list,
                                  external_data=color_ordered_amp.external_data,
                                  additional_options=color_ordered_amp.additional_options,
                                  masses=color_ordered_amp.masses,
                                  constants=color_ordered_amp.constants)
                        for col, dia_list in color_ordered_amps.items()]

        # no color-computation
        else:
            color_ordered_amps = {cc: [] for cc in set(self.color_strucs)}
            for i, cc in enumerate(self.color_strucs):
                color_ordered_amps[cc] = color_ordered_amps[cc] + \
                    [self.diagram_list[i]]
            amp_list = [Amplitude(name=self.name + '_color_'+str(i),
                                  diagram_list=dia_list,
                                  external_data=self.external_data,
                                  additional_options=self.additional_options,
                                  masses=self.masses,
                                  constants=self.constants)
                        for col, dia_list in color_ordered_amps.items()]

        return AmplitudeList(amp_list)


class FormProcessorAmp():
    """This class allows for the generation of c-code for amplitude integrand evaluation
    """

    def __init__(self, amplitude, alphaloop_path, form_options={}, form_wrk_space=None):

        plugin_path = alphaloop_path

        _FORM_OPTS_DEFAULT = {
            'FORM_path': str(Path(plugin_path).joinpath('libraries', 'form', 'sources', 'form').resolve()),
            'tFORM_path': str(Path(plugin_path).joinpath('libraries', 'form', 'sources', 'tform').resolve()),
            # Define the extra aguments for the compilation
            'cores': multiprocessing.cpu_count(),
            'extra-options': {'OPTIMITERATIONS': 1000},
            'OPTIMISATIONSTRATEGY': 'CSEgreedy',
            'DEBUGLVL':0,
            'OPTIMLVL': 3,
            'FORM_setup': {
                #   'MaxTermSize':'100K',
                #   'Workspace':'1G'
            }
        }

        self.FORM_workspace = form_wrk_space
        self.alphaloop_dir = alphaloop_path

        if isinstance(amplitude, Amplitude):
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
            self.form_options = copy.deepcopy(
                _FORM_OPTS_DEFAULT)
            self.form_options.update(copy.copy(form_options))
            
            # integrand list
            if self.additional_options.get('integrand_per_diag')==True:
                self.integrands = [
                    self.derive_integrand_from_diag(
                        diag).get('analytic_integrand')
                    for diag in (amplitude.diagram_list)]
                for diag in (amplitude.diagram_list):
                    self.FORM_variables['INDSHIFT_LIST'] = self.FORM_variables['INDSHIFT_LIST'] + [
                        diag.get('index_shift', 0)]
            else:
                self.integrands = ''
                ind_shift = 0
                for diag in (amplitude.diagram_list):
                    self.integrands += '+' + \
                        self.derive_integrand_from_diag(
                            diag).get('analytic_integrand')
                    # this would conflict for color. However, we initialize with 'integrand_per_diag':False when computing color
                    if diag.get('index_shift', 0) > ind_shift:
                        ind_shift = diag.get('index_shift', 0)
                self.integrands = [self.integrands]
                
                self.FORM_variables['INDSHIFT_LIST'] = [ind_shift]
            # constants
            self.additional_constants = copy.copy(amplitude.constants)
            self.additional_constants.update(copy.copy(amplitude.masses))

        else:
            sys.exit("not an amp")

    def derive_integrand_from_diag(self, diag):
        """Takes a diagram and computes the integrand

        Args:
            diag (dict): contains propagators and analytic numerator
        """
        my_diag = copy.deepcopy(diag)
        num = my_diag.get('analytic_num', '1')
        prop_list = []

        # determine propagators
        for prop in my_diag.get('propagators'):
            mom_strg = ''
            ext_sig = prop.get('incoming_signature') + \
                prop.get('outgoing_signature')

            for ii, sig in enumerate(prop.get('loop_signature')):
                if sig != 0:
                    mom_strg += '+('+str(sig)+')*'+'k'+str(ii+1)
            for ii, sig in enumerate(ext_sig):
                if sig != 0:
                    mom_strg += '+('+str(sig)+')*'+'p'+str(ii+1)
            prop_list += ['(sprop({},{}))^{}'.format(mom_strg,
                                                     prop.get('mass', '0'), prop.get('power', '1'))]

        integrand = '('+num+')*'+'{}'.format(format('*'.join(prop_list)))
        my_diag.update({'analytic_integrand': integrand})

        return my_diag

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
        eval_fct = """ void evaluate_PF_%(diagID)s%(mode)s(%(numbertype)s lm[], %(numbertype)s params[], int conf, %(numbertype)s* out) {
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
        form_output = form_output.replace('pi_', 'pi').replace('i_', 'I')
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
        f64_code = f64_code.replace('pi', 'M_PI')

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

# define  pi M_PI
# define  mUV params[0]
# define  mu params[1]
# define  small_mass_sq params[3]
"""
        for key, value in self.additional_constants.items():
            header_file += '\n'+'#define '+key+' '+str(float(value))
        header_file += '\n'+'# endif'
        out_file = pjoin(pjoin(self.FORM_workspace, '../'), 'numerator.h')
        with open(out_file, "w") as numfile:
            numfile.write(header_file)

    def generate_integrand_c_files(self):
        out_file = abspath(pjoin(self.FORM_workspace, "../"))
        out_file = pjoin(out_file, 'integrand.c')

        header = "# include <tgmath.h>\n# include <quadmath.h>\n# include <signal.h>\n\n // integrands \n"
        c_routines_f64 = '\n // f64 routines \n'
        c_routines_f128 = '\n// f128 routines \n'
        eval_f64 = "\nvoid evaluate_PF(double complex lm[], double complex params[], int diag, int conf, double complex* out) {\n\tswitch(diag) {"
        eval_f128 = "\nvoid evaluate_PF_f128(__complex128 lm[], __complex128 params[], int diag, int conf, __complex128* out) {\n\tswitch(diag) {"

        for i in range(len(self.integrands)):
            c_routines_f64 += "\nvoid evaluate_PF_{}(double complex[], double complex[], int conf, double complex* out);\n".format(
                i)
            c_routines_f128 += "\nvoid evaluate_PF_{}_f128(__complex128[], __complex128[], int conf, __complex128* out);\n".format(
                i)
            eval_f64 += "\n\t\tcase {}: evaluate_PF_{}(lm, params, conf, out); return;\n".format(
                i, i)
            eval_f128 += "\n\t\tcase {}: evaluate_PF_{}_f128(lm, params, conf, out); return;\n".format(i,
                                                                                                       i)

        eval_f64 += "\n\t\tdefault: raise(SIGABRT);\n\t}\n}\n"
        eval_f128 += "\n\t\tdefault: raise(SIGABRT);\n\t}\n}\n"

        dummy_fcts = """void evaluate_LTD_f128(__complex128 lm[], __complex128 params[], int diag, int conf, __complex128* out) {
        abort();
        }\n
        void evaluate_LTD(__complex128 lm[], __complex128 params[], int diag, int conf, __complex128* out) {
        abort();
        }\n"""

        full_integrand = header + c_routines_f64 + eval_f64 \
            + c_routines_f128 + eval_f128 + '\n\n\n//dummies \n'+dummy_fcts
        with open(out_file, "w") as c_code:
            c_code.write(full_integrand)

    def generate_form_integrand_files(self):
        """ generates input_`SGID'.h files used in integrand_amp.frm
        """

        # set additional symbols
        if len(self.additional_constants) > 0:
            additional_constants = format('\nS {};'.format(','.join(list(
                self.additional_constants.keys())) if len(self.additional_constants) > 0 else ''))
        else:
            additional_constants = ''

        integrand = additional_constants + "\n\nL F ="
        for SGID, inte in enumerate(self.integrands):
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
        form_vars['OPTIMISATIONSTRATEGY'] = self.form_options['OPTIMISATIONSTRATEGY']
        form_vars['OPTIMLVL'] = self.form_options['OPTIMLVL']
        form_vars['DEBUGLVL'] = self.form_options['DEBUGLVL']
        del form_vars['INDSHIFT_LIST']
        form_vars['INDSHIFT'] = self.FORM_variables['INDSHIFT_LIST'][diagID]
        form_vars['SGID'] = diagID

        workers = self.form_options['cores']
        if workers > 1:
            form_exec = self.form_options['tFORM_path']
            FORM_cmd = ' '.join([
                form_exec,
            ] +
                ['-D %s=%s' % (k, v) for k, v in form_vars.items()] +
                ['-w%s' % workers] +
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
        # r.returncode != 0 seems to happen at random.
        # sometimes there is an error even though success exists. I suppose thats due to the not having enough time between write and read
        if not os.path.isfile(pjoin(self.FORM_workspace, 'success_%s.proto_c' % diagID)):
            error_message = "FORM processing failed with error:\n\nFORM command to reproduce:\ncd %s; %s" % (
                self.FORM_workspace, FORM_cmd)
            sys.exit(error_message)

    def compile_integrand(self):
        """compiles the c-library for rust
        """

        root_output_path = abspath(pjoin(self.FORM_workspace, '../'))

        if os.path.isfile(pjoin(root_output_path, 'Makefile')):
            try:
                print("Now compiling FORM-generated integrands")
                t = time.time()
                cmd = 'make'
                nb_core = self.form_options["cores"]
                cwd = root_output_path
                try:
                    if nb_core > 1:
                        cmd + '-j%s' % nb_core

                    p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                         stderr=subprocess.STDOUT, cwd=cwd)
                    (out, err) = p.communicate()

                    if 'Error' in str(out) or not str(err) == 'None':
                        sys.exit(
                            'Could not compile integrand. \n make failed with\n: {} .'.format(
                                str(out))
                        )
                except OSError as error:
                    raise OSError(
                        'make all failed with\n: %s \n. Could not compile integrand.' % error)
                print("Compilation time: {:.2}s".format(time.time() - t))
            except OSError as e:
                print("Compilation of FORM-generated numerator failed:\n%s" % (str(e)))
        else:
            raise OSError(
                'Makefile in: \n %s \n does not exist. Compilation failure.' % root_output_path)

    def create_form_dir(self):
        """copies the neccessary files for form
        """
        Path(self.FORM_workspace).mkdir(parents=True, exist_ok=True)
        src_path = pjoin(self.alphaloop_dir, 'LTD',
                         'amplitudes_rewrite', 'form_codes_amplitudes')
        form_endings = [".h", ".frm", ".prc"]
        for ending in form_endings:
            flist = glob.glob(src_path + "/*"+ending)
            for ff in flist:
                shutil.copy(ff, self.FORM_workspace)
        src_path = pjoin(src_path, 'FORM_output_makefile')
        target = pjoin(pjoin(self.FORM_workspace, '../'), 'Makefile')
        shutil.copy2(src_path, target)

    def generate_output(self):
        """generates a compiled c-integrand
        """
        self.create_form_dir()
        self.generate_form_integrand_files()
        self.generate_c_header()
        with progressbar.ProgressBar(prefix='process integrand: ', max_value=len(self.integrands)) as bar:
            for int_id in range(len(self.integrands)):
                self.run_form(int_id)
                self.generate_c_files_from_proto_c(int_id)
                bar.update(int_id+1)
        self.generate_integrand_c_files()
        self.compile_integrand()

    def compute_color(self, color_option={'color_per_graph': False}):
        """Takes an integrand and performs color on it

        Args:
            color_option (dict, optional): Mode of color-treatment. Defaults to {'color_per_graph': False}.

        Returns:
            List: List of dictionaries. 
            Every element has the form: {color_struc1:numerator1,...,color_strucN:numeratorN}
            where color_struc denotes one of the independent color structures. The decomposition is trace based.
        """

        self.create_form_dir()
        self.generate_form_integrand_files()
        color_dicts = []
        with progressbar.ProgressBar(prefix='process color: ', max_value=len(self.integrands)) as bar:
            for diagID, inte in enumerate(self.integrands):
                self.run_color_form(diagID, color_option=color_option)
                with open(pjoin(self.FORM_workspace, 'SG_%d_color_decomp.txt' % diagID), 'r') as f:
                    color_dicts += [copy.deepcopy(
                        eval((f.read()).replace('\n', '')))]
                bar.update(diagID+1)
        return color_dicts

    def run_color_form(self, diagID, color_option={'color_per_graph': False}):
        """Runs form to compute the color_structure for the integrand of diagram diagID

        Args:
            diagID (int): identifier of diagram
        """

        form_vars = {
            'INDSHIFT': self.FORM_variables['INDSHIFT_LIST'][diagID],
            'SGID': diagID,
            'CPERGRAPH': 1 if color_option.get('color_per_graph') else 0}
        # verify file integrity
        FORM_source = pjoin(self.FORM_workspace, "color_basis.frm")
        FORM_target = pjoin(self.FORM_workspace, "input_"+str(diagID)+".h")
        if not os.path.isfile(FORM_source):
            self.create_form_dir()
        if not os.path.isfile(FORM_target):
            self.generate_form_integrand_files()

        form_exec = self.form_options['FORM_path']
        FORM_cmd = ' '.join([
            form_exec,
        ] +
            ['-D %s=%s' % (k, v) for k, v in form_vars.items()] +
            ['-M', '-l', '-C', 'color_computation_%s.log' % diagID] +
            [FORM_source, ]
        )

        with open(pjoin(self.FORM_workspace, 'FORM_run_cmd_%d.exe' % diagID), 'w') as f:
            f.write(FORM_cmd)

        r = subprocess.run(FORM_cmd,
                           shell=True,
                           cwd=self.FORM_workspace,
                           capture_output=True)
        # r.returncode != 0 seems to happen at random.
        if not os.path.isfile(pjoin(self.FORM_workspace, 'SG_%d_color_decomp.txt' % diagID)):
            error_message = "FORM processing failed with error:\n%s\nFORM command to reproduce:\ncd %s; %s" % (r.stdout.decode('UTF-8'),
                                                                                                               self.FORM_workspace, FORM_cmd)
            sys.exit(error_message)
