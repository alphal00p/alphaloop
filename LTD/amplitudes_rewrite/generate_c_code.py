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

plugin_path = "/home/armin/my_programs/MG5_py_3_alphaloop_master/PLUGIN/alphaloop"

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
mode = "f64"

# TODO Remove once FORM will have fixed its C output bug


def temporary_fix_FORM_output(FORM_output):

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


with open("./form_codes_amplitudes/out0.txt", 'r') as f:
    form_output = f.read()
form_output = temporary_fix_FORM_output(form_output)


# regular expressions to be replaced in proto_c
var_pattern = re.compile(r'Z\d*_')
var_pattern_rest = re.compile(r'\b\w+Z\d*_')
input_pattern = re.compile(r'lm(\d*)')
norm_pattern = re.compile(r'norm')
divide_pattern = re.compile(r'1/')
wrapper_function_pattern = re.compile(r'\((\d*)\)')
float_pattern = re.compile(r'((\d+\.\d*)|(\.\d+))')


temp_vars = list(set(var_pattern.findall(form_output)))
temp_vars += list(set(var_pattern_rest.findall(form_output)))
temp_vars = list(sorted(temp_vars))

# eval function (we will always have one config )
eval_fct = """ void evaluate_PF_%(diagID)s(%(numbertype)s lm[], %(numbertype)s params[], int conf, %(numbertype)s* out) {
   switch(conf) {
		case 0: evaluate_PF_%(diagID)s_0%(mode)s(lm, params, out); return;
		default: *out = 0.;
    }
}"""

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

diagID = 1;
# header ="""#include <tgmath.h>
# # include <quadmath.h>
# # include <signal.h>
# # include \"numerator.h\""""
header = ""
f64_code = form_output % {'mode':'',
                          'numbertype':'double complex', 'diagID':str(diagID),'header':header}
f128_code = form_output % {'mode': '_f128',
                           'numbertype': '__complex128', 'diagID': diagID,'header':header}
f128_code = float_pattern.sub(r'\1q', f128_code)
print(f128_code)
