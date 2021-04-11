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



protoc_file="./form_codes_amplitudes/out0.txt"
diagID=1

def generate_c_files_from_proto_c(protoc_file,diagID,out_dir="../"):
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
    header ="""#include <tgmath.h>
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
    form_output = temporary_fix_FORM_output(form_output)

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


    f64_code = form_output % {'mode':'',
                            'numbertype':'double complex', 'diagID':str(diagID),'header':header}
    # f128 needs function replacements
    f128_code = form_output % {'mode': '_f128',
                            'numbertype': '__complex128', 'diagID': diagID,'header':header}
    f128_code = float_pattern.sub(r'\1q', f128_code)
    f128_code = f128_code.replace('pow', 'cpowq').replace('sqrt', 'csqrtq').replace('log', 'clogq').replace('pi', 'M_PIq')

    #export

def generate_c_header(additional_variables:dict,out_dir="../"):
    """generates the c-header for amplitude integrands

    Args:
        additional_variables (dict): [contains the symbolic name and the numeric value of user-defined variables]
        out_dir (str, optional): [description]. Defaults to "../".
    """    
    header_file="""# ifndef NUM_H
# define NUM_H

#define  pi M_PI
#define  mUV params[0]
#define  mu params[1]
#define  small_mass_sq params[3]
"""
    for key,value in additional_variables.items():
        header_file+='\n'+'#define '+key+' '+str(value)
    header_file+='\n'+'# endif'
    # export


