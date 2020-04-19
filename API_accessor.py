#!/usr/bin/env python3
import argparse
import sys
import os
pjoin = os.path.join

root_path = os.path.dirname(os.path.realpath( __file__ ))
sys.path.insert(0, root_path)

from pprint import pprint

try:
    # Import the rust bindings
    from ltd import LTD, CrossSection
except ImportError:
    raise BaseException(
        "Could not import the rust back-end 'ltd' module. Compile it first with:"
        " ./make_lib from within the pyNLoop directory." )

_VALID_MODES = ['LTD', 'cross_section']

parser = argparse.ArgumentParser(description=
"""Wraps python bindings of Rust so as to expose them to Mathematica.
Available API functions exposed are:

    [mode: LTD and cross_section] exit
        --> Exits the wrapper

    [mode: LTD] get_deformation k1_E k1_x k1_y k1_z k2_E k2_x k2_y k2_z ...
        --> Returns: jac_re jac_im kappa1_E kappa1_x kappa1_y kappa1_z kappa2_E kappa2_x kappa2_y kappa2_z ...

    [mode: cross_section] get_deformation cut_ID k1_x k1_y k1_z k2_x k2_y k2_z ...
        --> Returns: kappa1_x_RE kappa1_y_RE kappa1_z_RE kappa2_x_RE kappa2_y_RE kappa2_z_RE ...
                      ... kappa1_x_IM kappa1_y_IM kappa1_z_IM kappa2_x_IM kappa2_y_IM kappa2_z_IM ...
        NOTE: The input momenta are in the loop momentum basis of the super graph while the deformed output
              momenta are in the cut basis of the select Cutkosky contribution

    [mode: cross_section] get_scaling cut_ID k1_x k1_y k1_z k2_x k2_y k2_z ...
        --> Returns: t_solution_A, t_jacobian_A, t_solution_B, t_jacobian_B
"""
)
requiredNamed = parser.add_argument_group('required named arguments')

# Required options
requiredNamed.add_argument('--name', '-n', type=str,
                    help='Name of the topology in LTD mode OR of the relative path to the squared topology file in cross_section mode.', required=True)

# Optional options
parser.add_argument('--mode', '-m', default=_VALID_MODES[0], choices=_VALID_MODES,
                    help='Select whether to access normal LTD bindings or Cross-section ones')
parser.add_argument('--hyperparameter-path', '--hpath', type=str, default=pjoin('LTD','hyperparameters.yaml'),
                    help='Path (relative to the pynloop root) of the hyperparameter yaml file.')
parser.add_argument('--topologies-path', '--tpath', type=str, default=pjoin('LTD','topologies.yaml'),
                    help='Path (relative to the pynloop root) of the topologies yaml file.')
parser.add_argument('--amplitudes-path', '--apath', type=str, default=pjoin('LTD','amplitudes.yaml'),
                    help='Path (relative to the pynloop root) of the amplitudes yaml file.')
parser.add_argument('--input', '-i', type=str, default='',
                    help='Specify an input string to process.')

args = parser.parse_args()


if args.mode=='LTD':
    rust_instance = LTD(
            settings_file = pjoin(root_path, args.hyperparameter_path),
            topology_file = pjoin(root_path, args.topologies_path),
            amplitude_file = pjoin(root_path, args.amplitudes_path),
            top_name = args.name,
            amp_name = '' 
        )
elif args.mode=='cross_section':
    rust_instance = CrossSection(
        pjoin(root_path,args.name),
        pjoin(root_path,args.hyperparameter_path)
    )

_SUPPORTED_API_FUNCTIONS = { 
    'LTD':
        ('exit', 'get_deformation'),
    'cross_section':
        ('exit', 'get_deformation', 'get_scaling'),
}

while True:
    if args.input == '':
        raw_input_str = input()
    else:
        raw_input_str = args.input
    raw_input_str = raw_input_str.split(' ')
    API_name = raw_input_str[0]
    raw_input_str = raw_input_str[1:]
    if API_name not in _SUPPORTED_API_FUNCTIONS[args.mode]:
        print("Unsupported API function: %s"%API_name)

    if API_name == 'exit':
        print("Exiting now.")
        break
    elif API_name == 'get_deformation':        
        if args.mode == 'LTD':
            momenta_input = [ [float(k) for k in raw_input_str[il*4:(il+1)*4]] for il in range(len(raw_input_str)//4) ]
            kappas, jac_re, jac_im = rust_instance.deform(momenta_input)
            print('%.16e '%jac_re+'%.16e '%jac_im+' '.join('%.16e'%(ke) for k in kappas for ke in k))
        elif args.mode == 'cross_section':
            cut_ID = int(raw_input_str[0])
            momenta_input = [ [float(k) for k in raw_input_str[1:][il*3:(il+1)*3]] for il in range(len(raw_input_str[1:])//3) ]
            deformed_momenta = rust_instance.get_cut_deformation(momenta_input, cut_ID)
            print(' '.join('%.16e'%(ke[0]) for k in deformed_momenta for ke in k)+' '+' '.join('%.16e'%(ke[1]) for k in deformed_momenta for ke in k))

    elif API_name == 'get_scaling':
        cut_ID = int(raw_input_str[0])
        momenta_input = [ [float(k) for k in raw_input_str[1:][il*3:(il+1)*3]] for il in range(len(raw_input_str[1:])//3) ]
        rescaling_solutions = rust_instance.get_scaling(momenta_input, cut_ID) 
        print(' '.join('%.16e'%(s) for ss in rescaling_solutions for s in ss))

    # Send out stdout
    sys.stdout.flush()
    if args.input != '':
        break
