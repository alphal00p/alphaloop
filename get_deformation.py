#!/usr/bin/env python2

import sys
import os
pjoin = os.path.join

root_path = os.path.dirname(os.path.realpath( __file__ ))
sys.path.insert(0, root_path)

import ltd
try:
    # Import the rust bindings
    from ltd import LTD
except ImportError:
    raise BaseException(
        "Could not import the rust back-end 'ltd' module. Compile it first with:"
        " ./make_lib from within the pyNLoop directory." )

args = sys.argv[1:]
if len(args)==0:
    print \
"""
Usage: ./get_deformation.py <topology_name>

Then, feed the desired loop momenta to deform in sys.stdin as follow: k1_E k1_x k1_y k1_z k2_E k2_x k2_y k2_z ...
and receive the deformations (kappas) on stdout
"""
    sys.exit(1)

studied_topology = args[0]
rust_instance = LTD(
        settings_file = pjoin(root_path, 'LTD','hyperparameters.yaml'),
        topology_file = pjoin(root_path, 'LTD','topologies.yaml'),
        amplitude_file = pjoin(root_path, 'LTD','amplitudes.yaml'),
        top_name = studied_topology,
        amp_name = '' 
    )

while True:
    momenta_to_deform = raw_input()
    if momenta_to_deform in ['','exit']:
        print "Exiting now."
        break
    try:
        momenta_to_deform = momenta_to_deform.split(' ')
        momenta_to_deform = [ [float(k) for k in momenta_to_deform[il*4:(il+1)*4]] for il in range(len(momenta_to_deform)/4) ]
    except Exception as e:
        print "Parsing error of: %s"%momenta_to_deform
        print "Exception: %s"%str(e)
        print "Exiting now."
        sys.exit(1)
    kappas, jac_re, jac_im = rust_instance.deform(momenta_to_deform)

    print ' '.join('%.16e'%(ke) for k in kappas for ke in k)
    sys.stdout.flush()
