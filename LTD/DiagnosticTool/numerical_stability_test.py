import os
import sys
import itertools
from pprint import pprint

pjoin = os.path.join

file_path = os.path.dirname(os.path.realpath( __file__ ))
sys.path.insert(0, pjoin(file_path,os.pardir))
sys.path.insert(0, pjoin(file_path,os.pardir,os.pardir))

import vectors
import ltd_commons 
hyperparameters = ltd_commons.hyperparameters
topology_collection = ltd_commons.hard_coded_topology_collection

try:
    # Import the rust bindings
    from ltd import LTD
except ImportError:
    raise BaseException(
        "Could not import the rust back-end 'ltd' module. Compile it first with:"
        " ./make_lib from within the pyNLoop directory." )

studied_topology = 'DoubleTriangle'
deformation_strategy = 'generic'
topology = topology_collection[studied_topology]

rust_instance = LTD(
        settings_file = pjoin(os.path.pardir,'hyperparameters.yaml'),
        topology_file = pjoin(os.path.pardir,'topologies.yaml'),
        name = studied_topology,
    )

test_point = [
    vectors.Vector([
        complex(7330.4365861740625, 1.2127674225949402e-13),
        complex(-11077.751601272048, -1.827989605688507e-13),
        complex(12526.631402205292, 2.0624688337851502e-13),
    ]),
    vectors.LorentzVector([
        complex(2706.9752464072853, -3.63830226778482e-13),
        complex(-4090.779453307183, 5.483968817065521e-13),
        complex(4625.820130630992, -6.187406501355451e-13),
    ]),
]

original_evaluation = {}
for ltd_cut_index, ltd_cut_structure in enumerate(topology.ltd_cut_structure):
    cut_propagator_indices = [[index*c for index in range(1,len(topology.loop_lines[i].propagators)+1)] if c!=0 else [0,] 
                                             for i, c in enumerate(ltd_cut_structure)]
    for cut_index, cut_structure  in enumerate(itertools.product( *cut_propagator_indices )):
        original_evaluation[((ltd_cut_index, ltd_cut_structure),(cut_index, cut_structure))] = rust_instance.evaluate_cut(
                        [ [ ( float(vi.real), float(vi.imag) ) for vi in v ] for v in test_point],
                        ltd_cut_index,  cut_index)

pprint(original_evaluation)
