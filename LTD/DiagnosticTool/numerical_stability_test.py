import os
import sys
import itertools
import math
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

test_pointA = [
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

test_pointB = [
    vectors.Vector([
        complex(7330.4365861740625, 1.2127674225949402e-13),
        complex(-11077.751601272048, -1.827989605688507e-13),
        complex(12526.631402205292, 2.0624688337851502e-13),
    ]),
    vectors.LorentzVector([
        complex(2706.9752464072853, -3.63830226778482e-13),
        complex(-4090.779453307183, 5.483968817065521e-13),
        complex(4625.820130630993, -6.187406501355451e-13),
    ]),
]

test_pointC = [
    vectors.Vector([
        complex(7330.4365861740625, 1.2127674225949402e-13),
        complex(-11077.751601272048, -1.827989605688507e-13),
        complex(12526.631402205292, 2.0624688337851502e-13),
    ]),
    vectors.LorentzVector([
        complex(2706.9752464072853, -3.63830226778482e-13),
        complex(-4090.779453307183, 5.483968817065521e-13),
        complex(4625.820130630994, -6.187406501355451e-13),
    ]),
]

problem_point = [
    vectors.Vector([
        complex(1.3383294957377334e-18, 5.113407100880837e-2),
        complex(2.1856579328334187e-2, -3.914697128714517e-2),
        complex(0.0, 2.3647879246892718e-2),
    ]),
    vectors.LorentzVector([
        complex(-4.0149884872132e-18, -6.006674816249325e-3),
        complex(-2.1856579328334187e-2, 1.570561593315422e-2),
        complex(0.0, -5.135935434433841e-2),
    ]),
]

problem_point_bis = [
    vectors.Vector([
        complex(1.3383294957377338e-18, 4.043743860338496e-2),
        complex(2.185657932833419e-2, -6.090935705221293e-2),
        complex(0.0, 2.3647879246892718e-2),
    ]),
    vectors.LorentzVector([
        complex(-4.0149884872132e-18, -2.3656298515124727e-2),
        complex(-2.1856579328334187e-2, -2.0202687258461784e-2),
        complex(0.0, -5.135935434433841e-2),
    ]),
]

def evaluate(point):
    evaluation = {}
    for ltd_cut_index, ltd_cut_structure in enumerate(topology.ltd_cut_structure):
        cut_propagator_indices = [[index*c for index in range(1,len(topology.loop_lines[i].propagators)+1)] if c!=0 else [0,] 
                                             for i, c in enumerate(ltd_cut_structure)]
        for cut_index, cut_structure  in enumerate(itertools.product( *cut_propagator_indices )):
            res = rust_instance.evaluate_cut(
                        [ [ ( float(vi.real), float(vi.imag) ) for vi in v ] for v in point],
                         ltd_cut_index,  cut_index)
            #res = (0., 0.)
            if not math.isnan(res[0]) and not math.isnan(res[1]):
                evaluation[((ltd_cut_index, ltd_cut_structure),(cut_index, cut_structure))] = complex(res[0],res[1])
            else:
                 print("Nan found for entry %s"%str(((ltd_cut_index, ltd_cut_structure),(cut_index, cut_structure))))

    return evaluation 

#evaluationA = evaluate(test_pointA)
#evaluationB = evaluate(test_pointB)
#evaluationC = evaluate(test_pointC)
problem_evaluation = evaluate(problem_point)
print("Orig eval")
pprint(problem_evaluation)
problem_evaluation_bis = evaluate(problem_point_bis)
print("Bis eval")
pprint(problem_evaluation_bis)

#print("evaluation A")
#pprint(evaluationA)
#print("evaluation B")
#pprint(evaluationB)
#print("evaluation C")
#pprint(evaluationC)

#print("Evaluation A total: %s"%str(sum(evaluationA.values())))
#print("Evaluation B total: %s"%str(sum(evaluationB.values())))
#print("Evaluation C total: %s"%str(sum(evaluationC.values())))

#pprint(rust_instance.evaluate([0.75, 0.875, 0.5, 0.75, 0.375, 0.5]))
print("--A--")
pprint(rust_instance.evaluate([0.022046342097730487, 0.25, 0.5, 0.022046342097730487, 0.75, 0.5]))
#print("--B--")
#pprint(rust_instance.evaluate([0.022046342097730488, 0.25, 0.5, 0.022046342097730487, 0.75, 0.5]))
print("--C--")
pprint(rust_instance.evaluate([0.022046342097730489, 0.25, 0.5, 0.022046342097730487, 0.75, 0.5]))
#print("--D--")
#pprint(rust_instance.evaluate([0.022046342097730490, 0.25, 0.5, 0.022046342097730487, 0.75, 0.5]))
