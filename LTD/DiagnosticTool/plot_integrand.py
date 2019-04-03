#!/usr/bin/env python

import os
import sys
import itertools
import math
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
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

N_points = 1000

studied_topology = 'Box_massless'
studied_topology = 'Decagon_P2_physical_massless'
studied_topology = 'Decagon_P2_one_ellipse_massless'
studied_topology = 'Tringigon_P1_one_ellipse'
studied_topology = 'Tringigon_P2_physical_few_ellipses'
studied_topology = 'Tringigon_P2_physical_many_ellipses'

studied_topology = 'DoubleBox'

studied_topology = 'Decagon_P2_physical_massless'

topology = topology_collection[studied_topology]


scan_args = eval(' '.join(sys.argv[1:]))
rust_instance = LTD(
        settings_file = pjoin(os.path.pardir,'hyperparameters.yaml'),
        topology_file = pjoin(os.path.pardir,'topologies.yaml'),
        name = studied_topology,
    )

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


plot_lines = {
    'integrand_re' : [],
    'integrand_im' : []
}

x_values = []
min_x = 0.0
max_x = 1.0
for t in range(1,N_points+1):
    if t%100==0:
        print "Currently at sample #%d..."%t
    x = min_x+(float(t)/float(N_points+1))*max_x
    x_values.append(x)
    random_variables = [(v if v>=0. else x) for v in scan_args]
    integrand = rust_instance.evaluate(random_variables)
    plot_lines['integrand_re'].append(integrand[0])
    plot_lines['integrand_im'].append(integrand[1])


lines = [(k, (x_values, [abs(vi) for vi in v])) for k,v in sorted(plot_lines.items(), key=lambda el: el[0]) if
        ((not k in []) and len(v)>0)]

NORMALISE = False
for line_name, (x_data, y_data) in lines:
    if NORMALISE:
        max_y_data = max(y_data)
        if max_y_data > 0.:
            y_data = [ y/float(max_y_data) for y in y_data ]

    plt.plot(x_data, y_data, label=line_name)
    
plt.yscale('log')
plt.legend(bbox_to_anchor=(0.75, 0.5))
plt.show()
