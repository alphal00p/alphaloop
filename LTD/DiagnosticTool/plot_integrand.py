#!/usr/bin/env python3

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
import ltd_utils
import ltd_commons 
hyperparameters = ltd_commons.hyperparameters
topology_collection = ltd_utils.TopologyCollection.import_from(pjoin(file_path,os.pardir,'topologies.yaml'))

try:
    # Import the rust bindings
    from ltd import LTD
except ImportError:
    raise BaseException(
        "Could not import the rust back-end 'ltd' module. Compile it first with:"
        " ./make_lib from within the pyNLoop directory." )

N_points = 100

studied_topology = 'Box_massless'
studied_topology = 'Decagon_P2_physical_massless'
studied_topology = 'Decagon_P2_one_ellipse_massless'
studied_topology = 'Tringigon_P1_one_ellipse'
studied_topology = 'Tringigon_P2_physical_few_ellipses'
studied_topology = 'Tringigon_P2_physical_many_ellipses'
studied_topology = 'DoubleBox'

studied_topology = 'Decagon_P1_one_ellipse_massless'
studied_topology = sys.argv[1]

topology = topology_collection[studied_topology]

scan_args = eval(' '.join(sys.argv[2:]))
rust_instance = LTD(
        settings_file = pjoin(file_path, os.path.pardir,'hyperparameters.yaml'),
        topology_file = pjoin(file_path, os.path.pardir,'topologies.yaml'),
        amplitude_file = pjoin(file_path, os.path.pardir,'amplitudes.yaml'),
        top_name = studied_topology,
        amp_name = ''
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
                evaluation[((ltd_cut_index, tuple(ltd_cut_structure)),(cut_index, cut_structure))] = complex(res[0],res[1])
            else:
                 print("Nan found for entry %s"%str(((ltd_cut_index, ltd_cut_structure),(cut_index, cut_structure))))

    return evaluation 


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        yield l[i:i + n]

plot_lines = {
    'integrand_re' : [],
    'integrand_im' : [],
    'deform_jac_re' : [],
    'deform_jac_im' : [],
    'param_jac' : []
}

x_values = []
#min_x = 0.15619610
#max_x = 0.15619612
min_x = 0.0
max_x = 1.0
# Problem at [0.4474308180164988, 0.6249338072526623, 0.32786022901064904]
#min_x = 0.4476
#max_x = 0.4478
#min_x = 0.3260626059739
#max_x = 0.3260626059740
for t in range(1,N_points+1):
    if t%100==0:
        print("Currently at sample #%d..."%t)
    x = min_x+(float(t)/float(N_points+1))*(max_x-min_x)
    x_values.append(x)
    random_variables = [(v if v>=0. else x) for v in scan_args]
    integrand = rust_instance.evaluate(random_variables)
    plot_lines['integrand_re'].append(integrand[0])
    plot_lines['integrand_im'].append(integrand[1])
    
    current_point = []
    parametrisation_jac = 1.
    for i, xis in enumerate(chunks(random_variables,3)):
        p = rust_instance.parameterize(xis, i)
        parametrisation_jac *= p[3]
        current_point.append(vectors.Vector([p[0],p[1],p[2]]))
    plot_lines['param_jac'].append(parametrisation_jac)

    kappas, jac_re, jac_im = rust_instance.deform([list(v) for v in current_point])
    deformed_point = [current_point[i]+vectors.Vector(kappas[i])*complex(0.,1.) for i in range(len(kappas))]

    plot_lines['deform_jac_re'].append(jac_re)
    plot_lines['deform_jac_im'].append(jac_im)

    duals = evaluate(deformed_point)
    for d, v in duals.items():
        if ('%d_%d_re'%(d[0][0],d[1][0]) in plot_lines) or ('%d_%d_im'%(d[0][0],d[1][0]) in plot_lines):
            plot_lines['%d_%d_re'%(d[0][0],d[1][0])].append(v.real)
            #plot_lines['%d_%d_im'%(d[0][0],d[1][0])].append(v.imag)            
        else:
            plot_lines['%d_%d_re'%(d[0][0],d[1][0])]=[v.real,]
            #plot_lines['%d_%d_im'%(d[0][0],d[1][0])]=[v.imag,]

selected = ['integrand_re', 'integrand_im', '0_0_re', '0_0_im', '0_1_re', '0_1_im','deform_jac_re', 'deform_jac_im','ALL']
veto_list=['param_jac','deform_jac_im']


#lines = [(k, (x_values, [abs(vi) for vi in v])) for k,v in sorted(plot_lines.items(), key=lambda el: el[0]) if
#        (((k in selected) or ('ALL' in selected)) and ((k not in veto_list) or 'NONE' in veto_list ) and len(v)>0)]

lines = [(k, (x_values, [vi for vi in v])) for k,v in sorted(plot_lines.items(), key=lambda el: el[0]) if
        (((k in selected) or ('ALL' in selected)) and ((k not in veto_list) or 'NONE' in veto_list ) and len(v)>0)]

NORMALISE = False
for line_name, (x_data, y_data) in lines:
    if NORMALISE:
        max_y_data = max(y_data)
        if max_y_data > 0.:
            y_data = [ y/float(max_y_data) for y in y_data ]
    if 'integrand' in line_name:
        plt.plot(x_data, y_data, label=line_name,linewidth=2)
    else:
        plt.plot(x_data, y_data, label=line_name)

plt.yscale('log')
plt.legend(bbox_to_anchor=(0.75, 0.5))
plt.show()
