#!/usr/bin/env python3

import os
import sys
import itertools
import math
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from pprint import pprint
import yaml

pjoin = os.path.join

file_path = os.path.dirname(os.path.realpath( __file__ ))
sys.path.insert(0, pjoin(file_path,os.pardir))
sys.path.insert(0, pjoin(file_path,os.pardir,os.pardir))

import vectors
import ltd_commons 
hyperparameters = ltd_commons.hyperparameters

try:
    # Import the rust bindings
    from ltd import LTD
except ImportError:
    raise BaseException(
        "Could not import the rust back-end 'ltd' module. Compile it first with:"
        " ./make_lib from within the pyNLoop directory." )

N_points = 100
studied_topology = 'manual_eeAA_amplitude_E'
studied_topology = sys.argv[1]
studied_amplitude = ''
if len(sys.argv) == 3:
    studied_amplitude = sys.argv[2]


with open("../topologies.yaml", 'r') as stream:
    try:
        topologies = yaml.load(stream)
    except yaml.YAMLError as exc:
        print(exc)
    for top in topologies:
        if top['name'] == studied_topology:
            topology = top
            break
    else:
        raise BaseException(
            "No topology known with this name %s" %studied_topology)

rust_instance = LTD(
        settings_file = pjoin(os.path.pardir,'hyperparameters.yaml'),
        topology_file = pjoin(os.path.pardir,'topologies.yaml'),
        top_name = studied_topology,
        amplitude_file = pjoin(os.path.pardir,'amplitudes.yaml'),
        amp_name = studied_amplitude, 
    )

def evaluate(point):
    evaluation = {}
    for ltd_cut_index, ltd_cut_structure in enumerate(topology['ltd_cut_structure']):
        cut_propagator_indices = [[index*c for index in range(1,len(topology['loop_lines'][i]['propagators'])+1)] if c!=0 else [0,] 
                                             for i, c in enumerate(ltd_cut_structure)]
        for cut_index, cut_structure  in enumerate(itertools.product( *cut_propagator_indices )):
            #res = rust_instance.evaluate_cut_ct_f128(
            res = rust_instance.evaluate_amplitude_cut_f128(
                        [ [ ( float(vi.real), float(vi.imag) ) for vi in v ] for v in point],
                         ltd_cut_index,  cut_index)
            
            if not math.isnan(res[0]) and not math.isnan(res[1]):
                evaluation[((ltd_cut_index, tuple(ltd_cut_structure)),(cut_index, cut_structure))] = complex(res[0],res[1])
            else:
                 print("Nan found for entry %s"%str(((ltd_cut_index, ltd_cut_structure),(cut_index, cut_structure))))

    return evaluation 

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        yield l[i:i + n]

def collinear_limit(mom, shift, delta):
    p1 = mom
    p2 = mom * [1,-1,-1,-1]
    
    q = [0., 1.2,12.,-1.45,0.03] #vectors.LorentzVector(np.random.rand(4))
    l_perp = -np.array([0.,p1[2]*q[3] - p1[3]*q[2],p1[3]*q[1] - p1[1]*q[3], p1[1]*q[2] - p1[2]*q[1]])

    c_mom = (shift + p1 + delta*p2 + np.sqrt(delta)*l_perp)
    return c_mom ,rust_instance.inv_parameterize_f128(c_mom[1:],0)[:3]

def soft_limit(mom, shift, delta):
    c_mom = (shift + delta * mom)
    return c_mom ,rust_instance.inv_parameterize_f128(c_mom[1:],0)[:3]

def uv_limit(mom, shift, delta):
    c_mom = shift + delta * mom
    return c_mom ,rust_instance.inv_parameterize_f128(c_mom[1:],0)[:3]

def line(x, a, b):
    return a*x+b 

import numpy as np
from scipy.optimize import curve_fit

min_logx = -15
max_logx = 0

#limits = ['Collinear p1', 'Collinear p2', 'Soft1','UV']
limits = ['UV']
#limits = ['Soft']
#limits = ['Collinear ' + p for p in ['p1','p2','p3','p4']] + ['Soft1', 'Soft3', 'UV']
#limits = ['Collinear ' + p for p in ['p1','p2']] + ['Soft1', 'UV']
#UV and SOFT random vector
np.random.seed(0)
l_UV = vectors.LorentzVector([3595.2511874644233,-6.148627708775928e-13,-3347.157026856964,-1312.391305414418])/1000. #vectors.LorentzVector(100*np.random.rand(4))
l_SOFT = vectors.LorentzVector(np.random.rand(4))

for limit in limits:
    plot_lines = {
        'integrand_re' : [],
        'integrand_im' : [],
        'deform_jac_re' : [],
        'deform_jac_im' : [],
        'param_jac' : [],
        'rescaled_integrand_re':[],
        'rescaled_integrand_im':[],
    }
    x_values = []
    plt.figure()
    for x in np.logspace(min_logx, max_logx, num=N_points):
        #Get limit
        p1= topology['external_kinematics'][0]
        p2= topology['external_kinematics'][2]
        p3= topology['external_kinematics'][3]
        p4= topology['external_kinematics'][4]
        scalar_scaling = 1.0
        if limit == 'Collinear p1':
            shift = -p1 ; y=0.3
            k ,rust_variables = collinear_limit(y*p1,shift, x)
            scalar_scaling = x**1 
        elif limit == 'Collinear p2': 
            shift = -p1 ; y=-0.3
            k ,rust_variables = collinear_limit(y*p2,shift, x)
            scalar_scaling = x**1 
        elif limit == 'Collinear p3': 
            shift = p4 ; y=0.3
            k ,rust_variables = collinear_limit(y*p3,shift, x)
            scalar_scaling = x**1 
        elif limit == 'Collinear p4': 
            shift = p4 ; y=-0.3
            k ,rust_variables = collinear_limit(y*p4,shift, x)
            scalar_scaling = x**1 
        elif limit == 'Soft1': 
            k ,rust_variables = soft_limit(l_SOFT,-p1, x)
            scalar_scaling = x**3
        elif limit == 'Soft3': 
            k ,rust_variables = soft_limit(l_SOFT,p4, x)
            scalar_scaling = x**3
        elif limit == 'UV': 
            x = 1./x
            k ,rust_variables = uv_limit(l_UV,[0,0,0,0], x)
            scalar_scaling = x**3

        print("Currently at log(x) = %f... x = %s"%(np.log(x),rust_variables))
        x_values.append(x) 

        integrand = rust_instance.evaluate(rust_variables)
        plot_lines['integrand_re'].append(integrand[0])
        plot_lines['integrand_im'].append(integrand[1])
        current_point = []
        parametrisation_jac = 1.
        for i, xis in enumerate(chunks(rust_variables,3)):
            p = rust_instance.parameterize_f128(xis, i)
            parametrisation_jac *= p[3]
            current_point.append(vectors.Vector([p[0],p[1],p[2]]))
        plot_lines['param_jac'].append(parametrisation_jac)
        #plot_lines['rescaled_integrand_re'].append(integrand[0] * scalar_scaling)
        #plot_lines['rescaled_integrand_im'].append(integrand[1] * scalar_scaling)
        
        kappas, jac_re, jac_im = rust_instance.deform([list(v) for v in current_point])
        deformed_point = [current_point[i]+vectors.Vector(kappas[i])*complex(0.,1.) for i in range(len(kappas))]

        #print jac_re, jac_im
        plot_lines['rescaled_integrand_re'].append(integrand[0] / parametrisation_jac * scalar_scaling * (2*math.pi)**3 )# / (jac_re**2-jac_im**2)*jac_re)
        plot_lines['rescaled_integrand_im'].append(integrand[1] / parametrisation_jac * scalar_scaling )# / (-jac_re**2+jac_im**2)*jac_im)
        
        plot_lines['deform_jac_re'].append(jac_re)
        plot_lines['deform_jac_im'].append(jac_im)
 
        duals = evaluate(deformed_point)
        for d, v in duals.items():
            if ('%d_%d_re'%(d[0][0],d[1][0]) in plot_lines) or ('%d_%d_im'%(d[0][0],d[1][0]) in plot_lines):
                plot_lines['%d_%d_re'%(d[0][0],d[1][0])].append(v.real * scalar_scaling)
                plot_lines['%d_%d_im'%(d[0][0],d[1][0])].append(v.imag * scalar_scaling)
            else:
                plot_lines['%d_%d_re'%(d[0][0],d[1][0])]=[v.real * scalar_scaling,]
                plot_lines['%d_%d_im'%(d[0][0],d[1][0])]=[v.imag * scalar_scaling,]

    #Sum of non-uv cuts
    plot_lines['sum_re'] = sum( [ np.array(plot_lines['0_%d_re'%i]) for i in range(len(topology['loop_lines'][0]['propagators'])-1)] )
    
    #popt, pcov = curve_fit(line, np.log(x_values),np.log(np.abs(plot_lines['0_0_re'])))
    #plt.plot(x_values,np.exp(popt[1])*x_values**popt[0], 'g--', label=r'fit subtracted: $\delta^{%2.1f}$' % popt[0], linewidth = 1.)
    
    popt, pcov = curve_fit(line, np.log(x_values[50:70]),np.log(np.abs(plot_lines['sum_re'][50:70])))
    plt.plot(x_values,np.exp(popt[1])*x_values**popt[0], 'g--', label=r'fit subtracted: $\delta^{%2.1f}$' % popt[0], linewidth = 1.)
    
    #Fit the integrand
    #popt, pcov = curve_fit(line, np.log(x_values[50:70]),np.log(np.abs(plot_lines['rescaled_integrand_re'][50:70])))
    #plt.plot(x_values,np.exp(popt[1])*x_values**popt[0], 'r--', label=r'fit subtracted: $\delta^{%2.1f}$' % popt[0], linewidth = 1.)
        
    selected = ['rescaled_integrand_re', '0_0_re','0_1_re', '0_2_re','0_3_re','0_4_re', '_sum_cut','sum_re','param_jac']
    veto_list= ['param_jac']
    #lines = [(k, (x_values, [abs(vi) for vi in v])) for k,v in sorted(plot_lines.items(), key=lambda el: el[0]) if
    #        (((k in selected) or ('ALL' in selected)) and ((k not in veto_list) or 'NONE' in veto_list ) and len(v)>0)]
    print(plot_lines['rescaled_integrand_re'])

    lines = [(k, (x_values, [abs(vi) for vi in v])) for k,v in sorted(plot_lines.items(), key=lambda el: el[0]) if
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
    plt.xscale('log')
    plt.title(limit)
    plt.legend(bbox_to_anchor=(0.75, 0.5))
plt.show()
