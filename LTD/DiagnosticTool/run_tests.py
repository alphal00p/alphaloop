#!/usr/bin/env python2

import sys

import os
pjoin = os.path.join

root_path = os.path.dirname(os.path.realpath( __file__ ))
sys.path.insert(0, pjoin(root_path,os.pardir,os.pardir))
sys.path.insert(0, pjoin(root_path,'rust_backend'))

import math
import time
import random
import yaml
from yaml import Loader, Dumper
noalias_dumper = Dumper
noalias_dumper.ignore_aliases = lambda self, data: True

from pprint import pprint
import itertools

#import matplotlib
#matplotlib.use('W')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

import sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath( __file__ )),os.path.pardir))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath( __file__ )),os.path.pardir,os.path.pardir))

from ltd_utils import LoopTopology, LoopLine, Propagator, TopologyCollection
from topologies import create_hard_coded_topology
import ltd_commons 
hyperparameters = ltd_commons.hyperparameters
topology_collection = ltd_commons.hard_coded_topology_collection
import vectors
import diagnostics

# Import below necessary for plotting the elliptic surface
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
import numpy as np

try:
    # Import the rust bindings
    from ltd import LTD
except ImportError:
    raise BaseException(
        "Could not import the rust back-end 'ltd' module. Compile it first with:"
        " ./make_lib from within the pyNLoop directory." )

def map_to_infinity(x):
    return ((1./(1.-x)) - 1./x)

def map_from_infinity(x):
    return -2./(-2.+x-math.sqrt(4.+x**2))

class DualCancellationScanner(object):

    def __init__(self,
                 log_stream         =   None,
                 topology           =   None,
                 direction_vectors  =   None,
                 parametrisation_vectors = None,
                 parametrisation_uv =   None,
                 t_values           =   None,
                 hyperparameters    =   None,
                 rust_instance      =   None,
                 rotated_instances  =   None,
                 diagnostic_tool    =   None,
                ):

        self.topology = topology 
        self.direction_vectors = direction_vectors
        self.parametrisation_vectors = parametrisation_vectors
        self.parametrisation_uv = parametrisation_uv
        self.t_values = t_values

        # Now save the default yaml configuration to the log stream if specified.
        self.configuration = {
            'topology': self.topology.name,
            'direction_vectors': [[float(v) for v in vec] for vec in self.direction_vectors],
            'parametrisation_vectors': [[float(v) for v in vec] for vec in self.parametrisation_vectors],
            'parametrisation_uv' : parametrisation_uv,
            'hyperparameters': hyperparameters,
        }

        self.log_stream = log_stream
        if self.log_stream is not None:
            log_stream.write(yaml.dump([('configuration', self.configuration), ], Dumper=noalias_dumper))

        self.rust_instance = rust_instance
        self.rotated_instances = rotated_instances
        self.diagnostic_tool = diagnostic_tool

    def scan(self, surface_ID, surface_characteristics, force=False):
        all_results = []

        n_cut = self.diagnostic_tool.cut.index(list(surface_characteristics[0]))
        
        # Choose a point to approach
        parametrisation_vectors = [ [ 0., ] + [v for v in vec]  for vec in self.parametrisation_vectors ]
        # u and v are defined in [-1.,1.]
        #u = self.parametrisation_uv['u']*2.-1.
        #v = self.parametrisation_uv['v']*2.-1.
        u = self.parametrisation_uv['u']
        v = self.parametrisation_uv['v']

        loop_mom = list(parametrisation_vectors)
        if len(loop_mom)==1:
            loop_mom.append([0.,0.,0.,0.])
        diagnostic_point = self.diagnostic_tool.get_parametrization(
            u=u, v=v, 
            loop_momenta=loop_mom,
            surface=surface_characteristics[1], n_cut=n_cut)
        
        if diagnostic_point is None:
            if not force:
                # Failed to construct a point for this surface, let's return immediately
                if self.log_stream is not None:
                    log_stream.write(yaml.dump([ ('surface_results', surface_ID, None), ], Dumper=noalias_dumper))

                return {
                    'configuration': self.configuration,
                    'surface'      : surface_ID,
                    'scan_results' : None
                }
            else:
                com_energy = topology.get_com_energy()
                while diagnostic_point is None:
                    parametrisation_vectors = [ [ 0., ] + list(vectors.Vector([(random.random() - 0.5) * 2 * com_energy for _ in range(3)]))
                                               for vec in parametrisation_vectors]
                    loop_mom = list(parametrisation_vectors)
                    if len(loop_mom)==1:
                        loop_mom.append(vectors.Vector([0.,0.,0.,0.]))
                    diagnostic_point = self.diagnostic_tool.get_parametrization(
                            u=u, v=v, 
                            loop_momenta=loop_mom,
                            surface=surface_characteristics[1], n_cut=n_cut)
        
        surface_point = [vectors.Vector(list(v)) for v in diagnostic_point]

        for t in self.t_values:

            one_result = {
                't_value': t,
            }
            approach_point = [ v*(1.+t) for v in surface_point ]
            one_result['parametric_equation'] = surface_characteristics[1].parametric_equation_eval(
                self.diagnostic_tool.propagators,
                approach_point,
                self.topology.n_loops
            ) 

            for i_rot, rotated_info in enumerate([None,]+self.rotated_instances):
                if rotated_info is not None:
                    point_to_test = [v.transform(rotated_info['rotation_matrix']) for v in approach_point]
                    one_result['rotation_%d'%i_rot] = {}
                    result_bucket = one_result['rotation_%d'%i_rot]
                    rust_worker = rotated_info['rust_instance']
                else:
                    point_to_test = approach_point
                    result_bucket = one_result
                    rust_worker = self.rust_instance
                
                print('i_rot=%d'%i_rot)
                if rotated_info:
                    print("rotation_matrix=%s"%str(rotated_info['rotation_matrix']))
                else:
                    print("rotation_matrix=identity")
                print('point=%s,%s'%(','.join('%.16f'%ve for ve in point_to_test[0]), ','.join('%.16f'%ve for ve in point_to_test[1]) ) )
                kappas, jac_re, jac_im = rust_worker.deform([list(v) for v in point_to_test])
                deformed_point = [point_to_test[i]+vectors.Vector(kappas[i])*complex(0.,1.) for i in range(len(kappas))]
            
                result_bucket['deformation_jacobian'] = complex(jac_re, jac_im)
                dual_integrands = {}
                point_weight = complex(0., 0.)
                for ltd_cut_index, ltd_cut_structure in enumerate(self.topology.ltd_cut_structure):
                    cut_propagator_indices = [[index*c for index in range(1,len(self.topology.loop_lines[i].propagators)+1)] if c!=0 else [0,] 
                                             for i, c in enumerate(ltd_cut_structure)]
                    dual_integrands[(ltd_cut_index, ltd_cut_structure)] = {}
                    for cut_index, cut_structure  in enumerate(itertools.product( *cut_propagator_indices )):
                        dual_integrand_real, dual_integrand_imag = rust_worker.evaluate_cut(
                            [ [ ( float(vi.real), float(vi.imag) ) for vi in v ] for v in deformed_point],
                            ltd_cut_index,  cut_index
                        )
                        if math.isnan(dual_integrand_real) or math.isnan(dual_integrand_imag):
                            print("WARNINNG : Rust returned NaN for the following dual integrand:")
                            print "ltd_cut_index, ltd_cut_structure = ",ltd_cut_index, ltd_cut_structure
                            print "cut_index, cut_structure         = ",cut_index, cut_structure
                            print "with input momenta = ",deformed_point
                            print("when approaching this surface:\n%s"%str(surface[1].__str__(cut_propagators=surface[0])))
                            print("Results now set to zero instead")
                            dual_integrand_real = 0.
                            dual_integrand_imag = 0.
                            sys.exit(1)

                        dual_integrands[(ltd_cut_index, ltd_cut_structure)][(cut_index, cut_structure)] = \
                                                                    complex(dual_integrand_real,dual_integrand_imag)
                        point_weight += complex(dual_integrand_real,dual_integrand_imag)

                result_bucket['dual_integrands'] = dual_integrands

                CT                                  = complex(0., 0.)
                result_bucket['ct']                 = CT
                result_bucket['point_weight_no_ct'] = point_weight
                result_bucket['point_weight']       = point_weight + CT

            # Regiter numerical stability as well
            target = one_result['point_weight']
            evaluations = []
            for k,v in one_result.items():
                if k.startswith('rotation_'):
                    evaluations.append(v['point_weight'])
            
            if len(evaluations)>0:
                one_result['accuracy'] = max(
                    max(abs(e.real-target.real) for e in evaluations),
                    max(abs(e.imag-target.imag) for e in evaluations)
                )
            else:
                one_result['accuracy'] = None

        ###     'point_weight': <point_weight_complex_value>,
        ###     'point_weight_no_ct': <point_weight_no_ct_complex_value>,
        ###     'ct' : <ct_value>,


            all_results.append(one_result)

        # Save the result into the logstream if specified
        if self.log_stream is not None:
            log_stream.write(yaml.dump([ ('surface_results', surface_ID, all_results), ], Dumper=noalias_dumper)) 

        return {
            'configuration': self.configuration,
            'surface'      : surface_ID,
            'scan_results' : all_results
        }

class PoleScanner(object):
    def __init__(self,
                 log_stream         =   None,
                 topology           =   None,
                 direction_vectors  =   None,
                 parametrisation_vectors = None,
                 parametrisation_uv =   None,
                 t_values           =   None,
                 hyperparameters    =   None,
                 rust_instance      =   None,
                 rotated_instances  =   None,
                 diagnostic_tool    =   None,
                 plot_surface       =   False,
                ):

        self.topology = topology 
        self.direction_vectors = direction_vectors
        self.parametrisation_vectors = parametrisation_vectors
        self.parametrisation_uv = parametrisation_uv
        self.t_values = t_values
        self.do_plot_surface = plot_surface
        if self.do_plot_surface:
            self.fig = plt.figure()
            self.ax = self.fig.add_subplot(111, projection='3d')
            self.ax.set_title('3D Surfaces')


        # Now save the default yaml configuration to the log stream if specified.
        self.configuration = {
            'topology': self.topology.name,
            'direction_vectors': [[float(v) for v in vec] for vec in self.direction_vectors],
            'parametrisation_vectors': [[float(v) for v in vec] for vec in self.parametrisation_vectors],
            'parametrisation_uv' : parametrisation_uv,
            'hyperparameters': hyperparameters,
        }

        self.log_stream = log_stream
        if self.log_stream is not None:
            log_stream.write(yaml.dump([('configuration', self.configuration), ], Dumper=noalias_dumper))

        self.rust_instance = rust_instance
        self.rotated_instances = rotated_instances
        self.diagnostic_tool = diagnostic_tool

        # Build a map between Zeno's conventions for propagators and the (loop_line_id, propagator_id) tuple one.
        self.propagators_map = {}
        self.inverse_propagators_map = {}
        propagator_id = 0
        for i_ll, ll in enumerate(self.topology.loop_lines):
            for i_prop, prop in enumerate(ll.propagators):
                self.propagators_map[propagator_id] = ((i_ll, i_prop), prop)
                self.inverse_propagators_map[(i_ll, i_prop)] = propagator_id
                propagator_id += 1

        # Build a map between Zeno's conventions for cuts and the (ltd_cut_struct_index, cut_index) tuple one.
        self.cut_indices_map = {}
        for ltd_cut_index, ltd_cut_structure in enumerate(self.topology.ltd_cut_structure):
            cut_propagator_indices = [[index*c for index in range(1,len(self.topology.loop_lines[i].propagators)+1)] if c!=0 else [0,] 
                                             for i, c in enumerate(ltd_cut_structure)]
            for cut_index, cut_structure  in enumerate(itertools.product( *cut_propagator_indices )):
                self.cut_indices_map[
                    tuple(sorted(list(set([self.inverse_propagators_map[(i_ll, abs(i_prop)-1)] for i_ll, i_prop in enumerate(cut_structure) if i_prop!=0 ]))))
                                                    ] = (ltd_cut_index, cut_index)

        # Turning off accuracy assessment below
        self.rotated_instances = []

    def plot_surface(self, surface, surface_ID, n_cut, loop_mom, ltd_cut_index, cut_index, onshell_propagator, onshell_propagator_id):
        """ Plot the elliptic surface specified, coloring it with the imaginary port of the corresponding onshell prop."""

        #first generate a linspace of u and v
        u_values = np.linspace(0.0001, 0.9999, len(self.t_values))
        v_values = np.linspace(0.0001, 0.9999, len(self.t_values))
        x, y, z, all_os_prop_evals = [], [], [], []
        for sheet in [0,1]:
            surface.sheet = sheet
            for u in u_values:
                xs, ys, zs, os_prop_evals = [], [], [], []
                for v in v_values:
                    ellipse_point = self.diagnostic_tool.get_parametrization(
                        u=u, v=v, 
                        loop_momenta=loop_mom,
                        surface=surface,
                        n_cut=n_cut
                    )
                    if ellipse_point is None:
                        xs.append(0.0)
                        ys.append(0.0)
                        zs.append(0.0)
                        os_prop_evals.append(-1.0)
                        continue
                    surface_point = ellipse_point[surface.param_variable]
                    xs.append(surface_point[0])
                    ys.append(surface_point[1])
                    zs.append(surface_point[2])

                    kappas, jac_re, jac_im = self.rust_instance.deform([list(v) for v in ellipse_point])
                    deformed_point = [ellipse_point[i]+vectors.Vector(kappas[i])*complex(0.,1.) for i in range(len(kappas))]

                    energies = self.rust_instance.get_loop_momentum_energies([ 
                            [(v_el.real, v_el.imag) for v_el in v] for v in deformed_point] ,ltd_cut_index, cut_index)

                    deformed_loop_four_momenta = []
                    for i_v, v in enumerate(deformed_point):
                        deformed_loop_four_momenta.append(vectors.LorentzVector(
                            [complex(energies[i_v][0],energies[i_v][1]),]+list(v)
                        ))

                    os_prop_evals.append(onshell_propagator.evaluate_inverse(deformed_loop_four_momenta).imag)

                if len(xs)>0:
                    x.append(xs)
                    y.append(ys)
                    z.append(zs)
                    all_os_prop_evals.append(os_prop_evals)
        
        # Convert to np.array
        x = np.array([np.array(x_i) for x_i in x])
        y = np.array([np.array(y_i) for y_i in y])
        z = np.array([np.array(z_i) for z_i in z])
        all_os_prop_evals = np.array([np.array(c_i) for c_i in all_os_prop_evals])
        
        color_specifier = np.array([ np.array([-1. if im > 0. else 1. for im in os_prop_evals]) for os_prop_evals in all_os_prop_evals])
        norm = Normalize()
        colors = norm(color_specifier)
        cmap = cm.get_cmap("coolwarm")

        #self.ax.set_title('Surface #%d with onshell prop %s'%(surface_ID,str(onshell_propagator_id)))        
        self.ax.plot_surface(x, y, z, linewidth=1, facecolors=cmap(colors), shade=True, alpha=0.35)
        

    def scan(self, surface_ID, surface_characteristics):
        
        all_results = []
        
        n_cut = self.diagnostic_tool.cut.index(list(surface_characteristics[0]))
        
        ltd_cut_index, cut_index = self.cut_indices_map[tuple(sorted(list(set([prop_index for prop_index in surface_characteristics[0]]))))]

        # Find the point for this t-value
        parametrisation_vectors = [ [ 0., ] + [v for v in vec]  for vec in self.parametrisation_vectors ]

        # u and v are defined in [-1.,1.]
        #u_base = self.parametrisation_uv['u']*2.-1.
        #v_base = self.parametrisation_uv['v']*2.-1.
        u_base = self.parametrisation_uv['u']
        v_base = self.parametrisation_uv['v']

        direction_u = self.direction_vectors[0][0]
        direction_v = self.direction_vectors[0][1]

        # Identify which is the onshell propagator
        onshell_propagator = self.propagators_map[surface_characteristics[1].n_surface][1] 
        onshell_propagator_id = self.propagators_map[surface_characteristics[1].n_surface][0]

        plot_surface_already_done = False
        for t in self.t_values:

            one_result = {
                't_value': t,
            }

            # Now build a point on the elliptic surface for this particular t-value
            loop_mom = list(parametrisation_vectors)
            if len(loop_mom)==1:
                loop_mom.append([0.,0.,0.,0.])
           
            u = u_base + direction_u*t
            v = v_base + direction_v*t
            if u < 0.:
                u = u%(-1.)
            else:
                u = u%(1.)
            if v < 0.:
                v = v%(-1.)
            else:
                v = v%(1.)
            
            diagnostic_point = self.diagnostic_tool.get_parametrization(
                u=u, v=v, 
                loop_momenta=loop_mom,
                surface=surface_characteristics[1],
                n_cut=n_cut
            )
            if diagnostic_point is None:
                #print("Point for surface #%d could not be generated."%surface_ID)
                one_result['accuracy'] = 0.
                one_result['loop_momenta'] = []
                one_result['generating_xs'] = []
                one_result['propagator_evaluation'] = complex(0., 0.)
                one_result['parametrisation_jacobian'] = 0.0
                one_result['deformation_jacobian'] = complex(0., 0.)
                all_results.append(one_result)
                continue
            
            surface_point = [vectors.Vector(list(v)) for i, v in enumerate(diagnostic_point) if i<self.topology.n_loops]
            
            if not plot_surface_already_done and self.do_plot_surface:
                self.plot_surface(surface_characteristics[1], surface_ID, n_cut, loop_mom, 
                                  ltd_cut_index, cut_index, onshell_propagator, onshell_propagator_id)
                plot_surface_already_done = True

            for i_rot, rotated_info in enumerate([None,]+self.rotated_instances):
                if rotated_info is not None:
                    point_to_test = [v.transform(rotated_info['rotation_matrix']) for v in surface_point]
                    one_result['rotation_%d'%i_rot] = {}
                    result_bucket = one_result['rotation_%d'%i_rot]
                    rust_worker = rotated_info['rust_instance']
                else:
                    point_to_test = surface_point
                    result_bucket = one_result
                    rust_worker = self.rust_instance
                
                if any(math.isnan(v_elem) for v in point_to_test for v_elem in v):
                    result_bucket['loop_momenta'] = []
                    result_bucket['generating_xs'] = []
                    result_bucket['parametrisation_jacobian'] = 0.0
                    result_bucket['deformation_jacobian'] = complex(0.0,0.0)
                    result_bucket['propagator_evaluation'] = complex(0.0,0.0)
                    continue

                result_bucket['loop_momenta'] = [vectors.Vector(v) for v in point_to_test]

                generating_variables = []
                parametrisation_jacobian = 1.0
                for i, v in enumerate(point_to_test):
                    res = rust_worker.inv_parameterize(list(v), i)
                    generating_variables.append(tuple(res[:3]))
                    parametrisation_jacobian *= res[3]
                result_bucket['generating_xs'] = tuple(generating_variables)

                # Now evaluate the energy component using the real momenta
                real_energies = rust_worker.get_loop_momentum_energies(
                        [[(v_elem,0.) for v_elem in v] for v in point_to_test], ltd_cut_index, cut_index)

                loop_four_momenta = []
                for i_v, v in enumerate(point_to_test):
                    assert(real_energies[i_v][1]==0.)
                    loop_four_momenta.append(vectors.LorentzVector(
                        [real_energies[i_v][0],]+list(v)
                    ))
                # Cross-check: verify that the onshell propagator is indeed onshell with real momenta
                os_prop_eval_with_real_momenta = onshell_propagator.evaluate_inverse(loop_four_momenta)
                #assert(os_prop_eval_with_real_momenta/(max(max(v) for v in loop_four_momenta)) < 1.0e-10)

                kappas, jac_re, jac_im = rust_worker.deform([list(v) for v in point_to_test])
                deformed_point = [point_to_test[i]+vectors.Vector(kappas[i])*complex(0.,1.) for i in range(len(kappas))]
                result_bucket['deformation_jacobian'] = complex(jac_re, jac_im)
                
                # Now evaluate the energy component using the imaginary momenta
                energies = rust_worker.get_loop_momentum_energies([ [(v_el.real, v_el.imag) for v_el in v] for v in deformed_point] ,ltd_cut_index, cut_index)

                deformed_loop_four_momenta = []
                for i_v, v in enumerate(deformed_point):
                    deformed_loop_four_momenta.append(vectors.LorentzVector(
                        [complex(energies[i_v][0],energies[i_v][1]),]+list(v)
                    ))

                result_bucket['propagator_evaluation'] = onshell_propagator.evaluate_inverse(deformed_loop_four_momenta)

            # Register numerical stability as well
            target = one_result['propagator_evaluation']
            evaluations = []
            for k,v in one_result.items():
                if k.startswith('rotation_'):
                    evaluations.append(v['propagator_evaluation'])
            
            if len(evaluations)>0:
                one_result['accuracy'] = max(
                    max(abs(e.real-target.real) for e in evaluations),
                    max(abs(e.imag-target.imag) for e in evaluations)
                )
            else:
                one_result['accuracy'] = None


            all_results.append(one_result)

        # Save the result into the logstream if specified
        if self.log_stream is not None:
            log_stream.write(yaml.dump([ ('surface_results', surface_ID, onshell_propagator_id, all_results), ], Dumper=noalias_dumper)) 
        
        return {
            'configuration': self.configuration,
            'surface'      : surface_ID,
            'onshell_prop_id' : onshell_propagator_id,
            'scan_results' : all_results
        }

class DualCancellationResultsAnalyser(object):
    def __init__(self, surfaces_repository, results):
        # Each entry of the surface_results has the following format:
        ###   [{'surface': surface_instance_approached, 'scan_results' : scan_results},]
        ###   with scan_results having the following format
        ###    {
        ###     't_value' : <parameter_t_float_value>,
        ###     'deformation_jacobian': <deformation_jacobian_complex_value>,
        ###     'dual_integrands': dictionary with keys (ltd_cut_index, cut_index),
        ###     'point_weight': <point_weight_complex_value>,
        ###     'point_weight_no_ct': <point_weight_no_ct_complex_value>,
        ###     'ct' : <ct_value>,
        ###     'accuracy': <float_accuracy>,
        ###    }
        ### And possibly copies of some of the above for each rotated PS point.
        self.surfaces = surfaces_repository
        self.surface_results   = results['surface_results']
        self.configuration  = results['configuration']
        self.direction_vectors = [ vectors.LorentzVector(lv) for lv in self.configuration['direction_vectors'] ]
        self.parametrisation_vectors    = [ vectors.LorentzVector(ov) for ov in self.configuration['parametrisation_vectors'] ]
        self.run_time       = results['run_time']


    def generate_plots(self, result, show_real=True, show_imag=False, normalise=True,
                       entries_to_plot = None,
                       log_y_scale = True,
                       log_x_scale = False):

        surface_ID = result['surface']
        cut_propagators = self.surfaces[surface_ID][0]
        surface = self.surfaces[surface_ID][1]
        scan_results = result['scan_results']
        first_result = scan_results[0]
        plt.title('Test approaching surface with id #%d'%surface_ID)
        plt.ylabel('Normalised weight')
        plt.xlabel('t')
        if log_y_scale:
            plt.yscale('log')
        if log_x_scale:
            plt.xscale('log')

        # Generate several lines to plot
        lines = []
        def special_abs(x):
            if log_y_scale:
                return abs(x)
            else:
                return x

        for ltd_cut_index, ltd_cut_structure in sorted(first_result['dual_integrands'].keys()):
            duals_for_ltd_cut = first_result['dual_integrands'][(ltd_cut_index, ltd_cut_structure)]
            for cut_index, cut_structure in sorted(duals_for_ltd_cut.keys()):
                dual_res = duals_for_ltd_cut[(cut_index, cut_structure)]
                
                dual_name = 'dual int. %-15s, %-15s' % ( 
                        '%d@%s'%(ltd_cut_index, '(%s)'%(','.join( '%s%d'%('+' if c>0 else 
                            (' ' if c==0 else ''),c) for c in ltd_cut_structure))), 
                        '%d@%s'%(cut_index, '(%s)'%(','.join( '%s%d'%('+' if c>0 else 
                            (' ' if c==0 else ''),c) for c in cut_structure)))
                )

        
                if show_real:
                    lines.append( ( 'Real %s'%dual_name , (
                        [ v['t_value'] for v in scan_results], 
                        [ special_abs(v['dual_integrands'][(ltd_cut_index, ltd_cut_structure)][(cut_index, cut_structure)].real)
                            for v in scan_results ] 
                    ) ) )
                if show_imag:
                    lines.append( ( 'Cmplx %s'%dual_name , (
                        [ v['t_value'] for v in scan_results],
                        [ special_abs(v['dual_integrands'][(ltd_cut_index, ltd_cut_structure)][(cut_index, cut_structure)].imag)
                            for v in scan_results ]
                    ) ) )

        if show_real:
            lines.append( ( 'Real weight', ( [ v['t_value'] for v in scan_results], 
                                             [ special_abs(v['point_weight'].real) for v in scan_results ] ) ) )
        if show_imag:
            lines.append( ( 'Cmplx weight', ( [ v['t_value'] for v in scan_results], 
                                              [ special_abs(v['point_weight'].imag) for v in scan_results ] ) ) )

        if scan_results[0]['accuracy'] is not None:
            lines.append( ( 'Accuracy', ( [ v['t_value'] for v in scan_results], 
                                             [ special_abs(v['accuracy']) for v in scan_results ] ) ) )

        lines.append( ( 'parametric_equation', (
            [ v['t_value'] for v in scan_results], [ special_abs(v['parametric_equation']) for v in scan_results ]
        ) ) )

        filtered_lines = []
        for ln, lv in lines:
            if all(lvel==0. for lvel in lv[1]):
                print('WARNING: All x-values for line "%s" are zero.'%ln)
            else:
                filtered_lines.append((ln, lv))
        lines = filtered_lines

        for line_name, (x_data, y_data) in lines:
            if normalise:
                max_y_data = max(y_data)
                if max_y_data > 0.:
                    y_data = [ y/float(max_y_data) for y in y_data ]

        for line_name, (x_data, y_data) in lines:
            plt.plot(x_data, y_data, label=line_name)
    
        plt.legend(bbox_to_anchor=(0.75, 0.5))
        plt.show()

    def analyse(self, **opts):

        print('='*80)
        print("run time: %.3e [h]"%(self.run_time/3600.0))
        for result in self.surface_results:
            print('-'*80)            
            surface_ID = result['surface']
            cut_propagators = self.surfaces[surface_ID][0]
            surface = self.surfaces[surface_ID][1]
            print('Results for the following surface #%d:\n%s'%(surface_ID,
                surface.__str__(cut_propagators=cut_propagators)))
            scan_results = result['scan_results']
            if scan_results is None:
                print('Surface did not exist for specified parametrisation vector.')
                continue
            print('Last 5 points of scan:')
            for r in scan_results[-5:]:
                print('>>> t_value = %-16e'%r['t_value'])
                print('Weight: (%-25.16e, %-25.16e) +/- %-25.16e'%(
                    r['point_weight'].real,r['point_weight'].imag,r['accuracy']))
                for ltd_cut_index, ltd_cut_structure in sorted(r['dual_integrands'].keys()):
                    duals_for_ltd_cut = r['dual_integrands'][(ltd_cut_index, ltd_cut_structure)]
                    for cut_index, cut_structure in sorted(duals_for_ltd_cut.keys()):
                        dual_res = duals_for_ltd_cut[(cut_index, cut_structure)]
                        print('Dual integrand %-15s, %-15s : (%-25s, %-25s) '%(
                            '%d@%s'%(ltd_cut_index, '(%s)'%(','.join( '%s%d'%('+' if c>0 else 
                                (' ' if c==0 else ''),c) for c in ltd_cut_structure))), 
                            '%d@%s'%(cut_index, '(%s)'%(','.join( '%s%d'%('+' if c>0 else 
                                (' ' if c==0 else ''),c) for c in cut_structure))), 
                            '%s%.16e' % ('+' if dual_res.real>0. else '', dual_res.real) if dual_res.real!=0. else ' 0.', 
                            '%s%.16e' % ('+' if dual_res.imag>0. else '', dual_res.imag) if dual_res.imag!=0. else ' 0.'))
        print('='*80)
        for result in self.surface_results:
            if result['scan_results'] is None:
                continue
            self.generate_plots(result, **opts)

class PoleResultsAnalyser(object): 
    def __init__(self, surfaces_repository, results):
        # Each entry of the surface_results has the following format:
        ###   [{'surface': surface_instance_approached, 'scan_results' : scan_results},]
        ###   with scan_results having the following format
        ###    {
        ###     't_value' : <parameter_t_float_value>,
        ###     'deformation_jacobian': <deformation_jacobian_complex_value>,
        ###     'propagator_evaluation': <complex_value>,
        ###     'accuracy': <float_accuracy>
        ###    }
        ### And possibly copies of some of the above for each rotated PS point.
        self.surfaces = surfaces_repository
        self.surface_results   = results['surface_results']
        self.configuration  = results['configuration']
        self.direction_vectors = [ vectors.LorentzVector(lv) for lv in self.configuration['direction_vectors'] ]
        self.parametrisation_vectors    = [ vectors.LorentzVector(ov) for ov in self.configuration['parametrisation_vectors'] ]
        self.run_time       = results['run_time']

    def test_imaginary_parts_in_scan(self):
        
        found_problem = False
        for result in self.surface_results:
            scan_results = result['scan_results'] 
            if scan_results is None:
                continue
            surface_ID = result['surface']                                    
            cut_propagators = self.surfaces[surface_ID][0]
            surface = self.surfaces[surface_ID][1]
            for r in scan_results:
                if r['propagator_evaluation'].imag < 0.:
                    found_problem += 1
                    print("==================================================================")                    
                    print(">>> INCORRECT ANALYTICAL CONTINUATION FOR THE FOLLOWING SETUP: <<<")
                    print("==================================================================")                    
                    print('Results for the following surface #%d:\n%s'%(surface_ID,
                        surface.__str__(cut_propagators=cut_propagators)))
                    print('Featuring the following onshell propagator: %s'%str(result['onshell_prop_id']))
                    print('>>> t_value = %.16e'%r['t_value'])
                    print('>>> loop_momenta = %s'%r['loop_momenta'])
                    print('Onshell propagator evaluation *violating causal prescription*:\n     %s%-25.16e %sI*%-25.16e'%(
                        '+' if r['propagator_evaluation'].real >= 0. else '-',
                        abs(r['propagator_evaluation'].real),
                        '+' if r['propagator_evaluation'].imag >= 0. else '-',
                        abs(r['propagator_evaluation'].imag),)
                    )
                    found_problem = True
                    return found_problem
            
        return found_problem

    def generate_plots(self, result, show_real=True, show_imag=False, normalise=True,
                       entries_to_plot = None,
                       log_y_scale = True,
                       log_x_scale = False):

        surface_ID = result['surface']
        scan_results = result['scan_results']
        first_result = scan_results[0]
        plt.title('Pole test for surface with id #%d (onshel prop: %s)'%
                        (surface_ID, str(result['onshell_prop_id']) ))
        plt.ylabel('Normalised weight')
        plt.xlabel('t')
        if log_y_scale:
            plt.yscale('log')
        if log_x_scale:
            plt.xscale('log')

        # Generate several lines to plot
        lines = []
        def special_abs(x):
            if log_y_scale:
                return abs(x)
            else:
                return x

        #if show_real:
        #    lines.append( ( 'Real eval prop.', ( [ v['t_value'] for v in scan_results], 
        #                                     [ special_abs(v['propagator_evaluation'].real) for v in scan_results ] ) ) )
        
        lines.append( ( 'Cmplx eval prop.', ( [ v['t_value'] for v in scan_results], 
                                              [ v['propagator_evaluation'].imag for v in scan_results ] ) ) )

        #if scan_results[0]['accuracy'] is not None:
        #    lines.append( ( 'Accuracy', ( [ v['t_value'] for v in scan_results], 
        #                                     [ special_abs(v['accuracy']) for v in scan_results ] ) ) )

        filtered_lines = []
        for ln, lv in lines:
            if all(lvel==0. for lvel in lv[1]):
                print('WARNING: All x-values for line "%s" are zero.'%ln)
            else:
                filtered_lines.append((ln, lv))
        lines = filtered_lines

        for line_name, (x_data, y_data) in lines:
            if normalise:
                max_y_data = max(y_data)
                if max_y_data > 0.:
                    y_data = [ y/float(max_y_data) for y in y_data ]

            plt.plot(x_data, y_data, label=line_name)
    
        plt.legend(bbox_to_anchor=(0.75, 0.5))
        plt.show()

    def analyse(self, n_points_to_show=5, **opts):

        print('='*80)
        print("run time: %.3e [h]"%(self.run_time/3600.0))
        for result in self.surface_results:
            print('-'*80)            
            surface_ID = result['surface']
            cut_propagators = self.surfaces[surface_ID][0]
            surface = self.surfaces[surface_ID][1]
            print('Results for the following surface #%d:\n%s'%(surface_ID,
                surface.__str__(cut_propagators=cut_propagators)))
            print('Featuring the following onshell propagator: %s'%str(result['onshell_prop_id']))
            scan_results = result['scan_results']
            if scan_results is None:
                print('Surface did not exist for specified parametrisation vector.')
                continue
            print('Last %d points of scan:'%n_points_to_show)
            for r in scan_results[-n_points_to_show:]:
                print('>>> t_value = %-16e'%r['t_value'])
                print('>>> loop_momenta = %s'%r['loop_momenta'])
                print('>>> generating_variables = %s'%r['generating_xs'])
                print('Onshell propagator evaluation: %s%-25.16e %sI*%-25.16e'%(
                    '+' if r['propagator_evaluation'].real >= 0. else '-',
                    abs(r['propagator_evaluation'].real),
                    '+' if r['propagator_evaluation'].imag >= 0. else '-',
                    abs(r['propagator_evaluation'].imag),)
                )
        print('='*80)
        for result in self.surface_results:
            if result['scan_results'] is None:
                continue
            self.generate_plots(result, **opts)

def load_results_from_yaml(log_file_path):
    """Load a full-fledged scan from a yaml dump"""

    raw_data = yaml.load(open(log_file_path,'r'), Loader=Loader)
    processed_data = {'surface_results': []}
    for entry in raw_data:
        entry_name = entry[0]
        if entry_name == 'surface_results':
            if len(entry)==3:
                processed_data['surface_results'].append(
                    { 'surface': entry[1],
                       'onshell_prop_id': entry[2],
                       'scan_results' : entry[3]
                    }
                )
            else:
                processed_data['surface_results'].append(
                    { 'surface': entry[1],
                       'scan_results' : entry[2]
                    }
                )
            continue
        processed_data[entry_name] = entry[1]
    
    return processed_data

if __name__ == '__main__':

    prog_options = ' '.join(sys.argv[1:]).split('--')

    topology            = 'DoubleTriangle'
    test_name           = 'cancellation_check'
    _supported_tests    = ['cancellation_check','pole_check']
    # Default values assume a two-loop topology
    # None means chosen at random
    direction_vectors   = [None,None]
    parametrisation_vectors      = [None,None]
    t_n_points          = 100
    t_range             = [0.,1.0]
    t_values_steps      = 'linear'
    # If None, it will be assigned to a default later on
    t_values            = None
    fixed_t_values      = [0.5, 0.5]
    # Specify None to run over all surfaces or a list of integer to select particular one. Specify [-1] to print all surfaces.
    surface_ids         = None
    # Offers the possibility of looping 'do_loop' times (or at infinity when set to -1) the test with different random
    # direction and offset vectors each time.
    do_loop             = None
    # For now the selected PS Point will have to be specified below
    PS_point            = None 
    # If specified to None, then a random PS Point is generated.
    #PS_point            = None #TOIMPLEMENT
    load_results_from   = None
    # Turn off the logging by default as it slows us down unnecessarily
    save_results_to     = None
    # save_results_to     = pjoin(os.getcwd(),'poles_imaginary_part_scan.yaml')
    plot_surface        = False

    # Number of points to be used for testing hyperboloid existence
    n_points_hyperboloid_test = 1000


    rotations_for_stability_check = [((1,0.3),(2,0.9),(3,1.7)),((3,1.2),)]

    # parametrisation_uv
    parametrisation_uv = { 'u' : None, 'v' : None }

    # Whether to show the imaginary and/or real part of the weights in the limit test
    analysis_options = {}

    for prog_option in prog_options:
        if prog_option.strip()=='':
            continue
        try:
            option_name, option_value = prog_option.strip().split('=')
        except ValueError:
            option_name, option_value = prog_option, None

        if option_name in ['name','n']:
            if option_value.lower == 'c':
                option_value = 'cancellation_check'
            elif option_value.lower == 'p':
                option_value = 'pole_check'
            if option_value.lower() not in _supported_tests:
                raise BaseException("Test name '%s' not recognized. Must be in %s"%(
                                              option_value.lower(), _supported_tests))
            test_name = option_value.lower()
        elif option_name in ['topology','t']:
            topology = option_value
        elif option_name in ['seed']:
            random.seed(eval(option_value))
        elif option_name in ['PS_point', 'ps']:
            option_eval = eval(option_value)
            PS_point = vectors.LorentzVectorList()
            for v in option_eval:
                PS_point.append(vectors.LorentzVector(v))
        elif option_name in ['direction_vectors','d']:
            # A vector set to None means it will be randomly chosen
            direction_vectors = [ (vectors.Vector(vec) if vec is not None else vec) for vec in eval(option_value) ]
        elif option_name in ['parametrisation_vectors','p']:
            # A vector set to None means it will be randomly chosen
            parametrisation_vectors = [ (vectors.Vector(vec) if vec is not None else vec) for vec in eval(option_value) ]
        elif option_name in ['t_values','tv']:
            # The ranges must be iterable of floats or None to select the default iterable.
            t_values = eval(option_value)
        elif option_name in ['t_values_steps','tvs']:
            if option_value not in ['linear','exponential']:
                raise BaseException("The only allowed values of the t values steps are 'linear' or 'exponential'.")
            t_values_steps = option_value
        elif option_name in ['t_range','tr']:
            t_range = eval(option_value)
        elif option_name in ['t_n_points','tn']:
            t_n_points = eval(option_value)
        elif option_name in ['fixed_t_values','ft']:
            # The ranges must be iterable of floats or None to select the default iterable.
            fixed_t_values = eval(option_value)
            if fixed_t_values is None:
                if len(direction_vectors)>0:
                    fixed_t_values = [0.5,]*len(direction_vectors)
                else:
                    # Assume two-loop momenta, as in the default fixed_t_values value set above
                    fixed_t_values = [0.5,0.5]
        elif option_name in ['surface_ids','s']:
            # The surface option must be an iterable yielding integers or [] to select them all, and [-1] to just list all surfaces.
            s_ids = eval(option_value)
            if len(s_ids)==0:
                surface_ids = None
            else:
                surface_ids = eval(option_value)
        elif option_name in ['load_results_from','load']:
            load_results_from = option_value
        elif option_name in ['plot_surface', 'psrfc']:
            if option_value is None:
                plot_surface = True
            else:
                plot_surface = option_name.lower() in ['on','true','t']
        elif option_name in ['save_results_to','save']:
            save_results_to = option_value
        elif option_name in ['n_points_hyperboloid_test','nh']:
            n_points_hyperboloid_test = int(option_value)
        elif option_name in ['loop','l']:
            do_loop = eval(option_value)
        elif option_name in ['analysis_options','ao']:
            analysis_options = eval(option_value)
        elif option_name in ['u','v']:
            parametrisation_uv[option_name] = float(option_value)
        elif option_name in ['uv']:
            parametrisation_uv['u'] = eval(option_value)[0]
            parametrisation_uv['v'] = eval(option_value)[1]
        else:
            raise BaseException("Option '%s' not reckognised."%option_name)
    
    # Adjust kinematics
    rotated_instances = []    
    if PS_point is not None:
        try:
            topology_collection[topology] = create_hard_coded_topology(topology, PS_point, name=topology) 

            # Now also create rotated copy of the topology at hand for numerical stability diagnostics
            for rotation_matrix_specifications in rotations_for_stability_check:
                rotation_matrix_4d = vectors.LorentzVector.rotation_matrix(*rotation_matrix_specifications[0])
                rotation_matrix_3d = vectors.Vector.rotation_matrix(*rotation_matrix_specifications[0])
                for rotation_matrix_specification in rotation_matrix_specifications[1:]:
                    rotation_matrix_4d = rotation_matrix_4d.dot(vectors.LorentzVector.rotation_matrix(*rotation_matrix_specification))
                    rotation_matrix_3d = rotation_matrix_3d.dot(vectors.Vector.rotation_matrix(*rotation_matrix_specification))                    
                rotated_topology_name = '%s_rotated_%s'%(topology,str(rotation_matrix_specifications))
                topology_collection[rotated_topology_name] = \
                            create_hard_coded_topology(topology, PS_point.transform(rotation_matrix_4d), name=rotated_topology_name)
                rotated_instances.append({
                    'rotation_matrix_angles':   rotation_matrix_specifications,
                    'rotation_matrix':          rotation_matrix_3d,
                    'rotated_topology_name':    rotated_topology_name,
                    'rust_instance':            None # Will be filled in later
                })

        except ZeroDivisionError:
            print("ERROR: topology '%s' cannot be generated with custom kinematics."%topology)
            sys.exit(1)
    try:
        # Make sure to synchronise the topology collection and hyperparameters with rust
        print("Synchronising ltd_commons.py with yaml configuration files for rust...")
        ltd_commons.synchronize('.')
    except Exception as e:
        print("ERROR: Could not synchronize yaml and python configuration files.")
        raise e

    # Generate the list of t_values if None
    if t_values is None:
        if test_name == 'pole_check':
            if t_values_steps == 'linear':
                t_values = [ ( t_range[0]+ ( ( (t_range[1]-t_range[0]) / float(t_n_points+1)) * i ) ) for i in range(1, t_n_points+1)]
            elif t_values_steps == 'exponential':
                t_values = [ ( t_range[1] - ( (t_range[1]-t_range[0]) * 10.**(-float(i)) ) ) for i in range(t_n_points)]
        elif test_name == 'cancellation_check':
            if t_values_steps == 'linear':
                t_values = [ ( t_range[1] - ( ( (t_range[1]-t_range[0]) / float(t_n_points+1)) * i ) ) for i in range(1, t_n_points+1)]
            elif t_values_steps == 'exponential':
                t_values = [ ( t_range[0] + ( (t_range[1]-t_range[0]) * 10.**(-float(i)) ) ) for i in range(t_n_points)]

    if topology not in topology_collection:
        raise BaseException('Topology "%s" not present in the hard-coded topology collection in ltd_commons.py'%topology)
    topology = topology_collection[topology]

    if topology.n_loops > 2:
        print("This testing tool only supports up to two-loop topologies, not %d."%topology.n_loops)
        sys.exit(1)

    # Analyse all the surfaces
    diagnostic_tool = diagnostics.Diagnostic_tool(topology)
    print('Now analysing surfaces ...')
    all_surfaces = diagnostic_tool.all_surfaces(n_points_hyperboloid_test)
    grouped_surfaces = diagnostic_tool.check_similarity(all_surfaces)
    if surface_ids is not None and any(s<0 for s in surface_ids):
        print("Listing of all surfaces:")
        print(diagnostics.print_all_surfaces(grouped_surfaces, show_group_members=True))
        sys.exit(1)
    else:
        print('\n'.join(diagnostics.print_all_surfaces(grouped_surfaces, show_group_members=False).split('\n')[-1:]))
    
    #print(diagnostics.print_all_surfaces(grouped_surfaces, show_group_members=False))
    #stop
    if surface_ids is not None and any(s>len(grouped_surfaces) for s in surface_ids):
        print("ERROR: Specified surface_ids '%s' exceed the number of detected group of identical surfaces (%d)."%(
            str(surface_ids), len(grouped_surfaces)) )
        sys.exit(1)

    # Now collect one representative of each elliptic and hyperboloid surface
    # The tuple contain (cut_propagators) , surface )
    elliptic_surfaces       = { i_surface+1 : (tuple(surface_group[0][0]), surface_group[0][1] ) 
        for i_surface, surface_group in enumerate(grouped_surfaces) if surface_group[0][1].is_ellipsoid()}
    hyperboloid_surfaces    = { i_surface+1 : (tuple(surface_group[0][0]), surface_group[0][1] )
        for i_surface, surface_group in enumerate(grouped_surfaces) if not surface_group[0][1].is_ellipsoid()}

    if test_name == 'pole_check' and len(elliptic_surfaces)==0:
        print("ERROR: No pole check to perform since no ellipsooid is detected for the topology '%s'."%topology.name)
        sys.exit(1)
    if test_name == 'cancellation_check' and len(hyperboloid_surfaces)==0:
        print("ERROR: No pole check to perform since no hyperboloid is detected for the topology '%s'."%topology.name)
        sys.exit(1)

    if test_name == 'cancellation_check':
        surfaces_to_analyse = hyperboloid_surfaces
    if test_name == 'pole_check':
        surfaces_to_analyse = elliptic_surfaces

    # Apply the filtering if specified
    if surface_ids is not None:
        surfaces_to_analyse = {s_id: s for s_id, s in surfaces_to_analyse.items() if s_id in surface_ids}
    
    # Start a rust instance to be used throughout the tests
    rust_instance = LTD(
        settings_file = pjoin(root_path,'hyperparameters.yaml'),
        topology_file = pjoin(root_path,'topologies.yaml'),
        name = topology.name,
    )
    for rotated_instance in rotated_instances:
        rotated_instance['rust_instance'] = LTD(
            settings_file = pjoin(root_path,'hyperparameters.yaml'),
            topology_file = pjoin(root_path,'topologies.yaml'),
            name = rotated_instance['rotated_topology_name'],
        ) 

    n_loops_performed = 0
    are_parametrisation_vectors_unspecified = any(v is None for v in parametrisation_vectors)
    com_energy = topology.get_com_energy()
    while True:
        try:
            # post-process the vector definitions by replacing the ones that are None and should be taken random.
            set_direction_vectors = [
                (vec if vec is not None else vectors.Vector([(random.random()-0.5)*2 for _ in range(3)]))
                for vec in direction_vectors]
            set_parametrisation_vectors = [
                (vec if vec is not None else vectors.Vector([(random.random() - 0.5) * 2 * com_energy for _ in range(3)]))
                for vec in parametrisation_vectors]

            set_parametrisation_uv = {}
            set_parametrisation_uv['u']  = (parametrisation_uv['u'] if parametrisation_uv['u'] is not None else
                                                                                                (random.random() - 0.5) * 2)
            set_parametrisation_uv['v']  = (parametrisation_uv['v'] if parametrisation_uv['v'] is not None else
                                                                                                (random.random() - 0.5) * 2)
            
            if do_loop is not None and n_loops_performed%1 == 0:
                print("\n")
                print("Now testing analytical continuation with the following setup:")
                print("  direction_vectors=\n     %s"%str(set_direction_vectors))
                print("  parametrisation_vectors=\n     %s"%str(set_parametrisation_vectors))
                print("  surface parametrisation:\n     u=%(u).16f, v=%(v).16f"%set_parametrisation_uv)
                print("\n")


            # Normalise the direction vectors if doing a pole check
            if test_name in ['pole_check','cancellation_check']:
                set_direction_vectors = [vec*(1.0/max(abs(v) for v in vec)) for vec in set_direction_vectors]
            if load_results_from is None:
                if save_results_to is not None:
                    log_stream = open(save_results_to,'w')
                else:
                    log_stream = None
                
                if test_name == 'pole_check':
                    scanner = PoleScanner(
                        log_stream              =   log_stream,
                        topology                =   topology,
                        direction_vectors       =   set_direction_vectors,
                        parametrisation_vectors =   set_parametrisation_vectors,
                        parametrisation_uv      =   set_parametrisation_uv,
                        t_values                =   t_values,
                        plot_surface            =   plot_surface,
                        diagnostic_tool         =   diagnostic_tool,
                        rust_instance           =   rust_instance,
                        rotated_instances       =   rotated_instances
                    )

                elif test_name == 'cancellation_check':
                    scanner = DualCancellationScanner(
                        log_stream              =   log_stream,
                        topology                =   topology,
                        direction_vectors       =   set_direction_vectors,
                        parametrisation_vectors =   set_parametrisation_vectors,
                        parametrisation_uv      =   set_parametrisation_uv,
                        t_values                =   t_values,
                        diagnostic_tool         =   diagnostic_tool,
                        rust_instance           =   rust_instance,
                        rotated_instances       =   rotated_instances
                    )
                start_time = time.time()
                all_test_results = {'surface_results':[]}
                for surface_ID in sorted(surfaces_to_analyse.keys()):
                    surface = surfaces_to_analyse[surface_ID]
                    #print('Now analysing surface with ID=%d :\n%s'%(
                    #    surface_ID, surface[1].__str__(cut_propagators = surface[0])
                    #))
                    if do_loop is None: print('Now analysing surface with ID #%d...'%surface_ID)
                    if test_name == 'pole_check':
                        test_results = scanner.scan(surface_ID, surface)
                    elif test_name == 'cancellation_check':
                        test_results = scanner.scan(surface_ID, surface,
                                                    force = are_parametrisation_vectors_unspecified)
                    try:
                        all_test_results['surface_results'].append(
                            {   'surface': test_results.pop('surface'),
                                'scan_results' : test_results.pop('scan_results'),
                                'onshell_prop_id' : test_results.pop('onshell_prop_id','N/A')
                            }
                        )
                    except KeyError:
                        pass
                    all_test_results.update(test_results)

                all_test_results['run_time'] = time.time()-start_time
                if log_stream is not None:
                    log_stream.write(yaml.dump([('run_time', all_test_results['run_time']), ], Dumper=noalias_dumper))
                    log_stream.close()

                if plot_surface:
                    plt.show()

            else:
                try:
                    print("Now loading results from '%s'..."%load_results_from)
                    all_test_results = load_results_from_yaml(load_results_from)
                    print('Done!')
                except Exception as e:
                    print("ERROR: Could not load scan results from file '%s'. Error: %s"%(
                        load_results_from, str(e) ))
                    raise e

            #pprint(test_results)
            if test_name == 'pole_check':
                analyser = PoleResultsAnalyser(surfaces_to_analyse, all_test_results)
            elif test_name == 'cancellation_check':
                analyser = DualCancellationResultsAnalyser(surfaces_to_analyse, all_test_results)

            n_loops_performed += 1
            if do_loop is not None:
                if analyser.test_imaginary_parts_in_scan():
                    print("The setup yielding wrong analytical continuation:")
                    print("  direction_vectors=\n     %s"%str(set_direction_vectors))
                    print("  parametrisation_vectors=\n     %s"%str(set_parametrisation_vectors))
                    print("  surface parametrisation:\n     u=%(u).16f, v=%(v).16f"%set_parametrisation_uv)
                    sys.exit(1)

            else:
                analyser.analyse(**analysis_options)

            if (do_loop is None) or (do_loop is not None and do_loop > 0 and n_loops_performed>=do_loop):
                break

            if n_loops_performed%1 == 0:
                print('Number of test scans performed: %d'%n_loops_performed)

        except KeyboardInterrupt:
            print('Run aborted by user. Number of test scans performed so far: %d' % n_loops_performed)
            break

    if do_loop is not None:
        print('%d test scans successfully performed.'%n_loops_performed)        

