#!/usr/bin/env python2

import logging
logging.basicConfig(level=logging.DEBUG)

import itertools
import math
import random
import os
import sys
import numpy as np

pjoin = os.path.join

root_path = os.path.dirname(os.path.realpath( __file__ ))
sys.path.insert(0, pjoin(root_path,os.pardir))
sys.path.insert(0, pjoin(root_path,os.pardir,os.pardir))

import ltd_utils
from ltd_utils import TopologyGenerator, TopologyCollection
import vectors

try:
    # Import the rust bindings
    from ltd import LTD
except ImportError:
    raise BaseException(
        "Could not import the rust back-end 'ltd' module. Compile it first with:"
        " ./make_lib from within the pyNLoop directory." )

def run_closure_test_for_master_topology(
        master_topology,
        external_momenta,
        loop_masses,
        loop_momenta,
        max_n_ltd_cut_structures_to_test=-1
):
    """ Run the consistency check of the evaluation of the residues with different choices for closing the contour.
    loop_momenta can be an integer which implies that it will randomly test on that many loop momentum kinematic configuration. """

    collection_of_topologies_to_test = TopologyCollection()

    loop_momenta_basis = master_topology.loop_momentum_bases()
    n_loops = master_topology.n_loops

    sorted_collection_names = []
    # Using the original momentum routing start by testing all possible ways of closing the contour
    for closing_option in itertools.product(*([(0,1),]*n_loops)):
        if max_n_ltd_cut_structures_to_test>0 and len(sorted_collection_names)>=max_n_ltd_cut_structures_to_test:
            break
        test_topology_name = "test_closure_%s" % (''.join('%d'%co for co in closing_option))
        sorted_collection_names.append(test_topology_name)
        logging.debug("Now generating topology for closing option: %s"%(''.join('%d'%co for co in closing_option)))
        collection_of_topologies_to_test.add_topology(
            master_topology.create_loop_topology(
                test_topology_name,
                ext_mom=external_momenta,
                mass_map=loop_masses,  # no masses
                loop_momenta_names=loop_momenta_basis[0], # Always take the first spanning tree as a defining one for this first test
                contour_closure=closing_option
            ),
            entry_name = test_topology_name
        )

    # Export this topology to a test yaml collection
    collection_of_topologies_to_test.export_to(os.path.join(root_path, 'test_topologies.yaml'))

    # Instantiate a generic multi-purpose rust worker (for access to parametrisationa and such)
    multipurpose_rust_worker = LTD(
            settings_file=pjoin(os.path.pardir, 'hyperparameters.yaml'),
            topology_file=pjoin('test_topologies.yaml'),
            name=sorted_collection_names[0],
        )

    # Now generate the loop momenta configurations to consider
    all_loop_momenta_configs = []
    if isinstance(loop_momenta, int):
        for i_loop_momenta in range(loop_momenta):
            this_loop_momenta_config = []
            for i_loop in range(n_loops):
                # Generate random unit hypercube input
                x_inputs = [random.random() for _ in range(3)]
                # Map the x inputs to an actual loop momenta configuration
                this_loop_momenta_config.append([float(ki) for ki in
                        multipurpose_rust_worker.parameterize(x_inputs, i_loop)][:-1])
            all_loop_momenta_configs.append(this_loop_momenta_config)
    else:
        all_loop_momenta_configs = loop_momenta

    # Loop over all topologies to test and collect the results
    for loop_momenta_config in all_loop_momenta_configs:
        logging.info("Now investigating the following loop momenta configuration:\n%s"%(
            '\n'.join('k%d = %s'%(i_loop,', '.join('%s%-20.16e'%('+' if ki >=0 else '-',abs(ki)) for ki in k)) for i_loop, k in enumerate(loop_momenta_config) )
        ))
        # Map back the loop momenta to x in hypercube
        x_inputs = []
        for i_loop, k in enumerate(loop_momenta_config):
            x_inputs.extend([float(xi) for xi in multipurpose_rust_worker.inv_parameterize(k, i_loop)][:-1])

        closing_test_results = {}
        for topology_name in sorted_collection_names:
            # Now start a rust worker
            rust_worker = LTD(
                settings_file=pjoin(os.path.pardir, 'hyperparameters.yaml'),
                topology_file=pjoin('test_topologies.yaml'),
                name=topology_name,
            )
            result_re, result_im = rust_worker.evaluate(x_inputs)
            logging.info('%-40s = %-20.16e %s%-20.16e*I'%(topology_name+'_f64', result_re, '+' if result_im>=0 else '-', abs(result_im)))
            closing_test_results[topology_name+'_f64'] = complex(result_re, result_im)
            result_re, result_im = rust_worker.evaluate_f128(x_inputs)
            logging.info('%-40s = %-20.16e %s%-20.16e*I'%(topology_name+'_f128', result_re, '+' if result_im>=0 else '-', abs(result_im)))
            closing_test_results[topology_name+'_f128'] = complex(result_re, result_im)


def find_defining_propagators(topology):
    """ find the defining edges, with the corresponding shifts in a given generated topology."""

    defining_propagators= []
    for i_loop in range(topology.n_loops):
        found_it = False
        target_signature = [0,]*topology.n_loops
        target_signature[i_loop] = 1
        target_signature = tuple(target_signature)
        for loop_line in topology.loop_lines:
            for prop in loop_line.propagators:
                if tuple(prop.signature) == target_signature and all(qi==0. for qi in prop.q):
                    found_it = True
                    defining_propagators.append((loop_line.start_node, loop_line.end_node, prop))
                    break
            if found_it:
                break
        if not found_it:
            logging.critical("Could not find defining propagators for loop momentum #%d in topology %s"%(i_loop, str(topology)))
            return None
    return defining_propagators

def find_equivalent_propagator(topology, target_propagator):
    """ Find the propagator with the same start and end node in specified topology.
    !!!
    WARNING: THIS EFFECTIVELY USES THE LOOP PROPAGATOR MASS AS A FLAG! It must be generalised using the edge names instead!
    !!!
    """

    for loop_line in topology.loop_lines:
        for prop in loop_line.propagators:
            if prop.m_squared != target_propagator[-1].m_squared:
                continue
            if loop_line.start_node != target_propagator[0] or loop_line.end_node != target_propagator[1]:
                logging.critical("ERROR the matching loop lines seem to have their edges flipped!")
                return None
            return prop

    logging.critical("ERROR Could not find a match for specified propagator: %s"%str(target_propagator))
    return None

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        yield l[i:i + n]

def manual_evaluation(rust_worker, topology, loop_momenta, f128=False):
    """ 'Manually' evalute the LTD expression by summing explicitely over the duals of the specified topology for the
    specified loop momentum."""
    evaluation = {}
    for ltd_cut_index, ltd_cut_structure in enumerate(topology.ltd_cut_structure):
        cut_propagator_indices = [[index*c for index in range(1,len(topology.loop_lines[i].propagators)+1)] if c!=0 else [0,]
                                             for i, c in enumerate(ltd_cut_structure)]
        for cut_index, cut_structure  in enumerate(itertools.product( *cut_propagator_indices )):
            if f128:
                res = rust_worker.evaluate_cut_f128(
                        [ [ ( float(vi.real), float(vi.imag) ) for vi in v ] for v in loop_momenta],
                         ltd_cut_index,  cut_index)
            else:
                res = rust_worker.evaluate_cut(
                        [ [ ( float(vi.real), float(vi.imag) ) for vi in v ] for v in loop_momenta],
                         ltd_cut_index,  cut_index)
            if not math.isnan(res[0]) and not math.isnan(res[1]):
                evaluation[((ltd_cut_index, tuple(ltd_cut_structure)),(cut_index, cut_structure))] = complex(res[0],res[1])
            else:
                 print("Nan found for entry %s"%str(((ltd_cut_index, ltd_cut_structure),(cut_index, cut_structure))))

    return evaluation

def run_parametrisation_test_for_master_topology(
        master_topology,
        external_momenta,
        loop_masses,
        loop_momenta,
        max_n_ltd_cut_structures_to_test=-1
):
    """ Run the consistency check of the evaluation of the residues with different way of building the momentum routing
    and closing the contour from below. loop_momenta can be an integer which implies that it will randomly test on that many
    loop momentum kinematic configuration. """

    collection_of_topologies_to_test = TopologyCollection()

    loop_momenta_basis = master_topology.loop_momentum_bases()
    n_loops = master_topology.n_loops

    sorted_collection_names = [ "test_parametrisation_%s"%('_'.join('%d'%i_line for i_line in loop_momenta_basis[0])) ]
    # First generate the "defining" topology with a momentum routing obtained by the first spanning tree.
    defining_topology = master_topology.create_loop_topology(
                sorted_collection_names[0],
                ext_mom=external_momenta,
                mass_map=loop_masses,  # no masses
                loop_momenta_names=loop_momenta_basis[0], # Always take the first spanning tree as a defining one
                contour_closure=[0,]*n_loops # Always close from below for this test
            )
    conversion_rules = { sorted_collection_names[0]: (np.diag([1,]*n_loops),tuple([vectors.LorentzVector(),]*n_loops) ) }

    collection_of_topologies_to_test.add_topology(defining_topology,entry_name = sorted_collection_names[0])

    # Sanity check
    all_masses = []
    for loop_line in defining_topology.loop_lines:
        for prop in loop_line.propagators:
            if prop.m_squared in all_masses:
                logging.exception("ERROR: This test currently only supports topologies with all loop masses differing!")
                return
            else:
                all_masses.append(prop.m_squared)

    # Find the defining edges for each loop in the defining topology:
    defining_propagators = find_defining_propagators(defining_topology)

    # Now generate all the others, while computing the conversion matrix
    for a_loop_momenta_basis in loop_momenta_basis[1:]:
        if max_n_ltd_cut_structures_to_test>0 and len(sorted_collection_names)>=max_n_ltd_cut_structures_to_test:
            break
        test_topology_name = "test_parametrisation_%s"%('_'.join('%d'%i_line for i_line in a_loop_momenta_basis))
        sorted_collection_names.append(test_topology_name)
        logging.debug("Now generating topology for mommetum routing option: %s"%('_'.join('%d'%co for co in a_loop_momenta_basis)))
        test_topology = master_topology.create_loop_topology(
                test_topology_name,
                ext_mom=external_momenta,
                mass_map=loop_masses,  # no masses
                loop_momenta_names=a_loop_momenta_basis, # Always take the first spanning tree as a defining one for this first test
                contour_closure=[0,]*n_loops
            )
        collection_of_topologies_to_test.add_topology(test_topology, entry_name = test_topology_name)

        # Now we must figure out the conversion matrix.
        defining_propagators_for_this_topology = find_defining_propagators(test_topology)

        conversion_matrix = []
        conversion_shifts = []
        for defining_propagator in defining_propagators_for_this_topology:
            matched_prop = find_equivalent_propagator(defining_topology, defining_propagator)
            conversion_matrix.append(matched_prop.signature)
            conversion_shifts.append(vectors.LorentzVector(matched_prop.q))
        conversion_rules[test_topology_name] = (np.array(conversion_matrix), tuple(conversion_shifts))

    # Export this topology to a test yaml collection
    collection_of_topologies_to_test.export_to(os.path.join(root_path, 'test_topologies.yaml'))

    # Instantiate a generic multi-purpose rust worker (for access to parametrisationa and such)
    multipurpose_rust_worker = LTD(
            settings_file=pjoin(os.path.pardir, 'hyperparameters.yaml'),
            topology_file=pjoin('test_topologies.yaml'),
            name=sorted_collection_names[0],
        )

    # Now generate the loop momenta configurations to consider
    all_loop_momenta_configs = []
    if isinstance(loop_momenta, int):
        for i_loop_momenta in range(loop_momenta):
            this_loop_momenta_config = []
            for i_loop in range(n_loops):
                # Generate random unit hypercube input
                x_inputs = [random.random() for _ in range(3)]
                # Map the x inputs to an actual loop momenta configuration, adding some random energy
                three_momentum = [float(ki) for ki in multipurpose_rust_worker.parameterize(x_inputs, i_loop)][:-1]
                this_loop_momenta_config.append(vectors.LorentzVector([random.random()*three_momentum[0]]+three_momentum))
            all_loop_momenta_configs.append(this_loop_momenta_config)
    else:
        all_loop_momenta_configs = loop_momenta

    # Loop over all topologies to test and collect the results
    for loop_momenta_config in all_loop_momenta_configs:
        logging.info("Now investigating the following loop momenta configuration:\n%s"%(
            '\n'.join('k%d = %s'%(i_loop,', '.join('%s%-20.16e'%('+' if ki >=0 else '-',abs(ki)) for ki in k)) for i_loop, k in enumerate(loop_momenta_config) )
        ))

        parametrisation_test_results = {}
        for topology_name in sorted_collection_names:

            # Generate the mapped kinematics
            conversion_matrix, conversion_shifts = conversion_rules[topology_name]
            mapped_momenta_config = []
            for i_loop in range(n_loops):
                mapped_momenta_config.append( sum(k*s for s,k in zip(conversion_matrix[i_loop],loop_momenta_config))
                                                                                            + conversion_shifts[i_loop] )

            # Printout the consistency check
            logging.info('%-60s = %-20.16e'%("%s : consistency 4D check"%topology_name,
                                       collection_of_topologies_to_test[topology_name].evaluate(mapped_momenta_config)))

            # Now we can project back onto 3D momenta
            mapped_momenta_config = [[ ki for ki in k[1:] ] for k in mapped_momenta_config]
            #logging.info(mapped_momenta_config)
            # And map back to x hypercube space
            x_inputs = []
            for i_loop, k in enumerate(mapped_momenta_config):
                x_inputs.extend([float(xi) for xi in multipurpose_rust_worker.inv_parameterize(k, i_loop)][:-1])
            # Sanity check on of the function inv_parametrize:
            three_momentum_reconstructed=[]
            for chunk in chunks(x_inputs,3):
                three_momentum_reconstructed.append([float(ki) for ki in multipurpose_rust_worker.parameterize(chunk, i_loop)][:-1])
            #print(mapped_momenta_config)
            #print(three_momentum_reconstructed)
            for i, k in enumerate(mapped_momenta_config):
                for j, ki in enumerate(k):
                    missmatch = abs(ki-three_momentum_reconstructed[i][j])/abs(ki)
                    if missmatch>1.0e-12:
                        logging.critical("Inverse parametrisation not working properly! relative mismatch = %e"%missmatch)
                        return


            # Now start a rust worker for that topology
            rust_worker = LTD(
                settings_file=pjoin(os.path.pardir, 'hyperparameters.yaml'),
                topology_file=pjoin('test_topologies.yaml'),
                name=topology_name,
            )
            result_re, result_im = rust_worker.evaluate(x_inputs)
            logging.info('%-60s = %-20.16e %s%-20.16e*I'%(topology_name+'_f64', result_re, '+' if result_im>=0 else '-', abs(result_im)))
            parametrisation_test_results[topology_name+'_f64'] = complex(result_re, result_im)
            result_re, result_im = rust_worker.evaluate_f128(x_inputs)
            logging.info('%-60s = %-20.16e %s%-20.16e*I'%(topology_name+'_f128', result_re, '+' if result_im>=0 else '-', abs(result_im)))
            parametrisation_test_results[topology_name+'_f128'] = complex(result_re, result_im)
            # Also perform a "manual" evaluation by summing explicitly over duals
            all_duals = manual_evaluation(rust_worker,collection_of_topologies_to_test[topology_name],mapped_momenta_config)
            manual_result = sum(all_duals.values())
            parametrisation_test_results[topology_name+'_manual'] = manual_result
            logging.info('%-60s = %-20.16e %s%-20.16e*I'%(topology_name+' : manual', manual_result.real,
                                                            '+' if manual_result.imag>=0 else '-', abs(manual_result.imag)))
            all_duals = manual_evaluation(rust_worker,collection_of_topologies_to_test[topology_name],mapped_momenta_config, f128=True)
            manual_result = sum(all_duals.values())
            parametrisation_test_results[topology_name+'_manual_f128'] = manual_result
            logging.info('%-60s = %-20.16e %s%-20.16e*I'%(topology_name+'_f128 : manual', manual_result.real,
                                                            '+' if manual_result.imag>=0 else '-', abs(manual_result.imag)))

if __name__ == '__main__':

    random.seed(1)

    bubble = TopologyGenerator([
        ('p1', 1, 2), ('p2', 2, 1),
        ('q1', 101, 1), ('q2', 102, 2),
    ])
    q1 = vectors.LorentzVector([0.1, 0.2, 0.5, 0.1])
    external_momenta = { 'q1': q1, 'q2': -q1 }
    loop_masses = {'p1': 1.1, 'p2': 1.2 }
    logging.info('Now testing closure independence of the one-loop bubble:')
    run_closure_test_for_master_topology(bubble, external_momenta, loop_masses, 1)
    logging.info('Now testing the parametrisation independence of the one-loop bubble:')
    run_parametrisation_test_for_master_topology(bubble, external_momenta, loop_masses, 1)
    logging.info('')

    box = TopologyGenerator([
        ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4), ('p4', 4, 1),
        ('q1', 101, 1), ('q2', 102, 2), ('q3', 103, 3), ('q4', 104, 4)
    ])
    q1 = vectors.LorentzVector([0.1, 0.2, 0.5, 0.1])
    q2 = vectors.LorentzVector([-0.3, 0.4, 0.1, 0.2])
    q3 = vectors.LorentzVector([0.1, 0.2, 0.5, 0.3])
    external_momenta = { 'q1': q1, 'q2': q2 , 'q3': q3, 'q4': -q1-q2-q3 }
    loop_masses = {'p1': 1.1, 'p2': 1.2, 'p3': 1.3, 'p4': 1.4}
    logging.info('Now testing closure independence of the one-loop box:')
    run_closure_test_for_master_topology(box, external_momenta, loop_masses, 1)
    logging.info('Now testing the parametrisation independence of the one-loop box:')
    run_parametrisation_test_for_master_topology(box, external_momenta, loop_masses, 1)
    logging.info('')

    doublebox = TopologyGenerator([
        ('q1', 101, 1), ('q2', 102, 2), ('q3', 103, 3), ('q4', 104, 4),
        ('p1', 1, 6), ('p2', 6, 7), ('p3', 7, 2), ('p4', 2, 1),
        ('p5', 6, 3), ('p6', 3, 4), ('p7', 4, 7)
    ])
    q1 = vectors.LorentzVector([1.2, 2.2, 1.0, 0.4])
    q2 = vectors.LorentzVector([2.0, -5.2, 2.1, 0.0])
    q3 = vectors.LorentzVector([-1.6, 0.1, -12.5, 2.4])
    q4 = -q1 - q2 - q3
    external_momenta = {'q1': q1, 'q2': q2, 'q3': q3, 'q4': q4}
    loop_masses = {'p1': 1.0, 'p2': 1.1, 'p3': 1.2, 'p4': 1.3, 'p5': 1.4, 'p6': 1.5, 'p7': 1.6}
    logging.info('Now testing closure independence of the double box:')
    run_closure_test_for_master_topology(doublebox, external_momenta, loop_masses, 1)
    logging.info('Now testing the parametrisation independence of the double box:')
    run_parametrisation_test_for_master_topology(doublebox, external_momenta, loop_masses, 1)
    logging.info('')

    # non planar four loop
    non_planar_four_loop = TopologyGenerator([
        ('q', 0, 1), ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),
        ('p4', 4, 1), ('p5', 5, 1), ('p6', 4, 6), ('p7', 5, 6),
        ('p8', 6, 2), ('p9', 5, 3), ('-q', 7, 3)
    ])
    q = vectors.LorentzVector([0.1, 0.2, 0.3, 0.4])
    external_momenta = {'q': q, '-q': -q}
    loop_masses = {'p1': 1.0, 'p2': 1.1, 'p3': 1.2, 'p4': 1.3, 'p5': 1.4, 'p6': 1.5, 'p7': 1.6,
                   'p8': 1.7, 'p9': 1.8}
    logging.info('Now testing closure independence of the non planar four loops:')
    run_closure_test_for_master_topology(non_planar_four_loop, external_momenta, loop_masses, 1, max_n_ltd_cut_structures_to_test=4)
    logging.info('Now testing the parametrisation independence of the non planar four loops:')
    run_parametrisation_test_for_master_topology(non_planar_four_loop, external_momenta, loop_masses, 1, max_n_ltd_cut_structures_to_test=4)

    logging.info("Uncomment the 'exit' statement below to perform tests in the more complicated topologies.")
    sys.exit(0)
