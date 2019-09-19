#!/usr/bin/env python2
import itertools

import PRL_cut_structure

from topologies import hard_coded_topology_collection

def compare_cut_structures(cs_A, cs_B):

    if len(cs_A)!=len(cs_B): return False
    for cs_A_line, cs_B_line in zip(cs_A, cs_B):
        if len(cs_A_line)!=len(cs_B_line): return False
        for cs_A_el, cs_B_el in zip(cs_A_line, cs_B_line):
            if cs_A_el != cs_B_el: return False
    return True

def specify_order(sig, order):
    
    return [sig[i] for i in order]

for topo_entry_name, topo in sorted(hard_coded_topology_collection.items()):
    if 'ala_weinzierl' in topo_entry_name.lower(): continue
    # Uncomment below to test a single topology
    #if topo_entry_name!='PRL_Mercedes': continue
    print "Now testing: %s"%topo_entry_name

    # The specify_order function below is really just to test independence upon relablling of the loop momenta
    # By default below it does not change the original order of the signature given the break at the very end
    # which can be commented for testing purposes.
    for order_choice in itertools.permutations(range(topo.n_loops)):
        print "Considering signature permutation: %s"%str(order_choice)
        loop_lines_signatures = [specify_order(ll.signature,order_choice) for ll in topo.loop_lines]
        cut_stucture_generator = PRL_cut_structure.CutStructureGenerator(loop_lines_signatures)
        
        if (topo_entry_name == 'PRL_Mercedes'
            or topo_entry_name == 'PRL_Mercedes_6p_massive'
            or topo_entry_name == 'TriangleBoxTriangle'):
            contour_closure = [PRL_cut_structure.CLOSE_BELOW,PRL_cut_structure.CLOSE_BELOW,PRL_cut_structure.CLOSE_ABOVE]
        elif topo_entry_name == 'manual_AltDoubleTriangle':
            contour_closure = [PRL_cut_structure.CLOSE_ABOVE,PRL_cut_structure.CLOSE_BELOW]
        elif (topo_entry_name == 'manual_TriangleBoxBox_alt'
            or topo_entry_name == 'manual_TriangleBoxBox_alt_ellipses'
            or topo_entry_name == 'manual_TriangleBoxTriangle_alt'
            or topo_entry_name == 'manual_TriangleBoxTriangle_alt_ellipses'):
            contour_closure = [PRL_cut_structure.CLOSE_ABOVE,PRL_cut_structure.CLOSE_BELOW,PRL_cut_structure.CLOSE_ABOVE]
        elif (topo_entry_name == 'manual_TriangleBoxBox'
            or topo_entry_name == 'manual_TriangleBoxBox_ellipses'
            or topo_entry_name == 'manual_TriangleBoxTriangle'
            or topo_entry_name == 'manual_TriangleBoxTriangle_ellipses'
            or topo_entry_name == 'manual_TripleBox'
            or topo_entry_name == 'manual_TripleBox_Weinzierl'
            or topo_entry_name == 'manual_TripleBox_no_ellipse'):
            contour_closure = [PRL_cut_structure.CLOSE_BELOW,PRL_cut_structure.CLOSE_BELOW,PRL_cut_structure.CLOSE_BELOW]
        else:
            contour_closure = [PRL_cut_structure.CLOSE_BELOW,]*topo.n_loops
        #print(contour_closure)
        cut_structure = sorted(cut_stucture_generator(contour_closure))
        if not compare_cut_structures(cut_structure, sorted(topo.ltd_cut_structure)):
            print 'From PRL: ', cut_structure
            print 'vs'
            print 'From  US: ', sorted(topo.ltd_cut_structure)
            print "ERROR, missmatch in cut structures for topology %s!"%topo_entry_name


        # Only consider the first original signature ordering
        break
