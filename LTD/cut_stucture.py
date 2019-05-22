# n loop cut structures
import numpy as numpy
import time
import numpy.linalg

import ltd_commons
import topologies
import itertools

import ltd_utils

def get_cut_stucture_string(cut_stucture):
	for i,cut_stuct in enumerate(cut_stucture):
		for j,cut in enumerate(cut_stuct):
			if cut == 0:
				cut_stucture[i][j] = 'LoopLine.NO_CUT'
			elif cut == -1:
				cut_stucture[i][j] = 'LoopLine.NEGATIVE_CUT'
			elif cut == 1:
				cut_stucture[i][j] = 'LoopLine.POSITIVE_CUT'
	cut_stucture_string = '(' + ',\n'.join([ '(' + ','.join([cut for cut in cut_stuct]) + ')' for cut_stuct in cut_stucture]) +')'
	return cut_stucture_string

topology = ltd_commons.hard_coded_topology_collection['non_planar_four_loop']
#print topology
graph = ltd_utils.TopologyGenerator([(None, loop_line.start_node, loop_line.end_node) for loop_line in topology.loop_lines])
momentum_bases = graph.loop_momentum_bases()
signature_matrix = [loop_line.signature for loop_line in topology.loop_lines]

allowed_systems = graph.find_allowed_systems(momentum_bases, signature_matrix)
residues = []
for contour_closure in itertools.product(*([[0,1]]*graph.n_loops)):
	#print contour_closure
	#if contour_closure == (1,1,1):
		residue = graph.evaluate_residues(allowed_systems, contour_closure)
		residue = graph.evaluate_thetas(residue)
		residue = graph.remove_cancelling_residues(residue)
		residues += residue
residues.sort()

fixed_sign = True
for residue in residues:
	#if not (residue[0][2] == -1.):
	#	fixed_sign = False
	print residue
#print 'fixed_sign: ', fixed_sign
"""
cut_stucture = []
for residue, momentum_basis in zip(residues, momentum_bases):
    cut_struct_iter = iter(residue[0][1])
    cut_stucture += [[next(cut_struct_iter)
                      if i in momentum_basis else 0 for i in range(len(topology.loop_lines))]]
print cut_stucture

print get_cut_stucture_string(cut_stucture)
"""



















