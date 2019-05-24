# n loop cut structures
import numpy as numpy
import time
import numpy.linalg

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

topology = topologies.hard_coded_topology_collection['mercedes']
graph = ltd_utils.TopologyGenerator([(None, loop_line.start_node, loop_line.end_node) for loop_line in topology.loop_lines])
momentum_bases = graph.loop_momentum_bases()
signature_matrix = [loop_line.signature for loop_line in topology.loop_lines]

allowed_systems = graph.find_allowed_systems(momentum_bases, signature_matrix)
residues = []
for contour_closure in itertools.product(*([[0,1]]*graph.n_loops)):
	#if contour_closure == (1,1,1,1):
		residues += graph.evaluate_residues(allowed_systems, contour_closure)
residues.sort()

residues = graph.evaluate_thetas(residues)
residues = graph.remove_cancelling_residues(residues)

for residue in residues:
	print residue

"""
cut_stucture = []
for residue, momentum_basis in zip(residues, momentum_bases):
    cut_struct_iter = iter(residue[0][1])
    cut_stucture += [[next(cut_struct_iter)
                      if i in momentum_basis else 0 for i in range(len(topology.loop_lines))]]
print cut_stucture

print get_cut_stucture_string(cut_stucture)
"""



















