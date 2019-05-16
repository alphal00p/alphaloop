# n loop cut structures
import numpy as numpy
import time
import numpy.linalg

import ltd_commons
import topologies
import itertools

def spanning_tree(tree, graph, accum, result):
    # find all edges that connect the tree to a new node
    edges = [(i, e) for i, e in enumerate(graph) if len(tree & set(e)) == 1]

    if len(edges) == 0:
        # no more new edges, so we are done
        s = list(sorted(accum))
        if s not in result:
            result.append(s)
    else:
        for i, e in edges:
            spanning_tree(tree | set(e), graph, accum + [i], result)

def get_thetas_for_residue(sigmas,allowed_system,close_contour):
	# close_contour is a list, 0 means contour closed on lower half plane, 1 upper half plane for every int variable
	n_down = sum([-1 for i in close_contour if i == 0])
	contour_sign = (-1)**n_down*(-1)**n_loops # normalize sign with n_loops
	sign = numpy.prod(sigmas)*numpy.linalg.det(allowed_system[2])
	residue = [[allowed_system[0],sigmas,sign*contour_sign]]
	thetas = []
	for n_iter in xrange(n_loops):
		theta = [0 for i in xrange(n_loops)]
		sub_matrix = numpy.array([[allowed_system[2][i][j] for j in xrange(n_iter+1)] for i in xrange(n_iter+1)])
		if n_iter == 0:
			theta[allowed_system[1][n_iter]] = 1./numpy.linalg.det(sub_matrix)*sigmas[allowed_system[1][n_iter]]
			theta[allowed_system[1][n_iter]] *= (-1.)**close_contour[n_iter]
		else:
			for r in xrange(n_iter+1):
				subsub_matrix = numpy.array([[allowed_system[2][i][j] for j in xrange(n_iter)] for i in xrange(n_iter+1) if i != r])
				theta[allowed_system[1][r]] = (-1.)**(n_iter+r+1)*numpy.linalg.det(subsub_matrix)/numpy.linalg.det(sub_matrix)
				theta[allowed_system[1][r]] *= sigmas[allowed_system[1][r]]
				theta[allowed_system[1][r]] *= (-1.)**close_contour[r]
		thetas += [theta]
	residue += [thetas]
	return residue

def evaluate_thetas(residues):
	evaluated_theta_resides = []
	for residue in residues:
		for i,theta in enumerate(residue[1]):
			if all(x <= 0. for x in theta):
				add = False
				break
			else:
				add = True
				residue[1] = [theta for theta in residue[1] if not all(x >= 0. for x in theta)]
		if add:
			evaluated_theta_resides += [residue]
	return evaluated_theta_resides

def find_allowed_systems(spanning_trees):
	# get all possible momentum bases from spanning trees
	momentum_bases = [[i for i in range(len(graph)) if i not in r] for r in spanning_trees]
	# get transformation matrices: cut_mometum basis <=> loop momentum basis
	s_matrices = []
	for i,momentum_basis in enumerate(momentum_bases):
		s_matrix = numpy.array([topology.loop_lines[i].signature for i in momentum_basis])#,dtype=float)
		s_matrices += [s_matrix]

	allowed_systems = []
	n_iter = 1
	for res_index,s_matrix in enumerate(s_matrices):
		permuted_s_matrices = numpy.array(list(itertools.permutations(s_matrix)))
		permuted_energies = numpy.array(list(itertools.permutations(range(n_loops))))
		for perm_index,permuted_s_matrix in enumerate(permuted_s_matrices):
			for n_iter in xrange(n_loops):
				sub_matrix = numpy.array([[permuted_s_matrix[i][j] for j in xrange(n_iter+1)] for i in xrange(n_iter+1)])
				if numpy.linalg.det(sub_matrix) == 0:
					allowed = False
					break
				else:
					allowed = True
			if allowed:
				allowed_systems += [[res_index,permuted_energies[perm_index],permuted_s_matrix]]
			else:
				continue
	return allowed_systems

def remove_cancelling_residues(residues):
	for index_i,residue_i in enumerate(list(residues)):
		cancelling_residue = [[residue_i[0][0],residue_i[0][1],-residue_i[0][2]],residue_i[1]]
		if cancelling_residue in residues:
			residues.remove(cancelling_residue)
			residues.remove(residue_i)
	if len(residues) != residues[-1][0][0] + 1:
		print('Check residues!')
	return residues

def evaluate_residues(allowed_systems,close_contour):
	residues = []
	for sigmas in itertools.product(*([[1,-1]]*n_loops)):
		for allowed_system in allowed_systems:
			residue = get_thetas_for_residue(sigmas,allowed_system,close_contour)
			residues += [residue]
	residues.sort()
	residues = evaluate_thetas(residues)
	residues = remove_cancelling_residues(residues)
	return residues

def generate_cut_structure(residues,spanning_trees):
	momentum_bases = [[i for i in range(len(graph)) if i not in r] for r in spanning_trees]
	assert(len(residues)==len(momentum_bases))
	cut_stucture = []
	for residue,momentum_basis in zip(residues,momentum_bases):
		cut_struct_iter = iter(residue[0][1])
		cut_stucture += [[next(cut_struct_iter) if i in momentum_basis else 0 for i in xrange(n_loop_lines)]]
	return cut_stucture

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

topology = ltd_commons.hard_coded_topology_collection['DoubleTriangle']
n_loops = topology.n_loops
n_loop_lines = len(topology.loop_lines)
graph = [(loop_line.start_node,loop_line.end_node) for loop_line in topology.loop_lines]
spanning_trees = []
spanning_tree({1}, graph, [], spanning_trees)

allowed_systems = find_allowed_systems(spanning_trees)
#for i,allowed_system in enumerate(allowed_systems):
#	print allowed_system

#close_contour = [0]*n_loops
#close_contour = [1]*n_loops
close_contour = [0,0]
residues = evaluate_residues(allowed_systems,close_contour)

for residue in residues:
	print residue

stop

sum_of_all_residues = []
for close_contour in itertools.product(*([[0,1]]*n_loops)):
	sum_of_all_residues += evaluate_residues(allowed_systems,close_contour)
sum_of_all_residues.sort()
sum_of_all_residues = remove_cancelling_residues(sum_of_all_residues)

for residue in sum_of_all_residues:
	print residue

stop

cut_stucture = generate_cut_structure(residues,spanning_trees)

#print get_cut_stucture_string(cut_stucture)




















