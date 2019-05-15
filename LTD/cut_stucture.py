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

topology = ltd_commons.hard_coded_topology_collection['Mercedes']

graph = [(loop_line.start_node,loop_line.end_node) for loop_line in topology.loop_lines]
spanning_trees = []
spanning_tree({1}, graph, [], spanning_trees)
# get all possible momentum bases
momentum_bases= [[i for i in range(len(graph)) if i not in r] for r in spanning_trees]

print momentum_bases

def evaluate_residue(sigmas,allowed_system):
	sign = numpy.prod(sigmas)*numpy.linalg.det(allowed_system[0])
	residue = [(allowed_system[2],sigmas,sign)]
	thetas = []
	for n_iter in xrange(n_loops):
		theta = [0 for i in xrange(n_loops)]
		sub_matrix = numpy.array([[allowed_system[0][i][j] for j in xrange(n_iter+1)] for i in xrange(n_iter+1)])
		if n_iter == 0:
			theta[allowed_system[1][n_iter]] = 1./numpy.linalg.det(sub_matrix)*sigmas[allowed_system[1][n_iter]]
		else:
			for r in xrange(n_iter+1):
				subsub_matrix = numpy.array([[allowed_system[0][i][j] for j in xrange(n_iter)] for i in xrange(n_iter+1) if i != r])
				theta[allowed_system[1][r]] = (-1.)**(n_iter+r+1)*numpy.linalg.det(subsub_matrix)/numpy.linalg.det(sub_matrix)
				theta[allowed_system[1][r]] *= sigmas[allowed_system[1][r]]
		thetas += [theta]
	residue += [thetas]
	return residue

# get transformation matrices: cut_mometum basis <=> loop momentum basis
s_matrices = []
for i,momentum_basis in enumerate(momentum_bases):
	s_matrix = numpy.array([topology.loop_lines[i].signature for i in momentum_basis])#,dtype=float)
	#n_loops = len(s_matrix)
	## get allowed orderings
	#s_matrix.sort(reverse=True,key = lambda s: numpy.linalg.det(s))
	#rows_to_be_permuted = [[i for i in xrange(n_loops) if s_matrix[i][j] != 0] for j in xrange(n_loops)]
	#print rows_to_be_permuted
	s_matrices += [s_matrix]

#print s_matrices
n_loops = len(s_matrices[0])

allowed_systems = []
n_iter = 1
for res_index,s_matrix in enumerate(s_matrices):
	permuted_s_matrices = numpy.array(list(itertools.permutations(s_matrix)))
	permuted_energies = numpy.array(list(itertools.permutations(range(n_loops))))
	#print permuted_s_matrices
	for perm_index,permuted_s_matrix in enumerate(permuted_s_matrices):
		#print "="*20
		#print permuted_s_matrix
		for n_iter in xrange(n_loops):
			sub_matrix = numpy.array([[permuted_s_matrix[i][j] for j in xrange(n_iter+1)] for i in xrange(n_iter+1)])
			#print sub_matrix
			if numpy.linalg.det(sub_matrix) == 0:
				allowed = False
				break
			else:
				allowed = True
		if allowed:
			allowed_systems += [(permuted_s_matrix,permuted_energies[perm_index],res_index)]
		else:
			continue
print "="*20

#for i,allowed_system in enumerate(allowed_systems):
#	print allowed_system

residues = []
for sigmas in itertools.product(*([[1,-1]]*n_loops)):
	for allowed_system in allowed_systems:
		residue = evaluate_residue(sigmas,allowed_system)
		residues += [residue]

residues.sort()#key= lambda row: row[0])
#for residue in residues:
#	print residue

print "="*10


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

for residue in evaluated_theta_resides:
	print residue

surviving_residues = []
for residue in evaluated_theta_resides:
	if residue[1] == []:
		surviving_residues += [residue[0]]

n_ll = len(topology.loop_lines)
cut_stuctures = []
for residue,momentum_basis in zip(surviving_residues,momentum_bases):
	cut_struct_iter = iter(residue[1])
	cut_stuctures += [[next(cut_struct_iter) if i in momentum_basis else 0 for i in xrange(n_ll)]]

topology_cut_structures = cut_stuctures
for i,cut_stucture in enumerate(topology_cut_structures):
	for j,cut in enumerate(cut_stucture):
		if cut == 0:
			topology_cut_structures[i][j] = 'LoopLine.NO_CUT'
		elif cut == -1:
			topology_cut_structures[i][j] = 'LoopLine.NEGATIVE_CUT'
		elif cut == 1:
			topology_cut_structures[i][j] = 'LoopLine.POSITIVE_CUT'

print('(' + ',\n'.join([ '(' + ','.join([cut for cut in cut_stucture]) + ')' for cut_stucture in topology_cut_structures]) +')') 





















