#!/usr/bin/env python2

#PRL_cut_structure.py
import itertools
import numpy

CLOSE_ABOVE = True
CLOSE_BELOW = False

class CutStructureGenerator(object):
	def __init__(self,loop_line_signatures):
		self.loop_line_signatures = loop_line_signatures
		self.n_loops = len(loop_line_signatures[0])
		self.n_loop_lines = len(loop_line_signatures)

	
	def __call__(self,contour_closure):
		assert(len(contour_closure) == self.n_loops)
		residue_elements = self.get_residues(contour_closure)
		cut_stucture = []
		for residue_element in residue_elements:
			cut_struct_iter = iter(residue_element['sigmas'])
			cut_stucture += [[next(cut_struct_iter) if i in residue_element['basis'] else 0 for i in range(self.n_loop_lines)]]
		return cut_stucture

	def get_residues(self,contour_closure,simplify=True):
		assert(len(contour_closure) == self.n_loops)
		residue_elements = []
		spanning_trees = self.get_spanning_tree_generators()
		for sigmas in itertools.product(*([[1, -1]]*self.n_loops)):
			for tree_nr, spanning_tree in enumerate(spanning_trees):
				residues_per_tree = spanning_tree.get_residues_per_tree(sigmas,contour_closure,simplify)
				for residue_element in residues_per_tree:
					residue_element.update({'nr': tree_nr})
					residue_elements += [residue_element]
		residue_elements.sort()
		return residue_elements

	def get_spanning_tree_generators(self):
		reference_signature_matrices = list(itertools.combinations(self.loop_line_signatures, self.n_loops))
		bases = list(itertools.combinations(range(len(self.loop_line_signatures)), self.n_loops))
		spanning_trees = []
		for basis,reference_signature_matrix in zip(bases,reference_signature_matrices):
			if numpy.linalg.det(reference_signature_matrix) == 0:
				continue
			else:
				spanning_trees += [SpanningTreeGenerator(reference_signature_matrix, basis)]
		return spanning_trees

class SpanningTreeGenerator(object):
	def __init__(self, reference_signature_matrix, basis):
		self.n_loops = len(reference_signature_matrix)
		self.reference_signature_matrix = reference_signature_matrix
		self.basis = basis

	def get_residues_per_tree(self,sigmas,contour_closure,simplify=True):
		residues_per_tree = []
		residues = self.get_residue_generators()
		for residue in residues:
			residue_element = residue.get_residue(sigmas,contour_closure,simplify)
			if simplify:
				if residue_element['sign'] != 0.:
					residues_per_tree += [residue_element]
			else:
				residues_per_tree += [residue_element]
		if simplify:
			for residue_element in list(residues_per_tree):
				cancelling_residue = {'thetas': list(residue_element['thetas']), 'sign': -residue_element['sign']}
				if cancelling_residue in list(residues_per_tree):
					residues_per_tree.remove(cancelling_residue)
					residues_per_tree.remove(residue_element)
			if not all(all(theta==[] for theta in residue_element['thetas']) for residue_element in residues_per_tree):
				print('Error: thetas left in residues: %s' % self.residues)		
		for residue_element in residues_per_tree:
			residue_element.update({'sigmas': sigmas, 'basis': self.basis}) 
		return residues_per_tree

	def get_residue_generators(self):
		residues = []
		permuted_signature_matrices = numpy.array(list(itertools.permutations(self.reference_signature_matrix)))
		permutations = numpy.array(list(itertools.permutations(range(self.n_loops))))
		for permuted_signature_matrix,permutation in zip(permuted_signature_matrices,permutations):
			for r in range(self.n_loops):
				sub_matrix = numpy.array([[permuted_signature_matrix[i][j] for j in range(r+1)] for i in range(r+1)])
				if numpy.linalg.det(sub_matrix) == 0:
					allowed = False
					break
				else:
					allowed = True
			if allowed:
				residues += [ResidueGenerator(permuted_signature_matrix,permutation)]
			else:
				continue
		return residues

class ResidueGenerator(object):
	def __init__(self,signature_matrix,permutation):
		self.n_loops = len(signature_matrix)
		self.signature_matrix = signature_matrix
		self.permutation = permutation

	def get_residue(self,sigmas,contour_closure,simplify = True):
		sign = self.get_sign(sigmas,contour_closure)
		thetas, theta_sign = self.get_thetas(sigmas,contour_closure,simplify)
		residue_element = {'sign': sign*theta_sign, 'thetas': thetas}
		return residue_element

	def get_sign(self,sigmas,contour_closure):
		contour_sign = (-1)**contour_closure.count(CLOSE_BELOW)
		cut_struct_sign = numpy.prod(sigmas)
		det_sign = numpy.linalg.det(self.signature_matrix)
		return contour_sign*cut_struct_sign*det_sign

	def get_thetas(self,sigmas,contour_closure,simplify = True):
		thetas = []
		for r in range(self.n_loops):
			theta = ThetaGenerator(self.n_loops)
			sub_matrix = numpy.array([[self.signature_matrix[i][j] for j in range(r+1)] for i in range(r+1)])
			if r == 0:
				basis_index = self.permutation[r]
				theta_value = 1./numpy.linalg.det(sub_matrix) * sigmas[basis_index]
				theta.add(basis_index,theta_value)
			else:
				for m in range(r+1):
					basis_index = self.permutation[m]
					# Leibniz formula: expansion along r th column
					subsub_matrix = numpy.array([[self.signature_matrix[i][j] for j in range(r)] for i in range(r+1) if i != m])
					theta_value = 1./numpy.linalg.det(sub_matrix)
					theta_value *= (-1)**(m+r)*numpy.linalg.det(subsub_matrix)
					theta_value *= sigmas[basis_index]
					theta.add(basis_index,theta_value)
			if contour_closure[r] == CLOSE_ABOVE:
				theta.invert_sign()
			if simplify:
				theta.simplify()
				if theta.sign == 0.:
					return [], 0.
			if theta.arg != []:
				thetas += [theta.arg]
		return thetas, 1.

class ThetaGenerator(object):
	def __init__(self,n_loops):
		self.arg = [None for i in range(n_loops)]
		self.sign = 1.

	def add(self,index,value):
		self.arg[index] = value

	def invert_sign(self):
		self.arg = [-th if th != None else None for th in self.arg]

	def simplify(self):
		if all(th is None for th in self.arg):
			self.arg = []
			self.sign = 0.
		elif all(th <= 0. or th is None for th in self.arg):
			self.arg = []
			self.sign = 0.
		elif all(th >= 0. or th is None for th in self.arg):
			self.arg = []
			self.sign = 1.
		return 

if __name__ == "__main__":
    loop_line_signatures = [(1,0,0),(0,1,0),(0,0,1),(1,-1,0),(-1,0,1),(0,1,-1)] # 3-loop mercedes
    #loop_line_signatures = [(1,0,0),(1,-1,0),(1,-1,-1),(1,-1,0),(0,1,0),(0,0,1)] # 3-loop ladder
    #loop_line_signatures = [(1,0),(0,1),(1,1)] # 2-loop
    #loop_line_signatures = [(1)] # 1-loop

    cut_stucture_generator = CutStructureGenerator(loop_line_signatures)

    for contour_closure in itertools.product(*([[CLOSE_BELOW,CLOSE_ABOVE]]*cut_stucture_generator.n_loops)):
        residues = cut_stucture_generator.get_residues(contour_closure,simplify=True)
        for residue in residues:
            print(residue)

    contour_closure = [CLOSE_ABOVE,CLOSE_BELOW,CLOSE_BELOW]
    cut_structure = cut_stucture_generator(contour_closure)
    print(cut_structure)
