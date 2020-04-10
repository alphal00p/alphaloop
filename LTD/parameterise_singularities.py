from ltd_utils import *
from squared_topologies import *
import warnings
import random
import os, sys
sys.path.insert(0,'../')
import ltd
import matplotlib.pyplot as plt
import itertools
import matplotlib.cm as cm

HEADER = '\033[95m'
OKBLUE = '\033[94m'
OKGREEN = '\033[92m'
WARNING = '\033[93m'
FAIL = '\033[91m'
ENDC = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'

# old: only for reference to cut group intersections
class OneLoopSingularityParameteriser(object):
	def __init__(self, squared_topology, cut_index, l_or_r, cut_momenta):
		self.external_momenta = [squared_topology.external_momenta['q%i'%(i+1)] for i in range(len(squared_topology.external_momenta))]
		self.squared_topology = squared_topology
		self.incoming_momenta_names = squared_topology.incoming_momenta
		self.cut_topology = self.squared_topology.cut_diagrams[cut_index]
		self.one_loop_topology = self.squared_topology.cut_diagrams[cut_index][l_or_r]['loop_topos']
		self.cut_legs = self.squared_topology.cuts[cut_index]
		self.cut_leg_signatures = self.squared_topology.cut_signatures[cut_index]
		self.tolerance = 1e-10
		self.cut_momenta = cut_momenta
		self.propagators = self.one_loop_topology.loop_lines[0].propagators
		self.shifts = [vectors.LorentzVector((numpy.array(propagator.parametric_shift[0]).dot(numpy.array(cut_momenta))
							+ numpy.array(propagator.parametric_shift[1]).dot(numpy.array(self.external_momenta)))) for propagator in self.propagators]
		self.existing_e_surfaces = self.get_existing_e_surfaces()
		self.cut_group_intersections = self.get_cut_group_intersections()
		self.pinched_surfaces = self.get_pinched_surfaces()
		self.e_cm_squared = sum(squared_topology.external_momenta[e][0]
								for e in squared_topology.incoming_momenta)**2 - sum(x*x
								for x in (sum(squared_topology.external_momenta[e][i]
								for e in squared_topology.incoming_momenta)
								for i in range(1, 4)))
		print('E-Surfaces: ',self.existing_e_surfaces)
		print('Cut Group Intersections:', self.cut_group_intersections)
		print('Pinched Surfaces: ', self.pinched_surfaces)

	def get_cut_basis(self,include_loop_momenta=True):
		cut_basis = []
		cut_signatures = []
		cut_shifts = []
		for cut_leg_index, name_and_sign in enumerate(self.cut_legs[:-1]):
			cut_basis += [name_and_sign]
			cut_signatures += [self.cut_leg_signatures[cut_leg_index][0]]
			cut_shifts += [self.cut_leg_signatures[cut_leg_index][1]]
		if include_loop_momenta:
			loop_momentum_index = 1
			for l_or_r in [0,1]:
				if self.cut_topology[l_or_r].n_loops > 0:
					for loop_index, signature in enumerate(self.cut_topology[l_or_r].loop_momentum_map):
						cut_basis += [('k%i'%loop_momentum_index,1)] #random cut sign, result is independent of it
						cut_signatures += [signature[0]]
						cut_shifts += [signature[1]]
						loop_momentum_index += 1
		print(cut_basis)
		trsf_loop_to_cut_basis = numpy.array(cut_signatures)
		if not include_loop_momenta:
			reduced_trsf_loop_to_cut_basis = []
			for column in trsf_loop_to_cut_basis.T:
				#print(column)
				if not all([c==0 for c in column]):
					reduced_trsf_loop_to_cut_basis += [column]
			#print(reduced_trsf_loop_to_cut_basis)
			trsf_loop_to_cut_basis = numpy.array(reduced_trsf_loop_to_cut_basis).T
		inv_trsf_loop_to_cut_basis = numpy.linalg.inv(trsf_loop_to_cut_basis)
		trsf_loop_to_cut_shift = numpy.array(cut_shifts)
		#spatial_external = numpy.array([external_momentum[1:] for external_momentum in self.external_momenta])
		#print(cut_basis)
		#print(trsf_loop_to_cut_shift)
		#print(trsf_loop_to_cut_basis)
		#supergraph_basis = inv_trsf_loop_to_cut_basis.dot(cut_basis-trsf_loop_to_cut_shift.dot(spatial_external))
		return cut_basis, trsf_loop_to_cut_basis, trsf_loop_to_cut_shift

	def get_lorentz_vector(self,time,space):
		return vectors.LorentzVector(numpy.append([time],space))

	def get_space_vector(self,r,unit_vector):
		return r*unit_vector

	def get_unit_vector(self,cos_theta,phi):
		sin_theta = numpy.sqrt(1-cos_theta**2)
		return numpy.array([sin_theta*numpy.cos(phi),sin_theta*numpy.sin(phi),cos_theta])

	def get_beta(self,p):
		assert(p.square()!=0)
		beta = p.space()/p[0]
		assert(beta.dot(beta)<1.)
		return beta

	def get_gamma(self,beta):
		assert(beta.dot(beta)<1.)
		return 1./numpy.sqrt(1.-beta.dot(beta))

	def kaellen_lambda(self,x,y,z):
		return x**2 + y**2 + z**2 - 2.*(x*y + x*z + y*z)

	def get_boost(self,beta):
		assert(beta.dot(beta)<1.)
		assert(beta.dot(beta)>=0.)
		if abs(beta.dot(beta)) < self.tolerance:
			boost = numpy.identity(4)
		else:
			gamma = self.get_gamma(beta)
			row0 = numpy.append([gamma],-gamma*beta)
			space_matrix = numpy.identity(3)+(gamma-1)/beta.dot(beta)*numpy.outer(beta,beta)
			rows1to3 = numpy.append([-gamma*beta],space_matrix,axis=0).T
			boost = numpy.append([row0],rows1to3,axis=0)
		return boost

	def get_cross_matrix(self,vector):
		cross_matrix = numpy.array([[0			, -vector[2], vector[1]	], 
									[vector[2]	, 0			, -vector[0]],
									[-vector[1]	, vector[0]	, 0			]])
		return cross_matrix

	def get_rotation(self,unit_vector):
		e_z = vectors.Vector([0,0,1])
		vec = numpy.cross(unit_vector,e_z)
		cos = unit_vector.dot(e_z)
		cross_matrix = self.get_cross_matrix(vec)
		rot = numpy.identity(3)+cross_matrix+cross_matrix.dot(cross_matrix)/(1.+cos)
		rotation = numpy.append([[1]+3*[0]],numpy.append([3*[0]],rot.T,axis=0).T,axis=0)
		return rotation

	def get_os_energy(self,os_energy_index,qi_radius_squared):
		mi_squared = self.propagators[os_energy_index].m_squared
		return numpy.sqrt(qi_radius_squared+mi_squared)

	def get_os_radius(self,os_energy_index,os_radius_index):
		assert(os_energy_index != os_radius_index)
		mi_squared = self.propagators[os_energy_index].m_squared
		pi = self.shifts[os_energy_index]
		mj_squared = self.propagators[os_radius_index].m_squared
		pj = self.shifts[os_radius_index]
		pji = pj-pi
		if (pji.square() - (numpy.sqrt(mi_squared)+numpy.sqrt(mj_squared))**2) > -self.tolerance and pji[0] < 0:
			# E-surface condition
			if abs(pji.square()) < self.tolerance: #pji.square()==0: # requires both masses to be zero
				#print("Pinched E-surface: Propagator %i of dual %i."%(os_radius_index,os_energy_index))
				qi_radius = -99
			else:
				qi_radius_squared = self.kaellen_lambda(pji.square(),mi_squared,mj_squared)/(4.*pji.square())
				assert(qi_radius_squared>=0.)
				qi_radius = numpy.sqrt(qi_radius_squared) 
		else:
			#print("Propagator %i of dual %i can't go on-shell."%(os_radius_index,os_energy_index))
			qi_radius = None
		return qi_radius

	def get_os_cos_theta(self,os_energy_index,os_radius_index,os_cos_theta_index):
		assert((os_energy_index,os_cos_theta_index) in self.existing_e_surfaces)
		assert(os_energy_index not in [os_radius_index,os_cos_theta_index] and os_radius_index != os_cos_theta_index)
		pi = self.shifts[os_energy_index]
		mi_squared = self.propagators[os_energy_index].m_squared
		pj = self.shifts[os_radius_index]
		mj_squared = self.propagators[os_radius_index].m_squared
		pji = pj - pi
		beta_ji = self.get_beta(pji)
		boost = self.get_boost(beta_ji)
		boost_inv = self.get_boost(-beta_ji)
		pk = self.shifts[os_cos_theta_index]
		mk_squared = self.propagators[os_cos_theta_index].m_squared
		pki = pk - pi
		assert((pji.square() - (numpy.sqrt(mi_squared)+numpy.sqrt(mj_squared))**2) > -self.tolerance and pji[0] < 0.)
		if abs(pki.square()) < self.tolerance:
			if mi_squared-mk_squared == 0.:
				raise ValueError("How can this happen? Pinched surfaces in the same cut group can't intersect.")
				cos_theta = -99
			else:
				#print("Kinematics don't allow propagators %i and %i of dual %i to go on-shell."%(os_energy_index,os_radius_index,os_energy_index))
				cos_theta = None
		else:
			beta_ki = self.get_beta(vectors.LorentzVector(boost.dot(pki)))
			gamma_ki = self.get_gamma(beta_ki)
			qi_prime_radius = self.get_os_radius(os_energy_index,os_radius_index)
			qi_prime_energy = self.get_os_energy(os_energy_index,qi_prime_radius**2)
			pki_squared = pki.square()
			assert(pki_squared>0)
			cos_theta = (mi_squared-mk_squared+pki_squared
						-2.*gamma_ki*qi_prime_energy*numpy.sqrt(pki_squared)
					)/(-2.*gamma_ki*numpy.sqrt(pki_squared)*qi_prime_radius*numpy.sqrt(beta_ki.dot(beta_ki)))
			if not abs(cos_theta) <= 1.:
				cos_theta = None
		return cos_theta

	def get_k_1c_os(self,os_energy_index,qi_space):
		qi_energy = self.get_os_energy(os_energy_index,qi_space.dot(qi_space))
		qi = self.get_lorentz_vector(qi_energy,qi_space)
		pi = self.shifts[os_energy_index]
		k = qi-pi
		return k

	def get_k_2c_os(self,os_energy_index,os_radius_index,qi_prime_unit):
		assert((os_energy_index,os_radius_index) in self.existing_e_surfaces)
		qi_prime_radius = self.get_os_radius(os_energy_index,os_radius_index)
		if qi_prime_radius == None:
			raise ValueError("Kinematics don't allow propagator %i of dual %i to go on-shell."%(os_radius_index,os_energy_index))
		assert(qi_prime_radius>=0.)
		qi_prime_energy = self.get_os_energy(os_energy_index,qi_prime_radius**2)
		qi_prime_space = self.get_space_vector(qi_prime_radius,qi_prime_unit)
		qi_prime = self.get_lorentz_vector(qi_prime_energy,qi_prime_space)
		pi = self.shifts[os_energy_index]
		pj = self.shifts[os_radius_index]
		pji = pj - pi
		assert(pji[0]<0.)
		beta = self.get_beta(pji)
		boost_inv = self.get_boost(-beta)
		k = boost_inv.dot(qi_prime) - pi
		return k

	def get_k_3c_os(self,os_energy_index,os_radius_index,os_cos_theta_index,phi_p2):
		assert((os_energy_index, os_radius_index) in self.existing_e_surfaces)
		assert((os_energy_index,os_cos_theta_index) in self.existing_e_surfaces)
		qi_p2_radius = self.get_os_radius(os_energy_index,os_radius_index)
		if qi_p2_radius == None:
			raise ValueError("Kinematics don't allow propagator %i of dual %i to go on-shell."%(os_radius_index,os_energy_index))
		assert(qi_p2_radius>=0.)
		qi_p2_energy = self.get_os_energy(os_energy_index,qi_p2_radius**2)
		cos_theta_p2 = self.get_os_cos_theta(os_energy_index,os_radius_index,os_cos_theta_index)
		if cos_theta_p2 == None:
			raise ValueError("Kinematics don't allow propagators %i and %i of dual %i to go on-shell."%(os_cos_theta_index,os_radius_index,os_energy_index))
		assert(abs(cos_theta_p2)<=1.)
		qi_p2_unit = self.get_unit_vector(cos_theta_p2,phi_p2)
		qi_p2_space = self.get_space_vector(qi_p2_radius,qi_p2_unit)
		qi_p2 = self.get_lorentz_vector(qi_p2_energy,qi_p2_space)
		pi = self.shifts[os_energy_index]
		pj = self.shifts[os_radius_index]
		pji = pj - pi
		beta_ji = self.get_beta(pji)
		boost = self.get_boost(beta_ji)
		boost_inv = self.get_boost(-beta_ji)
		pk = self.shifts[os_cos_theta_index]
		pki = pk - pi
		beta_ki = self.get_beta(vectors.LorentzVector(boost.dot(pki)))
		rotation = self.get_rotation(beta_ki/(numpy.sqrt(beta_ki.dot(beta_ki))))
		rotation_inv = rotation.T
		k = boost_inv.dot(rotation_inv.dot(qi_p2)) - pi
		return k

	def get_k_pinched(self,os_energy_index,os_radius_index,x):
		mi_squared = self.propagators[os_energy_index].m_squared
		pi = self.shifts[os_energy_index]
		mj_squared = self.propagators[os_radius_index].m_squared
		pj = self.shifts[os_radius_index]
		pji = pj-pi
		assert(abs(pji.square()) < self.tolerance and mi_squared == 0 and mj_squared == 0)
		qi_space = -x*pji.space()
		qi_energy = self.get_os_energy(os_energy_index,qi_space.dot(qi_space))
		qi = self.get_lorentz_vector(qi_energy,qi_space)
		k = qi-pi
		return k

	def get_spatial_loop_momenta(self,k):
		cut_basis, trsf_loop_to_cut_basis, trsf_loop_to_cut_shift = self.get_cut_basis()
		assert('k1' in [name for (name,energy_sign) in cut_basis])
		assert('k2' not in [name for (name,energy_sign) in cut_basis]) #not two loop ready for now
		cut_basis = [cut_momenta[name] for (name,energy_sign) in cut_basis if name in cut_momenta] + [k]
		#print(cut_momenta,k)
		#print(cut_basis)
		#stop
		#cut_basis = numpy.array([self.cut_momenta[0]]+[k])
		#trsf_loop_to_cut_basis = numpy.array([self.cut_leg_signatures[0][0]]+[self.one_loop_topology.loop_momentum_map[0][0]])
		inv_trsf_loop_to_cut_basis = numpy.linalg.inv(trsf_loop_to_cut_basis)
		#trsf_loop_to_cut_shift = numpy.array([self.cut_leg_signatures[0][1]]+[self.one_loop_topology.loop_momentum_map[0][1]])
		cut_basis_shift = -trsf_loop_to_cut_shift.dot(self.external_momenta)
		loop_basis = inv_trsf_loop_to_cut_basis.dot(cut_basis + cut_basis_shift)
		loop_basis = [vectors.LorentzVector(basis) for basis in loop_basis]
		spatial_loop_momenta = [basis.space() for basis in loop_basis]
		#print(spatial_loop_momenta)
		return spatial_loop_momenta

	def get_existing_e_surfaces(self):
		existing_e_surfaces = []
		for dual_index in range(len(self.propagators)): 
			existing_e_surfaces += self.get_existing_e_surfaces_in_dual(dual_index)
		return existing_e_surfaces

	def get_existing_e_surfaces_in_dual(self, dual_index):
		existing_e_surfaces = []
		for prop_index, propagator in enumerate(self.propagators):
			if prop_index != dual_index:
				pi = self.shifts[dual_index]
				pj = self.shifts[prop_index]
				os_radius = self.get_os_radius(dual_index,prop_index)
				if os_radius is not None:
					if os_radius != -99:
						existing_e_surfaces += [(dual_index,prop_index)]
						#print("Dual %i Propagator %i: E-surface"%(dual_index,prop_index))
					#else:
						#print("Dual %i Propagator %i: pinched E-surface"%(dual_index,prop_index))
				#else:
					#print("Dual %i Propagator %i: can't go on-shell"%(dual_index,prop_index))
		return existing_e_surfaces

	def get_pinched_surfaces(self):
		pinched_surfaces = []
		for dual_index in range(len(self.propagators)):
			pinched_surfaces += self.get_pinched_surfaces_in_dual(dual_index)
		return pinched_surfaces

	def get_pinched_surfaces_in_dual(self, dual_index):
		pinched_surfaces = []
		for prop_index, propagator in enumerate(self.propagators):
			if prop_index != dual_index:
				pi = self.shifts[dual_index]
				pj = self.shifts[prop_index]
				os_radius = self.get_os_radius(dual_index,prop_index)
				if os_radius == -99:
					pinched_surfaces += [(dual_index,prop_index)]
		return pinched_surfaces

	def get_cut_group_intersections(self):
		cut_group_intersections = []
		for (dual_index, e_surf_index) in self.existing_e_surfaces:
			cut_group_intersections += self.get_cut_group_intersections_with_e_surface(dual_index,e_surf_index)
		return cut_group_intersections

	def get_cut_group_intersections_with_e_surface(self,dual_index,e_surf_index):
		assert((dual_index,e_surf_index) in self.existing_e_surfaces)
		cut_group_intersections = []
		for prop_index, propagator in enumerate(self.propagators):
			if prop_index not in [dual_index,e_surf_index]:
				if (dual_index,prop_index) in self.existing_e_surfaces:
					os_cos_theta = self.get_os_cos_theta(dual_index,e_surf_index,prop_index)
					if os_cos_theta is not None and os_cos_theta != -99:
						cut_group_intersections += [(dual_index,e_surf_index,prop_index)]
					elif os_cos_theta == -99:
						print('pinched surface')
		return cut_group_intersections

	def evaluate_propagator(self,prop_index,k):
		return (k+self.shifts[prop_index]).square()+self.propagators[prop_index].m_squared

	def check_e_surface_parameterisations(self,unit_vector):
		for (dual_index,e_surf_index) in self.existing_e_surfaces:
			k = self.get_k_2c_os(dual_index,e_surf_index,unit_vector)
			for prop_index in (dual_index,e_surf_index):
				assert(abs(self.evaluate_propagator(prop_index,k)) < self.tolerance)
		return

	def check_cut_group_intersection_parameterisations(self,phi):
		for (dual_index,e_surf_index,intersection_index) in self.cut_group_intersections:
			k = self.get_k_3c_os(dual_index,e_surf_index,intersection_index,phi)
			for prop_index in (dual_index,e_surf_index,intersection_index):
				assert(abs(self.evaluate_propagator(prop_index,k)) < self.tolerance)
		return

	def check_pinched_surface_parameterisations(self,x):
		for (dual_index,pinched_surf_index) in self.pinched_surfaces:
			k = self.get_k_pinched(dual_index,pinched_surf_index,x)
			for prop_index in (dual_index,pinched_surf_index):
				assert(abs(self.evaluate_propagator(prop_index,k)) < self.tolerance)
		return

	def check_parameterisations(self):
		arbitrary_x = 0.5
		arbitrary_phi = 0.5
		arbitrary_unit_vector = self.get_unit_vector(0.7,arbitrary_phi)
		arbitrary_space_vector = self.get_space_vector(2,arbitrary_unit_vector)
		for dual_index in range(len(self.propagators)):
			k = self.get_k_1c_os(dual_index,arbitrary_space_vector)
			assert(abs(self.evaluate_propagator(dual_index,k)) < self.tolerance)
		self.check_e_surface_parameterisations(arbitrary_unit_vector)
		self.check_cut_group_intersection_parameterisations(arbitrary_phi)
		self.check_pinched_surface_parameterisations(arbitrary_x)
		return

class SquaredTopology(object):
	def __init__(self, info, hyperparameters_yaml):
		self.info = info
		info.export('squared_topology.yaml')
		self.cross_section = ltd.CrossSection('squared_topology.yaml', hyperparameters_yaml)
		self.n_cutkosky_cuts = len(self.info.cuts)
		self.n_loops = self.info.topo.n_loops
		self.name = self.info.name
		self.e_cm_squared = sum(self.info.external_momenta[e][0]
								for e in self.info.incoming_momenta)**2 - sum(x*x
								for x in (sum(self.info.external_momenta[e][i]
								for e in self.info.incoming_momenta)
								for i in range(1, 4)))
		self.external_momenta = ['q%i'%(i+1) for i in range(len(self.info.external_momenta))]
		self.external_momenta_num = [vectors.LorentzVector(self.info.external_momenta[name]) for name in self.external_momenta]
		self.incoming_momenta = self.info.incoming_momenta
		self.supergraph_integration_basis = [info.topo.edge_map_lin[e][0] for e in info.topo.loop_momenta]
		self.supergraph_integration_basis_signatures_and_shifts = self.info.topo.get_signature_map()
		self.cutkosky_cuts = [CutkoskyCut(self, index) for index in range(self.n_cutkosky_cuts)]
		self.pv_phase = 1 if self.n_loops % 2==0 else 0
		self.tolerance = 1e-10

		"""
		spatial_loop_momenta = numpy.array([[1.61569199, 0.53957407, 1.04804543], [ 8.04492459, 11.66211465, 14.11645388]])
		tot = 0
		for index in range(len(self.cutkosky_cuts)):
			c = self.cross_section.evaluate_cut_f128(spatial_loop_momenta,index,
						self.cross_section.get_scaling(spatial_loop_momenta,index)[1][0],
						self.cross_section.get_scaling(spatial_loop_momenta,index)[1][1])[0]
			#c2 = self.cross_section.evaluate_cut_f128(spatial_loop_momenta,index,
			#			self.cross_section.get_scaling(spatial_loop_momenta,index)[0][0],
			#			self.cross_section.get_scaling(spatial_loop_momenta,index)[0][1])
			print('cut %i:'%index, c)
			#print(c2)
			tot += c
		print('python sum:', tot)
		print('rust sum:', self.cross_section.evaluate_f128(spatial_loop_momenta)[0])
		stop
		"""

	def get_lorentz_vector(self,time,space):
		return vectors.LorentzVector(numpy.append([time],space))

	def get_space_vector(self,r,unit_vector):
		return r*unit_vector

	def get_unit_vector(self,cos_theta,phi):
		sin_theta = numpy.sqrt(1-cos_theta**2)
		return numpy.array([sin_theta*numpy.cos(phi),sin_theta*numpy.sin(phi),cos_theta])

	def get_beta(self,p):
		assert(p.square()!=0)
		beta = p.space()/p[0]
		assert(beta.dot(beta)<1.)
		return beta

	def get_gamma(self,beta):
		assert(beta.dot(beta)<1.)
		return 1./numpy.sqrt(1.-beta.dot(beta))

	def kaellen_lambda(self,x,y,z):
		return x**2 + y**2 + z**2 - 2.*(x*y + x*z + y*z)

	def get_boost(self,beta):
		assert(beta.dot(beta)<1.)
		assert(beta.dot(beta)>=0.)
		if abs(beta.dot(beta)) < self.tolerance:
			boost = numpy.identity(4)
		else:
			gamma = self.get_gamma(beta)
			row0 = numpy.append([gamma],-gamma*beta)
			space_matrix = numpy.identity(3)+(gamma-1)/beta.dot(beta)*numpy.outer(beta,beta)
			rows1to3 = numpy.append([-gamma*beta],space_matrix,axis=0).T
			boost = numpy.append([row0],rows1to3,axis=0)
		return boost

	def transform_to_integration_basis(self, spatial_loop_momenta, basis):
		trsf = numpy.array([self.supergraph_integration_basis_signatures_and_shifts[name][0] for name in basis])
		inv_trsf = numpy.linalg.inv(trsf)
		shift = numpy.array([self.supergraph_integration_basis_signatures_and_shifts[name][1] for name in basis])
		spatial_supergraph_integration_basis = inv_trsf.dot(spatial_loop_momenta - shift.dot([p.space() for p in self.external_momenta_num]))
		return spatial_supergraph_integration_basis

	def plot_singularity(self, spatial_loop_momenta, spatial_directions, title='', precision='f128', phase=None, simulate_jac_power=0):
		if phase is None:
			phase = self.pv_phase
		assert(len(spatial_loop_momenta)==self.n_loops)
		assert(len(spatial_directions)==self.n_loops)
		if precision == 'f128':
			evaluate = self.cross_section.evaluate_f128
			evaluate_cut = self.cross_section.evaluate_cut_f128
			inv_parameterize = self.cross_section.inv_parameterize_f128
		elif precision == 'f64':
			evaluate = self.cross_section.evaluate
			evaluate_cut = self.cross_section.evaluate_cut
			inv_parameterize = self.cross_section.inv_parameterize
		xx = numpy.linspace(-1e-0,1e-0,int(1e3))
		spatial_loop_momenta_list = [[p+x*spatial_direction*numpy.sqrt(self.e_cm_squared) for p,spatial_direction in zip(spatial_loop_momenta,spatial_directions)] for x in xx]
		yy = [[cutkosky_cut.evaluate_rust(spatial_loop_momenta)[phase] * x**simulate_jac_power 
					for spatial_loop_momenta,x in zip(spatial_loop_momenta_list,xx)] for cutkosky_cut in self.cutkosky_cuts]
		sum_manual = numpy.sum(yy,axis=0)
		sum_rust_f128 = [self.cross_section.evaluate_f128(spatial_loop_momenta)[phase] * x**simulate_jac_power
							for spatial_loop_momenta,x in zip(spatial_loop_momenta_list,xx)]
		sum_rust_f64 = [self.cross_section.evaluate(spatial_loop_momenta)[phase] * x**simulate_jac_power
							for spatial_loop_momenta,x in zip(spatial_loop_momenta_list,xx)]
		vegas_variables = [numpy.array([self.cross_section.inv_parameterize_f128(spatial_loop_momentum,loop_index,self.e_cm_squared)[:3]
							for loop_index, spatial_loop_momentum in enumerate(spatial_loop_momenta)]).flatten()
							for spatial_loop_momenta in spatial_loop_momenta_list]
		jacobians = [numpy.prod([self.cross_section.inv_parameterize_f128(spatial_loop_momentum,loop_index,self.e_cm_squared)[-1]
							for loop_index, spatial_loop_momentum in enumerate(spatial_loop_momenta)])
							for spatial_loop_momenta in spatial_loop_momenta_list]	
		sum_stabilised = [self.cross_section.evaluate_integrand(vegas_variable)[phase] * jacobian * x**simulate_jac_power
							for vegas_variable, jacobian, x in zip(vegas_variables,jacobians,xx)]
		markers = itertools.cycle(('v','<', '^', '>')) 
		colors = cm.rainbow(numpy.linspace(0, 1, self.n_cutkosky_cuts))
		plt.figure(figsize=(11,6))
		for cut_index, color in zip(range(self.n_cutkosky_cuts),colors):
			plt.plot(xx,yy[cut_index],marker=next(markers), color=color, ms=1,alpha=1,ls='None',
						label='cut %i: (%s)'%(cut_index,','.join([cut['edge']+('%+d'%cut['sign'])[0] for cut in self.info.cuts[cut_index]['cuts']])))
		plt.plot(xx,sum_manual,'b',marker='o',label='sum (python)',ms=0.5,alpha=1,ls='None')
		plt.plot(xx,sum_rust_f128,'g',marker='o',label='sum (rust f128)',ms=0.5,alpha=1,ls='None')
		plt.plot(xx,sum_rust_f64,'r',marker='o',label='sum (rust f64)',ms=0.5,alpha=1,ls='None')
		plt.plot(xx,sum_stabilised,'k',marker='o',label='sum (stablised)',ms=0.5,alpha=1,ls='None')
		plt.xlabel(r'$\mu$')
		plt.yscale('symlog')
		plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		momenta_str = ''.join(['\n'+r'$k_{'+'%i'%i+'}=$'
					+ '(%s)'%','.join(['%f'%k_i for k_i in spatial_loop_momentum])
					+ r'$+\mu E_{\mathrm{cm}}$'
					+ '(%s)'%','.join(['%f'%k_i for k_i in spatial_direction])
					for i,spatial_loop_momentum,spatial_direction in zip(range(self.n_loops),spatial_loop_momenta,spatial_directions)])
		simulate_jac_str = '\n multiplied by simulated jacobian factor ' + r'$\mu^%i$'%simulate_jac_power if simulate_jac_power !=0 else ''
		plt.title(self.name + '\n' + title + momenta_str + simulate_jac_str)
		plt.tight_layout()
		#plt.show()
		return

	def plot_uv(self, subgraph_basis, spatial_directions, title='', precision='f128', phase=None, simulate_jac_power=0):
		if phase is None:
			phase = self.pv_phase
		assert(len(spatial_directions)==self.n_loops)
		print(subgraph_basis,self.supergraph_integration_basis)
		print(set(subgraph_basis) <= set(self.supergraph_integration_basis))
		#xx = numpy.linspace(0.8,0.9999,int(1e3))
		#scaling_fct = lambda x: x/(1.-x)
		#scaling_fct_string = r'$\frac{\mu}{1-\mu}$'
		xx = numpy.logspace(1,10,int(1e3))
		scaling_fct = lambda x: x
		scaling_fct_string = r'$\mu$'
		if set(subgraph_basis) <= set(self.supergraph_integration_basis):
			spatial_loop_momenta_list = [[p * numpy.sqrt(self.e_cm_squared) * scaling_fct(x) if name in subgraph_basis
											else p*numpy.sqrt(self.e_cm_squared) 
											for p, name in zip(spatial_directions,self.supergraph_integration_basis)]
											for x in xx]
		else:
			print('Find possible bases and transform accordingly')
			raise
		#print(spatial_loop_momenta_list[0])
		#stop
		yy = [[cutkosky_cut.evaluate_rust(spatial_loop_momenta)[phase] * x**simulate_jac_power 
					for spatial_loop_momenta,x in zip(spatial_loop_momenta_list,xx)] for cutkosky_cut in self.cutkosky_cuts]
		sum_manual = numpy.sum(yy,axis=0)
		sum_rust_f128 = [self.cross_section.evaluate_f128(spatial_loop_momenta)[phase] * x**simulate_jac_power
							for spatial_loop_momenta,x in zip(spatial_loop_momenta_list,xx)]
		sum_rust_f64 = [self.cross_section.evaluate(spatial_loop_momenta)[phase] * x**simulate_jac_power
							for spatial_loop_momenta,x in zip(spatial_loop_momenta_list,xx)]
		vegas_variables = [numpy.array([self.cross_section.inv_parameterize_f128(spatial_loop_momentum,loop_index,self.e_cm_squared)[:3]
							for loop_index, spatial_loop_momentum in enumerate(spatial_loop_momenta)]).flatten()
							for spatial_loop_momenta in spatial_loop_momenta_list]
		jacobians = [numpy.prod([self.cross_section.inv_parameterize_f128(spatial_loop_momentum,loop_index,self.e_cm_squared)[-1]
							for loop_index, spatial_loop_momentum in enumerate(spatial_loop_momenta)])
							for spatial_loop_momenta in spatial_loop_momenta_list]	
		sum_stabilised = [self.cross_section.evaluate_integrand(vegas_variable)[phase] * jacobian * x**simulate_jac_power
							for vegas_variable, jacobian, x in zip(vegas_variables,jacobians,xx)]
		"""
		abs_yy = yy#abs(numpy.array(yy))
		abs_sum_manual = sum_manual#abs(numpy.array(sum_manual))
		abs_sum_rust_f128 = sum_rust_f128#abs_abs(numpy.array(sum_rust_f128))
		abs_sum_rust_f64 = sum_rust_f64#abs(numpy.array(sum_rust_f64))
		abs_sum_stabilised = sum_stabilised#abs(numpy.array(sum_stabilised))
		print('abs cuts', abs_yy)
		print('abs sum python', abs_sum_manual)
		print('abs sum rust 128',abs_sum_rust_f128)
		print('abs sum rust 64',abs_sum_rust_f64)
		print('abs sum rust stabilised',abs_sum_stabilised)
		"""
		markers = itertools.cycle(('v','<', '^', '>')) 
		colors = cm.rainbow(numpy.linspace(0, 1, self.n_cutkosky_cuts))
		plt.figure(figsize=(11,6))
		for cut_index, color in zip(range(self.n_cutkosky_cuts),colors):
			plt.plot(xx,yy[cut_index],marker=next(markers), color=color, ms=1,alpha=1,ls='None',
						label='cut %i: (%s)'%(cut_index,','.join([cut['edge']+('%+d'%cut['sign'])[0] for cut in self.info.cuts[cut_index]['cuts']])))
		plt.plot(xx,sum_manual,'b',marker='o',label='sum (python)',ms=0.5,alpha=1,ls='None')
		plt.plot(xx,sum_rust_f128,'g',marker='o',label='sum (rust f128)',ms=0.5,alpha=1,ls='None')
		plt.plot(xx,sum_rust_f64,'r',marker='o',label='sum (rust f64)',ms=0.5,alpha=1,ls='None')
		plt.plot(xx,sum_stabilised,'k',marker='o',label='sum (stablised)',ms=0.5,alpha=1,ls='None')
		plt.xlabel(r'$\mu$')
		plt.xscale('log')
		remove_nans_and_zeros = lambda x: (numpy.array(x)[~numpy.isnan(x)])[numpy.nonzero((numpy.array(x)[~numpy.isnan(x)]))]
		min_val = numpy.min([numpy.nanmin(numpy.abs(numpy.array(l)[numpy.nonzero(l)].flatten())) for l in [yy,sum_manual,sum_rust_f128,sum_rust_f64,sum_stabilised] if list(remove_nans_and_zeros(l).flatten()) != []])
		plt.yscale('symlog',linthreshy=float('1e%s'%('%e'%min_val)[-3:]))
		#plt.gca().set_yscale('symlog')
		#plt.semilogy()
		plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		momenta_str = ''.join(['\n'+r'$p_{'+'%i'%int(name[1:])+'}=$'
					+ scaling_fct_string
					+ r'$E_{\mathrm{cm}}$'
					+ '(%s)'%','.join(['%f'%k_i for k_i in spatial_direction])
					if name in subgraph_basis 
					else
					'\n'+r'$p_{'+'%i'%int(name[1:])+'}=$'
					+ r'$E_{\mathrm{cm}}$'
					+ '(%s)'%','.join(['%f'%k_i for k_i in spatial_direction])
					for name,spatial_direction in zip(self.supergraph_integration_basis,spatial_directions)])
		simulate_jac_str = '\n multiplied by simulated jacobian factor ' + r'$\mu^%i$'%simulate_jac_power if simulate_jac_power !=0 else ''
		plt.title(self.name + '\n' + title + momenta_str + simulate_jac_str)
		plt.tight_layout()
		# remove the ticks around zero because they overlap
		lables = [i.get_text() for i in plt.gca().get_yticklabels()]
		zero_index = lables.index('$\\mathdefault{0}$')
		plt.setp(plt.gca().get_yticklabels()[zero_index-1], visible=False)  
		plt.setp(plt.gca().get_yticklabels()[zero_index], visible=False)
	
	def probe_cancellation(self, spatial_loop_momenta, spatial_directions, phase=None, simulate_jac_power=0):
		if phase is None:
			phase = self.pv_phase
		max_ranges = numpy.logspace(-1,-5,5)
		curr_max_var = 1e99
		#max_factor = 1000
		passed = True
		for max_range in max_ranges:
			xx = numpy.linspace(-max_range,max_range,100)
			spatial_loop_momenta_list = [[p+x*spatial_direction*numpy.sqrt(self.e_cm_squared) for p,spatial_direction in zip(spatial_loop_momenta,spatial_directions)] for x in xx]
			vegas_variables = [numpy.array([self.cross_section.inv_parameterize_f128(spatial_loop_momentum,loop_index,self.e_cm_squared)[:3]
								for loop_index, spatial_loop_momentum in enumerate(spatial_loop_momenta)]).flatten()
								for spatial_loop_momenta in spatial_loop_momenta_list]
			jacobians = [numpy.prod([self.cross_section.inv_parameterize_f128(spatial_loop_momentum,loop_index,self.e_cm_squared)[-1]
								for loop_index, spatial_loop_momentum in enumerate(spatial_loop_momenta)])
								for spatial_loop_momenta in spatial_loop_momenta_list]	
			sum_stabilised = [self.cross_section.evaluate_integrand(vegas_variable)[phase] * jacobian * x**simulate_jac_power
								for vegas_variable, jacobian, x in zip(vegas_variables,jacobians,xx)]
			sum_stabilised_normalised = sum_stabilised/numpy.sqrt(self.e_cm_squared)
			#samples_without_zeros_and_outliers = [sample for sample in sum_stabilised_normalised if abs(sample) < max_factor]# and sample != 0.]
			mean = numpy.mean(sum_stabilised_normalised)
			var =  numpy.mean([(sample - mean)**2 for sample in sum_stabilised_normalised])
			print(var)
			if curr_max_var >= var:
				curr_max_var = var
			else:
				passed = False
				break
		text = 'Cancellation of singularities: '
		if passed:
			print(text + OKGREEN + u'\u2713' + ENDC)
		else:
			print(text + FAIL + u'\u2717' + ENDC)
		return

class CutkoskyCut(SquaredTopology):
	def __init__(self, squared_topology, index):
		self.squared_topology = squared_topology
		self.index = index
		self.cutkosky_cut_momenta = [cut['edge'] for cut in self.squared_topology.info.cuts[self.index]['cuts']]
		self.cutkosky_cut_energy_signs = [cut['sign'] for cut in self.squared_topology.info.cuts[self.index]['cuts']]
		self.cutkosky_cut_signatures = [cut['signature'][0] for cut in self.squared_topology.info.cuts[self.index]['cuts']]
		self.cutkosky_cut_shifts = [cut['signature'][1] for cut in self.squared_topology.info.cuts[self.index]['cuts']]
		self.subgraphs = [SubGraph(self, 0), SubGraph(self, 1)]
		self.tolerance = 1e-10

	def __str__(self):
		string = 'Cutkosky cut (%s)'%','.join([name+('%+d'%sign)[0] 
					for name, sign in zip(self.cutkosky_cut_momenta, self.cutkosky_cut_energy_signs)])
		return string

	def evaluate_rust(self, spatial_loop_momenta):
		#print(self.squared_topology.cross_section.get_scaling(spatial_loop_momenta,self.index))
		res =  self.squared_topology.cross_section.evaluate_cut_f128(spatial_loop_momenta,self.index,
						self.squared_topology.cross_section.get_scaling(spatial_loop_momenta,self.index)[1][0],
						self.squared_topology.cross_section.get_scaling(spatial_loop_momenta,self.index)[1][1])
		return res

	def get_singular_surfaces(self, cutkosky_cut_momenta_num):
		e_surfaces = []
		pinched_surfaces = []
		for subgraph in self.subgraphs:
			if subgraph.n_loops > 0:
				e_surfaces += subgraph.get_singular_surfaces(cutkosky_cut_momenta_num)['e_surfaces']
				pinched_surfaces += subgraph.get_singular_surfaces(cutkosky_cut_momenta_num)['pinched_surfaces']
		return {'e_surfaces': e_surfaces, 'pinched_surfaces': pinched_surfaces}

	def get_cut_basis(self,include_loop_momenta=True):
		cut_basis = []
		cut_signatures = []
		cut_shifts = []
		for cut_leg_index, name_and_sign in enumerate(zip(self.cutkosky_cut_momenta[:-1],self.cutkosky_cut_energy_signs[:-1])):
			cut_basis += [name_and_sign]
			cut_signatures += [self.cutkosky_cut_signatures[cut_leg_index]]
			cut_shifts += [self.cutkosky_cut_shifts[cut_leg_index]]
		if include_loop_momenta:
			loop_momentum_index = 1
			for subgraph in self.subgraphs:
				if subgraph.n_loops > 0:
					for loop_index, signature in enumerate(subgraph.info.loop_momentum_map):
						cut_basis += [('k%i'%loop_momentum_index,1)] #random cut sign, result is independent of it
						cut_signatures += [signature[0]]
						cut_shifts += [signature[1]]
						loop_momentum_index += 1
		trsf_loop_to_cut_basis = numpy.array(cut_signatures)
		if not include_loop_momenta:
			reduced_trsf_loop_to_cut_basis = []
			for column in trsf_loop_to_cut_basis.T:
				#print(column)
				if not all([c==0 for c in column]):
					reduced_trsf_loop_to_cut_basis += [column]
			#print(reduced_trsf_loop_to_cut_basis)
			trsf_loop_to_cut_basis = numpy.array(reduced_trsf_loop_to_cut_basis).T
		inv_trsf_loop_to_cut_basis = numpy.linalg.inv(trsf_loop_to_cut_basis)
		trsf_loop_to_cut_shift = numpy.array(cut_shifts)
		return cut_basis, trsf_loop_to_cut_basis, trsf_loop_to_cut_shift

	def get_PS_point(self,random_variables):
		counter = 0
		assert(len(random_variables)==3*(len(self.cutkosky_cut_momenta)-1)-1)
		# see last cut momentum as a dual propagator that is put on-shell by solving for a radius
		# express dual propagator in cut basis
		cut_basis, trsf_loop_to_cut_basis, trsf_loop_to_cut_shift = self.get_cut_basis(include_loop_momenta=True)
		prop_name = self.cutkosky_cut_momenta[-1]
		signature_loop_basis = numpy.array(self.cutkosky_cut_signatures[-1])
		shift_loop_basis = numpy.array(self.cutkosky_cut_shifts[-1])
		inv_trsf_loop_to_cut_basis = numpy.linalg.inv(trsf_loop_to_cut_basis)
		signature_cut_basis = signature_loop_basis.dot(inv_trsf_loop_to_cut_basis)
		shift_cut_basis = shift_loop_basis-signature_cut_basis.dot(trsf_loop_to_cut_shift)
		# count number of cut energies the propagator depends on
		n_independent_cuts = len([1 for sign in signature_cut_basis if sign != 0])
		assert(n_independent_cuts==len(self.cutkosky_cut_momenta)-1)
		os_cut_basis = [vectors.LorentzVector([0,0,0,0]) for i in range(len(cut_basis))]
		indep_counter = 0
		for index,((name,energy_sign),signature) in enumerate(zip(cut_basis,signature_cut_basis)):
			if signature != 0:
				indep_counter += 1
				mass_squared_basis = self.squared_topology.info.masses[name] if name in self.squared_topology.info.masses else 0.
				mass_prop = self.squared_topology.info.masses[prop_name] if prop_name in self.squared_topology.info.masses else 0.
				mass_prop += numpy.sum([self.squared_topology.info.masses
												if (name in self.squared_topology.info.masses and sign != 0) else 0.
												for sign,((name,_)) in zip(signature_cut_basis[:index],cut_basis[:index])])
				mass_squared_prop = mass_prop**2
				shift = vectors.LorentzVector(shift_cut_basis.dot(self.squared_topology.external_momenta_num))+numpy.sum([sign*q 
									for sign, q in zip(signature_cut_basis[:index],os_cut_basis[:index]) if sign != 0],axis=0)
				assert(shift.square() > (numpy.sqrt(mass_squared_prop)+numpy.sqrt(mass_squared_basis))**2)
				assert(shift[0]*energy_sign*signature < 0)
				radius_max = numpy.sqrt(self.kaellen_lambda(shift.square(),mass_squared_prop,mass_squared_basis)/(4.*shift.square()))
				if indep_counter == n_independent_cuts:
					radius = radius_max
				else:
					radius = random_variables[counter]*radius_max
					counter += 1
				cos_theta = 2*random_variables[counter]-1.
				counter += 1
				phi = 2*numpy.pi*random_variables[counter]
				counter += 1
				q_rest = self.get_lorentz_vector(energy_sign*numpy.sqrt(radius**2 + mass_squared_basis),
											self.get_space_vector(radius,
												self.get_unit_vector(cos_theta,phi)))
				q = vectors.LorentzVector(self.get_boost(-self.get_beta(shift)).dot(q_rest))
				os_cut_basis[index] = q
		os_cut_basis = numpy.array(os_cut_basis)
		os_prop_momentum = vectors.LorentzVector(signature_cut_basis.dot(os_cut_basis)+shift_cut_basis.dot(self.squared_topology.external_momenta_num))
		# assert momentum conservation
		incoming_momenta = [vectors.LorentzVector(momentum) for name,momentum in self.squared_topology.info.external_momenta.items() if name in self.squared_topology.info.incoming_momenta]
		all_cut_momenta = {name: vectors.LorentzVector(momentum) for (name,energy_sign),momentum in zip(cut_basis,os_cut_basis) if name in self.cutkosky_cut_momenta}
		all_cut_momenta.update({prop_name: os_prop_momentum})
		outgoing_momenta = [energy_sign*all_cut_momenta[name] for energy_sign, name in zip(self.cutkosky_cut_energy_signs, self.cutkosky_cut_momenta)]
		assert(all(numpy.zeros(4) == numpy.sum(incoming_momenta,axis=0)-numpy.sum(outgoing_momenta,axis=0)))
		return all_cut_momenta

class SubGraph(CutkoskyCut):
	def __init__(self, cutkosky_cut, left_or_right):
		self.cutkosky_cut = cutkosky_cut
		self.left_or_right = left_or_right
		self.info = self.cutkosky_cut.squared_topology.info.cut_diagrams[self.cutkosky_cut.index][0]['loop_topos'][self.left_or_right]
		self.n_loops = self.info.n_loops
		if self.n_loops > 0:
			self.ltd_cut_structures = [LTDCutStructure(self, ltd_cut_structure) for ltd_cut_structure in self.info.ltd_cut_structure]
			# supergraph integration basis to supergraph cut basis
			self.supergraph_basis_trsf = numpy.array([cutkosky_cut_signature for cutkosky_cut_signature in self.cutkosky_cut.cutkosky_cut_signatures[:-1]]
													+ [signature for (signature,shift) in self.info.loop_momentum_map])
			self.supergraph_shift_trsf = numpy.array([cutkosky_cut_shift for cutkosky_cut_shift in self.cutkosky_cut.cutkosky_cut_shifts[:-1]]
													+ [shift for (signature,shift) in self.info.loop_momentum_map])
			# if the matrix is not invertible it is because there's an intependent loop momentum in the other subgraph
			# although the basis trsf won't be dependent on it, it is necessary to invert the matrix
			# ideally the renduant loop momenta should be stripped away, leaving an invertible submatrix
			self.supergraph_inv_basis_trsf = numpy.linalg.inv(self.supergraph_basis_trsf)
			self.loop_basis = self.get_loop_basis()
		else:
			print('TODO: Complete the tree graph case!')
		
	def get_loop_basis(self):
		loop_basis = []
		index = 0
		for loop_line in self.info.loop_lines:
			for propagator in loop_line.propagators:
				basis_signature = numpy.zeros(self.n_loops)
				basis_signature[index] = 1
				if (all(a==0 for a in propagator.parametric_shift[0]) 
						and all(a==0 for a in propagator.parametric_shift[1])
						and all(s1==s2 for s1,s2 in zip(propagator.signature,basis_signature))):
					loop_basis += [propagator.name]
					index += 1
					if index == self.n_loops:
						break
			if index == self.n_loops:
				break
		return loop_basis

	def get_singular_surfaces(self,cutkosky_cut_momenta_num):
		e_surfaces = []
		pinched_surfaces = []
		for ltd_cut_structure in self.ltd_cut_structures:
			e_surfaces += ltd_cut_structure.get_singular_surfaces(cutkosky_cut_momenta_num)['e_surfaces']
			pinched_surfaces += ltd_cut_structure.get_singular_surfaces(cutkosky_cut_momenta_num)['pinched_surfaces']
		return {'e_surfaces': e_surfaces, 'pinched_surfaces': pinched_surfaces}

class LTDCutStructure(SubGraph):
	def __init__(self, subgraph, ltd_cut_structure):
		self.subgraph = subgraph
		self.ltd_cut_structure = ltd_cut_structure
		self.ltd_cut_signature = [line for line in self.ltd_cut_structure if abs(line) == 1]
		self.basis_trsf = numpy.array([loop_line.signature
									for loop_line,line in zip(self.subgraph.info.loop_lines,self.ltd_cut_structure)
									if abs(line)==1]) #cut line momenta from loop momenta
		self.inv_basis_trsf = numpy.linalg.inv(self.basis_trsf) #loop momenta from cut line momenta
		ltd_cut_propagator_indices = self.get_prop_cuts()
		self.ltd_cuts = [LTDCut(self, ltd_cut_indices) for ltd_cut_indices in ltd_cut_propagator_indices]
		
	def get_prop_cuts(self):
		prop_cuts = []
		cut_indices = [i for i,line in enumerate(self.ltd_cut_structure) if abs(line)==1]
		possible_cuts = [range(len(loop_line.propagators)) for loop_line,line in zip(self.subgraph.info.loop_lines,self.ltd_cut_structure) if abs(line)==1]
		for prop_cut in itertools.product(*possible_cuts):
			prop_cuts += [prop_cut]
		return prop_cuts

	def get_singular_surfaces(self,cutkosky_cut_momenta_num):
		e_surfaces = []
		pinched_surfaces = []
		for ltd_cut in self.ltd_cuts:
			e_surfaces += ltd_cut.get_singular_surfaces(cutkosky_cut_momenta_num)['e_surfaces']
			pinched_surfaces += ltd_cut.get_singular_surfaces(cutkosky_cut_momenta_num)['pinched_surfaces']
		return {'e_surfaces': e_surfaces, 'pinched_surfaces': pinched_surfaces}

class LTDCut(LTDCutStructure):
	def __init__(self, cut_structure, ltd_cut_indices):
		self.cut_structure = cut_structure
		self.ltd_cut_indices = ltd_cut_indices
		ltd_cut_structure = self.cut_structure.ltd_cut_structure
		cut_indices_iter = iter(self.ltd_cut_indices)
		self.ltd_cut_momenta = [loop_line.propagators[next(cut_indices_iter)].name for line, loop_line in zip(self.cut_structure.ltd_cut_structure,self.cut_structure.subgraph.info.loop_lines) if abs(line) != 0]
		# shifts of the LTD cut basis wrt to the cutkosky cuts and external momenta
		cut_indices_iter = iter(self.ltd_cut_indices)
		self.cutkosky_shift_trsf = numpy.array([numpy.array(loop_line.propagators[next(cut_indices_iter)].parametric_shift[0]) for line, loop_line in zip(self.cut_structure.ltd_cut_structure,self.cut_structure.subgraph.info.loop_lines) if abs(line) != 0])
		cut_indices_iter = iter(self.ltd_cut_indices)
		self.external_shift_trsf = numpy.array([numpy.array(loop_line.propagators[next(cut_indices_iter)].parametric_shift[1]) for line, loop_line in zip(self.cut_structure.ltd_cut_structure,self.cut_structure.subgraph.info.loop_lines) if abs(line) != 0])
		self.supergraph_cut_basis = [name for name in self.cut_structure.subgraph.cutkosky_cut.cutkosky_cut_momenta[:-1]]+[name for name in self.ltd_cut_momenta]
		self.ltd_loop_lines = [LTDLoopLine(self, i) for i,loop_line in enumerate(self.cut_structure.subgraph.info.loop_lines)]
		
	def __str__(self):
		string = ('LTD cut (%s)'%','.join([name+('%+d'%sign)[0] 
					for name, sign in zip(self.ltd_cut_momenta,self.cut_structure.ltd_cut_signature)])
					+ ' in ' + str(self.cut_structure.subgraph.cutkosky_cut))
		return string

	def get_subgraph_basis_shift(self, cutkosky_cut_momenta_num):
		return self.cutkosky_shift_trsf.dot(cutkosky_cut_momenta_num)+self.external_shift_trsf.dot(self.cut_structure.subgraph.cutkosky_cut.squared_topology.external_momenta_num)

	def get_singular_surfaces(self,cutkosky_cut_momenta_num):
		e_surfaces = []
		pinched_surfaces = []
		for ltd_loop_line in self.ltd_loop_lines:
			if not ltd_loop_line.is_tree_line:
				e_surfaces += ltd_loop_line.get_singular_surfaces(cutkosky_cut_momenta_num)['e_surfaces']
				pinched_surfaces += ltd_loop_line.get_singular_surfaces(cutkosky_cut_momenta_num)['pinched_surfaces']
		return {'e_surfaces': e_surfaces, 'pinched_surfaces': pinched_surfaces}

class LTDLoopLine(LTDCut):
	def __init__(self, ltd_cut, loop_line_index):
		self.ltd_cut = ltd_cut
		self.loop_line_index = loop_line_index
		this_loop_line = self.ltd_cut.cut_structure.subgraph.info.loop_lines[self.loop_line_index]
		self.is_tree_line = all([sign == 0 for sign in this_loop_line.signature])
		# signature in cut basis
		self.basis_signature = numpy.array(this_loop_line.signature).dot(self.ltd_cut.cut_structure.inv_basis_trsf)
		# signature when multiplied with cut structure signs
		self.energy_signature = self.basis_signature.dot(numpy.diag(self.ltd_cut.cut_structure.ltd_cut_signature))
		self.has_positive_signature = all(sign >= 0 for sign in self.energy_signature)
		self.has_negative_signature = all(sign <= 0 for sign in self.energy_signature)
		self.ltd_propagators = [LTDPropagator(self, a) for a,propagator in enumerate(this_loop_line.propagators)]

	def get_singular_surfaces(self,cutkosky_cut_momenta_num):
		e_surfaces = []
		pinched_surfaces = []
		for ltd_propagator in self.ltd_propagators:
			if ltd_propagator.is_e_surface(cutkosky_cut_momenta_num):
				e_surfaces += [ltd_propagator]
			if ltd_propagator.is_pinched_surface(cutkosky_cut_momenta_num):
				pinched_surfaces += [ltd_propagator]
		return {'e_surfaces': e_surfaces, 'pinched_surfaces': pinched_surfaces}

class LTDPropagator(LTDLoopLine):
	def __init__(self, ltd_loop_line, propagator_index):
		self.ltd_loop_line = ltd_loop_line
		self.propagator_index = propagator_index
		self.propagator_name = self.ltd_loop_line.ltd_cut.cut_structure.subgraph.info.loop_lines[self.ltd_loop_line.loop_line_index].propagators[self.propagator_index].name
		self.propagator_mass = self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.info.masses[self.propagator_name] if self.propagator_name in self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.info.masses else 0.
		cut_indices_iter = iter(self.ltd_loop_line.ltd_cut.ltd_cut_indices)
		full_cut_indices = [next(cut_indices_iter) if abs(line) != 0 else None for line in self.ltd_loop_line.ltd_cut.cut_structure.ltd_cut_structure]
		self.is_cut = (full_cut_indices[self.ltd_loop_line.loop_line_index] == self.propagator_index)
		self.e_surface_dimension = numpy.sum([1 for sign in self.ltd_loop_line.energy_signature if sign != 0])*3 - 1
		self.pinched_surface_dimension = numpy.sum([1 for sign in self.ltd_loop_line.energy_signature if sign != 0])
		#shift of propagator momentum with respect to the loop line momentum
		self.own_cutkosky_shift_trsf = self.ltd_loop_line.ltd_cut.cut_structure.subgraph.info.loop_lines[self.ltd_loop_line.loop_line_index].propagators[self.propagator_index].parametric_shift[0]
		self.own_external_shift_trsf = self.ltd_loop_line.ltd_cut.cut_structure.subgraph.info.loop_lines[self.ltd_loop_line.loop_line_index].propagators[self.propagator_index].parametric_shift[1]
		#total shift i.e. affine surface term
		self.total_cutkosky_shift_trsf = self.own_cutkosky_shift_trsf - self.ltd_loop_line.basis_signature.dot(self.ltd_loop_line.ltd_cut.cutkosky_shift_trsf)
		self.total_external_shift_trsf = self.own_external_shift_trsf - self.ltd_loop_line.basis_signature.dot(self.ltd_loop_line.ltd_cut.external_shift_trsf)
		self.tolerance = 1e-10
		self.supergraph_integration_basis_signature = numpy.array(self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.info.topo.get_signature_map()[self.propagator_name][0])
		self.supergraph_integration_basis_shift = numpy.array(self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.info.topo.get_signature_map()[self.propagator_name][1])

	def __str__(self):
		string = ('propagator %s'%self.propagator_name
				+ ' in ' + str(self.ltd_loop_line.ltd_cut))
		return string

	def get_total_shift(self, cutkosky_cut_momenta_num):
		total_cutkosky_shift = vectors.LorentzVector(self.total_cutkosky_shift_trsf.dot(cutkosky_cut_momenta_num))
		total_external_shift = vectors.LorentzVector(self.total_external_shift_trsf.dot(self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.external_momenta_num))
		return total_cutkosky_shift+total_external_shift
	
	def is_e_surface(self, cutkosky_cut_momenta_num):
		total_shift = self.get_total_shift(cutkosky_cut_momenta_num)
		s = total_shift.square()
		basis_masses = [abs(sign)*self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.info.masses[name] 
							if name in self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.info.masses else 0.
							for sign, name in zip(self.ltd_loop_line.basis_signature,self.ltd_loop_line.ltd_cut.ltd_cut_momenta)]
		tot_mass_squared = (numpy.sum(basis_masses)+self.propagator_mass)**2
		is_massless = (abs(tot_mass_squared) < self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.e_cm_squared*self.tolerance)
		if is_massless:
			is_causal = (abs(s) > self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.e_cm_squared*self.tolerance)
		else:
			is_causal = (s > tot_mass_squared)
		is_plus_e_surface = (self.ltd_loop_line.has_positive_signature and total_shift[0] < 0.) and is_causal
		is_minus_e_surface = (self.ltd_loop_line.has_negative_signature and total_shift[0] > 0.) and is_causal
		is_e_surface = (is_plus_e_surface or is_minus_e_surface)
		return is_e_surface

	def is_pinched_surface(self, cutkosky_cut_momenta_num):
		total_shift = self.get_total_shift(cutkosky_cut_momenta_num)
		s = total_shift.square()
		basis_masses = [abs(sign)*self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.info.masses[name] 
							if name in self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.info.masses else 0.
							for sign, name in zip(self.ltd_loop_line.basis_signature,self.ltd_loop_line.ltd_cut.ltd_cut_momenta)]
		tot_mass_squared = (numpy.sum(basis_masses)+self.propagator_mass)**2
		is_massless = abs(tot_mass_squared) < self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.e_cm_squared*self.tolerance
		is_plus_pinched_surface = (self.ltd_loop_line.has_positive_signature and total_shift[0] < 0.) and (abs(s)<self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.e_cm_squared*self.tolerance)
		is_minus_pinched_surface = (self.ltd_loop_line.has_negative_signature and total_shift[0] > 0.) and (abs(s)<self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.e_cm_squared*self.tolerance)
		is_pinched_surface = is_massless and (is_plus_pinched_surface or is_minus_pinched_surface)
		return is_pinched_surface

	def get_os_supergraph_integration_basis(self, cutkosky_cut_momenta_num, random_variables):
		assert(self.is_e_surface(cutkosky_cut_momenta_num) or self.is_pinched_surface(cutkosky_cut_momenta_num))
		os_subgraph_cut_basis = self.get_os_subgraph_cut_basis(cutkosky_cut_momenta_num, random_variables)
		os_subgraph_integration_basis = self.ltd_loop_line.ltd_cut.cut_structure.inv_basis_trsf.dot(os_subgraph_cut_basis - self.ltd_loop_line.ltd_cut.get_subgraph_basis_shift(cutkosky_cut_momenta_num))
		os_supergraph_cut_basis = numpy.array([momentum for momentum in cutkosky_cut_momenta_num[:-1]]+[momentum for momentum in os_subgraph_integration_basis])
		os_supergraph_integration_basis = self.ltd_loop_line.ltd_cut.cut_structure.subgraph.supergraph_inv_basis_trsf.dot(os_supergraph_cut_basis - self.ltd_loop_line.ltd_cut.cut_structure.subgraph.supergraph_shift_trsf.dot(self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.external_momenta_num))
		os_supergraph_integration_basis = [vectors.LorentzVector(p) for p in os_supergraph_integration_basis]
		assert(self.check_parameterisation(os_supergraph_integration_basis))
		return os_supergraph_integration_basis

	def get_os_subgraph_cut_basis(self, cutkosky_cut_momenta_num, random_variables):
		if self.is_e_surface(cutkosky_cut_momenta_num):
			os_subgraph_cut_basis = [vectors.LorentzVector([None,None,None,None]) for i in range(len(self.ltd_loop_line.ltd_cut.ltd_cut_momenta))]
			n_independent_cuts = len([1 for sign in self.ltd_loop_line.basis_signature if sign != 0])
			counter = 0
			indep_counter = 0
			for index,(name,energy_sign,signature) in enumerate(zip(self.ltd_loop_line.ltd_cut.ltd_cut_momenta, self.ltd_loop_line.ltd_cut.cut_structure.ltd_cut_signature, self.ltd_loop_line.basis_signature)):
				if signature != 0:
					indep_counter += 1
					mass_squared_basis = self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.info.masses[name] if name in self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.info.masses else 0.
					mass = self.propagator_mass
					mass += numpy.sum([self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.info.masses[name]
												if (name in self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.info.masses and sign != 0) else 0.
												for sign,((name,_)) in zip(self.ltd_loop_line.basis_signature[:index],self.ltd_loop_line.ltd_cut.ltd_cut_momenta[:index])])
					mass_squared = mass**2
					shift = self.get_total_shift(cutkosky_cut_momenta_num)
					shift += self.ltd_loop_line.basis_signature[:index].dot(os_subgraph_cut_basis[:index])
					assert(shift.square() > (numpy.sqrt(mass_squared)+numpy.sqrt(mass_squared_basis))**2)
					assert(shift[0]*energy_sign*signature < 0)
					radius_max = numpy.sqrt(self.kaellen_lambda(shift.square(),mass_squared,mass_squared_basis)/(4.*shift.square()))
					if indep_counter == n_independent_cuts:
						radius = radius_max
					else:
						radius = random_variables[counter]*radius_max
						counter += 1
					cos_theta = 2*random_variables[counter]-1.
					counter += 1
					phi = 2*numpy.pi*random_variables[counter]
					counter += 1
					basis_momentum_boosted = self.get_lorentz_vector(energy_sign*numpy.sqrt(radius**2 + mass_squared_basis),
												self.get_space_vector(radius, self.get_unit_vector(cos_theta,phi)))
					basis_momentum = vectors.LorentzVector(self.get_boost(-self.get_beta(shift)).dot(basis_momentum_boosted))
				else:
					mass_squared_basis = self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.info.masses[name] if name in self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.info.masses else 0.
					radius = numpy.sqrt(self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.e_cm_squared)*random_variables[counter]/(1.-random_variables[counter])
					counter += 1
					cos_theta = 2*random_variables[counter]-1.
					counter += 1
					phi = 2*numpy.pi*random_variables[counter]
					counter += 1
					basis_momentum = self.get_lorentz_vector(energy_sign*numpy.sqrt(radius**2 + mass_squared_basis),
												self.get_space_vector(radius, self.get_unit_vector(cos_theta,phi)))
				os_subgraph_cut_basis[index] = basis_momentum
		elif self.is_pinched_surface(cutkosky_cut_momenta_num):
			os_subgraph_cut_basis = [vectors.LorentzVector([None,None,None,None]) for i in range(len(self.ltd_loop_line.ltd_cut.ltd_cut_momenta))]
			n_independent_cuts = len([1 for sign in self.ltd_loop_line.basis_signature if sign != 0])
			counter = 0
			indep_counter = 0
			for index,(name,energy_sign,signature) in enumerate(zip(self.ltd_loop_line.ltd_cut.ltd_cut_momenta, self.ltd_loop_line.ltd_cut.cut_structure.ltd_cut_signature, self.ltd_loop_line.basis_signature)):
				if signature != 0:
					indep_counter += 1
					mass_squared_basis = self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.info.masses[name] if name in self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.info.masses else 0.
					mass = self.propagator_mass
					mass += numpy.sum([self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.info.masses[name]
												if (name in self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.info.masses and sign != 0) else 0.
												for sign,((name,_)) in zip(self.ltd_loop_line.basis_signature[:index],self.ltd_loop_line.ltd_cut.ltd_cut_momenta[:index])])
					mass_squared = mass**2
					assert(mass_squared == 0 and mass_squared_basis == 0)
					shift = self.get_total_shift(cutkosky_cut_momenta_num)
					shift += self.ltd_loop_line.basis_signature[:index].dot(os_subgraph_cut_basis[:index])
					assert(abs(shift.square()) < self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.e_cm_squared*self.tolerance)
					shift_energy_sign = 1 if shift[0] > 0 else -1
					spatial_basis_momentum = shift_energy_sign*energy_sign*random_variables[counter]*shift.space()
					counter += 1
					basis_momentum = self.get_lorentz_vector(energy_sign*numpy.sqrt(spatial_basis_momentum.square()),spatial_basis_momentum)
				else:
					mass_squared_basis = self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.info.masses[name] if name in self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.info.masses else 0.
					#if random_variables[counter] == 1:
					#	random = 0.5
					#	radius = numpy.sqrt(self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.e_cm_squared)*random/(1.-random)
					#else:
					radius = numpy.sqrt(self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.e_cm_squared)*random_variables[counter]/(1.-random_variables[counter])
					counter += 1
					cos_theta = 2*random_variables[counter]-1.
					counter += 1
					phi = 2*numpy.pi*random_variables[counter]
					counter += 1
					basis_momentum = self.get_lorentz_vector(energy_sign*numpy.sqrt(radius**2 + mass_squared_basis),
												self.get_space_vector(radius, self.get_unit_vector(cos_theta,phi)))
				os_subgraph_cut_basis[index] = basis_momentum
		else:
			raise ValueError("This propagator cannot go on-shell.")
		return os_subgraph_cut_basis	

	def check_parameterisation(self, os_supergraph_integration_basis):
		passed = True
		for this_propagator_name in (self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.cutkosky_cut_momenta+self.ltd_loop_line.ltd_cut.ltd_cut_momenta+[self.propagator_name]):
			signature = numpy.array(self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.supergraph_integration_basis_signatures_and_shifts[this_propagator_name][0])
			shift = numpy.array(self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.supergraph_integration_basis_signatures_and_shifts[this_propagator_name][1])
			momentum = vectors.LorentzVector(signature.dot(os_supergraph_integration_basis)+shift.dot(self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.external_momenta_num))
			mass = self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.info.masses[this_propagator_name] if this_propagator_name in self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.info.masses else 0.
			inv_prop = momentum.square() - mass**2
			if not (abs(inv_prop) < self.ltd_loop_line.ltd_cut.cut_structure.subgraph.cutkosky_cut.squared_topology.e_cm_squared*self.tolerance):
				passed = False
		return passed

if __name__ == "__main__":

	random.seed(1)
	#print(dir(ltd.CrossSection))
	#stop

	# the double triangle master graph
	t1 = SquaredTopologyGenerator([('q1', 0, 1), ('p1', 1, 2), ('p2', 2, 3), ('p3', 4, 3), ('p4', 4, 1), ('p5', 2, 4), ('q2', 3, 5)], "T", ['q1'], 2,
		{'q1': [1., 0., 0., 0.], 'q2': [1., 0., 0., 0.]},
		particle_ids={'p%s' % i: i for i in range(9)})
		#masses={'p1': 100, 'p2':100, 'p3': 100, 'p4': 100, 'p5': 100})

	var_two_to_two = SquaredTopologyGenerator([('q1', 0, 1), ('q2', 8, 2), ('q3', 5, 7), ('q4', 4, 9), ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4), ('p4', 4, 5),
		('p5', 5, 6), ('p6', 6, 1), ('p7', 6, 3), ], "var_two_to_two", ['q1', 'q2'], 2,
		{'q1': [1., 0., 0., .1], 'q2': [1., 0., 0., -.1], 'q3': [1., 0., 0., .1], 'q4': [1., 0., 0., -.1]})
		#masses={'p1': .1, 'p2':.1, 'p3': .1, 'p4': .1, 'p5': .1, 'p6': .1, 'p7': .1}, loop_momenta_names=('p1', 'p7'),)
	
	mercedes = SquaredTopologyGenerator([('q1', 0, 1), ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 6),
		('p4', 6, 5), ('p5', 5, 1), ('p6', 2, 4), ('p7', 3, 4), ('p8', 4, 5), ('q2', 6, 7)], "M", ['q1'], 2,
		{'q1': [1., 0., 0., 0.], 'q2': [1., 0., 0., 0.]},
		loop_momenta_names=('p1', 'p2', 'p3'),
		particle_ids={'p%s' % i: i for i in range(9)})

	ee_to_dd_2l_doubletriangle = SquaredTopologyGenerator([('q1', 101, 1), ('q2', 102, 1), ('q3', 6, 103), ('q4', 6, 104), ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 5), ('p4', 5, 6),
    ('p5',5, 4), ('p6', 4, 2), ('p7', 4, 3)],
        "ee_to_dd_2l_doubletriangle", ['q1', 'q2'], 2, {'q1': [1., 0., 0., 1.], 'q2': [1., 0., 0., -1.], 'q3': [1., 0., 0., 1.], 'q4': [1., 0., 0., -1.]},
        loop_momenta_names=('p2', 'p3'),
        particle_ids={'p%s' % i: i for i in range(8)},
        overall_numerator=1.0,
#        cut_filter={('p3', 'p5')},
        numerator_structure={('p2', 'p5', 'p7'):
            { (): # uv structure
            [[[0,4],[0.,+2.954161761482786e8]],
            [[0,4,4],[0.,+1.477080880741393e8]],
            [[0,7,7],[0.,-1.477080880741393e8]],
            [[1,5],[0.,-2.954161761482786e8]],
            [[1,4,5],[0.,-1.477080880741393e8]],
            [[2,6],[0.,-2.954161761482786e8]],
            [[2,4,6],[0.,-1.477080880741393e8]],
            [[3,7],[0.,-2.954161761482786e8]],
            [[0,0,4],[0.,-1.477080880741393e8]],
            [[0,0,4,4],[0.,-7.385404403706966e7]],
            [[0,0,5,5],[0.,+7.385404403706966e7]],
            [[0,0,6,6],[0.,+7.385404403706966e7]],
            [[0,0,7,7],[0.,+7.385404403706966e7]],
            [[0,1,5],[0.,+1.477080880741393e8]],
            [[0,1,4,5],[0.,+1.477080880741393e8]],
            [[0,2,6],[0.,+1.477080880741393e8]],
            [[0,2,4,6],[0.,+1.477080880741393e8]],
            [[1,1,5,5],[0.,-1.477080880741393e8]],
            [[1,2,5,6],[0.,-2.954161761482786e8]],
            [[1,3,5,7],[0.,-1.477080880741393e8]],
            [[2,2,6,6],[0.,-1.477080880741393e8]],
            [[2,3,6,7],[0.,-1.477080880741393e8]],
            [[3,3,4],[0.,+1.477080880741393e8]],
            [[3,3,4,4],[0.,+7.385404403706966e7]],
            [[3,3,5,5],[0.,-7.385404403706966e7]],
            [[3,3,6,6],[0.,-7.385404403706966e7]],
            [[3,3,7,7],[0.,-7.385404403706966e7]]]
            },
            ('p3', 'p6', 'p7'):
            { (): # uv structure
            [[[0,4],[0.,+2.954161761482786e8]],
            [[0,4,4],[0.,+1.477080880741393e8]],
            [[0,7,7],[0.,-1.477080880741393e8]],
            [[1,5],[0.,-2.954161761482786e8]],
            [[1,4,5],[0.,-1.477080880741393e8]],
            [[2,6],[0.,-2.954161761482786e8]],
            [[2,4,6],[0.,-1.477080880741393e8]],
            [[3,7],[0.,-2.954161761482786e8]],
            [[0,0,4],[0.,-1.477080880741393e8]],
            [[0,0,4,4],[0.,-7.385404403706966e7]],
            [[0,0,5,5],[0.,+7.385404403706966e7]],
            [[0,0,6,6],[0.,+7.385404403706966e7]],
            [[0,0,7,7],[0.,+7.385404403706966e7]],
            [[0,1,5],[0.,+1.477080880741393e8]],
            [[0,1,4,5],[0.,+1.477080880741393e8]],
            [[0,2,6],[0.,+1.477080880741393e8]],
            [[0,2,4,6],[0.,+1.477080880741393e8]],
            [[1,1,5,5],[0.,-1.477080880741393e8]],
            [[1,2,5,6],[0.,-2.954161761482786e8]],
            [[1,3,5,7],[0.,-1.477080880741393e8]],
            [[2,2,6,6],[0.,-1.477080880741393e8]],
            [[2,3,6,7],[0.,-1.477080880741393e8]],
            [[3,3,4],[0.,+1.477080880741393e8]],
            [[3,3,4,4],[0.,+7.385404403706966e7]],
            [[3,3,5,5],[0.,-7.385404403706966e7]],
            [[3,3,6,6],[0.,-7.385404403706966e7]],
            [[3,3,7,7],[0.,-7.385404403706966e7]]]
            },
            ('p2', 'p6'):
            {():
            [[[0,4],[0.,-2.954161761482786e8]],
            [[0,4,4],[0.,+1.477080880741393e8]],
            [[0,7,7],[0.,-1.477080880741393e8]],
            [[1,4,5],[0.,-1.477080880741393e8]],
            [[2,4,6],[0.,-1.477080880741393e8]],
            [[3,7],[0.,-2.954161761482786e8]],
            [[0,0,4],[0.,+1.477080880741393e8]],
            [[0,0,4,4],[0.,-7.385404403706966e7]],
            [[0,0,5,5],[0.,+7.385404403706966e7]],
            [[0,0,6,6],[0.,+7.385404403706966e7]],
            [[0,0,7,7],[0.,+7.385404403706966e7]],
            [[0,1,5],[0.,-1.477080880741393e8]],
            [[0,1,4,5],[0.,+1.477080880741393e8]],
            [[0,2,6],[0.,-1.477080880741393e8]],
            [[0,2,4,6],[0.,+1.477080880741393e8]],
            [[1,1,5,5],[0.,-1.477080880741393e8]],
            [[1,2,5,6],[0.,-2.954161761482786e8]],
            [[1,3,5,7],[0.,-1.477080880741393e8]],
            [[2,2,6,6],[0.,-1.477080880741393e8]],
            [[2,3,6,7],[0.,-1.477080880741393e8]],
            [[3,3,4],[0.,-1.477080880741393e8]],
            [[3,3,4,4],[0.,+7.385404403706966e7]],
            [[3,3,5,5],[0.,-7.385404403706966e7]],
            [[3,3,6,6],[0.,-7.385404403706966e7]],
            [[3,3,7,7],[0.,-7.385404403706966e7]]],
            ('p3',):
            [[[0],[0.,+4.726658818372458e9]],
            [[0,4,4],[0.,-1.477080880741393e8]],
            [[0,7,7],[0.,+1.477080880741393e8]],
            [[1,4,5],[0.,+1.477080880741393e8]],
            [[2,4,6],[0.,+1.477080880741393e8]],
            [[0,0],[0.,-2.363329409186229e9]],
            [[0,0,4,4],[0.,+7.385404403706966e7]],
            [[0,0,5,5],[0.,-7.385404403706966e7]],
            [[0,0,6,6],[0.,-7.385404403706966e7]],
            [[0,0,7,7],[0.,-7.385404403706966e7]],
            [[0,1,4,5],[0.,-1.477080880741393e8]],
            [[0,2,4,6],[0.,-1.477080880741393e8]],
            [[1,1,5,5],[0.,+1.477080880741393e8]],
            [[1,2,5,6],[0.,+2.954161761482786e8]],
            [[1,3,5,7],[0.,+1.477080880741393e8]],
            [[2,2,6,6],[0.,+1.477080880741393e8]],
            [[2,3,6,7],[0.,+1.477080880741393e8]],
            [[3,3],[0.,+2.363329409186229e9]],
            [[3,3,4,4],[0.,-7.385404403706966e7]],
            [[3,3,5,5],[0.,+7.385404403706966e7]],
            [[3,3,6,6],[0.,+7.385404403706966e7]],
            [[3,3,7,7],[0.,+7.385404403706966e7]]]
            },
            ('p3', 'p5'):
            {():
                [[[0,4],[0.,-2.954161761482786e8]],
                [[0,4,4],[0.,+1.477080880741393e8]],
                [[0,7,7],[0.,-1.477080880741393e8]],
                [[1,4,5],[0.,-1.477080880741393e8]],
                [[2,4,6],[0.,-1.477080880741393e8]],
                [[3,7],[0.,-2.954161761482786e8]],
                [[0,0,4],[0.,+1.477080880741393e8]],
                [[0,0,4,4],[0.,-7.385404403706966e7]],
                [[0,0,5,5],[0.,+7.385404403706966e7]],
                [[0,0,6,6],[0.,+7.385404403706966e7]],
                [[0,0,7,7],[0.,+7.385404403706966e7]],
                [[0,1,5],[0.,-1.477080880741393e8]],
                [[0,1,4,5],[0.,+1.477080880741393e8]],
                [[0,2,6],[0.,-1.477080880741393e8]],
                [[0,2,4,6],[0.,+1.477080880741393e8]],
                [[1,1,5,5],[0.,-1.477080880741393e8]],
                [[1,2,5,6],[0.,-2.954161761482786e8]],
                [[1,3,5,7],[0.,-1.477080880741393e8]],
                [[2,2,6,6],[0.,-1.477080880741393e8]],
                [[2,3,6,7],[0.,-1.477080880741393e8]],
                [[3,3,4],[0.,-1.477080880741393e8]],
                [[3,3,4,4],[0.,+7.385404403706966e7]],
                [[3,3,5,5],[0.,-7.385404403706966e7]],
                [[3,3,6,6],[0.,-7.385404403706966e7]],
                [[3,3,7,7],[0.,-7.385404403706966e7]]],
            ('p2',):
                [[[0],[0.,4.726658818372458e9]],
                [[0,4,4],[0.,-1.477080880741393e8]],
                [[0,7,7],[0.,1.477080880741393e8]],
                [[1,4,5],[0.,1.477080880741393e8]],
                [[2,4,6],[0.,1.477080880741393e8]],
                [[0,0],[0.,-2.363329409186229e9]],
                [[0,0,4,4],[0.,7.385404403706966e7]],
                [[0,0,5,5],[0.,-7.385404403706966e7]],
                [[0,0,6,6],[0.,-7.385404403706966e7]],
                [[0,0,7,7],[0.,-7.385404403706966e7]],
                [[0,1,4,5],[0.,-1.477080880741393e8]],
                [[0,2,4,6],[0.,-1.477080880741393e8]],
                [[1,1,5,5],[0.,1.477080880741393e8]],
                [[1,2,5,6],[0.,2.954161761482786e8]],
                [[1,3,5,7],[0.,1.477080880741393e8]],
                [[2,2,6,6],[0.,1.477080880741393e8]],
                [[2,3,6,7],[0.,1.477080880741393e8]],
                [[3,3],[0.,2.363329409186229e9]],
                [[3,3,4,4],[0.,-7.385404403706966e7]],
                [[3,3,5,5],[0.,7.385404403706966e7]],
                [[3,3,6,6],[0.,7.385404403706966e7]],
                [[3,3,7,7],[0.,7.385404403706966e7]]]
            },
            }
        )

	info = ee_to_dd_2l_doubletriangle

	squared_topology = SquaredTopology(info, "hyperparameters.yaml")

	CUT_INDEX = 1
	random_variables = [random.random() for i in range(3*len(info.cuts[CUT_INDEX]['cuts'])-4)]
	cut_momenta = squared_topology.cutkosky_cuts[CUT_INDEX].get_PS_point(random_variables)
	cutkosky_cut_momenta_num = [cut_momenta[cut['edge']] for cut in info.cuts[CUT_INDEX]['cuts']]
	singular_surfaces = squared_topology.cutkosky_cuts[CUT_INDEX].get_singular_surfaces(cutkosky_cut_momenta_num)
	
	for c,cutkosky_cut in enumerate(squared_topology.cutkosky_cuts):
		for i,subgraph in enumerate(cutkosky_cut.subgraphs):
			if subgraph.n_loops > 0:
				side = 'left' if i == 0 else 'right'
				title = 'UV region in ' + side + ' subgraph in ' + str(cutkosky_cut)
				random_directions = [numpy.array([random.random(),random.random(),random.random()]) for i in range(squared_topology.n_loops)]
				random_directions = [direction/numpy.linalg.norm(direction) for direction in random_directions]
				squared_topology.plot_uv(subgraph.loop_basis, random_directions, title=title, precision='f128', phase=0,simulate_jac_power=0)
	plt.show()
	stop

	for e_surface in singular_surfaces['e_surfaces']:
		title = 'E-surface (%i dim) '%e_surface.e_surface_dimension + str(e_surface)
		print(10*'='+'{:=<100}'.format(' '+title+' '))
		random_variables = [random.random() for i in range(3*e_surface.ltd_loop_line.ltd_cut.cut_structure.subgraph.n_loops-1)]
		os_supergraph_integration_basis = e_surface.get_os_supergraph_integration_basis(cutkosky_cut_momenta_num, random_variables)
		spatial_os_supergraph_integration_basis = [p.space() for p in os_supergraph_integration_basis]
		random_directions = [numpy.array([random.random(),random.random(),random.random()]) for i in range(squared_topology.n_loops)]
		random_directions = [direction/numpy.linalg.norm(direction) for direction in random_directions]
		squared_topology.plot_singularity(spatial_os_supergraph_integration_basis,random_directions, title=title, precision='f128', phase=0)
		#squared_topology.probe_cancellation(spatial_os_supergraph_integration_basis,random_directions)
	for pinched_surface in singular_surfaces['pinched_surfaces']:
		title = 'Pinched surface (%i dim) '%pinched_surface.pinched_surface_dimension + str(pinched_surface)
		print(10*'='+'{:=<100}'.format(' '+title+' '))
		random_variables = [random.random() for i in range(3*pinched_surface.ltd_loop_line.ltd_cut.cut_structure.subgraph.n_loops-1)]
		#random_variables = [1 for i in range(3*pinched_surface.ltd_loop_line.ltd_cut.cut_structure.subgraph.n_loops-1)]
		os_supergraph_integration_basis = pinched_surface.get_os_supergraph_integration_basis(cutkosky_cut_momenta_num, random_variables)
		spatial_os_supergraph_integration_basis = [p.space() for p in os_supergraph_integration_basis]
		random_directions = [numpy.array([random.random(),random.random(),random.random()]) for i in range(squared_topology.n_loops)]
		random_directions = [direction/numpy.linalg.norm(direction) for direction in random_directions]
		squared_topology.plot_singularity(spatial_os_supergraph_integration_basis, random_directions, title=title, precision='f128', phase=0)#,simulate_jac_power=2)
		#squared_topology.probe_cancellation(spatial_os_supergraph_integration_basis,random_directions)
	plt.show()








