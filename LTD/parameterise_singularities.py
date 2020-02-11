from ltd_utils import *


class OneLoopCutParameteriser(object):
	def __init__(self,one_loop_topology, cutkosky_cut_momenta, external_momenta):
		self.propagators = one_loop_topology.loop_lines[0].propagators
		self.shifts = [(numpy.array(propagator.parametric_shift[0]).dot(numpy.array(cutkosky_cut_momenta))
							+ numpy.array(propagator.parametric_shift[1]).dot(numpy.array(external_momenta))) for propagator in self.propagators]

	def square(self,q):
		return q[0]**2-q[1:].dot(q[1:])

	def get_lorentz_vector(self,time,space):
		return vectors.LorentzVector(numpy.append([time],space))

	def get_space_vector(self,r,unit_vector):
		return r*unit_vector

	def get_unit_vector(self,cos_theta,phi):
		sin_theta = numpy.sqrt(1-cos_theta**2)
		return numpy.array([sin_theta*numpy.cos(phi),sin_theta*numpy.sin(phi),cos_theta])

	def get_beta(self,p):
		return p[1:]/p[0]

	def get_gamma(self,beta):
		return 1./numpy.sqrt(1.-beta.dot(beta))

	def kaellen_lambda(self,x,y,z):
		return x**2 + y**2 + z**2 - 2.*(x*y + x*z + y*z)

	def get_boost(self,beta):
		assert(beta.dot(beta)<1.)
		assert(beta.dot(beta)>0.)
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
		e_z = numpy.array([0,0,1])
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
		pji = pj - pi
		qi_radius_squared = self.kaellen_lambda(self.square(pji),mi_squared,mj_squared)/(4.*self.square(pji))
		if not (qi_radius_squared>0):
			qi_radius_squared += 0j
		qi_radius = numpy.sqrt(qi_radius_squared)
		return qi_radius

	def get_os_cos_theta(self,os_energy_index,os_radius_index,os_cos_theta_index):
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
		beta_ki = self.get_beta(boost.dot(pki))
		gamma_ki = self.get_gamma(beta_ki)
		qi_prime_radius = self.get_os_radius(os_energy_index,os_radius_index)
		qi_prime_energy = self.get_os_energy(os_energy_index,qi_prime_radius**2)
		pki_squared = self.square(pki)
		assert(pki_squared>0)
		cos_theta = (mi_squared-mk_squared+pki_squared
						-2.*gamma_ki*qi_prime_energy*numpy.sqrt(pki_squared)
					)/(-2.*gamma_ki*numpy.sqrt(pki_squared)*qi_prime_radius*numpy.sqrt(beta_ki.dot(beta_ki)))
		return cos_theta

	def get_k_1c_os(self,os_energy_index,qi_space):
		qi_energy = self.get_os_energy(os_energy_index,qi_space.dot(qi_space))
		qi = self.get_lorentz_vector(qi_energy,qi_space)
		pi = self.shifts[os_energy_index]
		k = qi-pi
		return k

	def get_k_2c_os(self,os_energy_index,os_radius_index,qi_prime_unit):
		qi_prime_radius = self.get_os_radius(os_energy_index,os_radius_index)
		assert(numpy.real(qi_prime_radius)>=0. and numpy.imag(qi_prime_radius)==0)
		qi_prime_energy = self.get_os_energy(os_energy_index,qi_prime_radius**2)
		qi_prime_space = self.get_space_vector(qi_prime_radius,qi_prime_unit)
		qi_prime = self.get_lorentz_vector(qi_prime_energy,qi_prime_space)
		pi = self.shifts[os_energy_index]
		pj = self.shifts[os_radius_index]
		pji = pj - pi
		beta = self.get_beta(pji)
		boost_inv = self.get_boost(-beta)
		k = boost_inv.dot(qi_prime) - pi
		return k

	def get_k_3c_os(self,os_energy_index,os_radius_index,os_cos_theta_index,phi_p2):
		qi_p2_radius = self.get_os_radius(os_energy_index,os_radius_index)
		assert(numpy.real(qi_p2_radius)>=0. and numpy.imag(qi_p2_radius)==0)
		qi_p2_energy = self.get_os_energy(os_energy_index,qi_p2_radius**2)
		cos_theta_p2 = self.get_os_cos_theta(os_energy_index,os_radius_index,os_cos_theta_index)
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
		beta_ki = self.get_beta(boost.dot(pki))
		rotation = self.get_rotation(beta_ki/(numpy.sqrt(beta_ki.dot(beta_ki))))
		rotation_inv = rotation.T
		k = boost_inv.dot(rotation_inv.dot(qi_p2)) - pi
		return k

	def get_existing_e_surfaces(self, dual_index):
		existing_esurfaces = []
		for prop_index, propagator in enumerate(self.propagators):
			if prop_index != dual_index:
				pi = self.shifts[dual_index]
				pj = self.shifts[prop_index]
				os_radius = self.get_os_radius(dual_index,prop_index)
				if numpy.imag(os_radius) == 0 and numpy.real(os_radius) >= 0 and (pj-pi)[0]<0:
					existing_esurfaces += [prop_index]
		return existing_esurfaces

	def get_existing_cut_group_intersections(self,dual_index,e_surf_index):
		assert(e_surf_index in self.get_existing_e_surfaces(dual_index))
		existing_cut_group_intersections = []
		for prop_index, propagator in enumerate(self.propagators):
			if prop_index not in [dual_index,e_surf_index]:
				if prop_index in self.get_existing_e_surfaces(dual_index):
					os_cos_theta = self.get_os_cos_theta(dual_index,e_surf_index,prop_index)
					if abs(os_cos_theta) <= 1.:
						existing_cut_group_intersections += [prop_index]
		return existing_cut_group_intersections

	def evaluate_propagator(self,prop_index,k):
		return self.square(k+self.shifts[prop_index])+self.propagators[prop_index].m_squared

# the double triangle master graph
t1 = SquaredTopologyGenerator([('q1', 0, 1), ('p1', 1, 2), ('p2', 2, 3), ('p3', 4, 3), ('p4', 4, 1), ('p5', 2, 4), ('q2', 3, 5)], "T", ['q1'], 2,
	{'q1': [1., 0., 0., 0.], 'q2': [1., 0., 0., 0.]},
	particle_ids={'p%s' % i: i for i in range(9)})
		#masses={'p1': 100, 'p2':100, 'p3': 100, 'p4': 100, 'p5': 100})

# select the cut with left tree and right loop
#print(t1.loop_topologies[1])
l_tree = t1.loop_topologies[1][0]
r_loop = t1.loop_topologies[1][1]
#print(l_tree)
#print(r_loop)

#print(t1.cuts[1])
p1 = -vectors.LorentzVector([5, -1, 2, -3])
p4 = vectors.LorentzVector([-9, 4, -4, 4])
cuts = [p1,p4]

#print(t1.external_momenta)
#q1 = t1.external_momenta['q1']
#q2 = t1.external_momenta['q2']
q1 = -vectors.LorentzVector([-4, 3, -2, 1])
q2 = vectors.LorentzVector([-4, 3, -2, 1])
ext = [q1,q2]

prop_momenta = {}
for propagator in r_loop.loop_lines[0].propagators:
	prop_momenta[propagator.name] = (numpy.array(propagator.parametric_shift[0]).dot(numpy.array(cuts))
									+ numpy.array(propagator.parametric_shift[1]).dot(numpy.array(ext)))
#print(prop_momenta)

myParam = OneLoopCutParameteriser(r_loop,cuts,ext)

for dual_index in range(len(myParam.propagators)):
	existing_esurfaces = myParam.get_existing_e_surfaces(dual_index)
	for e_surf_index in existing_esurfaces:
		k = myParam.get_k_2c_os(dual_index,e_surf_index,myParam.get_unit_vector(0.5,1))
		print([dual_index,e_surf_index],myParam.evaluate_propagator(e_surf_index,k))
		existing_cut_group_intersections = myParam.get_existing_cut_group_intersections(dual_index,e_surf_index)
		for prop_index in existing_cut_group_intersections:
			k = myParam.get_k_3c_os(dual_index,e_surf_index,prop_index,0.5)
			print([dual_index,e_surf_index,prop_index],myParam.evaluate_propagator(prop_index,k))





