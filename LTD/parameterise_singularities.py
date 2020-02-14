from ltd_utils import *
import warnings
import os, sys
sys.path.insert(0,'../')
import ltd
import matplotlib.pyplot as plt

class CrossSectionSingularityProber(object):
	def __init__(self, squared_topology, hyperparameters_yaml):
		squared_topology.export('squared_topology.yaml')
		self.cross_section = ltd.CrossSection('squared_topology.yaml', hyperparameters_yaml)
		self.squared_topology = squared_topology
		self.external_momenta = squared_topology.external_momenta
		self.external_momenta = [self.external_momenta['q1'],self.external_momenta['q2']]
		self.loop_momentum_basis = []
		self.n_cuts = len(self.squared_topology.cuts)

	def probe_singularities_in_cut_topology(self,index,cut_momenta,plot=True,precision='f128'):
		cut_topology = self.squared_topology.loop_topologies[index]
		cut_legs = self.squared_topology.cuts[index]
		cut_leg_signatures = self.squared_topology.cut_signatures[index]
		assert(set(cut_momenta.keys()) == {p_name for (p_name,p_energy_sign) in cut_legs})
		cut_momenta = [cut_momenta[p_name] for (p_name,p_energy_sign) in cut_legs]
		for l_or_r in [0,1]:
			if cut_topology[l_or_r].n_loops == 0:
				continue
			elif cut_topology[l_or_r].n_loops == 1:
				parameteriser = OneLoopSingularityParameteriser(self.squared_topology,index,l_or_r,cut_momenta,self.external_momenta)
				parameteriser.check_parameterisations()
				for (dual_index,e_surf_index) in parameteriser.existing_e_surfaces:
					print(10*'='+' probing E-surface (%i,%i) '%(dual_index,e_surf_index)+10*'=')
					k = parameteriser.get_k_2c_os(dual_index,e_surf_index,parameteriser.get_unit_vector(0.7,0.5))
					spatial_loop_momenta = parameteriser.get_spatial_loop_momenta(k)
					self.probe_singular_surface(spatial_loop_momenta,2*[numpy.ones(3)],plot=plot,precision=precision)
					self.evaluate_cuts(spatial_loop_momenta)
				for (dual_index,e_surf_index,intersection_index) in parameteriser.cut_group_intersections:
					print(10*'='+' probing E-surface (%i,%i) '%(dual_index,e_surf_index)+10*'=')
					k = parameteriser.get_k_3c_os(dual_index,e_surf_index,intersection_index,0.5)
					spatial_loop_momenta = parameteriser.get_spatial_loop_momenta(k)
					self.probe_singular_surface(spatial_loop_momenta,2*[numpy.ones(3)],plot=plot,precision=precision)
					self.evaluate_cuts(spatial_loop_momenta)
				for (dual_index,pinched_surf_index) in parameteriser.pinched_surfaces:
					print(10*'='+' probing pinch surface (%i,%i) '%(dual_index,e_surf_index)+10*'=')
					k = parameteriser.get_k_pinched(dual_index,pinched_surf_index,0.5)
					spatial_loop_momenta = parameteriser.get_spatial_loop_momenta(k)
					self.probe_singular_surface(spatial_loop_momenta,2*[numpy.ones(3)],plot=plot,precision=precision)
					self.evaluate_cuts(spatial_loop_momenta)
			else:
				print('Not implemented yet.')

	def probe_singular_surface(self,os_spatial_loop_momenta, spatial_directions,plot=True,precision='f128'):
		if precision == 'f128':
			evaluate = self.cross_section.evaluate_f128
			evaluate_cut = self.cross_section.evaluate_cut_f128
		elif precision == 'f64':
			evaluate = self.cross_section.evaluate
			evaluate_cut = self.cross_section.evaluate_cut
		mu_list = numpy.logspace(0,-4,5)
		results = []
		for mu in mu_list:
			spatial_loop_momenta = [p_space+mu*dir_space for p_space,dir_space in zip(os_spatial_loop_momenta,spatial_directions)]
			results += [evaluate(spatial_loop_momenta)]
		for result in results:
			print(result)
		if plot == True:
			plt.figure()
			xx = numpy.linspace(-0.000002,0.000002,int(1e4))
			spatial_loop_momenta_list = [[p+x*spatial_direction for p,spatial_direction in zip(os_spatial_loop_momenta,spatial_directions)] for x in xx]
			yy = [[evaluate_cut(spatial_loop_momenta,c)[1] for spatial_loop_momenta in spatial_loop_momenta_list] for c in range(self.n_cuts)]
			tot = numpy.sum(yy,axis=0)
			for c in range(self.n_cuts):
				plt.plot(xx,yy[c],'o',label='cut%i'%c,markersize=0.5)
			plt.plot(xx,tot,'ko',label='sum',markersize=0.5)
			plt.xlabel(r'$\mu$')
			plt.yscale('symlog')
			plt.legend(loc='lower right')
			plt.title('t1_squared'
				+ ''.join(['\n'+r'$k_{'+'%i'%c+'}=$'+str(spatial_loop_momentum)+r'$+\mu$'+str(spatial_direction)
						for c,spatial_loop_momentum,spatial_direction in zip(range(self.n_cuts),os_spatial_loop_momenta,spatial_directions)]))
			plt.show()

	def evaluate_cuts(self,spatial_loop_momenta,mu=1e-6,spatial_dir=vectors.Vector([1.,1.,1.])):
		#mu = 0.54*1e-8
		#spatial_loop_momenta = [p_space+mu*spatial_dir for p_space in spatial_loop_momenta]
		print('Total', self.cross_section.evaluate_f128(spatial_loop_momenta))

		tot = 0
		for x in range(self.n_cuts):
		    cut = self.cross_section.evaluate_cut_f128(spatial_loop_momenta, x)
		    print('cut%d: (%e, %e)' % (x, cut[0], cut[1]))
		    tot += numpy.array(cut)

		print('Sum cuts: (%e, %e)' % (tot[0],tot[1]))
		return

class OneLoopSingularityParameteriser(object):
	def __init__(self,squared_topology, cut_index, l_or_r, cut_momenta, external_momenta):
		self.external_momenta = external_momenta
		self.squared_topology = squared_topology
		self.one_loop_topology = self.squared_topology.loop_topologies[cut_index][l_or_r]
		self.cut_legs = self.squared_topology.cuts[cut_index]
		self.cut_leg_signatures = self.squared_topology.cut_signatures[cut_index]
		self.tolerance = 1e-10
		self.cut_momenta = cut_momenta
		self.propagators = self.one_loop_topology.loop_lines[0].propagators
		self.shifts = [vectors.LorentzVector((numpy.array(propagator.parametric_shift[0]).dot(numpy.array(cut_momenta))
							+ numpy.array(propagator.parametric_shift[1]).dot(numpy.array(external_momenta)))) for propagator in self.propagators]
		self.existing_e_surfaces = self.get_existing_e_surfaces()
		self.cut_group_intersections = self.get_cut_group_intersections()
		self.pinched_surfaces = self.get_pinched_surfaces()
		print('E-Surfaces: ',self.existing_e_surfaces)
		print('Cut Group Intersections:', self.cut_group_intersections)
		print('Pinched Surfaces: ', self.pinched_surfaces)

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
		cut_basis = numpy.array([self.cut_momenta[0]]+[k])
		trsf_loop_to_cut_basis = numpy.array([self.cut_leg_signatures[0][0]]+[self.one_loop_topology.loop_momentum_map[0][0]])
		inv_trsf_loop_to_cut_basis = numpy.linalg.inv(trsf_loop_to_cut_basis)
		trsf_loop_to_cut_shift = numpy.array([self.cut_leg_signatures[0][1]]+[self.one_loop_topology.loop_momentum_map[0][1]])
		cut_basis_shift = -trsf_loop_to_cut_shift.dot(self.external_momenta)
		loop_basis = inv_trsf_loop_to_cut_basis.dot(cut_basis + cut_basis_shift)
		loop_basis = [vectors.LorentzVector(basis) for basis in loop_basis]
		spatial_loop_momenta = [basis.space() for basis in loop_basis]
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

# the double triangle master graph
t1 = SquaredTopologyGenerator([('q1', 0, 1), ('p1', 1, 2), ('p2', 2, 3), ('p3', 4, 3), ('p4', 4, 1), ('p5', 2, 4), ('q2', 3, 5)], "T", ['q1'], 2,
	{'q1': [1., 0., 0., 0.], 'q2': [1., 0., 0., 0.]},
	particle_ids={'p%s' % i: i for i in range(9)})
		#masses={'p1': 100, 'p2':100, 'p3': 100, 'p4': 100, 'p5': 100})

prober = CrossSectionSingularityProber(t1,"hyperparameters.yaml")

q1 = t1.external_momenta['q1']
q2 = t1.external_momenta['q2']
ext = [q1,q2]

p1_space = vectors.Vector([1.,2.,3.])
p4_space = p1_space
cuts_space = [p1_space,p4_space]

cut_index = 1
# set arbitrary outgoing
cuts_lorentz_vector = []
for cut_p, p_space in zip(t1.cuts[cut_index],cuts_space):
	cut_p_name = cut_p[0]
	cut_p_sign = cut_p[1]
	m_squared = t1.masses[cut_p_name]**2 if cut_p_name in t1.masses else 0.
	p_energy  = cut_p_sign * numpy.sqrt(p_space.dot(p_space)+m_squared)
	p = vectors.LorentzVector(numpy.append([p_energy],p_space))
	cuts_lorentz_vector += [p]
#print(cuts_lorentz_vector)
# set momentum conservation with incomming and rescale
energy_sum = 0.
for cut_p, p in zip(t1.cuts[cut_index],cuts_lorentz_vector):
	cut_p_sign = cut_p[1]
	energy_sum += cut_p_sign*p[0]
t = ext[0][0]/energy_sum

#print(t)
cuts = [t*p for p in cuts_lorentz_vector]
#print(cuts)
cut_momenta = {'p1': cuts[0],'p4':cuts[1]}

# check momentum conservation
print(ext[0]-cuts[0]+cuts[1])

prober.probe_singularities_in_cut_topology(1,cut_momenta,plot=True,precision='f128')

#myParam = OneLoopSingularityParameteriser(t1,1,1,cuts,ext)
