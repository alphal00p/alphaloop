from ltd_utils import *
import warnings
import random
import os, sys
sys.path.insert(0,'../')
import ltd
import matplotlib.pyplot as plt
import itertools

class CrossSectionSingularityProber(object):
	def __init__(self, squared_topology, hyperparameters_yaml):
		squared_topology.export('squared_topology.yaml')
		self.cross_section = ltd.CrossSection('squared_topology.yaml', hyperparameters_yaml)
		self.squared_topology = squared_topology
		self.incoming_momenta_names = squared_topology.incoming_momenta
		self.external_momenta = [squared_topology.external_momenta['q%i'%(i+1)] for i in range(len(squared_topology.external_momenta))]
		self.n_cuts = len(self.squared_topology.cuts)
		self.n_loops = squared_topology.topo.n_loops
		self.name = squared_topology.name
		self.e_cm_squared = sum(squared_topology.external_momenta[e][0]
								for e in squared_topology.incoming_momenta)**2 - sum(x*x
								for x in (sum(squared_topology.external_momenta[e][i]
								for e in squared_topology.incoming_momenta)
								for i in range(1, 4)))
		#print(self.cross_section.evaluate_f128([[0.5, 1.5, 1.5],[5.5, 1.5, 0.5]]))
		#print(self.cross_section.get_scaling([[0.5, 0.5, 0.5],[0.5, 0.5, 0.5]],2))
		"""
		test_loop_momenta = [[1.,2.1,3.],[1.,2.,3.]]
		print('tot = ', self.cross_section.evaluate(test_loop_momenta))
		tot = numpy.array([0.,0.])
		for c in range(self.n_cuts):
			scaling_and_jac = self.cross_section.get_scaling(test_loop_momenta,c)
			scaling_and_jac = scaling_and_jac[0]
			cut = self.cross_section.evaluate_cut(test_loop_momenta,c,scaling_and_jac[0],scaling_and_jac[1])
			tot += numpy.array(cut)
			print('cut %i =' %c, cut)
		print('sum = ', tuple(tot))
		"""

	def rescale_cut_momenta(self,index,spatial_cut_momenta_dict):
		assert(len(spatial_cut_momenta_dict) == len(self.squared_topology.cuts[index])-1)
		cut_topology = self.squared_topology.loop_topologies[index]
		cut_legs = self.squared_topology.cuts[index]
		cut_leg_signatures = self.squared_topology.cut_signatures[index]
		cut_basis = []
		cut_signatures = []
		cut_shifts = []
		for (name, cut_momentum) in spatial_cut_momenta_dict.items():
			cut_basis += [cut_momentum]
			for cut_leg_index, cut_leg in enumerate(cut_legs):
				if name == cut_leg[0]:
					cut_signatures += [cut_leg_signatures[cut_leg_index][0]]
					cut_shifts += [cut_leg_signatures[cut_leg_index][1]]
		for l_or_r in [0,1]:
			if cut_topology[l_or_r].n_loops > 0:
				for loop_index, signature in enumerate(cut_topology[l_or_r].loop_momentum_map):
					cut_basis += [vectors.Vector([0.,0.,0.])] #random cut basis vector, result is independent of it
					cut_signatures += [signature[0]]
					cut_shifts += [signature[1]]
		cut_basis = numpy.array(cut_basis)
		trsf_loop_to_cut_basis = numpy.array(cut_signatures)
		inv_trsf_loop_to_cut_basis = numpy.linalg.inv(trsf_loop_to_cut_basis)
		trsf_loop_to_cut_shift = numpy.array(cut_shifts)
		spatial_external = numpy.array([external_momentum[1:] for external_momentum in self.external_momenta])
		supergraph_basis = inv_trsf_loop_to_cut_basis.dot(cut_basis-trsf_loop_to_cut_shift.dot(spatial_external))
		scalings_and_jac = self.cross_section.get_scaling(supergraph_basis,index)
		positive_scalings_and_jac = [scaling_and_jac for scaling_and_jac in scalings_and_jac if scaling_and_jac[0]>0]
		if len(positive_scalings_and_jac) == 1:
			scaled_supergraph_basis = positive_scalings_and_jac[0][0]*supergraph_basis
		else:
			print('Select manually which scaling to use.')
		all_spatial_cut_momenta_scaled = [numpy.array(cut_leg_signatures[cut_leg_index][0]).dot(scaled_supergraph_basis)
											+ numpy.array(cut_leg_signatures[cut_leg_index][1]).dot(spatial_external)
											for cut_leg_index, cut_leg in enumerate(cut_legs)]
		all_cut_momenta_scaled = []
		for cut_p, p_space in zip(cut_legs,all_spatial_cut_momenta_scaled):
			cut_p_name = cut_p[0]
			cut_p_sign = cut_p[1]
			m_squared = self.squared_topology.masses[cut_p_name]**2 if cut_p_name in self.squared_topology.masses else 0.
			p_energy  = cut_p_sign * numpy.sqrt(p_space.dot(p_space)+m_squared)
			p = vectors.LorentzVector(numpy.append([p_energy],p_space))
			all_cut_momenta_scaled += [p]
		# check momentum conservation
		incoming_momenta = [momentum for name,momentum in self.squared_topology.external_momenta.items() if name in self.incoming_momenta_names]
		assert(all(numpy.zeros(4) == numpy.sum(incoming_momenta,axis=0)-numpy.sum([name_and_sign[1]*p for name_and_sign,p in zip(cut_legs,all_cut_momenta_scaled)],axis=0)))
		all_cut_momenta_scaled_dict = {name_and_sign[0]: p for name_and_sign,p in zip(cut_legs,all_cut_momenta_scaled)}
		return all_cut_momenta_scaled_dict

	def probe_singularities_in_cut_topology(self,index,cut_momenta,directions,plot=True,precision='f128'):
		cut_topology = self.squared_topology.loop_topologies[index]
		cut_legs = self.squared_topology.cuts[index]
		cut_leg_signatures = self.squared_topology.cut_signatures[index]
		print('Cut considered: '+str(cut_legs))
		#os_cut_momenta = self.rescale_cut_momenta(index,cut_momenta)
		os_cut_momenta = cut_momenta
		assert(set(os_cut_momenta.keys()) == {p_name for (p_name,p_energy_sign) in cut_legs})
		os_cut_momenta = [os_cut_momenta[p_name] for (p_name,p_energy_sign) in cut_legs]
		#print(os_cut_momenta)
		for l_or_r in [0,1]:
			if cut_topology[l_or_r].n_loops == 0:
				continue
			elif cut_topology[l_or_r].n_loops == 1:
				parameteriser = OneLoopSingularityParameteriser(self.squared_topology,index,l_or_r,os_cut_momenta)
				parameteriser.check_parameterisations()
				title_suffix = ' in Cutkosky cut %i: (%s)'%(index,','.join([name_and_sign[0]+('%+d'%name_and_sign[1])[0] for name_and_sign in self.squared_topology.cuts[index]]))
				for (dual_index,e_surf_index) in parameteriser.existing_e_surfaces:
					title = 'E-surface (%i,%i), i.e. (%s+,%s-) os'%(dual_index,e_surf_index,
								cut_topology[l_or_r].loop_lines[0].propagators[dual_index].name,
								cut_topology[l_or_r].loop_lines[0].propagators[e_surf_index].name)+title_suffix
					print(10*'='+' probing ' + title + ' ' +10*'=')
					k = parameteriser.get_k_2c_os(dual_index,e_surf_index,parameteriser.get_unit_vector(0.7,0.5))
					spatial_loop_momenta = parameteriser.get_spatial_loop_momenta(k)
					self.probe_singular_surface(spatial_loop_momenta,directions,plot=plot,plot_title=title,precision=precision)
					self.evaluate_cuts(spatial_loop_momenta)
				for (dual_index,e_surf_index,intersection_index) in parameteriser.cut_group_intersections:
					title = 'E-surface intersection (%i+,%i-,%i-), i.e. (%s+,%s,%s) os'%(dual_index,e_surf_index,intersection_index,
								cut_topology[l_or_r].loop_lines[0].propagators[dual_index].name,
								cut_topology[l_or_r].loop_lines[0].propagators[e_surf_index].name,
								cut_topology[l_or_r].loop_lines[0].propagators[intersection_index].name)+title_suffix
					print(10*'='+' probing ' + title + ' ' +10*'=')
					k = parameteriser.get_k_3c_os(dual_index,e_surf_index,intersection_index,0.5)
					spatial_loop_momenta = parameteriser.get_spatial_loop_momenta(k)
					self.probe_singular_surface(spatial_loop_momenta,directions,plot=plot,plot_title=title,precision=precision)
					self.evaluate_cuts(spatial_loop_momenta)
				for (dual_index,pinched_surf_index) in parameteriser.pinched_surfaces:
					title = 'pinch surface (%i,%i), i.e. (%s+,%s-) os'%(dual_index,pinched_surf_index,
								cut_topology[l_or_r].loop_lines[0].propagators[dual_index].name,
								cut_topology[l_or_r].loop_lines[0].propagators[pinched_surf_index].name)+title_suffix
					print(10*'='+' probing ' + title + ' ' +10*'=')
					k = parameteriser.get_k_pinched(dual_index,pinched_surf_index,1.-1e-10)
					spatial_loop_momenta = parameteriser.get_spatial_loop_momenta(k)
					self.probe_singular_surface(spatial_loop_momenta,directions,plot=plot,plot_title=title,precision=precision)
					self.evaluate_cuts(spatial_loop_momenta)
			else:
				print('Not implemented yet.')
		if plot:
			plt.show()
		return

	def probe_singular_surface(self,os_spatial_loop_momenta, spatial_directions,plot=True,plot_title='',precision='f128',phase=0):
		if precision == 'f128':
			evaluate = self.cross_section.evaluate_f128
			evaluate_cut = self.cross_section.evaluate_cut_f128
		elif precision == 'f64':
			evaluate = self.cross_section.evaluate
			evaluate_cut = self.cross_section.evaluate_cut
		mu_list = numpy.logspace(-1,-5,5)
		results = []
		for mu in mu_list:
			spatial_loop_momenta = [p_space+mu*dir_space*numpy.sqrt(self.e_cm_squared) for p_space,dir_space in zip(os_spatial_loop_momenta,spatial_directions)]
			results += [evaluate(spatial_loop_momenta)]
		print('Approaching the singularity:')
		for result in results:
			print(result)
		if plot == True:
			plt.figure()
			xx = numpy.linspace(-1e-1,1e-1,int(1e3))
			spatial_loop_momenta_list = [[p+x*spatial_direction*numpy.sqrt(self.e_cm_squared) for p,spatial_direction in zip(os_spatial_loop_momenta,spatial_directions)] for x in xx]
			yy = [[evaluate_cut(spatial_loop_momenta,c,
					self.cross_section.get_scaling(spatial_loop_momenta,c)[1][0],self.cross_section.get_scaling(spatial_loop_momenta,c)[1][1])[phase]
					for spatial_loop_momenta in spatial_loop_momenta_list] for c in range(self.n_cuts)]
			tot = numpy.sum(yy,axis=0)
			markers = itertools.cycle(('v','<', '^', '>')) 
			for c in range(self.n_cuts):
				plt.plot(xx,yy[c],marker=markers.next(),label='cut %i: (%s)'%(c,','.join([name_and_sign[0]+('%+d'%name_and_sign[1])[0] for name_and_sign in self.squared_topology.cuts[c]])),
							ms=1,alpha=1,ls='None')
			plt.plot(xx,tot,'k',marker='o',label='sum',ms=0.5,alpha=1,ls='None')
			plt.xlabel(r'$\mu$')
			plt.yscale('symlog')
			plt.legend(loc='lower right')
			plt.title(self.name + '\n' + plot_title
				+ ''.join(['\n'+r'$k_{'+'%i'%c+'}=$'
							+'(%s)'%','.join(['%f'%k_i for k_i in spatial_loop_momentum])
							+r'$+\mu E_{\mathrm{cm}}$'
							+'(%s)'%','.join(['%f'%k_i for k_i in spatial_direction])
						for c,spatial_loop_momentum,spatial_direction in zip(range(self.n_cuts),os_spatial_loop_momenta,spatial_directions)]))
			plt.tight_layout()
			#plt.show()

	def evaluate_cuts(self,spatial_loop_momenta,mu=1e-6,spatial_dir=vectors.Vector([1.,1.,1.])):
		mu = 1e-3
		spatial_loop_momenta = [p_space+mu*spatial_dir*numpy.sqrt(self.e_cm_squared) for p_space in spatial_loop_momenta]
		print('On the singularity:')
		print('Total: (%e, %e)'%self.cross_section.evaluate_f128(spatial_loop_momenta))

		tot = 0
		for c in range(self.n_cuts):
			#print(self.cross_section.get_scaling(spatial_loop_momenta,c))
			cut = self.cross_section.evaluate_cut_f128(spatial_loop_momenta, c,
				self.cross_section.get_scaling(spatial_loop_momenta,c)[1][0],self.cross_section.get_scaling(spatial_loop_momenta,c)[1][1])
			print('cut%d: (%e, %e)' % (c, cut[0], cut[1]))
			tot += numpy.array(cut)

		print('Sum cuts: (%e, %e)' % (tot[0],tot[1]))
		return

class OneLoopSingularityParameteriser(object):
	def __init__(self, squared_topology, cut_index, l_or_r, cut_momenta):
		self.external_momenta = [squared_topology.external_momenta['q%i'%(i+1)] for i in range(len(squared_topology.external_momenta))]
		self.squared_topology = squared_topology
		self.incoming_momenta_names = squared_topology.incoming_momenta
		self.cut_topology = self.squared_topology.loop_topologies[cut_index]
		self.one_loop_topology = self.squared_topology.loop_topologies[cut_index][l_or_r]
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

class CutkoskyCutParametiser(object):
	def __init__(self,squared_topology,cut_index):
		self.external_momenta = [squared_topology.external_momenta['q%i'%(i+1)] for i in range(len(squared_topology.external_momenta))]
		self.squared_topology = squared_topology
		self.incoming_momenta_names = squared_topology.incoming_momenta
		self.cut_topology = self.squared_topology.loop_topologies[cut_index]
		self.cut_legs = self.squared_topology.cuts[cut_index]
		self.cut_leg_signatures = self.squared_topology.cut_signatures[cut_index]
		self.tolerance = 1e-10

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
	
	def get_PS_point(self,random_variables):
		counter = 0
		assert(len(random_variables)==3*(len(self.cut_legs)-1)-1)
		# see last cut momentum as a dual propagator that is put on-shell by solving for a radius
		# express dual propagator in cut basis
		# check that the number of cut energies (# square roots -1 in e-surface) indeed is len(self.cut_legs)-1
		cut_basis, trsf_loop_to_cut_basis, trsf_loop_to_cut_shift = self.get_cut_basis(include_loop_momenta=True)
		prop_name, prop_sign = self.cut_legs[-1]
		signature_loop_basis = numpy.array(self.cut_leg_signatures[-1][0])
		shift_loop_basis = numpy.array(self.cut_leg_signatures[-1][1])
		inv_trsf_loop_to_cut_basis = numpy.linalg.inv(trsf_loop_to_cut_basis)
		signature_cut_basis = signature_loop_basis.dot(inv_trsf_loop_to_cut_basis)
		shift_cut_basis = shift_loop_basis-signature_cut_basis.dot(trsf_loop_to_cut_shift)
		#print(prop_name,prop_sign)
		#print('sig_l',signature_loop_basis)
		#print('sh_l',shift_loop_basis)
		#print('inv_trsf_l',inv_trsf_loop_to_cut_basis)
		#print('trsf_s',trsf_loop_to_cut_shift)
		#print('sig_c',signature_cut_basis)
		#print('sh_c',shift_cut_basis)
		# count number of cut energies the propagator depends on
		n_independent_cuts = len([1 for sign in signature_cut_basis if sign != 0])
		assert(n_independent_cuts==len(self.cut_legs)-1)
		#os_momenta_flow = [(name,energy_sign,signature) for (name,energy_sign),signature in zip(cut_basis,signature_cut_basis) if signature != 0]
		#print(os_momenta_flow)
		os_cut_basis = [vectors.LorentzVector([0,0,0,0]) for i in range(len(cut_basis))]
		#print(os_cut_basis)
		indep_counter = 0
		for index,((name,energy_sign),signature) in enumerate(zip(cut_basis,signature_cut_basis)):
			if signature != 0:
				#print(name,energy_sign)
				indep_counter += 1
				mass_squared_basis = self.squared_topology.masses[name] if name in self.squared_topology.masses else 0.
				mass_prop = self.squared_topology.masses[prop_name] if prop_name in self.squared_topology.masses else 0.
				mass_prop += numpy.sum([self.squared_topology.masses[name]
												if (name in self.squared_topology.masses and sign != 0) else 0.
												for sign,((name,_)) in zip(signature_cut_basis[:index],cut_basis[:index])])
				mass_squared_prop = mass_prop**2
				shift = vectors.LorentzVector(shift_cut_basis.dot(self.external_momenta))+numpy.sum([sign*q 
									for sign, q in zip(signature_cut_basis[:index],os_cut_basis[:index]) if sign != 0],axis=0)
				#print(mass_squared_basis,mass_squared_prop,shift)
				#print('signature',signature_cut_basis[:index])
				#print(shift,vectors.LorentzVector(shift_cut_basis.dot(self.external_momenta)))
				#print('shift.sq',shift.square())
				#print('sum mass sq',(numpy.sqrt(mass_squared_prop)+numpy.sqrt(mass_squared_basis))**2)
				assert(shift.square() > (numpy.sqrt(mass_squared_prop)+numpy.sqrt(mass_squared_basis))**2)
				#print('0 comp sign', shift[0],energy_sign,signature)
				assert(shift[0]*energy_sign*signature < 0)
				radius_max = numpy.sqrt(self.kaellen_lambda(shift.square(),mass_squared_prop,mass_squared_basis)/(4.*shift.square()))
				#print(indep_counter,n_independent_cuts)
				if indep_counter == n_independent_cuts:
					radius = radius_max
				else:
					#radius = numpy.sqrt(self.e_cm_squared)*random_variables[counter]/(1.-random_variables[counter])+radius_min
					radius = random_variables[counter]*radius_max
					#print(radius,radius_max)
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
				#print(os_cut_basis)
		os_cut_basis = numpy.array(os_cut_basis)
		os_prop_momentum = vectors.LorentzVector(signature_cut_basis.dot(os_cut_basis)+shift_cut_basis.dot(self.external_momenta))
		#print(os_prop_momentum)
		# assert momentum conservation
		incoming_momenta = [vectors.LorentzVector(momentum) for name,momentum in self.squared_topology.external_momenta.items() if name in self.incoming_momenta_names]
		all_cut_momenta = {name: vectors.LorentzVector(momentum) for (name,energy_sign),momentum in zip(cut_basis,os_cut_basis) if (name,energy_sign) in self.cut_legs}
		all_cut_momenta.update({prop_name: os_prop_momentum})
		outgoing_momenta = [energy_sign*all_cut_momenta[name] for (name,energy_sign) in self.cut_legs]
		assert(all(numpy.zeros(4) == numpy.sum(incoming_momenta,axis=0)-numpy.sum(outgoing_momenta,axis=0)))
		for name,p in all_cut_momenta.items():
			print(name,p)
		print(incoming_momenta)
		return all_cut_momenta

class SquaredTopology(object):
	def __init__(self, **opts):
		self.opts = opts
		self.squared_topology_dump = self.opts['squared_topology_dump']
		self.n_cutkosky_cuts = len(self.squared_topology_dump.cuts)
		self.n_loops = self.squared_topology_dump.topo.n_loops
		self.name = self.squared_topology_dump.name
		self.e_cm_squared = sum(self.squared_topology_dump.external_momenta[e][0]
								for e in self.squared_topology_dump.incoming_momenta)**2 - sum(x*x
								for x in (sum(self.squared_topology_dump.external_momenta[e][i]
								for e in self.squared_topology_dump.incoming_momenta)
								for i in range(1, 4)))
		self.external_momenta = ['q%i'%(i+1) for i in range(len(self.squared_topology_dump.external_momenta))]
		self.external_momenta_num = [squared_topology.external_momenta[name] for name in self.external_momenta]
		self.incoming_momenta = self.squared_topology_dump.incoming_momenta
		self.supergraph_loop_basis = ['p%i'%i for i in self.squared_topology_dump.topo.loop_momenta]

	def generate(self):
		self.cutkosky_cuts = [CutkoskyCut(cutkosky_cut_index=cutkosky_cut_index, **self.opts) for cutkosky_cut_index in range(self.n_cutkosky_cuts)]
		for cutkosky_cut in self.cutkosky_cuts:
			cutkosky_cut.generate()

class CutkoskyCut(SquaredTopology):
	def __init__(self, **opts):
		super(CutkoskyCut, self).__init__(**opts)
		self.opts = opts
		self.cutkosky_cut_index = self.opts['cutkosky_cut_index']
		self.cutkosky_cut_momenta = [momentum for (momentum, energy_sign) in self.squared_topology_dump.cuts[self.cutkosky_cut_index]]
		self.cutkosky_cut_energy_signs = [energy_sign for (momentum, energy_sign) in self.squared_topology_dump.cuts[self.cutkosky_cut_index]]
		self.cutkosky_cut_signatures = [signature for (signature, shift) in self.squared_topology_dump.cut_signatures[self.cutkosky_cut_index]]
		self.cutkosky_cut_shifts = [shift for (signature, shift) in self.squared_topology_dump.cut_signatures[self.cutkosky_cut_index]]

	def generate(self):
		self.subgraphs = [SubGraph(left_or_right=0,**self.opts), SubGraph(left_or_right=1,**self.opts)]
		for subgraph in self.subgraphs:
			subgraph.generate()

	def get_e_surfaces(self, cutkosky_cut_momenta_num):
		e_surfaces = []
		#for subgraph in [self.subgraph_left,self.subgraph_right]:
		#	if subgraph.n_loops_subgraph > 0:
		#		for e_surface in subgraph.get_e_surfaces(cutkosky_cut_momenta_num):
		#			if e_surface != []:
		#				e_surface['subgraph'] = subgraph.left_or_right
		#				e_surfaces += [e_surface]
		for subgraph in self.subgraphs:
			if subgraph.n_loops_subgraph > 0:
				e_surfaces += subgraph.get_e_surfaces(cutkosky_cut_momenta_num)
		return e_surfaces

	def get_PS_point(self, random_variables):
		""" put a function here"""

class SubGraph(CutkoskyCut):
	def __init__(self, **opts):
		super(SubGraph, self).__init__(**opts)
		self.opts = opts
		self.left_or_right = self.opts['left_or_right']
		self.subgraph_dump = self.squared_topology_dump.loop_topologies[self.cutkosky_cut_index][self.left_or_right]
		self.n_loops_subgraph = self.subgraph_dump.n_loops

	def generate(self):
		if self.n_loops_subgraph > 0:
			self.ltd_cut_structures = [LTDCutStructure(ltd_cut_structure=ltd_cut_structure,**self.opts) for ltd_cut_structure in self.subgraph_dump.ltd_cut_structure]
			for ltd_cut_structure in self.ltd_cut_structures:
				ltd_cut_structure.generate()
		else:
			print('TODO: Complete the tree graph case!')

	def get_e_surfaces(self,cutkosky_cut_momenta_num):
		e_surfaces = []
		#for ltd_cut_structure in self.ltd_cut_structures:
		#	for e_surface in ltd_cut_structure.get_e_surfaces(cutkosky_cut_momenta_num):
		#		if e_surface != []:
		#			e_surface['ltd_cut_structure'] = ltd_cut_structure.ltd_cut_structure_index
		#			e_surfaces += [e_surface]
		for ltd_cut_structure in self.ltd_cut_structures:
			e_surfaces += ltd_cut_structure.get_e_surfaces(cutkosky_cut_momenta_num)
		return e_surfaces

class LTDCutStructure(SubGraph):
	def __init__(self, **opts):
		super(LTDCutStructure, self).__init__(**opts)
		self.opts = opts
		self.ltd_cut_structure = self.opts['ltd_cut_structure']
		self.ltd_cut_signature = [line for line in self.ltd_cut_structure if abs(line) == 1]
		self.basis_trsf = numpy.array([loop_line.signature
									for loop_line,line in zip(self.subgraph_dump.loop_lines,self.ltd_cut_structure)
									if abs(line)==1]) #cut line momenta from loop momenta
		self.inv_basis_trsf = numpy.linalg.inv(self.basis_trsf) #loop momenta from cut line momenta

	def generate(self):
		ltd_cut_propagator_indices = self.get_prop_cuts()
		self.ltd_cuts = [LTDCut(ltd_cut_indices=ltd_cut_indices, **self.opts) for ltd_cut_indices in ltd_cut_propagator_indices]
		for ltd_cut in self.ltd_cuts:
			ltd_cut.generate()

	def get_prop_cuts(self):
		prop_cuts = []
		cut_indices = [i for i,line in enumerate(self.ltd_cut_structure) if abs(line)==1]
		possible_cuts = [range(len(loop_line.propagators)) for loop_line,line in zip(self.subgraph_dump.loop_lines,self.ltd_cut_structure) if abs(line)==1]
		for prop_cut in itertools.product(*possible_cuts):
			prop_cuts += [prop_cut]
		return prop_cuts

	def get_e_surfaces(self,cutkosky_cut_momenta_num):
		e_surfaces = []
		#for ltd_cut in self.ltd_cuts:
		#	for e_surface in ltd_cut.get_e_surfaces(cutkosky_cut_momenta_num):
		#		if e_surface != []:
		#			e_surface['ltd_cut'] = ltd_cut.ltd_cut_indices
		#			e_surfaces += [e_surface]
		for ltd_cut in self.ltd_cuts:
			e_surfaces += ltd_cut.get_e_surfaces(cutkosky_cut_momenta_num)
		return e_surfaces

class LTDCut(LTDCutStructure):
	def __init__(self, **opts):
		super(LTDCut, self).__init__(**opts)
		self.opts = opts
		self.ltd_cut_indices = self.opts['ltd_cut_indices']
		cut_indices_iter = iter(self.ltd_cut_indices)
		self.ltd_cut_momenta = [loop_line.propagators[next(cut_indices_iter)].name for line, loop_line in zip(self.ltd_cut_structure,self.subgraph_dump.loop_lines) if abs(line) != 0]
		# shifts of the LTD cut basis wrt to the cutkosky cuts and external momenta
		cut_indices_iter = iter(self.ltd_cut_indices)
		self.cutkosky_shift_trsf = numpy.array([numpy.array(loop_line.propagators[next(cut_indices_iter)].parametric_shift[0]) for line, loop_line in zip(self.ltd_cut_structure,self.subgraph_dump.loop_lines) if abs(line) != 0])
		cut_indices_iter = iter(self.ltd_cut_indices)
		self.external_shift_trsf = numpy.array([numpy.array(loop_line.propagators[next(cut_indices_iter)].parametric_shift[1]) for line, loop_line in zip(self.ltd_cut_structure,self.subgraph_dump.loop_lines) if abs(line) != 0])
		self.supergraph_cut_basis = [name for name in self.cutkosky_cut_momenta[:-1]]+[name for name in self.ltd_cut_momenta]
		# supergraph loop basis to supergraph cut basis
		self.supergraph_basis_trsf = numpy.array([cutkosky_cut_signature for cutkosky_cut_signature in self.cutkosky_cut_signatures[:-1]]
												+ [signature for (signature,shift) in self.subgraph_dump.loop_momentum_map])
		self.supergraph_shift_trsf = numpy.array([cutkosky_cut_shift for cutkosky_cut_shift in self.cutkosky_cut_shifts[:-1]]
												+ [shift for (signature,shift) in self.subgraph_dump.loop_momentum_map])

		self.supergraph_inv_basis_trsf = numpy.linalg.inv(self.supergraph_basis_trsf)

	def get_subgraph_basis_shift(self, cutkosky_cut_momenta_num):
		return self.cutkosky_shift_trsf.dot(cutkosky_cut_momenta_num)+self.external_shift_trsf.dot(self.external_momenta_num)

	def generate(self):
		self.ltd_loop_lines = [LTDLoopLine(loop_line_index=i,**self.opts) for i,loop_line in enumerate(self.subgraph_dump.loop_lines)]
		for ltd_loop_line in self.ltd_loop_lines:
			ltd_loop_line.generate()

	def get_e_surfaces(self,cutkosky_cut_momenta_num):
		e_surfaces = []
		#for ltd_loop_line in self.ltd_loop_lines:
		#	for e_surface in ltd_loop_line.get_e_surfaces(cutkosky_cut_momenta_num):
		#		if e_surface != []:
		#			e_surface['loop_line'] = ltd_loop_line.loop_line_index
		#			e_surfaces += [e_surface]
		for ltd_loop_line in self.ltd_loop_lines:
			e_surfaces += ltd_loop_line.get_e_surfaces(cutkosky_cut_momenta_num)
		return e_surfaces

class LTDLoopLine(LTDCut):
	def __init__(self,**opts):
		super(LTDLoopLine, self).__init__(**opts)
		self.opts = opts
		self.loop_line_index = self.opts['loop_line_index']
		this_loop_line = self.subgraph_dump.loop_lines[self.loop_line_index]
		# signature in cut basis
		self.basis_signature = numpy.array(this_loop_line.signature).dot(self.inv_basis_trsf)
		# signature when multiplied with cut structure signs
		self.energy_signature = self.basis_signature.dot(numpy.diag(self.ltd_cut_signature))
		self.has_positive_signature = all(sign >= 0 for sign in self.energy_signature)
		self.has_negative_signature = all(sign <= 0 for sign in self.energy_signature)

	def generate(self):
		this_loop_line = self.subgraph_dump.loop_lines[self.loop_line_index]
		self.ltd_propagators = [LTDPropagator(propagator_index = a, **self.opts) for a,propagator in enumerate(this_loop_line.propagators)]

	def get_e_surfaces(self,cutkosky_cut_momenta_num):
		e_surfaces = []
		for ltd_propagator in self.ltd_propagators:
		#	if ltd_propagator.is_ellipsoid(cutkosky_cut_momenta_num):
		#		e_surfaces += [{'propagator': ltd_propagator.propagator_index}]
			if ltd_propagator.is_ellipsoid(cutkosky_cut_momenta_num):
				e_surfaces += [ltd_propagator]
		return e_surfaces

class LTDPropagator(LTDLoopLine):
	def __init__(self,**opts):
		super(LTDPropagator, self).__init__(**opts)
		self.opts = opts
		self.propagator_index = self.opts['propagator_index']
		self.propagator_name = self.subgraph_dump.loop_lines[self.loop_line_index].propagators[self.propagator_index].name
		cut_indices_iter = iter(self.ltd_cut_indices)
		full_cut_indices = [next(cut_indices_iter) if abs(line) != 0 else None for line in self.ltd_cut_structure]
		self.is_cut = (full_cut_indices[self.loop_line_index] == self.propagator_index)
		self.dimension = numpy.sum([1 for sign in self.energy_signature if sign != 0])*3 - 1
		#shift of propagator momentum with respect to the loop line momentum
		self.own_cutkosky_shift_trsf = self.subgraph_dump.loop_lines[self.loop_line_index].propagators[self.propagator_index].parametric_shift[0]
		self.own_external_shift_trsf = self.subgraph_dump.loop_lines[self.loop_line_index].propagators[self.propagator_index].parametric_shift[1]
		#total shift i.e. affine surface term
		self.total_cutkosky_shift_trsf = self.own_cutkosky_shift_trsf - self.basis_signature.dot(self.cutkosky_shift_trsf)
		self.total_external_shift_trsf = self.own_external_shift_trsf - self.basis_signature.dot(self.external_shift_trsf)
		self.tolerance = 1e-10

	def get_total_shift(self, cutkosky_cut_momenta_num):
		total_cutkoksy_shift = vectors.LorentzVector(self.total_cutkosky_shift_trsf.dot(cutkosky_cut_momenta_num))
		total_external_shift = vectors.LorentzVector(self.total_external_shift_trsf.dot(self.external_momenta_num))
		return total_cutkoksy_shift+total_external_shift
	
	def is_ellipsoid(self, cutkosky_cut_momenta_num):
		total_shift = self.get_total_shift(cutkosky_cut_momenta_num)
		s = total_shift.square()
		is_plus_ellipsoid = (self.has_positive_signature and total_shift[0] < 0.) and (s > 0.)
		is_minus_ellipsoid = (self.has_negative_signature and total_shift[0] > 0.) and (s > 0.)
		is_ellipsoid = (is_plus_ellipsoid or is_minus_ellipsoid)
		return is_ellipsoid

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

	def probe_e_surface(self, cutkosky_cut_momenta_num, random_variables):
		assert(self.is_ellipsoid(cutkosky_cut_momenta_num))
		assert(len(random_variables)==self.dimension)
		title = ('E_surface %s'%self.propagator_name + ' (%i dim.)'%(self.dimension)
				+ ' in LTD cut (%s)'%','.join([name+('%+d'%sign)[0] 
					for name, sign in zip(self.ltd_cut_momenta,self.ltd_cut_signature)])
				+ ' in Cutkosky cut (%s)'%','.join([name+('%+d'%sign)[0] 
					for name, sign in zip(self.cutkosky_cut_momenta, self.cutkosky_cut_energy_signs)]))
		print(10*'='+' probing ' + title + ' ' +10*'=')
		os_subgraph_cut_basis = self.get_k_non_pinched(cutkosky_cut_momenta_num, random_variables)
		print(os_subgraph_cut_basis)
		os_subgraph_loop_basis = self.inv_basis_trsf.dot(os_subgraph_cut_basis - self.get_subgraph_basis_shift(cutkosky_cut_momenta_num))
		print(os_subgraph_loop_basis)
		os_supergraph_cut_basis = numpy.array([momentum for momentum in cutkosky_cut_momenta_num[:-1]]+[momentum for momentum in os_subgraph_loop_basis])
		print(os_supergraph_cut_basis)
		os_supergraph_loop_basis = self.supergraph_inv_basis_trsf.dot(os_supergraph_cut_basis - self.supergraph_shift_trsf.dot(self.external_momenta_num))
		os_supergraph_loop_basis = [vectors.LorentzVector(p) for p in os_supergraph_loop_basis]
		print(os_supergraph_loop_basis)
		spatial_os_supergraph_loop_basis = [p.space() for p in os_supergraph_loop_basis]
		prober = CrossSectionSingularityProber(squared_topology,"hyperparameters.yaml")
		#directions = [numpy.array([5.,-2.,4.])]*squared_topology.topo.n_loops
		directions = [numpy.array([3.,-2.,1.]),numpy.array([-2.,-5.,6.]),numpy.array([1.,2.,3.])]
		directions = [direction/numpy.linalg.norm(direction) for direction in directions]
		prober.probe_singular_surface(spatial_os_supergraph_loop_basis,directions,plot=True,plot_title=title,precision='f128')
		#prober.evaluate_cuts(spatial_os_supergraph_loop_basis)


	def get_k_non_pinched(self, cutkosky_cut_momenta_num, random_variables):
		assert(len(random_variables)==self.dimension)
		os_subgraph_cut_basis = [vectors.LorentzVector([0,0,0,0]) for i in range(len(self.ltd_cut_momenta))]
		prop_name = self.propagator_name
		n_independent_cuts = len([1 for sign in self.basis_signature if sign != 0])
		counter = 0
		indep_counter = 0
		for index,(name,energy_sign,signature) in enumerate(zip(self.ltd_cut_momenta, self.ltd_cut_signature, self.basis_signature)):
			if signature != 0:
				indep_counter += 1
				mass_squared_basis = self.squared_topology_dump.masses[name] if name in self.squared_topology_dump.masses else 0.
				mass_prop = self.squared_topology_dump.masses[prop_name] if prop_name in self.squared_topology_dump.masses else 0.
				mass_prop += numpy.sum([self.squared_topology_dump.masses[name]
											if (name in self.squared_topology_dump.masses and sign != 0) else 0.
											for sign,((name,_)) in zip(self.basis_signature[:index],self.ltd_cut_momenta[:index])])
				mass_squared_prop = mass_prop**2
				shift = self.get_total_shift(cutkosky_cut_momenta_num)
				shift += self.basis_signature[:index].dot(os_subgraph_cut_basis[:index])
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
				basis_momentum_boosted = self.get_lorentz_vector(energy_sign*numpy.sqrt(radius**2 + mass_squared_basis),
											self.get_space_vector(radius,
												self.get_unit_vector(cos_theta,phi)))
				basis_momentum = vectors.LorentzVector(self.get_boost(-self.get_beta(shift)).dot(basis_momentum_boosted))
				os_subgraph_cut_basis[index] = basis_momentum
		return os_subgraph_cut_basis

if __name__ == "__main__":

	#print(dir(ltd.CrossSection))

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

	squared_topology = mercedes

	CUT_INDEX = 2
	print(squared_topology.cuts)
	sqtopo = SquaredTopology(squared_topology_dump=squared_topology)
	sqtopo.generate()

	parameteriser = CutkoskyCutParametiser(squared_topology,CUT_INDEX)
	parameteriser.get_cut_basis()
	cut_momenta = parameteriser.get_PS_point([random.random() for i in range(3*len(squared_topology.cuts[CUT_INDEX])-4)])
	#print(cut_momenta)
	cutkosky_cut_momenta_num = [cut_momenta[name] for (name,_) in squared_topology.cuts[CUT_INDEX]]
	#print(cutkosky_cut_momenta_num)
	# BUILD: sqtopo.cutkosky_cuts[CUT_INDEX].get_PS_point([random.random() for i in range(3*len(squared_topology.cuts[CUT_INDEX])-4)])
	e_surfaces = sqtopo.cutkosky_cuts[CUT_INDEX].get_e_surfaces(cutkosky_cut_momenta_num)
	#print(e_surfaces[0].dimension)
	for e_surface in e_surfaces:
		print(e_surface)
		print('dimension: ', e_surface.dimension)
		random_variables = [random.random() for i in range(e_surface.dimension)]
		e_surface.probe_e_surface(cutkosky_cut_momenta_num, random_variables)
	#print(sqtopo.get_e_surfaces())
	plt.show()
	stop

	#for cut_index, cut in enumerate(squared_topology.cuts):
	#	print(cut_index)
	#	parameteriser = CutkoskyCutParametiser(squared_topology,cut_index)
	#	#parameteriser.get_cut_basis()
	#	parameteriser.get_PS_point([random.random() for i in range(3*len(squared_topology.cuts[cut_index])-4)])

	print(squared_topology.cuts[1])
	print(squared_topology.cut_signatures[1])
	stop
	#print(squared_topology.external_momenta)
	#print(squared_topology.topo.loop_momenta)
	CUT_INDEX = 4
	#print(squared_topology.cuts[CUT_INDEX])
	#print(squared_topology.loop_topologies[CUT_INDEX][1])
	#print(squared_topology.cut_signatures[CUT_INDEX])
	
	#print(squared_topology.masses)

	parameteriser = CutkoskyCutParametiser(squared_topology,CUT_INDEX)
	parameteriser.get_cut_basis()
	cut_momenta = parameteriser.get_PS_point([random.random() for i in range(3*len(squared_topology.cuts[CUT_INDEX])-4)])



	# set len(cuts[CUT_INDEX])-1 spatial cut momenta (randomly!)
	#arbitrary_p = vectors.Vector([1.,3.,6.])
	#cut_momenta = {}
	#for name, sign in squared_topology.cuts[CUT_INDEX][:-1]:
	#	cut_momenta[name] = arbitrary_p
	#print(cut_momenta)
	
	prober = CrossSectionSingularityProber(squared_topology,"hyperparameters.yaml")
	#directions = [numpy.array([5.,-2.,4.])]*squared_topology.topo.n_loops
	directions = [numpy.array([3.,-2.,1.]),numpy.array([-2.,-5.,6.]),numpy.array([1.,2.,3.])]
	normalised_directions = [direction/numpy.linalg.norm(direction) for direction in directions]
	prober.probe_singularities_in_cut_topology(CUT_INDEX,cut_momenta,normalised_directions,plot=True,precision='f128')








