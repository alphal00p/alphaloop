# n loop LTD
import sys
import numpy as numpy
import vegas
from adipy import adipy
import time
import itertools
import mpmath
import math
import ltd_commons
import topologies
#import manually_crafted_topologies as topologies
import random
from numpy.linalg import inv
import pycuba
from scipy import integrate

#import warnings
#warnings.simplefilter("error")
#warnings.simplefilter("ignore", DeprecationWarning)

class LTDnLoop:

	def __init__(self,topology,hyperparameters):
		self.ltd_cut_structure  = topology.ltd_cut_structure
		self.ltd_loop_lines		= topology.loop_lines
		self.n_loops            = topology.n_loops 
		self.name               = topology.name
		self.external_kinematics= topology.external_kinematics
		self.scale				= self.get_alt_scale()
		if topology.analytic_result is None:
			self.analytical_result 	= 0.
		else:
			self.analytical_result 	= topology.analytic_result
		self.cut_structures 	= [CutStructure(cut_structure=ltd_cut_structure, ltd_loop_lines = self.ltd_loop_lines) for ltd_cut_structure in self.ltd_cut_structure]
		self.lambda_ij			= hyperparameters['Deformation']['scaling']['lambda']
		self.a_ij				= hyperparameters['Deformation']['additive']['a_ij']
		self.surfaces 			= self.get_surfaces()
		self.abs_min_lambda		= 1.
		self.abs_min_distance	= 1e10
		self.n_channels			= numpy.sum([len(cut_structure.cuts) for cut_structure in self.cut_structures])

		print('{:=^90}'.format(' '+ self.name + ' '))
		print('{:=^90}'.format(' '+ 'scale = %e' %self.scale + ' '))
		print('{:=^90}'.format(' analytical result = {num.real} + {num.imag}j '.format(num=self.analytical_result)))

	def get_surfaces(self):
		surfaces = []
		for cut_structure in self.cut_structures:
			for cut in cut_structure.cuts:
				for loop_line in cut.loop_lines:
					for propagator in loop_line.propagators:
						new_surfaces = propagator.get_surfaces()
						for surface in new_surfaces:
							surface = self.normalize_surface(surface)
							surfaces += [surface]
		surfaces = sorted(surfaces, key = lambda surf: (surf['is_ellipsoid'],len(surf['deltas']),[(delta['loop_line_index'],delta['propagator_index']) for delta in surf['deltas']]))
		surfaces = self.group_surfaces(surfaces)
		return surfaces

	def normalize_surface(self,surface):
		sorted_deltas = sorted(surface['deltas'], key=lambda x: x['loop_line_index'], reverse=False)
		#python 2: sorted_deltas = sorted(surface['deltas'])#, key=lambda delta: (delta['loop_line_index'],delta['propagator_index']))
		ll_i = sorted_deltas[-1]['loop_line_index']
		prop_i = sorted_deltas[-1]['propagator_index']
		signatures = numpy.array([self.ltd_loop_lines[delta['loop_line_index']].signature for delta in sorted_deltas])
		inv_trsf = numpy.linalg.pinv(signatures[:-1])
		shift_space = [self.ltd_loop_lines[delta['loop_line_index']].propagators[delta['propagator_index']].q[1:] for delta in sorted_deltas]
		p_shift_space = -signatures[-1].dot(inv_trsf).dot(shift_space[:-1]) + shift_space[-1]
		surface['p_shift'][1:] = p_shift_space
		if surface['p_shift'][0] > 0.:
			for i,sorted_delta in enumerate(sorted_deltas):
				sorted_deltas[i]['sign'] *= -1
			surface['p_shift'][0] = -surface['p_shift'][0]
		elif surface['p_shift'][0] == 0.:
			if sorted_deltas[0]['sign'] == -1:
				for i,sorted_delta in enumerate(sorted_deltas):
					sorted_deltas[i]['sign'] *= -1
		surface['deltas'] = sorted_deltas
		return surface

	def group_surfaces(self,surfaces):
		""" groups identical surfaces, does not yet group linearly dependent hyperboloids"""
		grouped_surfaces = []
		for i,surface in enumerate(surfaces):
			add = True
			for j,grouped_surface in enumerate(grouped_surfaces):
				if surface['deltas'] == grouped_surface['deltas']:
					if surface['p_shift'] == grouped_surface['p_shift']:
						grouped_surfaces[j]['multiplicity'] += 1
						surfaces[i]['multiplicity'] += 1
						add = False
			if add:
				grouped_surfaces += [surface]
		return grouped_surfaces

	def display_surfaces(self,show_ellipsoids=True,show_hyperboloids=True):
		surfaces = self.get_surfaces()
		ellipsoid_multi = [surface['multiplicity'] for surface in surfaces if surface['is_ellipsoid'] == True]
		hyperboloid_multi = [surface['multiplicity'] for surface in surfaces if surface['is_ellipsoid'] == False]
		n_ellipsoids = len(ellipsoid_multi)
		n_hyperboloids = len(hyperboloid_multi)
		if show_ellipsoids and show_hyperboloids:
			title = '%i surfaces (%i ellipsoids, %i hyperboloids)' %(n_ellipsoids+n_hyperboloids,n_ellipsoids,n_hyperboloids)
		elif show_ellipsoids:
			title = '%i ellipsoids / %i surfaces' %(n_ellipsoids,n_ellipsoids+n_hyperboloids)
		elif show_hyperboloids:
			title = '%i hyperboloids/ %i surfaces' %(n_hyperboloids,n_ellipsoids+n_hyperboloids)
		print('{:=^90}'.format(' '+title+' '))
		digits = 7
		title_element_width = 24
		title_prop_str = '{:^{width}}'.format('loopline & propagator',width=title_element_width)
		title_prop_str += '{:^{width}}'.format('signature',width=int(title_element_width/2))
		title_prop_str += '{:^{width}}'.format('type',width=int(title_element_width/(self.n_loops+1)))
		title_prop_str += '{:^{width}}'.format('multiplicity',width=int(title_element_width/2))
		title_prop_str += '{:^{width}}'.format('p_shift',width=4*digits+5)
		if (show_ellipsoids and n_ellipsoids>0) or (show_hyperboloids and n_hyperboloids>0):
			print(title_prop_str)
		for surface in surfaces:
			if (surface['is_ellipsoid'] and not show_ellipsoids) or (not surface['is_ellipsoid'] and not show_hyperboloids):
					continue
			if surface['is_ellipsoid']:
				type_str = '{:^{width}}'.format('E', width = int(title_element_width/(self.n_loops+1)))
			else:
				type_str = '{:^{width}}'.format('H', width = int(title_element_width/(self.n_loops+1)))
			prop_arr = ['{:^{width}}'.format('({} {})'.format(	delta['loop_line_index'],delta['propagator_index']),
																width = int(title_element_width/len(surface['deltas']))) for delta in surface['deltas']]
			sign_arr = ['+' if delta['sign'] == 1. else '-' for delta in surface['deltas']]
			sign_arr = ['{:^{width}}'.format(sign, width = int(title_element_width/2/len(sign_arr))) for sign in sign_arr]
			multiplicity = '{:^{width}}'.format(float(surface['multiplicity']), width = int(title_element_width/2))
			p_shift_str = '[{0:{digits}},{1:{digits}},{2:{digits}},{3:{digits}}]'.format(*surface['p_shift'],digits=digits)
			print(''.join(prop_arr) + ''.join(sign_arr) + type_str + multiplicity + p_shift_str)

	def delta(self,q_space,m):
		on_f_shell = adipy.sqrt(self.norm(q_space)**2 + m**2)
		return on_f_shell
	
	def get_rust_scale(self):
		if len(self.external_kinematics) > 2:
			p = self.external_kinematics[0] + self.external_kinematics[1]
		else:
			p = self.external_kinematics[0]
		scale = numpy.sqrt(abs(p[0]**2 - self.norm(p[1:])**2))
		return scale

	def get_alt_scale(self):
		p_posE = [p for p in self.external_kinematics if p[0] > 0.]
		if -numpy.sum(self.external_kinematics,axis=0)[0] > 0.:
			p_posE += [-sum(self.external_kinematics)]
		if (len(p_posE) == 0):
			print("All incoming momentum energies are zero!! Improvising scale...")
			scale = self.norm(numpy.sum(self.external_kinematics[1:],axis=0))**2
		else:
			p_sum = sum(p_posE)
			s = p_sum[0]**2 - self.norm(p_sum[1:])**2
			scale = numpy.sqrt(abs(s))
		return scale
	
	def norm(self,q):
		"""for complex and real q arrays
		not the same as numpy.linalg.norm for complex numbers!!"""
		return adipy.sqrt(numpy.sum([q_i*q_i for q_i in q]))
		
	def parametrize_analytic(self,random):
		wgt = 1.
		radius = self.scale*random[0]/(1.-random[0])
		wgt *= (self.scale+radius)**2/self.scale
		assert(len(random) > 1)
		phi = 2*numpy.pi*random[1]
		wgt *= 2*numpy.pi
		if len(random)==3:	
			cos_theta = -1.+2.*random[2]
			wgt *= 2.
			sin_theta = numpy.sqrt(1.-cos_theta**2)
			l_space = radius*numpy.array([sin_theta*numpy.cos(phi), sin_theta*numpy.sin(phi),cos_theta])
			wgt *= radius**2 #spherical coord
		elif len(random)==2:
			l_space = radius*numpy.array([numpy.cos(phi), numpy.sin(phi)])
			wgt *= radius
		return l_space, wgt

	def get_deltas(self,loop_momenta):
		""" computes the sqrt[(l_space+p_space)**2 + m**2] of all looplines
			looks like looplines_prop = [[q_i0p1],[q_i0p2],[q_i0p3]]"""
		deltas = []
		for loop_line in self.ltd_loop_lines:
			loop_momentum = numpy.sum([l*sign for l,sign in zip(loop_momenta,loop_line.signature)],axis=0)
			deltas += [[self.delta(loop_momentum + prop.q[1:],numpy.sqrt(prop.m_squared)) for prop in loop_line.propagators]]
		return deltas
		
	def dual(self,loop_momenta):
		vec = numpy.array(loop_momenta).flatten()
		dim = len(vec)
		dual_vec = [adipy.ad(vec[i], numpy.array([0. if j !=i  else 1. for j in range(dim)])) for i in range(dim)]
		dual_momenta = [dual_vec[i*3:(i+1)*3] for i in range(self.n_loops)]
		return numpy.array(dual_momenta)

	def rescale_with_expansion_parameter(self,kappas,loop_momenta,real_deltas):
		self.curr_min_lambda = adipy.ad(1.,numpy.array([0. for i in range(self.n_loops*3)]))
		all_prop_deltas = [{'loop_line_index': loop_line_index,'propagator_index': propagator_index}
									for loop_line_index,loop_line in enumerate(self.ltd_loop_lines)
									for propagator_index,propagator in enumerate(loop_line.propagators)]
		#all_prop_deltas = []
		#for cut_structure in self.cut_structures:
		#	for cut in cut_structure.cuts:
		#		for loop_line in cut.loop_lines:
		#			for propagator in loop_line.propagators:
		#				delta = {'loop_line_index': loop_line.loop_line_index, 'propagator_index': propagator.propagator_index}
		#				all_prop_deltas += [delta]
		#for i,surface in enumerate(self.surfaces):
		#	for delta in surface['deltas']:
		for delta in all_prop_deltas:
				shift_space = self.ltd_loop_lines[delta['loop_line_index']].propagators[delta['propagator_index']].q[1:]
				signature = numpy.array(self.ltd_loop_lines[delta['loop_line_index']].signature)
				prop_momentum = numpy.array(signature.dot(loop_momenta) + shift_space)
				prop_kappa = signature.dot(kappas)
				prop_delta = real_deltas[delta['loop_line_index']][delta['propagator_index']]
				#lambda_j= abs(expansion_parameter*prop_delta/prop_kappa.dot(prop_momentum))
				if prop_kappa.dot(prop_kappa) > 0.:
					lambda_j_sq = 2*prop_delta**4/(prop_kappa.dot(prop_kappa)*prop_delta**2-prop_kappa.dot(prop_momentum)**2)
					if lambda_j_sq > 0.:
						lambda_j = adipy.sqrt(lambda_j_sq)/10.
					else:
						lambda_j = self.curr_min_lambda
				else:
					lambda_j = self.curr_min_lambda
				if lambda_j.nom < self.curr_min_lambda.nom:
					self.curr_min_lambda = lambda_j
		if self.curr_min_lambda.nom < self.abs_min_lambda:
			self.abs_min_lambda = self.curr_min_lambda.nom
			print(self.abs_min_lambda)
		kappas = numpy.array([[self.curr_min_lambda*k for k in kappa] for kappa in kappas])
		return kappas

	def rescale_with_branch_cut_condition(self,kappas,loop_momenta):
		self.curr_min_lambda = adipy.ad(1.,numpy.array([0. for i in range(self.n_loops*3)]))
		all_prop_deltas = [{'loop_line_index': loop_line_index,'propagator_index': propagator_index}
									for loop_line_index,loop_line in enumerate(self.ltd_loop_lines)
									for propagator_index,propagator in enumerate(loop_line.propagators)]
		for delta in all_prop_deltas:
			shift_space = self.ltd_loop_lines[delta['loop_line_index']].propagators[delta['propagator_index']].q[1:]
			prop_m_squared = self.ltd_loop_lines[delta['loop_line_index']].propagators[delta['propagator_index']].m_squared
			signature = numpy.array(self.ltd_loop_lines[delta['loop_line_index']].signature)
			prop_momentum = numpy.array(signature.dot(loop_momenta) + shift_space)
			prop_kappa = signature.dot(kappas)
			if prop_kappa.dot(prop_kappa) > 0.:
				lambda_j_sq = (prop_momentum.dot(prop_momentum)+prop_m_squared)/prop_kappa.dot(prop_kappa)
				if lambda_j_sq > 0:
					lambda_j = adipy.sqrt(lambda_j_sq)/10.
				else:
					lambda_j = self.curr_min_lambda
			else:
				lambda_j = self.curr_min_lambda
			if lambda_j.nom < self.curr_min_lambda.nom:
				self.curr_min_lambda = lambda_j
		if self.curr_min_lambda.nom < self.abs_min_lambda:
			self.abs_min_lambda = self.curr_min_lambda.nom
			print(self.abs_min_lambda)
		kappas = numpy.array([[self.curr_min_lambda*k for k in kappa] for kappa in kappas])
		return kappas

	def rescale_from_complex_zeros(self,kappas,loop_momenta,real_deltas):
		curr_min = adipy.ad(1.,numpy.array([0. for i in range(self.n_loops*3)]))
		for surface in self.surfaces:
			# solve (linearized_surface(lambda) = 0) <=> A*lambda^2 + i*B*lambda + C = 0
			A = B = C = 0.
			for delta in surface['deltas']:
				shift_space = self.ltd_loop_lines[delta['loop_line_index']].propagators[delta['propagator_index']].q[1:]
				signature = numpy.array(self.ltd_loop_lines[delta['loop_line_index']].signature)
				prop_momentum = numpy.array(signature.dot(loop_momenta) + shift_space)
				prop_kappa = signature.dot(kappas)
				prop_delta = real_deltas[delta['loop_line_index']][delta['propagator_index']]
				C += prop_delta*delta['sign']
				B += prop_momentum.dot(prop_kappa)/prop_delta*delta['sign']
				A += 0.5/prop_delta*(-prop_kappa.dot(prop_kappa)+0.5*(prop_momentum.dot(prop_kappa))**2/prop_delta**2)
			C += surface['p_shift'][0]
			if A*A > 1e-20:
				sqrt_X = -B/(2.*A)
				X = sqrt_X**2
				Y = -C/A
				scaling_param = self.scaling_condition(X,Y)
				if scaling_param.nom < curr_min.nom:
					curr_min = scaling_param
		if curr_min.nom < self.abs_min_lambda:
			self.abs_min_lambda = curr_min.nom
			print('abs_min_lambda: ', self.abs_min_lambda)
		kappas = numpy.array([[curr_min*k for k in kappa] for kappa in kappas])
		return kappas

	def estimate_distance_to_ellipsoid(self,this_surface,kappas,loop_momenta,real_deltas):
		curr_min = adipy.ad(10,numpy.array([0. for i in range(self.n_loops*3)]))
		for i,kappa_i in enumerate(kappas):
			if kappa_i.dot(kappa_i) > 0:
				for surface in self.surfaces:
					A = B = C = 0.
					if surface['is_ellipsoid'] and this_surface != surface:
						for delta in surface['deltas']:
							shift_space = self.ltd_loop_lines[delta['loop_line_index']].propagators[delta['propagator_index']].q[1:]
							signature = numpy.array(self.ltd_loop_lines[delta['loop_line_index']].signature)
							prop_momentum = numpy.array(signature.dot(loop_momenta) + shift_space)
							prop_delta = real_deltas[delta['loop_line_index']][delta['propagator_index']]
							if signature[i] != 0:
								A += 0.5/prop_delta*(kappa_i.dot(kappa_i)-0.5*(kappa_i.dot(prop_momentum)/prop_delta)**2)
								B += signature[i]*kappa_i.dot(prop_momentum)/prop_delta
							C += prop_delta #add for every prop delta, independent of signature[i] (outside loop!)
						C += surface['p_shift'][0]
					determinant = B**2-4.*A*C
					if A*A > 1e-20 and determinant.nom > 0:
						sqrt_det = adipy.sqrt(B**2-4.*A*C)
						line_parameter_p = (-B+sqrt_det)/(2.*A)
						line_parameter_m = (-B-sqrt_det)/(2.*A)
						for line_parameter in [line_parameter_p,line_parameter_m]:
							if line_parameter.nom > 0. and line_parameter.nom < curr_min.nom:
								curr_min = line_parameter
		distance_vector = numpy.array([[curr_min*k for k in kappa] for kappa in kappas])
		min_distance = self.norm(distance_vector)
		if min_distance.nom < self.abs_min_distance:
			self.abs_min_distance = min_distance.nom
			print('abs_min_distance: ', self.abs_min_distance)
		return min_distance

	def scaling_condition(self,X,Y):
		if Y > 2.*X:
			scaling_param_ij_sq = .25*Y
		elif Y < 0.:
			scaling_param_ij_sq = X - .5*Y
		else:
			scaling_param_ij_sq = X - .25*Y
		scaling_param_ij = adipy.sqrt(scaling_param_ij_sq)
		return scaling_param_ij

	def get_kappa(self,loop_momenta,real_deltas,surface):
		kappas = numpy.array([[coord*0. for coord in loop_momentum] for loop_momentum in loop_momenta])
		if surface['is_ellipsoid']:
			my_cuts = [(delta['loop_line_index'], delta['propagator_index']) for delta in surface['deltas']]
			all_signatures = numpy.array([self.ltd_loop_lines[delta['loop_line_index']].signature for delta in surface['deltas']])
			shift_space = [self.ltd_loop_lines[delta['loop_line_index']].propagators[delta['propagator_index']].q[1:] for delta in surface['deltas']]
			surf_momenta = all_signatures.dot(loop_momenta) + shift_space
			surf_deltas = [real_deltas[surf_delta['loop_line_index']][surf_delta['propagator_index']] for surf_delta in surface['deltas']]
			for index in range(self.n_loops):
				kappa_index = numpy.sum([signature[index]*surf_momenta[i]/surf_deltas[i] for i,signature in enumerate(all_signatures)],axis=0)
				#if kappa_index.dot(kappa_index) > 0.:
				#	kappa_index *= 1./self.norm(kappa_index)
				kappas[index] = kappa_index
			min_distance = self.estimate_distance_to_ellipsoid(surface,kappas,loop_momenta,real_deltas)
			surface_eq = numpy.sum([surf_deltas[i]*delta['sign'] for i,delta in enumerate(surface['deltas'])])+surface['p_shift'][0]
			#interpolation = adipy.exp(-surface_eq**2/self.scale**2/self.a_ij)
			interpolation = adipy.exp(-surface_eq**2/min_distance**2)
			#other_interpolations = []
			#for other_surface in self.surfaces:
			#	if surface != other_surface:
			#		if other_surface['is_ellipsoid']:
			#			other_surf_deltas = [real_deltas[surf_delta['loop_line_index']][surf_delta['propagator_index']] for surf_delta in other_surface['deltas']]
			#			other_surface_eq = numpy.sum([other_surf_deltas[i]*delta['sign'] for i,delta in enumerate(other_surface['deltas'])])+other_surface['p_shift'][0]
			#			other_interpolations += [1.-adipy.exp(-other_surface_eq**2/self.scale**2/self.a_ij)]
			#kappas *= numpy.prod(other_interpolations)
			kappas *= -interpolation#*surface['multiplicity']
		return kappas

	def deform(self,loop_momenta,):
		loop_momenta = self.dual(loop_momenta)
		real_deltas = self.get_deltas(loop_momenta)
		kappas = numpy.array([[coord*0. for coord in loop_momentum] for loop_momentum in loop_momenta])
		for i,surface in enumerate(self.surfaces):
			kappas += self.get_kappa(loop_momenta,real_deltas,surface)
		kappas *= self.scale
		#kappas = numpy.array([[k*self.norm(loop_momentum) for k in kappa] for kappa,loop_momentum in zip(kappas,loop_momenta)])
		if self.lambda_ij > 0:
			#kappas = self.rescale_with_branch_cut_condition(kappas,loop_momenta)
			#kappas = self.rescale_with_expansion_parameter(kappas,loop_momenta,real_deltas)
			kappas = self.rescale_from_complex_zeros(kappas,loop_momenta,real_deltas)
		else:
			kappas *= abs(self.lambda_ij)
		full_jac = adipy.jacobian(loop_momenta.flatten() + 1j*kappas.flatten())
		det_jac = numpy.linalg.det(full_jac)
		kappas = [numpy.array([kapp.nom for kapp in kappa]) for kappa in kappas]
		return kappas, det_jac

	def deform_on_prop(self,loop_momenta):
		loop_momenta = self.dual(loop_momenta)
		real_deltas = self.get_deltas(loop_momenta)
		kappas = numpy.array([[coord*0. for coord in loop_momentum] for loop_momentum in loop_momenta])
		opts = {'scale': self.scale, 'lambda_ij': self.lambda_ij, 'a_ij': self.a_ij}
		for cut_structure in self.cut_structures:
			for cut in cut_structure.cuts:
				for loop_line in cut.loop_lines:
					for propagator in loop_line.propagators:
						kappas += propagator.deform(loop_momenta,real_deltas,**opts)
		full_jac = adipy.jacobian(loop_momenta.flatten() + 1j*kappas.flatten())
		det_jac = numpy.linalg.det(full_jac)
		kappas = [numpy.array([kapp.nom for kapp in kappa]) for kappa in kappas]
		return kappas, det_jac

	def ltd_integrand(self,x,multi_channel_index=None):
		assert(len(x) == self.n_loops*3)
		loop_momenta = [numpy.zeros(3) for i in range(self.n_loops)]
		wgt = numpy.zeros(self.n_loops)
		for i in range(self.n_loops):
			loop_momenta[i], wgt[i] = self.parametrize_analytic(x[(i*3):(i+1)*3])
		if abs(numpy.prod(wgt)) < self.scale*1e-20:
			return 0.
		if multi_channel_index is not None:
			# basis transformation to cut basis st jacobian cancels integrable singularities ~1/E
			channel_cut_structure = self.cut_structures[multi_channel_index[0]]
			channel_cut = channel_cut_structure.cuts[multi_channel_index[1]]
			channel_basis_trsf = channel_cut_structure.inv_basis_trsf
			channel_basis_shift = numpy.array([q[1:] for q in channel_cut.basis_shift])
			loop_momenta = channel_basis_trsf.dot(loop_momenta - channel_basis_shift)

		kappas,jac = self.deform_on_prop(loop_momenta)
		#stop
		loop_momenta = [loop_momentum+1j*kappa for loop_momentum,kappa in zip(loop_momenta,kappas)]
		deltas = self.get_deltas(loop_momenta)
	
		if multi_channel_index is not None:
			multi_channel_factors = [cut_structure.get_multi_channel_factors(deltas) for cut_structure in self.cut_structures]
			channel_factor = multi_channel_factors[multi_channel_index[0]][multi_channel_index[1]]/numpy.sum(numpy.sum(multi_channel_factors))
			if abs(channel_factor) < self.scale*1e-20:
				return 0.
		else:
			channel_factor = 1.

		full_residue = self.evaluate(deltas)
		full_residue *= (-2*numpy.pi*1j)**self.n_loops
		full_residue *= (1./(2.*numpy.pi)**4)**self.n_loops

		return full_residue*numpy.prod(wgt)*jac*channel_factor

	def evaluate(self,deltas):
		res = [cut_structure.evaluate(deltas) for cut_structure in self.cut_structures]
		return numpy.sum(res)

	def vegas_integrate(self,integrand, N_train=1000, itn_train = 7, N_refine=1000, itn_refine = 10):
		print('{:=^90}'.format(' integrating ... '))
		print('{:=^90}'.format(' N_train = %d, itn_train = %d, N_refine = %d, itn_refine = %d ' %(N_train,itn_train,N_refine,itn_refine)))
		N_tot = N_train*itn_train+N_refine*itn_refine
		t0 = time.time()
		numpy.random.seed(0)
		batch_size = 1000
		real_integr = lambda x: integrand(x).real
		imag_integr = lambda x: integrand(x).imag
		real_vegas3_integrator = vegas.Integrator(3 * self.n_loops * [[0., 1.]])#, nhcube_batch=batch_size)
		imag_vegas3_integrator = vegas.Integrator(3 * self.n_loops * [[0., 1.]])#, nhcube_batch=batch_size)
		# train
		real_vegas3_integrator(real_integr,nitn=itn_train, neval=N_train)
		imag_vegas3_integrator(imag_integr,nitn=itn_train, neval=N_train)
		# refine
		real_result = real_vegas3_integrator(real_integr,nitn=itn_refine, neval=N_refine)
		imag_result = imag_vegas3_integrator(imag_integr,nitn=itn_refine, neval=N_refine)
		print('\n'+real_result.summary()+'\n'+imag_result.summary())
		real_vegas3_integrator.map.show_grid()
		imag_vegas3_integrator.map.show_grid()
		result = [real_result,imag_result]
		print('{:=^90}'.format(' I = ' + '{0} + {1}j '.format(result[0],result[1])))
		print('{:=^90}'.format(' Time: %.2f s ' % (time.time()-t0)))
		return result

	def cuba_integrate(self,integrand,integrator_name='vegas',neval=10000,**opts):
		self.n_bad_points = 0
		def cuba_integrand(ndim,xx,ncomp,ff,userdata):
			random_variables = [xx[i] for i in range(ndim.contents.value)]
			result = integrand(random_variables)
			if numpy.isfinite(result):
				ff[0] = result.real
				ff[1] = result.imag
			else:
				self.n_bad_points += 1
				print('bad point = ', self.n_bad_points)
				ff[0] = 0.
				ff[1] = 0.
			return 0
		
		opts = {'ndim': self.n_loops*3,
				'ncomp': 2,
				'maxeval': neval,
				'mineval': neval,
				'seed': 1,
				'verbose': 1,
				#'nincrease': 100, #vegas only
				#'nstart': neval,
				}
		

  		#NNEW = hyperparameters['Integrator']['n_new']
  		#NMIN = hyperparameters['Integrator']['n_start']
  		#FLATNESS = hyperparameters['Integrator']['flatness']
		#MINEVAL = hyperparameters['Integrator']['n_min']
		#MAXEVAL = hyperparameters['Integrator']['n_max']

		if integrator_name=='vegas':
			return pycuba.Vegas(cuba_integrand, **opts)
		elif integrator_name=='cuhre':
			return pycuba.Cuhre(cuba_integrand, **opts)

	def cuba_mc_integrate(self,integrator_name='vegas',neval=10000,**opts):
		self.n_bad_points = 0
		def cuba_integrand(ndim,xx,ncomp,ff,userdata):
			random_variables = [xx[i] for i in range(ndim.contents.value)]
			single_channel_res = self.ltd_integrand(random_variables)
			multi_channel_res = [self.ltd_integrand(random_variables,(cut_structure_index,cut_index))
											for cut_structure_index, cut_structure in enumerate(self.cut_structures)
											for cut_index,  cut in enumerate(cut_structure.cuts)]
			if numpy.isfinite(single_channel_res):
				ff[0] = single_channel_res.real
				ff[1] = single_channel_res.imag
			else:
				self.n_bad_points += 1
				print('bad point = ', self.n_bad_points)
				ff[0] = 0.
				ff[1] = 0.

			for channel_nr, channel_result in enumerate(multi_channel_res):
				if numpy.isfinite(channel_result):
					ff[2*(channel_nr+1)] = channel_result.real
					ff[2*(channel_nr+1)+1] = channel_result.imag
				else:
					self.n_bad_points += 1
					print('bad point = ', self.n_bad_points)
					ff[2*(channel_nr+1)] = 0.
					ff[2*(channel_nr+1)+1] = 0.
			return 0
		
		opts = {'ndim': int(self.n_loops*3),
				'ncomp': int((self.n_channels+1)*2),
				'maxeval': neval,
				'mineval': neval,
				'seed': 1,
				'verbose': 1,
				#'nincrease': 100, #vegas only
				#'nstart': neval, #vegas only
				}
		

  		#NNEW = hyperparameters['Integrator']['n_new']
  		#NMIN = hyperparameters['Integrator']['n_start']
  		#FLATNESS = hyperparameters['Integrator']['flatness']
		#MINEVAL = hyperparameters['Integrator']['n_min']
		#MAXEVAL = hyperparameters['Integrator']['n_max']

		if integrator_name=='vegas':
			return pycuba.Vegas(cuba_integrand, **opts)
		elif integrator_name=='cuhre':
			return pycuba.Cuhre(cuba_integrand, **opts)

	def cuba_mc_random_integrate(self,integrator_name='vegas',neval=10000,**opts):
		def cuba_integrand(ndim,xx,ncomp,ff,userdata):
			random_variables = [xx[i] for i in range(ndim.contents.value)]
			channel_index = [(cut_structure_index,cut_index)
											for cut_structure_index, cut_structure in enumerate(self.cut_structures)
											for cut_index,  cut in enumerate(cut_structure.cuts)]
			single_channel_res = self.ltd_integrand(random_variables)
			i = random.randint(0,self.n_channels-1)
			multi_channel_res = self.n_channels*self.ltd_integrand(random_variables,channel_index[i])
			ff[0] = single_channel_res.real
			ff[1] = single_channel_res.imag
			ff[2] = multi_channel_res.real
			ff[3] = multi_channel_res.imag
			return 0
		
		opts = {'ndim': int(self.n_loops*3),
				'ncomp': 4,
				'maxeval': neval,
				'mineval': neval,
				'seed': 1,
				'verbose': 1,
				#'nincrease': 100, #vegas only
				#'nstart': neval, #vegas only
				}
		

  		#NNEW = hyperparameters['Integrator']['n_new']
  		#NMIN = hyperparameters['Integrator']['n_start']
  		#FLATNESS = hyperparameters['Integrator']['flatness']
		#MINEVAL = hyperparameters['Integrator']['n_min']
		#MAXEVAL = hyperparameters['Integrator']['n_max']

		if integrator_name=='vegas':
			return pycuba.Vegas(cuba_integrand, **opts)
		elif integrator_name=='cuhre':
			return pycuba.Cuhre(cuba_integrand, **opts)


	def print_mc_results(self,result):
		for channel_nr in range(self.n_channels):
			channel_re = result['results'][2*channel_nr+2]
			channel_im = result['results'][2*channel_nr+3]
			if channel_nr != self.n_channels:
				print('Multi-Channel %i: 	Re %e +- %e 	Im %e +- %e'%(channel_nr,channel_re['integral'],channel_re['error'],
																					channel_im['integral'],channel_im['error']))
		channel_sum_int_re = numpy.sum([result['results'][2*channel_nr+2]['integral'] for channel_nr in range(self.n_channels)])
		channel_sum_int_im = numpy.sum([result['results'][2*channel_nr+3]['integral'] for channel_nr in range(self.n_channels)])
		channel_sum_err_re = numpy.sqrt(numpy.sum([result['results'][2*channel_nr+2]['error']**2 for channel_nr in range(self.n_channels)]))
		channel_sum_err_im = numpy.sqrt(numpy.sum([result['results'][2*channel_nr+3]['error']**2 for channel_nr in range(self.n_channels)]))
		print('Multi-Channel sum: 	Re %e +- %e 	Im %e +- %e'%(channel_sum_int_re,channel_sum_err_re,
																					channel_sum_int_im,channel_sum_err_im))
		single_channel_re = result['results'][0]
		single_channel_im = result['results'][1]
		print('Single-Channel: 	Re %e +- %e 	Im %e +- %e'%(single_channel_re['integral'],single_channel_re['error'],
																					single_channel_im['integral'],single_channel_im['error']))
		return

	def nested_cuhre_integrate(self,integrand,**opts):
		def inner_integrand(ndim,xx,ncomp,ff,userdata):
			random_variables = [xx[i] for i in range(ndim.contents.value)]
			all_variables = random_variables+self.userdata
			result = integrand(all_variables)
			if numpy.isfinite(result):
				ff[0] = result.real
				ff[1] = result.imag
			else:
				ff[0] = 0.
				ff[1] = 0.
			return 0

		def outer_integrand(ndim,xx,ncomp,ff,userdata):
			random_variables = [xx[i] for i in range(ndim.contents.value)]
			self.userdata = random_variables
			res_dict = pycuba.Cuhre(inner_integrand, ndim=3, ncomp= 2, maxeval = 1000, seed = 1, verbose=0,key=11)
			if numpy.isfinite(res_dict['results'][0]['integral']) and numpy.isfinite(res_dict['results'][1]['integral']):
				ff[0] = res_dict['results'][0]['integral']
				ff[1] = res_dict['results'][1]['integral']
			else:
				ff[0] = 0.
				ff[1] = 0.
			return 0
		
		print(pycuba.Cuhre(outer_integrand, ndim=3, ncomp= 2, maxeval = 1000, seed = 0, verbose=2,key=11))

	def nested_integrate(self,integrand,**opts):
		from scipy import integrate
		real_integrand = lambda x1,x2,x3: integrand([x1,x2,x3]).real
		print(integrate.nquad(real_integrand,[[0,1],[0,1],[0,1]],opts={'epsrel': 1.49e-01,'epsabs': 1.49e-01},full_output=True))

	def analytic_three_point_ladder(self,s1,s2,s3,l):
		aa = (1j/(16*numpy.pi**2*s3))**l
		lambd = lambda x,y: numpy.sqrt((1.-x-y)**2-4*x*y+0j)
		rho = lambda x,y: 2./(1.-x-y+lambd(x,y))
		def phi(x,y,l):
			bb = -1./(math.factorial(l)*lambd(x,y))
			summand = 0.
			for j in range(l,2*l+1):
				cc = (-1.)**j*math.factorial(j)*numpy.log(y/x)**(2*l-j)
				dd = math.factorial(j-l)*math.factorial(2*l-j)
				ee = mpmath.polylog(j,-1./(x*rho(x,y)))-mpmath.polylog(j,-y*rho(x,y))
				summand += cc*ee/dd
			return bb*summand
		x = s1/s3
		y = s2/s3
		return aa*phi(x,y,l)

class CutObject(object):

	def __init__(self,**opts):
		if 'cut_structure' in opts:
			self.cut_structure = opts.pop('cut_structure')
		if 'ltd_loop_lines' in opts:
			self.ltd_loop_lines = opts.pop('ltd_loop_lines')
		if 'cut_signature' in opts:
			self.cut_signature = opts.pop('cut_signature')
		if 'basis_trsf' in opts:
			self.basis_trsf = opts.pop('basis_trsf')
		if 'basis_shift' in opts:
			self.basis_shift = opts.pop('basis_shift')
		if 'inv_basis_trsf' in opts:
			self.inv_basis_trsf = opts.pop('inv_basis_trsf') 
		if 'cut_indices' in opts:
			self.cut_indices = opts.pop('cut_indices')
		if 'loop_line_index' in opts:
			self.loop_line_index = opts.pop('loop_line_index')
		if 'propagator_index' in opts:
			self.propagator_index = opts.pop('propagator_index')
		if 'basis_signature' in opts:
			self.basis_signature = opts.pop('basis_signature')
		if 'delta_signature' in opts:
			self.delta_signature = opts.pop('delta_signature')
		if 'has_positive_signature' in opts:
			self.has_positive_signature = opts.pop('has_positive_signature')
		if 'has_negative_signature' in opts:
			self.has_negative_signature = opts.pop('has_negative_signature')
		if opts != {}:
			print('Warning: Not all options used!')

class CutStructure(CutObject):

	def __init__(self,**opts):
		super(CutStructure, self).__init__(**opts)
		self.cut_signature = self.get_cut_signature()
		self.basis_trsf = self.get_basis_trsf() #cut line momenta from loop momenta
		self.inv_basis_trsf = numpy.linalg.inv(self.basis_trsf) #loop momenta from cut line momenta
		opts.update({'cut_signature': self.cut_signature, 'basis_trsf': self.basis_trsf, 'inv_basis_trsf': self.inv_basis_trsf})
		cut_propagator_indices = self.get_prop_cuts()
		self.cuts = [Cut(cut_indices=cut_indices,**opts) for cut_indices in cut_propagator_indices]

	def get_basis_trsf(self):
		trsf = []
		for loop_line,line in zip(self.ltd_loop_lines,self.cut_structure):
			if abs(line) == 1:
				trsf += [loop_line.signature]
		return numpy.array(trsf)

	def get_prop_cuts(self):
		prop_cuts = []
		cut_indices = [i for i,line in enumerate(self.cut_structure) if abs(line)==1]
		possible_cuts = [range(len(loop_line.propagators)) for loop_line,line in zip(self.ltd_loop_lines,self.cut_structure) if abs(line)==1]
		for prop_cut in itertools.product(*possible_cuts):
			prop_cuts += [prop_cut]
		return prop_cuts

	def get_cut_signature(self):
		signature = []
		for line in self.cut_structure:
			if abs(line) == 1:
				signature += [line]
		return signature

	def evaluate(self,deltas):
		cut_residues = [cut.evaluate(deltas) for cut in self.cuts]
		return numpy.sum(cut_residues)

	def get_multi_channel_factors(self,deltas):
		return [cut.get_multi_channel_factor(deltas) for cut in self.cuts]

class Cut(CutObject):

	def __init__(self,**opts):
		super(Cut, self).__init__(**opts)
		assert('cut_indices' in opts)
		self.basis_shift = self.get_basis_shift()
		opts.update({'basis_shift': self.basis_shift})
		self.loop_lines = [CutLoopLine(loop_line_index=i,**opts) for i,loop_line in enumerate(self.ltd_loop_lines)]

	def get_basis_shift(self):
		cut_indices_iter = iter(self.cut_indices)
		basis_shift = [loop_line.propagators[next(cut_indices_iter)].q for line, loop_line in zip(self.cut_structure,self.ltd_loop_lines) if abs(line) != 0]
		return basis_shift
	
	def get_multi_channel_factor(self,deltas):
		channel_factor = 1.
		for ll_deltas, loop_line in zip(deltas,self.loop_lines):
			for prop_delta, propagator  in zip(ll_deltas,loop_line.propagators):
				if not propagator.is_cut:
					channel_factor *= prop_delta
		return channel_factor

	def evaluate(self,deltas):
		ll_residues = [loop_line.evaluate(deltas) for loop_line in self.loop_lines]
		return numpy.prod(ll_residues)

class CutLoopLine(CutObject):

	def __init__(self,**opts):
		super(CutLoopLine, self).__init__(**opts)
		assert('loop_line_index' in opts)
		this_loop_line = self.ltd_loop_lines[self.loop_line_index]
		self.basis_signature = numpy.array(this_loop_line.signature).dot(self.inv_basis_trsf)
		self.delta_signature = self.basis_signature.dot(numpy.diag(self.cut_signature))
		self.has_positive_signature = all(sign >= 0 for sign in self.delta_signature)
		self.has_negative_signature = all(sign <= 0 for sign in self.delta_signature)
		opts.update({'delta_signature': self.delta_signature, 'basis_signature': self.basis_signature,
					'has_positive_signature': self.has_positive_signature, 'has_negative_signature': self.has_negative_signature})
		self.propagators = [CutPropagator(propagator_index = a, **opts) for a,propagator in enumerate(this_loop_line.propagators)]

	def evaluate(self,deltas):
		propagator_res = [propagator.evaluate(deltas) for propagator in self.propagators]
		return numpy.prod(propagator_res)

class CutPropagator(CutObject):

	def __init__(self,**opts):
		super(CutPropagator, self).__init__(**opts)
		assert('propagator_index' in opts)
		cut_indices_iter = iter(self.cut_indices)
		full_cut_indices = [next(cut_indices_iter) if abs(line) != 0 else None for line in self.cut_structure]
		self.is_cut = (full_cut_indices[self.loop_line_index] == self.propagator_index)
		#shift of propagator momentum with respect to the loop line momentum
		self.own_shift = self.ltd_loop_lines[self.loop_line_index].propagators[self.propagator_index].q
		# shift of the propagator momentum with respect to the cut propagator momentum basis
		self.cut_shift = -self.basis_signature.dot(self.basis_shift)
		#total shift i.e. affine surface term
		self.total_shift = self.own_shift + self.cut_shift
	
		s = self.total_shift[0]**2 - self.norm(self.total_shift[1:])**2
		self.is_plus_ellipsoid = (self.has_positive_signature and self.total_shift[0] < 0.) and (s > 0.)
		self.is_minus_ellipsoid = (self.has_negative_signature and self.total_shift[0] > 0.) and (s > 0.)
		self.is_ellipsoid = (self.is_plus_ellipsoid or self.is_minus_ellipsoid)
		#self.surfaces = self.get_surfaces()

	def get_surfaces(self):
		if not self.is_cut:
			cut_loop_line_indices = [i for i,line in enumerate(self.cut_structure) if line != 0]
			cut_deltas = [{	'loop_line_index'	: ll_i,
							'propagator_index'	: prop_i,
							'sign'				: int(sign)}
							for ll_i, prop_i, sign in zip(cut_loop_line_indices,self.cut_indices,self.delta_signature) if sign != 0]
			is_one_surface = True
			if self.is_plus_ellipsoid:
				sign = 1
			elif self.is_minus_ellipsoid:
				sign = -1
			elif self.has_positive_signature:
				sign = -1
			elif self.has_negative_signature:
				sign = 1
			else:
				is_one_surface = False
				sign = 1
			prop_delta = [{	'loop_line_index'	: self.loop_line_index,
							'propagator_index'	: self.propagator_index,
							'sign'				: sign}]
			surface = {'deltas': cut_deltas+prop_delta, 'p_shift': self.total_shift.copy(), 'is_ellipsoid': self.is_ellipsoid,'multiplicity': 1}
			if is_one_surface:
				return [surface]
			else:
				another_prop_delta = [{	'loop_line_index'	: self.loop_line_index,
										'propagator_index'	: self.propagator_index,
										'sign'				: -sign}]
				another_surface = {'deltas': cut_deltas+another_prop_delta, 'p_shift': self.total_shift.copy(), 'is_ellipsoid': self.is_ellipsoid,'multiplicity': 1}
				return [surface,another_surface]
		else:
			return []

	def deform(self,loop_momenta,deltas,**opts):
		scale = opts.pop('scale')
		a_ij = opts.pop('a_ij')
		lambda_ij = opts.pop('lambda_ij')
		kappas = [numpy.zeros(3) for loop_momentum in loop_momenta]
		if not self.is_cut:
			if self.is_ellipsoid:
				cut_indices_iter = iter(self.cut_indices)
				my_cuts =  [(i,next(cut_indices_iter)) for i, line in enumerate(self.cut_structure) if abs(line) != 0]
				my_cuts += [(self.loop_line_index,self.propagator_index)]
				basis_shift_space = [shift[1:] for shift in self.basis_shift]
				cut_prop_momenta = self.basis_trsf.dot(loop_momenta) + basis_shift_space
				cut_prop_kappas = numpy.array([[coord*0. for coord in loop_momentum] for loop_momentum in loop_momenta])
				cut_deltas = self.get_cut_deltas(deltas)
				delta = deltas[self.loop_line_index][self.propagator_index]
				v1 = numpy.array(self.ltd_loop_lines[self.loop_line_index].signature).dot(loop_momenta) + self.own_shift[1:]
				for i,sign in enumerate(self.delta_signature):
					if sign != 0:
						v2 = cut_prop_momenta[i]
						cut_prop_kappas[i] = self.basis_signature[i]*v1/delta
						cut_prop_kappas[i] += v2/cut_deltas[i]
				energy = self.delta_signature.dot(cut_deltas) + self.total_shift[0]
				interpolation = adipy.exp(-(energy-delta*numpy.sign(self.total_shift[0]))**2/scale**2/a_ij)
				#interpolation = adipy.exp(-(energy**2-delta**2)**2/scale**4/a_ij)
				"""
				interpol = 1.
				# this is just for 2-loop!!
				if abs(numpy.sum(self.delta_signature)) == len(cut_deltas):
					for i,sign in enumerate(self.delta_signature):
						new_signature = numpy.append(self.delta_signature[:i],self.delta_signature[i+1:])
						new_deltas = numpy.append(cut_deltas[:i],cut_deltas[i+1:])
						energy_i = new_signature.dot(new_deltas)
						#interpol *= (1.-adipy.exp(-(energy_i**2-delta**2)**2/scale**4/a_ij))
						interpol *= (1.-adipy.exp(-(energy_i-delta**numpy.sign(self.total_shift[0]))**2/scale**2/a_ij))
					#interpol *= (1. -adipy.exp(-(energy**2)**2/scale**4/a_ij))
					interpol *= (1. -adipy.exp(-(energy)**2/scale**2/a_ij))
				interpolation *= interpol
				"""
				kappas = self.inv_basis_trsf.dot(cut_prop_kappas)
				#lambda_factor = -0.01*scale
				lambda_factor = -abs(lambda_ij)
				kappas *= lambda_factor*scale*interpolation
				#if my_cuts == [(0, 0), (2, 0), (0, 1)]:
					#print my_cuts
					#print interpolation
					#print energy
					#print cut_prop_kappas
					#print kappas
					#print '='*100
				#print kappas
		return kappas

	def norm(self,q):
		"""for complex and real q arrays, not the same as numpy.linalg.norm for complex numbers!!"""
		return adipy.sqrt(numpy.sum([q_i*q_i for q_i in q]))

	def inv_G(self,x,y):
		""" the dual propagator """
		return x**2 - y**2

	def get_cut_deltas(self,deltas):
		cut_indices_iter = iter(self.cut_indices)
		cut_deltas = [delta[next(cut_indices_iter)] for line,delta in zip(self.cut_structure,deltas) if abs(line) != 0]
		return cut_deltas

	def evaluate(self,deltas):
		cut_deltas = self.get_cut_deltas(deltas)
		delta = deltas[self.loop_line_index][self.propagator_index]
		if not self.is_cut:
			energy = self.delta_signature.dot(cut_deltas) + self.total_shift[0]
			inv_prop = self.inv_G(energy,delta)
			if abs(inv_prop) < 1e-20:
				if ((all(sign >= 0. for sign in self.delta_signature) and energy < 0.)
					or (all(sign <= 0. for sign in self.delta_signature) and energy > 0.)):
					raise ValueError('Hit ellipsoid singularity')
				#print('a non cut propagator is zero')
				#print(cut_deltas,self.cut_structure,self.cut_indices,energy,delta,self.delta_signature)
				#raise ValueError('a non cut propagator is zero')
				return 0.
		else:
			inv_prop = 2*delta
			if inv_prop == 0.:
				print('a cut propagator is zero')
				print(cut_deltas,self.cut_structure,self.cut_indices,energy,delta,self.delta_signature)
				#raise ValueError('a cut propagator is zero')
		return 1./inv_prop

if __name__ == "__main__":	
	
	my_topology = topologies.hard_coded_topology_collection['Triangle_no_ellipse']
	hyperparameters = ltd_commons.hyperparameters
	my_LTD = LTDnLoop(my_topology,hyperparameters)
	my_LTD.display_surfaces(show_ellipsoids=True,show_hyperboloids=False)

	#print my_LTD.ltd_integrand4([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
	
	#external_momenta = my_LTD.external_kinematics
	#s = [p[0]**2-numpy.sum(p[1:]**2) for p in external_momenta]
	#print my_LTD.analytic_three_point_ladder(s[0],s[1],s[2],3)

	#result = my_LTD.cuba_integrate(my_LTD.ltd_integrand,integrator_name='cuhre',neval=800000)
	#result = my_LTD.cuba_integrate(lambda x: my_LTD.ltd_integrand(x,(0,0)),integrator_name='vegas',neval=80000)

	print(my_LTD.ltd_integrand([0.5,0.4,0.3],multi_channel_index=(0,0)))

	#result = my_LTD.cuba_mc_integrate(integrator_name='vegas',neval=100000)
	#my_LTD.print_mc_results(result)

	result = my_LTD.cuba_mc_random_integrate(integrator_name='vegas',neval=10000)

	#result = my_LTD.vegas_integrate(my_LTD.ltd_integrand, N_refine=1000)#, N_train = 500)


