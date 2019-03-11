# n loop LTD
import numpy as numpy
import vegas
import adipy as adipy
import time
import itertools
import mpmath
import math
import ltd_commons
import random
from numpy.linalg import inv
#import warnings
#warnings.simplefilter("error")
#warnings.simplefilter("ignore", DeprecationWarning)
import sys
sys.setrecursionlimit(100)

class LTDnLoop:

	def __init__(self,topology):
		self.loop_lines         = topology.loop_lines
		self.ltd_cut_structure  = topology.ltd_cut_structure
		self.n_loops            = topology.n_loops 
		self.name               = topology.name
		self.external_kinematics= topology.external_kinematics
		self.scale				= self.get_scale()
		self.old_possible_cuts	= self.old_get_possible_cuts()
		self.analytical_result 	= topology.analytical_result

		self.cut_structures		= [CutStructure(ltd_cut_structure,self.loop_lines) for ltd_cut_structure in self.ltd_cut_structure]
		for cut_structure in self.cut_structures:
			for cut_indices in cut_structure.cut_propagator_indices:
				cut_structure.ltd_cuts += [Cut(cut_indices,cut_structure)]

		for cut_structure in self.cut_structures:
			for cut in cut_structure.ltd_cuts:
				for i,loop_line in enumerate(cut.loop_lines):
					ll_propagators = []
					for a in xrange(len(loop_line.propagators)): 
						ll_propagators += [Propagator(a,i,cut)]
					cut.ltd_propagators += [ll_propagators]

	def delta(self,q_space,m):
		on_f_shell = adipy.sqrt(self.norm(q_space)**2 + m**2)
		return on_f_shell
	
	def get_scale(self):
		p_posE = [p for p in self.external_kinematics if p[0] > 0.]
		if -numpy.sum(self.external_kinematics,axis=0)[0] > 0.:
			p_posE += [-sum(self.external_kinematics)]
		if (len(p_posE) == 0):
			print "All incoming momentum energies are zero!! Improvising scale..."
			scale = self.norm(numpy.sum(self.external_kinematics[1:],axis=0))**2
		else:
			p_sum = sum(p_posE)
			s = p_sum[0]**2 - self.norm(p_sum[1:])**2
			scale = numpy.sqrt(abs(s))
		print 'scale = ', scale
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
	
	def old_get_parents(self,ll_nr):
		""" find parent looplines of a given loopline """
		which_node = bool(random.getrandbits(1))
		if which_node:
			node = self.loop_lines[ll_nr].start_node
		else:
			node = self.loop_lines[ll_nr].end_node
		parents = []
		for i,loop_line in enumerate(self.loop_lines):
			if i != ll_nr:
				if node == loop_line.start_node:
					parents += [[i,False]]
				elif node == loop_line.end_node:
					parents += [[i,True]]
		for i,parent in enumerate(parents):
			if not which_node:
				parents[i][1] = not parent[1]
		return parents

	def get_deltas(self,loop_momenta):
		""" computes the sqrt[(l_space+p_space)**2 + m**2] of all looplines
			looks like looplines_prop = [[q_i0p1],[q_i0p2],[q_i0p3]]"""
		deltas = []
		for loop_line in self.loop_lines:
			loop_momentum = numpy.sum([l*sign for l,sign in zip(loop_momenta,loop_line.signature)],axis=0)
			deltas += [[self.delta(loop_momentum + prop.q[1:],numpy.sqrt(prop.m_squared)) for prop in loop_line.propagators]]
		return deltas
	
	def old_inv_G(self,x,y):
		""" the dual propagator """
		return x**2 - y**2
		
	def old_get_four_momenta(self,line_energies,loop_momenta):
		four_momenta = [numpy.zeros(4) for i in xrange(self.n_loops)]
		for loop_line in self.loop_lines:
			aa = numpy.array(loop_line.signature)
			if numpy.sum(aa**2) == 1:
				index = numpy.where(numpy.absolute(aa)==1)[0][0]
				four_momenta[index][1:] = aa[index]*loop_momenta[index]
				four_momenta[index][0] = aa[index]*line_energies[index]
		return four_momenta
	
	def old_evaluate_ll(self,loop_line,c,line_energy,delta):
		""" ll: loopline, c: positive (cut index+1), line_energy: energy comp of line momentum, ll_prop: Deltas of this line """
		""" evaluates a loop line at a fixed cut momentum """
		ll_residue = 1.
		for i,prop in enumerate(loop_line.propagators):
			if i != (c-1):
				ll_residue *= self.old_inv_G(line_energy + prop.q[0],delta[i])
			else:
				ll_residue *= 2*delta[i]
		ll_residue = 1./ll_residue	
		return ll_residue
	
	def old_evaluate_cut(self,deltas,cut):
		""" cut looks like cut = [1,-2,0], 0 is no cut, the abs(c) is the (propagator index +1), - means negative cut
			computes the residue of cuts of specific propagators """
		line_energies = self.old_get_all_line_energies(deltas,cut)
		ll_residues = []
		for loop_line,c,line_energy,delta in zip(self.loop_lines,cut,line_energies,deltas):
			ll_residues += [self.old_evaluate_ll(loop_line,abs(c),line_energy,delta)]
		cut_residue = numpy.prod(ll_residues)
		return cut_residue

	def old_get_all_line_energies(self,deltas,cut):
		line_energies = [numpy.sign(c)*delta[abs(c)-1] - loop_line.propagators[abs(c)-1].q[0] if c != 0
						else None
						for loop_line,delta,c in zip(self.loop_lines,deltas,cut)]
		for ll_nr,line_energy in enumerate(line_energies):
			if line_energy == None:
				line_energies[ll_nr] = self.old_get_line_energy(ll_nr,line_energies)
		return line_energies
	
	def old_get_line_energy(self,ll_nr,line_energies):
		"""a recursive algorythm that computes the energy of the corresponding loop line,
			should work for loops larger than 2 as well"""
		new_line_energy = 0.
		parents = self.old_get_parents(ll_nr)
		for nr,sign in parents:
			if line_energies[nr] == None:
				line_energies[nr] = self.get_line_energy(nr,line_energies)
			new_line_energy += (-1.)**(sign+1)*line_energies[nr]
		return new_line_energy
				
	def old_evaluate_cut_structure(self,deltas,cut_structure):
		""" computes all residues of a cut structure """
		cut_indices = [i for i,line in enumerate(cut_structure) if abs(line)==1]
		cut_residues = []
		for ll_cuts in itertools.product(*[self.old_possible_cuts[i] for i in cut_indices]):
			cut = [0 for i in xrange(len(self.loop_lines))]
			for i,ll_cut in zip(cut_indices,ll_cuts):
				cut[i] = ll_cut*cut_structure[i]
			cut_residues += [self.old_evaluate_cut(deltas,cut)]
		cut_structure_residue = numpy.sum(cut_residues)
		return cut_structure_residue
				
	def old_get_possible_cuts(self):
		""" returns a list with (indices + 1) of all possible cuts per loop line """
		possible_cuts = [range(1,len(loop_line.propagators)+1) for loop_line in self.loop_lines]
		return possible_cuts
	
	def old_get_dual_vec(self,k_space,l_space):
		vec = numpy.append(k_space,l_space)
		dim = 3*self.n_loops
		vec_dual = [0.]*dim
		for i in xrange(dim):
			vec_dual[i] = adipy.ad(vec[i], numpy.array([0. if j !=i  else 1. for j in xrange(dim)]))
		return vec_dual
	
	def dual(self,loop_momenta):
		vec = numpy.array(loop_momenta).flatten()
		dim = len(vec)
		dual_vec = [0.]*dim
		for i in xrange(dim):
			dual_vec[i] = adipy.ad(vec[i], numpy.array([0. if j !=i  else 1. for j in xrange(dim)]))
		dual_momenta = [dual_vec[i*3:(i+1)*3] for i in xrange(self.n_loops)]
		return numpy.array(dual_momenta)

	def old_hardcoded_dt_def(self,k_space,l_space,p):
		# according to dario convention of loop line and momentum flow assignments
		vec_dual = self.old_get_dual_vec(k_space,l_space)
		k_space = vec_dual[:3]
		l_space = vec_dual[3:]
		p_space = p[1:]
		p0 = p[0]
		
		k_minus_l = numpy.array([k_sp-l_sp for k_sp,l_sp in zip(k_space,l_space)])
		
		sigma_squared = 0.1*abs(p0**2-numpy.sum(p_space**2))
		k_p = self.norm(k_space+p_space)
		l_p = self.norm(l_space+p_space)
		k_l = self.norm(k_minus_l)
		l = self.norm(l_space)
		k = self.norm(k_space)
		
		inv_G = lambda x,y: x**2-y**2
		exp_tri_k = k_p+k_l+l-p0
		exp_tri_l = l_p+k_l-p0+k
		exp_dua_k = k_p-p0+k
		exp_dua_l = l_p-p0+l
		
		f_tri_k = adipy.exp(-(exp_tri_k**2)/sigma_squared)
		f_tri_l = adipy.exp(-(exp_tri_l**2)/sigma_squared)
		f_tri_k2 = f_tri_l
		f_tri_l2 = f_tri_k
 				
		f_dua_k = adipy.exp(-(exp_dua_k**2)/sigma_squared)
		f_dua_l = adipy.exp(-(exp_dua_l**2)/sigma_squared)
		
		dir_kp = numpy.array([(k_sp+p_sp)/k_p for k_sp,p_sp in zip(k_space,p_space)])
		dir_lp = numpy.array([(l_sp+p_sp)/l_p for l_sp,p_sp in zip(l_space,p_space)])
		dir_kl = numpy.array([k_minus_l_sp/k_l for k_minus_l_sp in k_minus_l])
		dir_k = numpy.array([k_sp/k for k_sp in k_space])
		dir_l = numpy.array([l_sp/l for l_sp in l_space])
		
		dir_tri_k = dir_kp + dir_kl
		dir_tri_l = dir_lp - dir_kl
		dir_tri_k2 = dir_k + dir_kl
		dir_tri_l2 = dir_l - dir_kl
		
		dir_dua_k = dir_k + dir_kp
		dir_dua_l = dir_l + dir_lp

		kapp_tri_k = numpy.array([dir*f_tri_k for dir in dir_tri_k])
		kapp_tri_l = numpy.array([dir*f_tri_l for dir in dir_tri_l])
		kapp_tri_k2 = numpy.array([dir*f_tri_k2 for dir in dir_tri_k2])
		kapp_tri_l2 = numpy.array([dir*f_tri_l2 for dir in dir_tri_l2])
		kapp_dua_k = numpy.array([dir*f_dua_k for dir in dir_dua_k])
		kapp_dua_l = numpy.array([dir*f_dua_l for dir in dir_dua_l])
				
		lambda_k = lambda_l = .1
		kappa_k = - lambda_k * ( kapp_tri_k + kapp_tri_k2 + kapp_dua_k)
		kappa_l = - lambda_l * ( kapp_tri_l + kapp_tri_l2 + kapp_dua_l)
		
		"""
		scaling_param_k = self.get_scaling_param(k_space,kappa_k)
		scaling_param_l = self.get_scaling_param(l_space,kappa_l)
		for i in xrange(self.dim):
			kappa_k[i] *= scaling_param_k
			kappa_l[i] *= scaling_param_l
			if scaling_param_k < self.curr_min_scaling:
				self.curr_min_scaling = scaling_param_k
				print 'Current min. scaling', self.curr_min_scaling
			elif scaling_param_l < self.curr_min_scaling:
				self.curr_min_scaling = scaling_param_l
				print 'Current min. scaling', self.curr_min_scaling
		"""
		
		full_jac = adipy.jacobian(vec_dual + 1j*numpy.append(kappa_k,kappa_l))
		det_jac = numpy.linalg.det(full_jac)
		kappa_k = numpy.array([kapp.nom for kapp in kappa_k])
		kappa_l = numpy.array([kapp.nom for kapp in kappa_l])	
		
		return kappa_k,kappa_l,det_jac
	
	def old_deform(self,loop_momenta):
		""" doesn't deform yet """
		kappas = [numpy.zeros(3) for i in xrange(self.n_loops)]
		jac = 1.
		
		#double triangle (with dario convention of loop line assignments!!)
		#kappa_k,kappa_l,jac = self.old_hardcoded_dt_def(loop_momenta[0],loop_momenta[1],self.loop_lines[0].propagators[1].q)
		#kappas = [kappa_k,kappa_l]
		
		return kappas, jac
	
	def deform(self,loop_momenta):
		loop_momenta = self.dual(loop_momenta)
		real_deltas = self.get_deltas(loop_momenta)
		kappas = numpy.array([[coord*0. for coord in loop_momentum] for loop_momentum in loop_momenta])
		for cut_structure in self.cut_structures:
			for cut in cut_structure.ltd_cuts:
				for ll_nr,loop_line in enumerate(self.loop_lines):
					for propagator in cut.ltd_propagators[ll_nr]:
						kappas += propagator.deform(loop_momenta,real_deltas,self.scale)
		full_jac = adipy.jacobian(loop_momenta.flatten() + 1j*kappas.flatten())
		det_jac = numpy.linalg.det(full_jac)
		kappas = [numpy.array([kapp.nom for kapp in kappa]) for kappa in kappas]
		return kappas, det_jac

	def old_ltd_integrand(self,x):
		loop_momenta = [numpy.zeros(3) for i in xrange(self.n_loops)]
		wgt = numpy.zeros(self.n_loops)
		for i in xrange(self.n_loops):
			loop_momenta[i], wgt[i] = self.parametrize_analytic(x[(i*3):(i+1)*3])
		kappas,jac = self.old_deform(loop_momenta)
		loop_momenta = [loop_momentum+1j*kappa for loop_momentum,kappa in zip(loop_momenta,kappas)]
		deltas = self.get_deltas(loop_momenta)
		full_residue = numpy.sum([self.old_evaluate_cut_structure(deltas,cut_structure) for cut_structure in self.ltd_cut_structure])
		ltd_factor = (-2*numpy.pi*1j)**self.n_loops
		standard_factor = (1./(2.*numpy.pi)**4)**self.n_loops
		return full_residue*numpy.prod(wgt)*ltd_factor*standard_factor*jac
	
	def ltd_integrand(self,x):
		loop_momenta = [numpy.zeros(3) for i in xrange(self.n_loops)]
		wgt = numpy.zeros(self.n_loops)
		for i in xrange(self.n_loops):
			loop_momenta[i], wgt[i] = self.parametrize_analytic(x[(i*3):(i+1)*3])
		kappas,jac = self.deform(loop_momenta)
		loop_momenta = [loop_momentum+1j*kappa for loop_momentum,kappa in zip(loop_momenta,kappas)]
		deltas = self.get_deltas(loop_momenta)
		full_residue = self.evaluate(deltas)
		ltd_factor = (-2*numpy.pi*1j)**self.n_loops
		standard_factor = (1./(2.*numpy.pi)**4)**self.n_loops
		return full_residue*numpy.prod(wgt)*ltd_factor*standard_factor*jac

	def evaluate(self,deltas):
		res = 0.
		for cut_structure in self.cut_structures:
			res += cut_structure.evaluate(deltas)
		return res

	def ltd_integrate(self,integrand, N_train=1000, itn_train = 7, N_refine=1000, itn_refine = 10):
		print '{:=^90}'.format(' '+ self.name + ' ')
		print '{:=^90}'.format(' analytical result = {num.real} + {num.imag}j'.format(num =self.analytical_result))
		print '{:=^90}'.format(' integrating ... ')
		print '{:=^90}'.format(' N_train = %d, itn_train = %d, N_refine = %d, itn_refine = %d ' %(N_train,itn_train,N_refine,itn_refine))
		N_tot = N_train*itn_train+N_refine*itn_refine
		t0 = time.time()
		numpy.random.seed(0)
		batch_size = 600
		real_integr = lambda x: integrand(x).real
		imag_integr = lambda x: integrand(x).imag
		real_vegas3_integrator = vegas.Integrator(3 * self.n_loops * [[0., 1.]], nhcube_batch=batch_size)
		imag_vegas3_integrator = vegas.Integrator(3 * self.n_loops * [[0., 1.]], nhcube_batch=batch_size)
		# train
		real_vegas3_integrator(real_integr,nitn=itn_train, neval=N_train)
		imag_vegas3_integrator(imag_integr,nitn=itn_train, neval=N_train)
		# refine
		real_result = real_vegas3_integrator(real_integr,nitn=itn_refine, neval=N_refine)
		imag_result = imag_vegas3_integrator(imag_integr,nitn=itn_refine, neval=N_refine)
		print('\n'+real_result.summary()+'\n'+imag_result.summary())
		#print('\n'+real_result.summary())
		real_vegas3_integrator.map.show_grid()
		imag_vegas3_integrator.map.show_grid()
		result = [real_result,imag_result]
		print '{:=^90}'.format(' I = ' + '{0} + {1}j '.format(result[0],result[1]))
		print '{:=^90}'.format(' Time: %.2f s ' % (time.time()-t0))
		return result

	def analytic_three_point_ladder(self,s1,s2,s3,l):
		aa = (1j/(16*numpy.pi**2*s3))**l
		lambd = lambda x,y: numpy.sqrt((1.-x-y)**2-4*x*y+0j)
		rho = lambda x,y: 2./(1.-x-y+lambd(x,y))
		def phi(x,y,l):
			bb = -1./(math.factorial(l)*lambd(x,y))
			summand = 0.
			for j in xrange(l,2*l+1):
				cc = (-1.)**j*math.factorial(j)*numpy.log(y/x)**(2*l-j)
				dd = math.factorial(j-l)*math.factorial(2*l-j)
				ee = mpmath.polylog(j,-1./(x*rho(x,y)))-mpmath.polylog(j,-y*rho(x,y))
				summand += cc*ee/dd
			return bb*summand
		x = s1/s3
		y = s2/s3
		return aa*phi(x,y,l)

class CutStructure(LTDnLoop):

	def __init__(self,cut_structure,loop_lines):
		self.loop_lines = loop_lines
		self.cut_structure = cut_structure
		self.n_loop_lines = len(self.loop_lines)
		self.cut_signature = self.get_cut_signature()
		self.basis_trsf = self.get_basis_trsf() #cut line momenta from loop momenta
		self.inv_basis_trsf = numpy.linalg.inv(self.basis_trsf) #loop momenta from cut line momenta
		self.cut_propagator_indices = self.get_prop_cuts()
		self.ltd_cuts = []

	def get_energy_signature(self,ll_nr):
		loop_line = self.loop_lines[ll_nr]
		sign = numpy.array(loop_line.signature)
		cut_sign = numpy.diag(self.cut_signature)
		energy_signature = sign.dot(self.inv_basis_trsf).dot(cut_sign)
		return energy_signature

	def get_basis_trsf(self):
		trsf = []
		for i,line in enumerate(self.cut_structure):
			if abs(line) == 1:
				trsf += [self.loop_lines[i].signature]
		return numpy.array(trsf)

	def get_prop_cuts(self):
		prop_cuts = []
		cut_indices = [i for i,line in enumerate(self.cut_structure) if abs(line)==1]
		possible_cuts = [range(len(loop_line.propagators)) for loop_line,line in zip(self.loop_lines,self.cut_structure) if abs(line)==1]
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
		cut_residues = []
		for cut in self.ltd_cuts:
			cut_residues += [cut.evaluate_cut(deltas)]
		return numpy.sum(cut_residues)

class Cut(CutStructure):

	def __init__(self, cut, cut_structure):
		assert(isinstance(cut_structure,CutStructure))
		self.cut_structure = cut_structure
		self.loop_lines = cut_structure.loop_lines
		self.cut = cut
		self.basis_shift = self.get_basis_shift() #new basis CutPropagatorBasis = CutLineBasis + basis_shift
		self.ltd_propagators = []

	def get_basis_shift(self):
		basis_shift = []
		i = 0
		for loop_line,line in zip(self.loop_lines,self.cut_structure.cut_structure):
			if line != 0:
				basis_shift += [loop_line.propagators[self.cut[i]].q]
				i += 1
		return basis_shift

	def evaluate_ll(self,ll_nr,cut_deltas,ll_deltas):
		ll_propagators = self.ltd_propagators[ll_nr]
		ll_residue = 1.
		for a,prop in enumerate(self.loop_lines[ll_nr].propagators):
			ll_residue *= ll_propagators[a].evaluate(cut_deltas,ll_deltas[a])
		return ll_residue

	def get_cut_deltas(self,deltas):
		cut_deltas = []
		i = 0
		for line,loop_line,delta in zip(self.cut_structure.cut_structure,self.loop_lines,deltas):
			if line != 0:
				cut_deltas += [delta[self.cut[i]]]
				i += 1
		return cut_deltas

	def evaluate_cut(self,deltas):
		cut_deltas = self.get_cut_deltas(deltas)
		ll_residues = []
		for ll_nr,loop_lines in enumerate(self.loop_lines):
			ll_residues += [self.evaluate_ll(ll_nr,cut_deltas,deltas[ll_nr])]
		return numpy.prod(ll_residues)

class Propagator(Cut):

	def __init__(self, a, i, cut):
		""" propagator a of loop line i"""
		assert(isinstance(cut,Cut))
		self.cut = cut
		self.cut_structure = cut.cut_structure
		self.loop_lines = cut.loop_lines
		self.basis_trsf = self.cut_structure.basis_trsf
		self.inv_basis_trsf = self.cut_structure.inv_basis_trsf
		self.basis_shift = cut.basis_shift
		self.a = a
		self.i = i
		self.loop_line = self.loop_lines[i]

		# find out if this particular propagator is cut
		cut_indices = [i for i,line in enumerate(self.cut_structure.cut_structure) if abs(line)==1]
		ll_cut = [None for i in xrange(len(self.loop_lines))]
		for i,c in zip(cut_indices,self.cut.cut):
			ll_cut[i] = c
		self.is_cut = False
		for i,c in enumerate(ll_cut):
			if i == self.i and c == self.a:
				self.is_cut = True

		self.own_shift = self.loop_line.propagators[a].q #shift of propagator momentum with respect to the loop line momentum
		self.cut_shift = self.get_cut_shift() # shift of the propagator momentum with respect to the cut propagator momentum basis
		self.total_shift = self.own_shift + self.cut_shift #total shift i.e. affine surface term

	def deform(self,loop_momenta,deltas,scale):
		kappas = [numpy.zeros(3) for loop_momentum in loop_momenta]
		if not self.is_cut:
			energy_signature = self.cut_structure.get_energy_signature(self.i)
			plus_ellipsoid = (all(sign >= 0 for sign in energy_signature) and self.total_shift[0] < 0.)
			minus_ellipsoid = (all(sign <= 0 for sign in energy_signature) and self.total_shift[0] > 0.)
			if (minus_ellipsoid or plus_ellipsoid):
				s = self.total_shift[0]**2-self.total_shift[1]**2-self.total_shift[2]**2 -self.total_shift[3]**2
				#assert(s > 0.) #this is not a complete condition, sometimes ellipsoids don't exist for certain loop_momenta
				basis_shift_space = [shift[1:] for shift in self.basis_shift]
				cut_prop_momenta = self.basis_trsf.dot(loop_momenta) + basis_shift_space
				cut_prop_signature = numpy.array(self.loop_line.signature).dot(self.inv_basis_trsf)
				cut_prop_kappas = numpy.array([[coord*0. for coord in loop_momentum] for loop_momentum in loop_momenta])
				for i,sign in enumerate(energy_signature):
					if sign != 0:
						v1 = cut_prop_momenta[i]
						v2 = numpy.array(self.loop_line.signature).dot(loop_momenta) + self.own_shift[1:]
						cut_prop_kappas[i] = v1/self.norm(v1)
						cut_prop_kappas[i] += cut_prop_signature[i]*v2/self.norm(v2)
				cut_deltas = self.cut.get_cut_deltas(deltas)
				delta = deltas[self.i][self.a]
				energy = energy_signature.dot(cut_deltas) + self.total_shift[0]
				A_ij = 10.
				if minus_ellipsoid:
					interpolation = adipy.exp(-(energy-delta)**2/scale**2/A_ij) 
				else:
					interpolation = adipy.exp(-(energy+delta)**2/scale**2/A_ij)
				kappas = self.inv_basis_trsf.dot(cut_prop_kappas)
				lambda_factor = -0.01*scale
				#lambda_factor = -10.
				kappas *= lambda_factor*interpolation
		return kappas

	def norm(self,q):
		"""for complex and real q arrays
		not the same as numpy.linalg.norm for complex numbers!!"""
		return adipy.sqrt(numpy.sum([q_i*q_i for q_i in q]))

	def get_cut_shift(self):
		cut_shift = -numpy.array(self.loop_line.signature).dot(self.inv_basis_trsf.dot(self.basis_shift))
		return cut_shift

	def inv_G(self,x,y):
		""" the dual propagator """
		return x**2 - y**2

	def evaluate(self,cut_deltas,delta):
		if not self.is_cut:
			sign = numpy.array(self.loop_line.signature)
			cut_sign = numpy.diag(self.cut_structure.cut_signature)
			energy = sign.dot(self.inv_basis_trsf).dot(cut_sign).dot(cut_deltas) + self.total_shift[0]
			inv_prop = self.inv_G(energy,delta)
		else:
			inv_prop = 2*delta
		return 1./inv_prop

if __name__ == "__main__":	
	
	def example():
		# example to evaluate specific residues at given loop momenta configuration
		my_topology = ltd_commons.hard_coded_topology_collection['DoubleTriangle']
		my_LTD = LTDnLoop(my_topology)
		loop_momenta = [numpy.array([1.1,2.2,3.3]),numpy.array([4.4,5.5,6.6])]
		deltas = my_LTD.get_deltas(loop_momenta)
		cut = [1,-2,0] # propagator numbering starts with 1, 0 means no cut, - negative cut.
		cut_residue = my_LTD.evaluate_cut(deltas,cut)
		cut_structure = [-1,0,1] # the usual convention
		cut_structure_residue = my_LTD.evaluate_cut_structure(deltas,cut_structure)
		full_residue = numpy.sum([my_LTD.evaluate_cut_structure(deltas,cut_structure) for cut_structure in my_LTD.ltd_cut_structure])

		# example to evaluate propagator i of loop line j for a given cut:
		i = 0
		j = 1
		line_energies = my_LTD.get_all_line_energies(deltas,cut)
		inv_propagator = my_LTD.inv_G(line_energies[j] + my_LTD.loop_lines[j].propagators[i].q[0],deltas[j][i])
		# or alternatively build it yourself from the four momenta
		four_momenta = my_LTD.get_four_momenta(line_energies,loop_momenta)
		four_loop_momentum_j = numpy.sum([l*sign for l,sign in zip(four_momenta,my_LTD.loop_lines[j].signature)],axis=0) #loop momentum in line j
		four_ext_momentum_i = my_LTD.loop_lines[j].propagators[i].q #(sum of) external momenta that flow into propagator i
		m_squared_i = my_LTD.loop_lines[j].propagators[i].m_squared # mass of propagator i
		feyn_prop = lambda k,m_squared: k[0]**2-my_LTD.norm(k[1:])**2 - m_squared
		alt_inv_propagator = feyn_prop(four_loop_momentum_j+four_ext_momentum_i,m_squared_i)
		print('Convince yourself:', inv_propagator, ' = ', alt_inv_propagator)
		return
	
	#example()
	
	
	random.seed(0)
	my_topology = ltd_commons.hard_coded_topology_collection['Triangle']
	my_LTD = LTDnLoop(my_topology)
		
	#print my_LTD.ltd_integrand4([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
	
	external_momenta = my_LTD.external_kinematics
	s = [p[0]**2-numpy.sum(p[1:]**2) for p in external_momenta]
	print my_LTD.analytic_three_point_ladder(s[0],s[1],s[2],3)

	result = my_LTD.ltd_integrate(my_LTD.ltd_integrand, N_refine=1000)


