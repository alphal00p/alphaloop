# n loop LTD

class LTDnLoop:

	def __init__(self,topology):
		self.loop_lines         = topology.loop_lines
		self.ltd_cut_structure  = topology.ltd_cut_structure
		self.n_loops            = topology.n_loops 
		self.name               = topology.name
		self.external_kinematics= topology.external_kinematics
		self.scale				= self.get_scale()
		self.possible_cuts		= self.get_possible_cuts()
	
	def delta(self,q_space,m):
		on_f_shell = adipy.sqrt(self.norm(q_space)**2 + m**2)
		return on_f_shell
	
	def get_scale(self):
		p_posE = [p for p in self.external_kinematics if p[0] > 0.]
		if -numpy.sum(self.external_kinematics,axis=0)[0] > 0.:
			p_posE += [-sum(self.external_kinematics)]
		assert(len(p_posE)!=0)
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
	
	def get_parents(self,ll_nr):
		""" find parent looplines of a given loopline """
		start = self.loop_lines[ll_nr].start_node
		parents = []
		for i,loop_line in enumerate(self.loop_lines):
			if i != ll_nr:
				if start == loop_line.start_node:
					parents += [(i,False)]
				elif start == loop_line.end_node:
					parents += [(i,True)]
		return parents
	
	def get_deltas(self,loop_momenta):
		""" computes the sqrt[(l_space+p_space)**2 + m**2] of all looplines
			looks like looplines_prop = [[q_i0p1],[q_i0p2],[q_i0p3]]"""
		deltas = []
		for loop_line in self.loop_lines:
			loop_momentum = numpy.sum([l*sign for l,sign in zip(loop_momenta,loop_line.signature)],axis=0)
			deltas += [[self.delta(loop_momentum + prop.q[1:],numpy.sqrt(prop.m_squared)) for prop in loop_line.propagators]]
		return deltas
	
	def inv_G(self,x,y):
		return x**2 - y**2
		
	def get_four_momenta(self,line_energies,loop_momenta):
		four_momenta = [numpy.zeros(4) for i in xrange(self.n_loops)]
		for loop_line in self.loop_lines:
			aa = numpy.array(loop_line.signature)
			if numpy.sum(aa**2) == 1:
				index = numpy.where(numpy.absolute(aa)==1)[0][0]
				four_momenta[index][1:] = aa[index]*loop_momenta[index]
				four_momenta[index][0] = aa[index]*line_energies[index]
		return four_momenta
	
	def evaluate_ll(self,loop_line,c,line_energy,delta):
		""" ll: loopline, c: positive (cut index+1), line_energy: energy comp of line momentum, ll_prop: Deltas of this line """
		""" evaluates a loop line at a fixed cut momentum """
		ll_residue = 1.
		for i,prop in enumerate(loop_line.propagators):
			if i != (c-1):
				ll_residue *= self.inv_G(line_energy + prop.q[0],delta[i])
			else:
				ll_residue *= 2*delta[i]
		return 1./ll_residue	
	
	def evaluate_cut(self,deltas,cut):
		""" cut looks like cut = [1,-2,0], 0 is no cut, the abs(c) is the (propagator index +1), - means negative cut
			computes the residue of cuts of specific propagators """
		line_energies = self.get_all_line_energies(deltas,cut)
		ll_residues = []
		for loop_line,c,line_energy,delta in zip(self.loop_lines,cut,line_energies,deltas):
			ll_residues += [self.evaluate_ll(loop_line,abs(c),line_energy,delta)]
		cut_residue = numpy.prod(ll_residues)
		return cut_residue
	
	def get_all_line_energies(self,deltas,cut):
		line_energies = [numpy.sign(c)*delta[abs(c)-1] - loop_line.propagators[abs(c)-1].q[0] if c != 0
						else None
						for loop_line,delta,c in zip(self.loop_lines,deltas,cut)]
		for ll_nr,line_energy in enumerate(line_energies):
			if line_energy == None:
				line_energies[ll_nr] = self.get_line_energy(ll_nr,line_energies)
		return line_energies
	
	def get_line_energy(self,ll_nr,line_energies):
		"""a recursive algorythm that computes the energy of the corresponding loop line,
			should work for loops larger than 2 as well"""
		new_line_energy = 0.
		parents = self.get_parents(ll_nr)
		for nr,sign in parents:
			if line_energies[nr] == None:
				line_energies[nr] = self.get_line_energy(nr,line_energies)
			new_line_energy += (-1.)**(sign+1)*line_energies[nr]
		return new_line_energy
				
	def evaluate_cut_structure(self,deltas,cut_structure):
		""" computes all residues of a cut structure """
		cut_indices = [i for i,line in enumerate(cut_structure) if abs(line)==1]
		cut_residues = []
		for ll_cuts in itertools.product(*[self.possible_cuts[i] for i in cut_indices]):
			cut = [0 for i in xrange(len(self.loop_lines))] #True/False arbitrary here
			for i,ll_cut in zip(cut_indices,ll_cuts):
				cut[i] = ll_cut*cut_structure[i]
			cut_residues += [self.evaluate_cut(deltas,cut)]
		cut_structure_residue = numpy.sum(cut_residues)
		return cut_structure_residue
				
	def get_possible_cuts(self):
		""" returns a list with (indices + 1) of all possible cuts per loop line """
		possible_cuts = [range(1,len(loop_line.propagators)+1) for loop_line in self.loop_lines]
		return possible_cuts
	
	def get_dual_vec(self,k_space,l_space):
		vec = numpy.append(k_space,l_space)
		dim = 3*self.n_loops
		vec_dual = [0.]*dim
		for i in xrange(dim):
			vec_dual[i] = adipy.ad(vec[i], numpy.array([0. if j !=i  else 1. for j in xrange(dim)]))
		return vec_dual
	
	def hardcoded_dt_def(self,k_space,l_space,p):
		# according to dario convention of loop line and momentum flow assignments
		vec_dual = self.get_dual_vec(k_space,l_space)
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
	
	def deform(self,loop_momenta):
		""" doesn't deform yet """
		kappas = [numpy.zeros(3) for i in xrange(self.n_loops)]
		jac = 1.
		
		#double triangle (with dario convention of loop line assignments!!)
		#kappa_k,kappa_l,jac = self.hardcoded_dt_def(loop_momenta[0],loop_momenta[1],self.loop_lines[0].propagators[1].q)
		#kappas = [kappa_k,kappa_l]
		
		return kappas, jac
		
	def ltd_integrand(self,x):
		loop_momenta = [numpy.zeros(3) for i in xrange(self.n_loops)]
		wgt = numpy.zeros(self.n_loops)
		for i in xrange(self.n_loops):
			loop_momenta[i], wgt[i] = self.parametrize_analytic(x[(i*3):(i+1)*3])
		kappas,jac = self.deform(loop_momenta)
		loop_momenta = [loop_momentum+1j*kappa for loop_momentum,kappa in zip(loop_momenta,kappas)]
		deltas = self.get_deltas(loop_momenta)
		full_residue = numpy.sum([self.evaluate_cut_structure(deltas,cut_structure) for cut_structure in self.ltd_cut_structure])
		ltd_factor = (-2*numpy.pi*1j)**self.n_loops
		standard_factor = (-1j/numpy.pi**2)**self.n_loops
		#alt_factor = (1./(2.*numpy.pi)**4)**self.n_loops
		return full_residue*numpy.prod(wgt)*ltd_factor*standard_factor*jac
	
	def ltd_integrate(self,integrand, N_train=1000, N_refine=1000):
		numpy.random.seed(0)
		batch_size = 600
		real_integr = lambda x: integrand(x).real
		imag_integr = lambda x: integrand(x).imag
		real_vegas3_integrator = vegas.Integrator(3 * self.n_loops * [[0., 1.]], nhcube_batch=batch_size)
		imag_vegas3_integrator = vegas.Integrator(3 * self.n_loops * [[0., 1.]], nhcube_batch=batch_size)
		# train
		real_vegas3_integrator(real_integr,nitn=7, neval=N_train)
		imag_vegas3_integrator(imag_integr,nitn=7, neval=N_train)
		# refine
		real_result = real_vegas3_integrator(real_integr,nitn=10, neval=N_refine)
		imag_result = imag_vegas3_integrator(imag_integr,nitn=10, neval=N_refine)
		print('\n'+real_result.summary()+'\n'+imag_result.summary())
		#print('\n'+real_result.summary())
		real_vegas3_integrator.map.show_grid()
		imag_vegas3_integrator.map.show_grid()
		return [real_result,imag_result]

if __name__ == "__main__":
	
	import numpy as numpy
	import scipy.special as sc
	import vegas
	import adipy as adipy
	import time
	import itertools
	import mpmath
	import math
	import ltd_commons
	#import warnings
	#warnings.simplefilter("error")
	#warnings.simplefilter("ignore", DeprecationWarning)
	
	t0 = time.time()

	print '='*(2*36+7) + '\n' + '='*36+' hello '+'='*36 + '\n' + '='*(2*36+7)
	
	def example():
		# example to evaluate specific residues at given loop momenta configuration
		my_topology = ltd_commons.hard_coded_topology_collection['DoubleTriange']
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
	
	my_topology = ltd_commons.hard_coded_topology_collection['DoubleTriange']
	my_LTD = LTDnLoop(my_topology)
	result = my_LTD.ltd_integrate(my_LTD.ltd_integrand,N_refine=1000)
	
	print '='*(2*36+7)
	print 'I = ', result[0], '+ i',result[1]
	print '='*(2*36+7)

	print 'Time: %.2f s' % (time.time()-t0)