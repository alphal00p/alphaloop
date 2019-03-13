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

	def __init__(self,topology,hyperparameters):
		self.ltd_cut_structure  = topology.ltd_cut_structure
		self.ltd_loop_lines		= topology.loop_lines
		self.n_loops            = topology.n_loops 
		self.name               = topology.name
		self.external_kinematics= topology.external_kinematics
		self.scale				= self.get_scale()
		self.analytical_result 	= topology.analytical_result
		self.cut_structures 	= [CutStructure(cut_structure=ltd_cut_structure, ltd_loop_lines = self.ltd_loop_lines) for ltd_cut_structure in self.ltd_cut_structure]
		self.lambda_ij			= hyperparameters['Deformation']['lambda']
		self.a_ij				= hyperparameters['Deformation']['additive']['a_ij']

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
		dual_vec = [adipy.ad(vec[i], numpy.array([0. if j !=i  else 1. for j in xrange(dim)])) for i in xrange(dim)]
		dual_momenta = [dual_vec[i*3:(i+1)*3] for i in xrange(self.n_loops)]
		return numpy.array(dual_momenta)
	
	def deform(self,loop_momenta):
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

	def ltd_integrand(self,x):
		assert(len(x) == self.n_loops*3)
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
		print '{:=^90}'.format(' '+ 'scale = %f' %self.scale + ' ')
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
		if opts != {}:
			print 'Warning: Not all options used!'

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
		cut_residues = []
		for cut in self.cuts:
			cut_residues += [cut.evaluate(deltas)]
		return numpy.sum(cut_residues)

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
	
	def evaluate(self,deltas):
		ll_residues = []
		for ll_nr,loop_line in enumerate(self.loop_lines):
			ll_residues += [loop_line.evaluate(deltas)]
		return numpy.prod(ll_residues)

class CutLoopLine(CutObject):

	def __init__(self,**opts):
		super(CutLoopLine, self).__init__(**opts)
		assert('loop_line_index' in opts)
		this_loop_line = self.ltd_loop_lines[self.loop_line_index]
		self.basis_signature = numpy.array(this_loop_line.signature).dot(self.inv_basis_trsf)
		self.delta_signature = self.basis_signature.dot(numpy.diag(self.cut_signature))
		opts.update({'delta_signature': self.delta_signature, 'basis_signature': self.basis_signature})
		self.propagators = [CutPropagator(propagator_index = a, **opts) for a,propagator in enumerate(this_loop_line.propagators)]

	def evaluate(self,deltas):
		ll_residue = 1.
		for a,propagator in enumerate(self.propagators):
			ll_residue *= propagator.evaluate(deltas)
		return ll_residue

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
	
		plus_ellipsoid = (all(sign >= 0 for sign in self.delta_signature) and self.total_shift[0] < 0.)
		minus_ellipsoid = (all(sign <= 0 for sign in self.delta_signature) and self.total_shift[0] > 0.)
		s = self.total_shift[0]**2 - self.norm(self.total_shift[1:])
		self.is_ellipsoid = (plus_ellipsoid or minus_ellipsoid) and (s > 0.)

	def deform(self,loop_momenta,deltas,**opts):
		scale = opts.pop('scale')
		a_ij = opts.pop('a_ij')
		lambda_ij = opts.pop('lambda_ij')
		kappas = [numpy.zeros(3) for loop_momentum in loop_momenta]
		if not self.is_cut:
			if self.is_ellipsoid:
				basis_shift_space = [shift[1:] for shift in self.basis_shift]
				cut_prop_momenta = self.basis_trsf.dot(loop_momenta) + basis_shift_space
				cut_prop_kappas = numpy.array([[coord*0. for coord in loop_momentum] for loop_momentum in loop_momenta])
				for i,sign in enumerate(self.delta_signature):
					if sign != 0:
						v1 = cut_prop_momenta[i]
						v2 = numpy.array(self.ltd_loop_lines[self.loop_line_index].signature).dot(loop_momenta) + self.own_shift[1:]
						cut_prop_kappas[i] = v1/self.norm(v1)
						cut_prop_kappas[i] += self.basis_signature[i]*v2/self.norm(v2)
				cut_deltas = self.get_cut_deltas(deltas)
				delta = deltas[self.loop_line_index][self.propagator_index]
				energy = self.delta_signature.dot(cut_deltas) + self.total_shift[0]
				interpolation = adipy.exp(-(energy-delta*numpy.sign(self.total_shift[0]))**2/scale**2/a_ij) 
				#interpolation = adipy.exp(-(energy**2-delta**2)**2/scale**4/a_ij) 
				kappas = self.inv_basis_trsf.dot(cut_prop_kappas)
				#lambda_factor = -0.01*scale
				lambda_factor = -abs(lambda_ij)
				kappas *= lambda_factor*interpolation
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
	
	my_topology = ltd_commons.hard_coded_topology_collection['AltDoubleTriangle']
	hyperparameters = ltd_commons.hyperparameters
	my_LTD = LTDnLoop(my_topology,hyperparameters)
		
	#print my_LTD.ltd_integrand4([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
	
	#external_momenta = my_LTD.external_kinematics
	#s = [p[0]**2-numpy.sum(p[1:]**2) for p in external_momenta]
	#print my_LTD.analytic_three_point_ladder(s[0],s[1],s[2],3)

	result = my_LTD.ltd_integrate(my_LTD.ltd_integrand, N_refine=1000)


