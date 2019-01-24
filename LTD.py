# One Loop LTD

class LTD:
	
	def __init__(self,config_string):
		self.p_i, self.m_i = self.__select_config(config_string)
		self.scale = self.__get_scale()
		self.tolerance = 1e-7
		self.k_i = self.__get_k_i()
		self.n = len(self.k_i)
		self.dim = len(self.k_i[0]) - 1
		self.sing_matrix = self.__find_singularity_matrix()
		self.lambda_ij,self.A_ij = self.__set_deformation(1.,1e-1,show=True)
		self.overall_lambda = 0.1
		print 'My overall scaling = ', self.overall_lambda
		self.curr_max_lambda = self.overall_lambda
	
	def __select_config(self,s):
		if s == '3d_bub':
			m1 = 0.
			m2 = 0.
			p1 = numpy.array([1,0,0])
			p_i = [p1]
			m_i = [m1,m2]
			s = p1[0]**2-self.norm(p1[1:])**2
			if s < 0:
				sqrt_s = numpy.sqrt(-s)
			elif s > 0:
				sqrt_s = -1j*numpy.sqrt(s)
			factor = numpy.pi**(3./2.)/(2*numpy.pi)**len(p_i[0])
			analytic = lambda s: factor*numpy.pi**(3./2.)/sqrt_s
			print 'analytic: ', analytic(s)
		elif s == '3d_test':
			m1 = 10.
			m2 = 50.
			p1 = numpy.array([1,3,7])
			p_i = [p1]
			m_i = [m1,m2]
		elif s == 'test':
			m1 = m2 = 0.
			m3 = 1.
			p1 = numpy.array([-2,0.1,0,.1])
			p2 = numpy.array([3,0,0,4])
			p_i = [p1,p2]
			m_i = [m1,m2,m3]
		elif s == 'neg_s_triangle':
			E_cm = 5000 #GeV
			q = E_cm/2.
			p1 = q*numpy.array([1.,0,0,1.])
			p2 = q*numpy.array([-1.,0,0.,1.])
			m = m1 = m2 = m3 = 100 #GeV
			p_i = [p1,p2]
			m_i = [m1,m2,m3]
			# factor to match with the definition in http://qcdloop.fnal.gov/bkvs_a11.pdf
			factor = (2*numpy.pi)**len(p_i[0])/(numpy.pi**2)
			analytic = lambda s,m: 1./factor*1./(2*s)*(numpy.log(-(1/2.*(1-numpy.sqrt(1.-4*m**2/s)))/(1/2.*(1+numpy.sqrt(1.-4*m**2/s)))))**2
			p3 = -p1 - p2
			s = p3[0]**2 - self.norm(p3[1:])**2
			print 'analytic: ', analytic(s,m)
		elif s == 'P1':
			#P1 in https://arxiv.org/pdf/1510.00187.pdf
			p1 = numpy.array([5.23923,-4.18858,0.74966,-3.05669])
			p2 = numpy.array([6.99881,-2.93659,5.03338,3.87619])
			m1 = m2 = m3 = 7.73358
			p_i = [p1,p2]
			m_i = [m1,m2,m3]
		elif s == 'P2':
			#P2 in https://arxiv.org/pdf/1510.00187.pdf
			p1 = numpy.array([13.42254,58.79478,-73.11858,-91.95015])
			p2 = numpy.array([81.65928,-68.52173,8.75578,-95.05353])
			m1 = 49.97454
			m2 = 86.92490
			m3 = 80.22567
			p_i = [p1,p2]
			m_i = [m1,m2,m3]
		elif s == 'P3':
			#P3 in https://arxiv.org/pdf/1510.00187.pdf
			p1 = numpy.array([10.51284,6.89159,-7.40660,-2.85795])
			p2 = numpy.array([6.45709,2.46635,5.84093,1.22257])
			m1 = m2 = m3 = 0.52559
			p_i = [p1,p2]
			m_i = [m1,m2,m3]
		elif s == 'P4':
			#P4 in https://arxiv.org/pdf/1510.00187.pdf
			p1 = numpy.array([95.77004,31.32025,-34.08106,-9.38565])
			p2 = numpy.array([94.54738,-53.84229,67.11107,45.56763])
			m1 = 83.02643
			m2 = 76.12873
			m3 = 55.00359
			p_i = [p1,p2]
			m_i = [m1,m2,m3]
		elif s == 'P5':
			p1 = numpy.array([31.54872,-322.40325,300.53015,-385.58013])
			p2 = numpy.array([103.90430,202.00974,-451.27794,-435.12848])
			p3 = numpy.array([294.76653,252.88958,447.09194,311.71630])
			m1 = m2 = m3 = m4 = 4.68481
			p_i = [p1,p2,p3]
			m_i = [m1,m2,m3,m4]
		elif s == 'P6':
			p1 = numpy.array([50.85428,-55.74613,11.69987,94.92591])
			p2 = numpy.array([0.69914,67.19262,-5.78627,91.52776])
			p3 = numpy.array([52.35768,76.32258,43.82222,13.05874])
			m1 = 54.29650
			m2 = 53.54058
			m3 = 55.96814
			m4 = 51.74438
			p_i = [p1,p2,p3]
			m_i = [m1,m2,m3,m4]
		elif s == 'P7':
			p1 = numpy.array([62.80274,-49.71968,-5.53340,-79.44048])
			p2 = numpy.array([48.59375,-1.65847,34.91140,71.89564])
			p3 = numpy.array([76.75934,-19.14334,-17.10279,30.22959])
			m1 = m2 = m3 = m4 = 9.82998
			p_i = [p1,p2,p3]
			m_i = [m1,m2,m3,m4]
		elif s == 'P8':
			p1 = numpy.array([98.04093,77.37405,30.53434,-81.88155])
			p2 = numpy.array([73.67657,-53.78754,13.69987,14.20439])
			p3 = numpy.array([68.14197,-36.48119,59.89499,-81.79030])
			m1 = 81.44869
			m2 = 94.39003
			m3 = 57.53145
			m4 = 0.40190
			p_i = [p1,p2,p3]
			m_i = [m1,m2,m3,m4]
		elif s == 'P9':
			p1 = numpy.array([76.50219,-72.36197,10.95225,-99.79612])
			p2 = numpy.array([99.02723,27.27133,-25.11907,86.10825])
			p3 = numpy.array([64.19420,13.10011,18.37737,-29.16095])
			m1 = m2 = 37.77809
			m3 = m4 = 36.84323
			p_i = [p1,p2,p3]
			m_i = [m1,m2,m3,m4]
		elif s == 'P10':
			p1 = numpy.array([13.62303,-64.20757,-17.59085,-8.81785])
			p2 = numpy.array([96.67650,89.65623,-18.47276,40.73203])
			p3 = numpy.array([66.21913,-39.49917,3.640139,-82.31669])
			m1 = m3 = 64.67282
			m2 = m4 = 51.13181
			p_i = [p1,p2,p3]
			m_i = [m1,m2,m3,m4]	
		else:
			print "No configuration with name ", s, " found." 
			return None, None
	
		propagator_string = '-i / (2*pi)^%d  x  ' %(len(p_i[0])) +'\int d^%dl / ' %(len(p_i[0]))
		for i in xrange(len(p_i)):
			propagator_string += '[(l+p%d)^2-m%d^2] ' %(i+1,i+1)
		propagator_string += '[l^2-m%d^2]' %(len(p_i)+1)
		print propagator_string
		for i,p in enumerate(p_i):
			print 'p%d = '%(i+1), p, 'p%d^2 = '%(i+1), p[0]**2-self.norm(p[1:])**2
		#p_n = numpy.sum(p_i,axis=0)
		#print 'p%d = '%(len(p_i)+1), p_n, 'p%d^2 = '%(len(p_i)+1), p_n[0]**2-self.norm(p_n[1:])**2
		m_string = {'m%d'%(i+1): m for i,m in enumerate(m_i)}
		print m_string
		return p_i, m_i
	
	def __get_k_i(self):
		k_i = [numpy.sum(self.p_i[:(i+1)],axis=0) for i in xrange(len(self.p_i))]
		k_i += [numpy.zeros(len(k_i[0]))]
		return k_i
	
	def __get_scale(self):
		p_posE = [p for p in self.p_i if p[0] > 0.]
		assert(len(p_posE)!=0)
		if -sum(self.p_i)[0] > 0.:
			p_posE += [-sum(self.p_i)]
		p_sum = sum(p_posE)
		s = p_sum[0]**2 - self.norm(p_sum[1:])**2
		scale = numpy.sqrt(abs(s))
		print 'scale = ', scale
		return scale
		
	def __find_singularity_matrix(self):
		sing_matrix = numpy.empty([self.n,self.n],dtype=str)
		for i in xrange(self.n):
			for j in xrange(self.n):
				if i != j:
					k_ji = self.k_i[j]-self.k_i[i]
					if k_ji[0]**2 - self.norm(k_ji[1:])**2 - (self.m_i[j]+self.m_i[i])**2 > 0. and k_ji[0] < 0.:
						sing_matrix[i,j] = 'E'
					elif k_ji[0]**2 - self.norm(k_ji[1:])**2 - (self.m_i[j]+self.m_i[i])**2 < 0.:
						sing_matrix[i,j] = 'H'
					elif k_ji[0]**2 - self.norm(k_ji[1:])**2 - (self.m_i[j]+self.m_i[i])**2 == 0.:
						sing_matrix[i,j] = 'c' #collinear
					else:
						sing_matrix[i,j] = '0'
				else:
					sing_matrix[i,i] = '0'
		print 'Singularity Matrix = \n', sing_matrix
		return sing_matrix
		
	def __set_deformation(self,l,A,show=True):
		assert(A>0.)
		lambda_ij = -1.*(numpy.ones([self.n,self.n]) - numpy.diag(numpy.ones([self.n])))
		A_ij = numpy.ones([self.n,self.n])
		lambda_ij *= l
		A_ij *= A
		if show:
			print 'lambda_ij ~ ', l
			print 'A_ij ~ ', A
		return lambda_ij, A_ij
	
	def norm(self,q):
		"""for complex and real q arrays
		not the same as numpy.linalg.norm for complex numbers!!"""
		return adipy.sqrt(numpy.sum(q**2))

	def q_0p(self,q_space,m):
		on_f_shell = adipy.sqrt(self.norm(q_space)**2 + m**2)
		return on_f_shell

	def inv_G_D(self,q_i0p,q_j0p,k_ji0):
		denominator = (q_i0p + k_ji0)**2-q_j0p**2
		return denominator

	def G(self,q_j,m_j):
		denominator = (q_j[0]**2-self.norm(q_j[1:])**2-m_j**2)
		return 1./denominator

	def dual_function(self,which,l_space):
		# Resiude at q_i0p[which]
		q_ispace = [l_space+k[1:] for k in self.k_i]
		q_i0p = [self.q_0p(q_space,self.m_i[i]) for i,q_space in enumerate(q_ispace)]
	 
		residue_factor = -1j*2.*numpy.pi
		delta_factor = 1./(2.*q_i0p[which])
		propagators = [1./self.inv_G_D(q_i0p[which],q_j0p,self.k_i[j][0]-self.k_i[which][0])
									for j,q_j0p in enumerate(q_i0p) if j != which]

		return residue_factor*delta_factor*numpy.prod(propagators)
	
	def parametrize_numerical(self,random):
		random = adipy.ad(random)
		radius = self.scale*random[0]/(1.-random[0])
		phi = 2*numpy.pi*random[1]
		if len(random)==3:	
			cos_theta = -1.+2.*random[2]
			sin_theta = adipy.sqrt(1.-cos_theta**2)
			l_space = numpy.array([radius*sin_theta*adipy.cos(phi), radius*sin_theta*adipy.sin(phi),radius*cos_theta])
		elif len(random)==2:
			l_space = numpy.array([radius*adipy.cos(phi), radius*adipy.sin(phi)])
		wgt = numpy.linalg.det(adipy.jacobian(l_space))
		l_space = numpy.array([l.nom for l in l_space])
		return l_space,wgt
			
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
	
	def deform_contour_numerical(self,l_space):
		""" Pick  lambda_ij < 0 and 0 on diag, A_ij > 0 """
		assert(all(i == 0. for i in numpy.diag(self.lambda_ij)))
		if all(i == 0. for i in self.lambda_ij.flatten()):
			return l_space, 1.
			
		l_space_dual = [0.]*self.dim
		for i in xrange(self.dim):
			l_space_dual[i] = adipy.ad(l_space[i], numpy.array([0. if j !=i  else 1. for j in xrange(self.dim)]))
		
		q_ispace = [l_space_dual+k[1:] for k in self.k_i]
		q_i0p = [self.q_0p(q_space,self.m_i[i]) for i,q_space in enumerate(q_ispace)]
		
		inv_G_D_ij = [[0. for i in xrange(self.n)] for i in xrange(self.n)]
		kappa_ij = [[[0. for i in xrange(self.dim)] for i in xrange(self.n)] for i in xrange(self.n)]
		surpr_ij = [[0. for i in xrange(self.n)] for i in xrange(self.n)]
		dir_ij = [[[0. for i in xrange(self.dim)] for i in xrange(self.n)] for i in xrange(self.n)]
		
		dir_i = numpy.array([q_space/self.norm(q_space) for q_space in q_ispace])

		for i in xrange(self.n):
			for j in xrange(self.n):
				if i != j:
					inv_G_D_ij[i][j]  = self.inv_G_D(q_i0p[i],q_i0p[j],self.k_i[j][0] - self.k_i[i][0])
					surpr_ij[i][j] = adipy.exp(-inv_G_D_ij[i][j]**2/(self.A_ij[i,j]*self.scale**4))
					dir_ij[i][j] = dir_i[i] + dir_i[j]
					kappa_ij[i][j] =  self.lambda_ij[i][j]*dir_ij[i][j]*surpr_ij[i][j]	
		
		kappa = numpy.sum(numpy.sum(kappa_ij,axis=0),axis=0)
		l_contour = l_space_dual + 1j*self.overall_lambda*kappa
		jac = adipy.jacobian(l_contour)
		det_jac = numpy.linalg.det(jac)
		
		l_contour = numpy.array([l.nom for l in l_contour])
		kappa = numpy.array([kapp.nom for kapp in kappa])
		
		self.find_curr_max_scaling(l_space,kappa)
		self.find_root_exp_max_scaling(l_space,kappa,exp_factor=0.1)
		
		return l_contour, det_jac
			
	def deform_contour_analytic(self,l_space):
		""" Pick  lambda_ij < 0 and 0 on diag, A_ij > 0 """
		assert(all(i == 0. for i in numpy.diag(self.lambda_ij)))
		if all(i == 0. for i in self.lambda_ij.flatten()):
			return l_space, 1.
		
		q_ispace = [l_space+k[1:] for k in self.k_i]
		q_i0p = [self.q_0p(q_space,self.m_i[i]) for i,q_space in enumerate(q_ispace)]
		
		inv_G_D_ij = numpy.zeros([self.n,self.n])
		kappa_ij = numpy.zeros([self.n,self.n,self.dim])
		surpr_ij = numpy.zeros([self.n,self.n])
		dir_ij = numpy.zeros([self.n,self.n,self.dim])
		jac_surpr_ij = numpy.zeros([self.n,self.n,self.dim])
		jac_ij = numpy.zeros([self.n,self.n,self.dim,self.dim])
		
		
		dir_i = numpy.array([q_space/self.norm(q_space) for q_space in q_ispace])
		
		jac_dir_i = numpy.array([(-numpy.dot(dir_i[i].reshape(self.dim,1),dir_i[i].reshape(1,self.dim))
										+numpy.diag(numpy.ones(self.dim)))/self.norm(q_space) for i,q_space in enumerate(q_ispace)])
		
		for i in xrange(self.n):
			for j in xrange(self.n):
				if i != j:
					k_ji0 = self.k_i[j][0] - self.k_i[i][0]
					inv_G_D_ij[i,j]  = self.inv_G_D(q_i0p[i],q_i0p[j],k_ji0)
					surpr_ij[i,j] = numpy.exp(-inv_G_D_ij[i,j]**2/(self.A_ij[i,j]*self.scale**4))
					dir_ij[i,j] = dir_i[i] + dir_i[j]
					kappa_ij[i,j] =  self.lambda_ij[i,j]*dir_ij[i,j]*surpr_ij[i,j]
					jac_surpr_ij[i,j] = 2.*((q_i0p[i] + k_ji0)*q_ispace[i]/q_i0p[i]-q_ispace[j])
					jac_surpr_ij[i,j] *= -2.*inv_G_D_ij[i,j]*surpr_ij[i,j]/(self.A_ij[i,j]*self.scale**4)
					jac_ij[i,j] = self.lambda_ij[i,j]*((jac_dir_i[i]+jac_dir_i[j])*surpr_ij[i,j]
									+ numpy.dot(dir_ij[i,j].reshape(self.dim,1),jac_surpr_ij[i,j].reshape(1,self.dim)))
		
		kappa = numpy.sum(numpy.sum(kappa_ij,axis=0),axis=0)
		l_contour = l_space + 1j*self.overall_lambda*kappa
		
		self.find_curr_max_scaling(l_space,kappa)
		self.find_root_exp_max_scaling(l_space,kappa,exp_factor=0.1)
		
		jac = numpy.sum(numpy.sum(jac_ij,axis=0),axis=0)
		full_jac = numpy.diag(numpy.ones(self.dim))+1j*self.overall_lambda*jac
		det_jac = numpy.linalg.det(full_jac)
	
		return l_contour, det_jac

	def test_function(self,l_space,mu,D):
		if self.dim == 3:
			normalize = 1./(-(1.-1j)*numpy.exp(-0.5*1j*D*numpy.pi)*numpy.pi**(3./2)*sc.gamma(D-3./2)/(numpy.sqrt(2)*sc.gamma(D)))
		elif self.dim == 2:
			normalize = 1./(1j**(1.-D)*numpy.pi/(D-1.))
		return mu**(2.*D-self.dim)/(self.norm(l_space)**2 + mu**2*1j)**D*normalize

	def test_integrand(self,x,mu):
		l_space,wgt = self.get_integration_point(x)
		l_space,wgt = self.deform_contour(l_space,wgt)
		result = wgt*self.test_function(l_space,mu,D=3.)
		return result
	
	def gaussian_integrand(self,x):
		l_space,wgt = self.get_integration_point(x)
		l_space,wgt = self.deform_contour(l_space,wgt)
		gaussian = numpy.exp(-self.norm(l_space)**2)
		normalize = 1./numpy.pi**(self.dim/2.)
		return wgt*normalize*gaussian
	
	def dual_integrand(self,x,which_duals):
		l_space,wgt = self.parametrize_analytic(x)
		l_contour,jac_wgt = self.deform_contour_analytic(l_space)
		residues = [self.dual_function(this_dual,l_contour) for this_dual in which_duals]
		factor = -1j/(2*numpy.pi)**(self.dim+1)
		result = wgt*jac_wgt*sum(residues)*factor
		return result

	def integrate(self,integrand, N_train=1000, N_refine=1000,share_grid=False):
		batch_size = 600
		if share_grid:
			integr = lambda x: [integrand(x).real,integrand(x).imag]
			vegas3_integrator = vegas.Integrator(self.dim * [[0., 1.]], nhcube_batch=batch_size)
			# train
			vegas3_integrator(integr,nitn=7, neval=N_train)
			# refine
			real_result,imag_result = vegas3_integrator(integrand,nitn=10, neval=N_refine)
			vegas3_integrator.map.show_grid()
		else:
			real_integr = lambda x: integrand(x).real
			imag_integr = lambda x: integrand(x).imag
			real_vegas3_integrator = vegas.Integrator(self.dim * [[0., 1.]], nhcube_batch=batch_size)
			imag_vegas3_integrator = vegas.Integrator(self.dim * [[0., 1.]], nhcube_batch=batch_size)
			# train
			real_vegas3_integrator(real_integr,nitn=7, neval=N_train)
			imag_vegas3_integrator(imag_integr,nitn=7, neval=N_train)
			# refine
			real_result = real_vegas3_integrator(real_integr,nitn=10, neval=N_refine)
			imag_result = imag_vegas3_integrator(imag_integr,nitn=10, neval=N_refine)
			print('\n'+real_result.summary()+'\n'+imag_result.summary())
			real_vegas3_integrator.map.show_grid()
			imag_vegas3_integrator.map.show_grid()
		return [real_result,imag_result]

	def rot(self,k,n,cos_theta):
		l = n*(numpy.dot(n,k)) + cos_theta*numpy.cross(numpy.cross(n,k),n)+numpy.sqrt(1.-cos_theta**2)*numpy.cross(n,k)
		return l
	
	def generate_point_on_ff_intersect(self,k_a,m_a,k_b,m_b):
		"""	generate random point on ff intersection (hyperboloid) """
		# REMEMBER FOCI ARE AT -k_a and -k_b respectively
		# is SYMMETRIC in k_a,m_a,k_b,m_b = k_b,m_b,k_a,m_a
		
		# reducing the general case to to the case when k_ba0 > 0
		if (k_b[0]-k_a[0]) < 0:
			k_a,m_a,k_b,m_b = k_b,m_b,k_a,m_a
		k_ba = k_b-k_a
		k_ba0 = k_ba[0]
		k_baspace = k_ba[1:]
		
		cos_theta = 2*(numpy.random.rand(1)[0])-1.
	
		# solve A l^2 + B l + C = 0
		C = k_ba0**4 + (self.norm(k_baspace)**2 -m_a**2+m_b**2)**2 - 2*k_ba0**2*(self.norm(k_baspace)**2 +m_a**2+m_b**2)
		B = 4.*cos_theta*self.norm(k_baspace)*(self.norm(k_baspace)**2 - k_ba0**2 - m_a**2 + m_b**2)
		A = -4*(k_ba0**2-cos_theta**2*self.norm(k_baspace)**2)
		
		absq_aplus = (-B + numpy.sqrt(B**2-4*A*C))/(2*A)
		absq_aminus = (-B - numpy.sqrt(B**2-4*A*C))/(2*A)
		
		q_b0p = lambda q_a: numpy.sqrt(q_a**2+ self.norm(k_baspace)**2 + 2*self.norm(k_baspace)*q_a*cos_theta+m_b**2)
		q_a0p = lambda q_a: numpy.sqrt(q_a**2+m_a**2)
		denominator = lambda q_a: q_a0p(q_a) - q_b0p(q_a) + k_ba0
		
		ness_cond = lambda absq_a: absq_a >= 0. #magnitude obviously > 0
		suff_cond = lambda absq_a: abs(denominator(absq_a)) < self.tolerance*self.scale
		is_solution = lambda absq_a: ness_cond(absq_a) and suff_cond(absq_a)
		
		if is_solution(absq_aplus) and is_solution(absq_aminus) and absq_aplus != absq_aminus:
			print 'Two solutions found, this is suspicious. We pick the first one.'
		if is_solution(absq_aplus):
			absq_a = absq_aplus
		elif is_solution(absq_aminus):
			absq_a = absq_aminus
		else:
			print 'No forward-forward intersection found with random cos_theta = ',cos_theta
			return numpy.array([None,None,None])
	
		q_aspace = self.generate_fixed_SP_vector(absq_a,k_baspace,cos_theta)
		l = q_aspace - k_a[1:]
	
		# check total
		total_cond = self.q_0p(l+k_a[1:],m_a) - self.q_0p(l+k_b[1:],m_b)+k_ba0 < self.tolerance*self.scale
		assert(total_cond)
		return l
	
	def generate_fixed_SP_vector(self,magnitude,dir,cos_theta):
		""" generates a (random in dim>2) vector q with length magnitude, so that
			the scalar product q*dir = magnitude*abs(dir)*cos_theta """
		if self.norm(dir) != 0.:
			n_dir = dir/self.norm(dir)
		else:
			n_dir = numpy.random.rand(self.dim)
			n_dir = n_dir/self.norm(n_dir)
		
		if self.dim == 3:
			# n is perpendicular to k
			n_a, n_b = numpy.random.rand(2)
			c = numpy.where(n_dir != 0.)[0][0]
			a,b = [i for i in xrange(self.dim) if i != c]
			n = numpy.empty(3)
			n[a] = n_a
			n[b] = n_b
			n[c] = -(n_dir[a]*n_a+n_dir[b]*n_b)/n_dir[c]
			n = n/self.norm(n)
			q = self.rot(n_dir,n,cos_theta)*magnitude
		elif self.dim == 2:
			q = numpy.array([n_dir[0]*cos_theta - n_dir[1]*numpy.sqrt(1-cos_theta**2),
					n_dir[0]*numpy.sqrt(1-cos_theta**2) + n_dir[1]*cos_theta])
			q *= magnitude

		# check rotation
		rot_cond = numpy.abs(numpy.dot(q,dir) - cos_theta*self.norm(dir)*magnitude) < self.tolerance*self.scale**2
		assert(rot_cond)
		return q
	
	def generate_point_on_fb_intersect(self,k_a,m_a,k_b,m_b):
		""" generate random point on fb intersection (ellipsoid) """
		# REMEMBER FOCI ARE AT -k_a and -k_b respectively
		# solve q_0p(q_a,m_a) + q_0p(q_b,m_b) + k_ba0 = 0 for abs(q_a) depending on cos_theta, the angle between q_a,q_b
		# NOT SYMMETRIC in k_a,m_a,k_b,m_b = k_b,m_b,k_a,m_a
		k_ba = k_b-k_a
		if k_ba[0] > 0.:
			print 'No forward-backward intersection for k_ba0 =',k_ba[0]
			return
		k_ba0 = k_ba[0]
		k_baspace = k_ba[1:]

		cos_theta = 2*(numpy.random.rand(1)[0])-1.
	
		# solve A l^2 + B l + C = 0
		C = k_ba0**4 + (self.norm(k_baspace)**2 -m_a**2+m_b**2)**2 - 2*k_ba0**2*(self.norm(k_baspace)**2 +m_a**2+m_b**2)
		B = 4.*cos_theta*self.norm(k_baspace)*(self.norm(k_baspace)**2 - k_ba0**2 - m_a**2 + m_b**2)
		A = -4*(k_ba0**2-cos_theta**2*self.norm(k_baspace)**2)
	
		absq_aplus = (-B + numpy.sqrt(B**2-4*A*C))/(2*A)
		absq_aminus = (-B - numpy.sqrt(B**2-4*A*C))/(2*A)
	
		q_b0p = lambda q_a: numpy.sqrt(q_a**2+ self.norm(k_baspace)**2 + 2*self.norm(k_baspace)*q_a*cos_theta+m_b**2)
		q_a0p = lambda q_a: numpy.sqrt(q_a**2+m_a**2)
		denominator = lambda q_a: q_a0p(q_a) + q_b0p(q_a) + k_ba0
		
		ness_cond = lambda absq_a: absq_a >= 0. # necessary condition magnitude obviously > 0
		suff_cond = lambda absq_a: abs(denominator(absq_a)) < self.tolerance*self.scale
		is_solution = lambda absq_a: ness_cond(absq_a) and suff_cond(absq_a)
		
		if is_solution(absq_aplus) and is_solution(absq_aminus) and absq_aplus != absq_aminus:
			print 'Two solutions found, this is suspicious. We pick the first one.'
		if is_solution(absq_aplus):
			absq_a = absq_aplus
		elif is_solution(absq_aminus):
			absq_a = absq_aminus
		else:
			print 'No forward-backward intersection found.'
			return numpy.array([None,None,None])
		
		q_aspace = self.generate_fixed_SP_vector(absq_a,k_baspace,cos_theta)
		l = q_aspace - k_a[1:]
		
		# check total
		total_cond = self.q_0p(l+k_a[1:],m_a) + self.q_0p(l+k_b[1:],m_b)+k_ba0 < self.tolerance*self.scale
		assert(total_cond)
		return l
	
	def test_sign_on_fb_intersect(self):
		"""test sign on all ellipsoid intersections"""
		eps_sing = []
		for i in xrange(self.n):
			for j in xrange(self.n):
				if j != i:
					if self.sing_matrix[i,j] == 'E':
						eps_sing += [[i,j]]
	
		for tuple in eps_sing:	
			which_f = tuple[0]
			which_b = tuple[1]
			k_f = self.k_i[which_f]
			m_f = self.m_i[which_f]
			k_b = self.k_i[which_b]
			m_b = self.m_i[which_b]
			k_bf0 = k_b[0]-k_f[0]
	
			for i in xrange(1000):
				l_space = self.generate_point_on_fb_intersect(k_f,m_f,k_b,m_b)
				if any(x == None for x in l_space):
					print 'Sign test failed'
					return
				l_contour, wgt = self.deform_contour_analytic(l_space)
				q_fspace = l_contour + k_f[1:]
				q_bspace = l_contour + k_b[1:]
				dual_prop = self.inv_G_D(self.q_0p(q_fspace,m_f),self.q_0p(q_bspace,m_b),k_bf0)
				# ASSERT
				if not (numpy.sign(dual_prop.imag) == -1*numpy.sign(k_bf0)):
					print 'dual_prop ', dual_prop
		print 'Sign test successful'
		return
	
	def update_curr_max_scaling(self,X_j,Y_j):
		if 2*X_j.real < Y_j.real:
			lambda_j_squared = Y_j/4.
		elif 0 < Y_j.real < 2*X_j.real:
			lambda_j_squared = X_j - Y_j/4.
		elif Y_j.real < 0:
			lambda_j_squared = X_j - Y_j/2.
		elif Y_j.real == 0:
			return
		lambda_j = numpy.sqrt(lambda_j_squared)
		assert(lambda_j.imag == 0.)
		lambda_j = lambda_j.real
		if lambda_j < self.curr_max_lambda:
			self.curr_max_lambda = lambda_j
			print 'Required max. scaling currently', self.curr_max_lambda
		return
	
	def solve_quadr_eq(self,A,B,C):
		sqrt_X_j = 1j*B/(2.*A) #we need this for the sign
		X_j = sqrt_X_j**2
		Y_j = -C/A
		lambda_p = 1j*sqrt_X_j + numpy.sqrt(Y_j-X_j)
		lambda_m = 1j*sqrt_X_j - numpy.sqrt(Y_j-X_j)
		return lambda_p, lambda_m, X_j, Y_j
	
	def find_curr_max_scaling(self,l_space,kappa):
		for i in xrange(self.n):
			for j in xrange(self.n):
				if j == i: # "propagator" coming from delta function cut
					q_a = self.k_i[i][1:] + l_space
					m_a = self.m_i[i]
					A = -self.norm(kappa)**2 #REAL
					B = 2j*numpy.dot(q_a,kappa) #IMAGINARY
					C = self.norm(q_a)**2 + m_a**2 #REAL
					prop = lambda x: self.q_0p(q_a+1j*x*kappa,m_a)	
					cond = lambda x: abs(prop(x).real) < self.tolerance*self.scale and abs(prop(x).imag) < self.tolerance*self.scale
				else: # Dual propagators
					q_a = self.k_i[i][1:] + l_space
					q_b = self.k_i[j][1:] + l_space
					k_ba0 = self.k_i[j][0]-self.k_i[i][0]
					m_a = self.m_i[i]
					m_b = self.m_i[j]
					C1 = self.norm(q_a)**2 + m_a**2 + k_ba0**2 - self.norm(q_b)**2 - m_b**2
					B1 = 2*1j*numpy.dot(q_a-q_b,kappa)
					A = B1**2 + 4.*k_ba0**2*self.norm(kappa)**2 #symmetric, REAL
					B = 2*B1*C1 - 4.*k_ba0**2*2j*numpy.dot(q_a,kappa) #not symmetric but seems to be if kappa small, IMAGINARY
					C = C1**2 - 4.*k_ba0**2*(self.norm(q_a)**2 + m_a**2) #not symmetric but seems to be if kappa small, REAL
					prop = lambda x: (self.q_0p(q_a+1j*x*kappa,m_a) + k_ba0)**2 - (self.q_0p(q_b+1j*x*kappa,m_b))**2
					alt_prop = lambda x: (self.q_0p(q_a+1j*x*kappa,m_a) - k_ba0)**2 - (self.q_0p(q_b+1j*x*kappa,m_b))**2
					cond = lambda x: abs(prop(x).real) < self.tolerance*self.scale**2 and abs(prop(x).imag) < self.tolerance*self.scale**2
					alt_cond = lambda x: abs(alt_prop(x).real) < self.tolerance*self.scale**2 and abs(alt_prop(x).imag) < self.tolerance*self.scale**2
				
				assert(A.imag == 0)
				assert(B.real == 0)
				assert(C.imag == 0)
					
				if abs(A.real) > 1e-50*self.scale**2:
					lambda_p, lambda_m, X_j, Y_j = self.solve_quadr_eq(A,B,C)
				else:
					continue
				
				is_large = lambda x: abs(x) > self.curr_max_lambda				
				if is_large(lambda_p) and is_large(lambda_m):
					continue
				else:
					if cond(lambda_p) or cond(lambda_m):
						self.update_curr_max_scaling(X_j,Y_j)
					else:
						continue
		return
	
	def find_root_exp_max_scaling(self,l_space,kappa,exp_factor=0.1):
		for exp_f in [exp_factor,-exp_factor]: # expand factor negative gives different solutions but is also valid
			for i in xrange(self.n):
				q_a = self.k_i[i][1:] + l_space
				m_a = self.m_i[i]
				A = -self.norm(kappa)**2
				B = -2j*numpy.dot(q_a,kappa)/exp_f
				C = self.norm(q_a)**2 + m_a**2
			
				if abs(A.real) > 1e-50*self.scale**2:
					lambda_p, lambda_m, X_j, Y_j = self.solve_quadr_eq(A,B,C)
				else:
					continue
				
				is_large = lambda x: abs(x) > self.curr_max_lambda
				prop = lambda x: self.norm(q_a)**2 + m_a**2 - x**2*self.norm(kappa)**2 - x*2j*numpy.dot(q_a,kappa)/exp_f
				cond = lambda x: abs(prop(x).real) < self.tolerance*self.scale**2 and abs(prop(x).imag) < self.tolerance*self.scale**2

				assert(cond(lambda_p) and cond(lambda_m))
					
				if is_large(lambda_p) and is_large(lambda_m):
					continue
				else:
					self.update_curr_max_scaling(X_j,Y_j)
		return
		
			
if __name__ == "__main__":
	
	import numpy as numpy
	import scipy.special as sc
	import vegas
	import adipy as adipy
	
	print '='*(2*36+7) + '\n' + '='*36+' hello '+'='*36 + '\n' + '='*(2*36+7)
	
	my_LTD = LTD('P3')
	my_LTD.test_sign_on_fb_intersect()

	#pair = [0,1]
	#k_a, m_a, k_b, m_b = my_LTD.k_i[pair[0]],my_LTD.m_i[pair[0]],my_LTD.k_i[pair[1]],my_LTD.m_i[pair[1]]
	#l = my_LTD.generate_point_on_fb_intersect(k_a,m_a,k_b,m_b)

	all_duals = range(my_LTD.n)
	integr = lambda x: my_LTD.dual_integrand(x,which_duals=all_duals)
	
	result = my_LTD.integrate(integr,N_refine=1000,share_grid=False)
	
	print '='*(2*36+7)
	print 'I = ', result[0], '+ i',result[1] 
	print '='*(2*36+7)