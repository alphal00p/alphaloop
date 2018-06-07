#####################################################
#                                                   #
#  Source file of the pyNLoop GGVV MG5aMC plugin.   #
#  Use only with consent of its authors.            #
#                                                   #
#  author: Valentin Hirschi, Ben Ruij               #
#                                                   #
#####################################################

import os
import logging

import madgraph.integrator.integrators as integrators
from madgraph import InvalidCmd, MadGraph5Error, MG5DIR, ReadWrite, MPI_RANK, MPI_SIZE, MPI_ACTIVE
import madgraph.various.misc as misc
from madgraph.iolibs.files import cp, ln, mv

import nloop_integrands

from pySecDec.integral_interface import IntegralLibrary as pySecDecIntegralLibrary
import sympy as sp

pjoin = os.path.join

logger = logging.getLogger('pyNLoop.pySecDecIntegrator')
class pySecDecIntegratorError(MadGraph5Error):
    """ Error for pySecDecIntegrator."""

class pySecDecIntegrator(integrators.VirtualIntegrator):
    """ Wrapping class around pySecDec integration facilities."""

    def __init__(self, integrands, pySecDec_path=None, **opts):

        super(pySecDecIntegrator, self).__init__(integrands, **opts)

        if len(self.integrands)>1:
            raise pySecDecIntegratorError('pySecDecIntegrator only supports a single loop integration at once for now.')
        else:
            self.integrand=self.integrands[0]

        if not isinstance(self.integrand, nloop_integrands.NLoopIntegrand):
            raise pySecDecIntegratorError('pySecDecIntegrator can only integrate an integrand instance of NLoopIntegrand.')
        
        self.pySecDec_path = pySecDec_path
        self.options       = dict(opts)

        # On initialization, compile already
        logger.info("Compiling pySecDec output of integrand '%s' ..."%self.integrand.nice_string())
        self.compile_pySecDec_package()
    
    def compile_pySecDec_package(self):
        """ Compile the pySecDec package of the integrand specified."""
        
        # For some reason the default makefile does not correctly specify the path to its
        # include and libraries, so we soft-link them here locally for now
        necessary_lib     = ['libcuba.a','libgsl.a','libgslcblas.a']
        necessary_include = ['cuba.h','gsl','secdecutil']
        
        for lib in necessary_lib:
            ln(pjoin(self.pySecDec_path,'lib',lib),
               pjoin(self.integrand.output_folder,self.integrand._pySecDecOutputName))
        for include in necessary_include:
            ln(pjoin(self.pySecDec_path,'include',include),
               pjoin(self.integrand.output_folder,self.integrand._pySecDecOutputName))
        
        # Now we can hopefully successfully compile
        misc.compile(arg=[], 
                     cwd=pjoin(self.integrand.output_folder,self.integrand._pySecDecOutputName), 
                     mode = 'cpp', nb_core=1)

        if not os.path.isfile(pjoin(self.integrand.output_folder,
            self.integrand._pySecDecOutputName,'%s_pylink.so'%self.integrand._pySecDecOutputName)):
            raise pySecDecIntegratorError('Could not successfully compile pySecDec output.')
    
    def integrate(self):
        """ Steer the integration performed by pySecDec output code. """

        # Load c++ library
        pySecDec_integral = pySecDecIntegralLibrary(pjoin(self.integrand.output_folder,
            self.integrand._pySecDecOutputName,'%s_pylink.so'%self.integrand._pySecDecOutputName))
        
        # Choose an integrator
        pySecDec_integral.use_Vegas(
                flags   = 2,  # ``flags=2``: verbose --> see Cuba manual
                epsrel  = self.options['target_accuracy'], 
                # epsabs  = 1e-07, 
                maxeval = int(1e7)
        )
    
        # Integrate
        #logger.info("Starting integration of loop integrand '%s' with pySecDec ..."%
        #                                                      self.integrand.nice_string())
        with misc.Silence(active = self.options['verbosity']==0):
            str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = \
              pySecDec_integral(real_parameters=self.integrand.integrand_parameters_values)
        
        # Convert complex numbers from c++ to sympy notation
        str_integral_with_prefactor = str_integral_with_prefactor.replace(',','+I*')
        str_prefactor = str_prefactor.replace(',','+I*')
        str_integral_without_prefactor = str_integral_without_prefactor.replace(',','+I*')
        
        # Convert result to sympy expressions
        integral_with_prefactor = sp.sympify(
            str_integral_with_prefactor.replace('+/-','*value+error*')).expand().coeff('value')
        integral_with_prefactor_err = sp.sympify(
            str_integral_with_prefactor.replace('+/-','*value+error*')).expand().coeff('error')
        prefactor = sp.sympify(str_prefactor)
        integral_without_prefactor = sp.sympify(
            str_integral_without_prefactor.replace('+/-','*value+error*')).expand().coeff('value')
        integral_without_prefactor_err = sp.sympify(
            str_integral_without_prefactor.replace('+/-','*value+error*')).expand().coeff('error')
        
        return integral_with_prefactor, integral_with_prefactor_err
