################################################################################
#
# Copyright (c) 2009 The MadGraph5_aMC@NLO Development team and Contributors
#
# This file is a part of the MadGraph5_aMC@NLO project, an application which
# automatically generates Feynman diagrams and matrix elements for arbitrary
# high-energy processes in the Standard Model and beyond.
#
# It is subject to the MadGraph5_aMC@NLO license which should accompany this
# distribution.
#
# For more information, visit madgraph.phys.ucl.ac.be and amcatnlo.web.cern.ch
#
################################################################################

import traceback
import os
import logging
import math
import shutil
import numpy as np
import sys

if __name__ == '__main__':
    sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.path.pardir, os.path.pardir))

import madgraph

MADEVENT= False
import madgraph.various.misc as misc
import madgraph.iolibs.files as files
import madgraph.various.cluster as cluster
import madgraph.various.lhe_parser as lhe_parser
import alpha_loop.integrator.integrands as integrands
import alpha_loop.integrator.functions as functions

logger = logging.getLogger('madgraph.integrator')
pjoin = os.path.join

class IntegratorError(Exception):
    """Exception raised if an exception is triggered in integrators.""" 

class VirtualIntegrator(object):
    """A mother base class that specifies the feature that any integrator should implement."""
    
    def __init__(self, integrands, **opts):
        if isinstance(integrands, list):
            self.integrands = integrands
        else:
            self.integrands = [integrands]

    def integrate(self):
        """ Return the final integral and error estimators. """
        raise NotImplementedError
    
    def get_name(self):
        """ Returns the integrator name."""
        return self.__class__.__name__

    def observables_normalization(self, n_integrand_calls):
        """ Given the number of integrand calls, return the appropriate overall normalization
        to apply to the observables."""
        return 1.0/float(n_integrand_calls)

class SimpleMonteCarloIntegrator(VirtualIntegrator):
    """ Implements the simplest Monte-Carlo integrator. """
    
    def __init__(self, integrands, 
                 accuracy_target=0.01,
                 n_iterations=None,
                 n_points_per_iterations=100,
                 verbosity = 2,
                 save_points_to_file = None,
                 **opts):
        """ Initialize the simplest MC integrator."""
        
        self.accuracy_target = accuracy_target
        self.n_iterations = n_iterations
        self.n_points_per_iterations = n_points_per_iterations
        self.verbosity = verbosity
        self.save_points_to_file = save_points_to_file
        self.tot_func_evals = 0

        super(SimpleMonteCarloIntegrator, self).__init__(integrands, **opts)
        #misc.sprint(self.integrands)
        
    def integrate(self):
        """ Return the final integral and error estimates."""
        
        iteration_number   = 0
        sum_int            = 0.0
        sum_squared        = 0.0
        n_points           = 0
        error_estimate     = sys.maxsize
        integral_estimate  = 0.0
        
        for i, integrand in enumerate(self.integrands):
            integrand.counter = 0

        phase_space_volumes = [integrand.get_dimensions().volume() for integrand in self.integrands]
    
        if self.save_points_to_file is not None:
            logger.info("Saving all integration sample points to file '%s'."%self.save_points_to_file)
            out_stream = open(self.save_points_to_file, 'w')
            out_stream.write('IntegrandNumber, continuous_dimensions, discrete_dimension, complete_weight\n')
        else:
            out_stream = None

        points_to_write_to_file = []
        while (self.n_iterations is None or iteration_number < self.n_iterations ) and \
              (self.accuracy_target is None or error_estimate/(1e-99+integral_estimate) > self.accuracy_target):   

            iteration_number += 1
            n_curr_points = 0
            while n_curr_points < self.n_points_per_iterations:
                n_points += 1
                n_curr_points += 1
                # Compute phase-space volume
                new_wgt = 0.0
                for i, integrand in enumerate(self.integrands):
                    discrete_dimensions = integrand.discrete_dimensions.random_sample()
                    continuous_dimensions = integrand.continuous_dimensions.random_sample()
                    try:
                        new_wgt += phase_space_volumes[i]*integrand(continuous_dimensions,discrete_dimensions)
                        self.tot_func_evals += 1
                    except AssertionError as err:
                        traceback.print_tb(sys.exc_info()[-1])
                        logger.warning('Assertion error encountered.')
                        pass
                    if out_stream is not None:
                        points_to_write_to_file.append(', '.join([str(i+1),str(list(continuous_dimensions)),str(list(discrete_dimensions)),'%.16e'%new_wgt]))
                    if len(points_to_write_to_file) > 1 and len(points_to_write_to_file)%1000==1:
                        out_stream.write('\n'.join(points_to_write_to_file)+'\n')
                        points_to_write_to_file = []

                sum_int += new_wgt
                sum_squared += new_wgt**2

            integral_estimate = sum_int / n_points
            error_estimate =  math.sqrt( ((sum_squared / n_points) - integral_estimate**2)/n_points)
            msg = '%s :: iteration # %d / %s :: point #%d :: %.4e +/- %.2e'%(
                self.__class__.__name__, iteration_number, 
                '%d'%self.n_iterations if self.n_iterations else 'inf' ,n_points, 
                integral_estimate, error_estimate)
            if self.verbosity > 0:
                logger.info(msg)
        
        if out_stream is not None:
            out_stream.write('\n'.join(points_to_write_to_file))
            points_to_write_to_file = []
            out_stream.close()
            logger.info("Saved all integration sample points to file '%s'."%self.save_points_to_file)

        return integral_estimate, error_estimate
   
                
if __name__ == "__main__":

    # Example usage of the new integrator framework
    
    # First build the integrand
    my_integrand = integrands.VirtualIntegrand()
    # Define its dimensions
    my_integrand.set_dimensions( integrands.DimensionList([
            integrands.ContinuousDimension('x',lower_bound=2.0, upper_bound=5.0),
            integrands.ContinuousDimension('y' ,lower_bound=1.0, upper_bound=3.0)
             ]) )
                    
    # Define its constituting functions
    my_integrand.function_list = functions.FunctionList([
                functions.FunctionFromPythonExpression('math.sin(x/y)', dimensions=my_integrand.get_dimensions()),
                functions.FunctionFromPythonExpression('math.cos(x*y)', dimensions=my_integrand.get_dimensions())
            ])

    # Then the integrator
    my_integrator = SimpleMonteCarloIntegrator(my_integrand, 
        n_iterations=50, n_points_per_iterations=10000, accuracy_target=None)
    
    # Finally integrate
    print('\nFinal result: %.4e +/- %.2e'%my_integrator.integrate())
