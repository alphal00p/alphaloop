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

import sys
import os
import logging
import threading
import signal
import shutil
import traceback
from ctypes import c_bool
from multiprocessing import Value, Array, Event

pjoin = os.path.join

logger = logging.getLogger('pyCubaIntegrator')

import alpha_loop.integrator.cuba_interface as pycuba

try:
    import numpy as np
except ImportError:
    logger.critical("ERROR :: The pyCubaIntegrator requires the numpy python module. Install it first.")
    sys.exit(1)

MADEVENT= False
import madgraph.various.misc as misc
import madgraph.iolibs.files as files
import alpha_loop.integrator.integrators as integrators
from alpha_loop.integrator.integrators import IntegratorError
import alpha_loop.integrator.integrands as integrands
import alpha_loop.integrator.functions as functions
import alpha_loop.utils as utils

class tmpGlobal(object):
    """Designed so as to temporary add a global variable
    """

    def __init__(self, name, value):
        self.global_name  = name
        self.global_value = value

    def __enter__(self):
        globals()[self.global_name] = self.global_value

    def __exit__(self, ctype, value, traceback ):
        try:
            del globals()[self.global_name]
        except KeyError:
            pass

class SignalHandler:
    """
    The object that will handle signals and stop the worker threads.
    """

    #: The stop event that's shared by this handler and threads.
    stopper = None

    def __init__(self, stopper):
        self.stopper = stopper

    def __call__(self, signum, frame):
        """
        This will be called by the python signal module

        https://docs.python.org/3/library/signal.html#signal.signal
        """
        logger.critical('Abort signal detected.')
        with self.stopper.get_lock():
            self.stopper.value = True

class pyCubaIntegrator(integrators.VirtualIntegrator):

    def __init__(self, integrands, **opts):
        """ Initialize the pyCuba integrator."""

        # General steering parameters
        default_opts= {
           'algorithm':'Vegas',
           'verbosity':2,
           'target_accuracy':1.0e-3,
           'state_file_folder' : '.'
        }

        # Number of cores to be used in Cuba.
        # Recover nb_core in the option'cluster'
        if 'cluster' in opts and hasattr(opts['cluster'], 'nb_core'):
            default_opts['nb_core'] = opts['cluster'].nb_core
        else:
            # By default turn off parallelization as it is potentially problematic
            # (see warning message later).
            default_opts['nb_core'] = 0

        # Maximum number of evaluation before returning an answer
        default_opts['max_eval'] = int(1e10)
        # Minimum number of evaluation before returning an answer
        default_opts['min_eval'] = 5000

        # Parameter relevant for Vegas integration method
        # ------------------------------------------------
        
        # the number of integrand evaluations per iteration to start with. 
        default_opts['n_start'] = 1000
        # the increase in the number of integrand evaluations per iteration. 
        default_opts['n_increase'] = 500
        
        # the number of integrand evaluations per iteration to start with during the survey
        default_opts['n_start_survey'] = 1000
        # the increase in the number of integrand evaluations per iteration. 
        default_opts['n_increase_survey'] = 500
        default_opts['target_accuracy_survey'] = 5.0e-2
        # Maximum number of evaluation before returning an answer during the survey
        # Set to 0 if you want to skip the survey
        default_opts['max_eval_survey'] = 10000

        default_opts['save_grids'] = None
        default_opts['load_grids'] = None
        default_opts['keep_last_run_state_file'] = True

        # the batch size for sampling.
        # Vegas samples points not all at once, but in batches of size nbatch, to avoid exces-
        # sive memory consumption. 1000 is a reasonable value, though it should not affect
        # performance too much. 
        default_opts['n_batch'] = 100
        # Random seed
        default_opts['seed'] = None

        # Parameter relevant for Suave integration method
        # ------------------------------------------------
        
        # the number of new integrand evaluations in each subdivision.
        default_opts['n_new'] = 1000
        # the minimum number of samples a former pass must contribute
        # to a subregion to be considered in that region's compound integral value. Increasing
        # nmin may reduce jumps in the chi^2 value.
        default_opts['n_min'] = 2

        # the parameter p in Eq. (1), i.e. the type of norm
        # used to compute the fluctuation of a sample. This determines how prominently 'out-liers'
        # . i.e. individual samples with a large fluctuation, figure in the total fluctuation,
        # which in turn determines how a region is split up. As suggested by its name, flatness
        # should be chosen large for 'flat' integrands and small for 'volatile' integrands with
        # high peaks. Note that since flatness appears in the exponent, one should not use
        # too large values (say, no more than a few hundred) lest terms be truncated internally
        # to prevent overflow.
        default_opts['flatness'] = 50.

        # Parameter relevant for Divonne integration method
        # ------------------------------------------------

        default_opts['key1'] = 47
        # *key1*: determines sampling in the partitioning phase:
        # key1 = 7, 9, 11, 13 selects the cubature rule of degree key1. Note that the degree-11
        # rule is available only in 3 dimensions, the degree-13 rule only in 2 dimensions.
        # For other values of key1, a quasi-random sample of n1 = |key1| points is used, where
        # the sign of key1 determines the type of sample,
        # - key1 > 0, use a Korobov quasi-random sample,
        # - key1 < 0, use a "standard" sample (a Sobol quasi-random sample if seed = 0,
        # otherwise a pseudo-random sample).

        default_opts['key2'] = 1
        # *key2*: determines sampling in the final integration phase:
        # key2 = 7, 9, 11, 13 selects the cubature rule of degree key2. Note that the degree-11
        # rule is available only in 3 dimensions, the degree-13 rule only in 2 dimensions.
        # For other values of key2, a quasi-random sample is used, where the sign of key2
        # determines the type of sample,
        # - key2 > 0, use a Korobov quasi-random sample,
        # - key2 < 0, use a 'standard' sample (see description of key1 above),
        # and n2 = |key2| determines the number of points,
        # - n2 >= 40, sample n2 points,
        # - n2 < 40, sample n2 need points, where nneed is the number of points needed to
        # reach the prescribed accuracy, as estimated by Divonne from the results of the
        # partitioning phase.

        default_opts['key3'] = 1
        # *key3* sets the strategy for the refinement phase:
        # key3 = 0, do not treat the subregion any further.
        # key3 = 1, split the subregion up once more.
        # Otherwise, the subregion is sampled a third time with key3 specifying the sampling
        # parameters exactly as key2 above.

        default_opts['max_pass'] = 5
        # *maxpass*: controls the thoroughness of the partitioning phase: The
        # partitioning phase terminates when the estimated total number of integrand evalu-
        # ations (partitioning plus final integration) does not decrease for maxpass successive
        # iterations.
        # A decrease in points generally indicates that Divonne discovered new structures of
        # the integrand and was able to find a more effective partitioning. maxpass can be
        # understood as the number of 'safety' iterations that are performed before the par-
        # tition is accepted as final and counting consequently restarts at zero whenever new
        # structures are found.

        default_opts['border'] = 1.e-5
        # *border*: the width of the border of the integration region.
        # Points falling into this border region will not be sampled directly, but will be extrap-
        # olated from two samples from the interior. Use a non-zero border if the integrand
        # subroutine cannot produce values directly on the integration boundary.
        
        default_opts['maxchisq'] = 10.
        # *maxchisq*: the maximum chisquare value a single subregion is al-
        # lowed to have in the final integration phase. Regions which fail this chisquare test and whose
        # sample averages differ by more than mindeviation move on to the refinement phase.
 
        default_opts['min_deviation'] = .25 
        # *mindeviation*: a bound, given as the fraction of the re-
        #  quested error of the entire integral, which determines whether it is worthwhile fur-
        # ther examining a region that failed the chisquare test. Only if the two sampling averages
        # obtained for the region differ by more than this bound is the region further treated.

    
        # *ldxgiven*: the leading dimension of xgiven, i.e. the offset between one
        # point and the next in memory.
        default_opts['ldxgiven'] = None
        default_opts['x_given'] = None
        # *xgiven(ldxgiven,n_given)*: a list of points where the inte-
        # grand might have peaks. Divonne will consider these points when partitioning the
        # integration region. The idea here is to help the integrator find the extrema of the in-
        # tegrand in the presence of very narrow peaks. Even if only the approximate location
        # of such peaks is known, this can considerably speed up convergence.

        default_opts['n_extra'] = 0
        # *nextra*: the maximum number of extra points the peak-finder subrou-
        # tine will return. If nextra is zero, peakfinder is not called and an arbitrary object
        # may be passed in its place, e.g. just 0.

        default_opts['peak_finder'] = None
        # *peakfinder*: the peak-finder subroutine. This subroutine is called
        # whenever a region is up for subdivision and is supposed to point out possible peaks
        # lying in the region, thus acting as the dynamic counterpart of the static list of points
        # supplied in xgiven. It is expected to be declared as
        # def peakfinder(ndim, b, n, x):
        # The bounds of the subregion are passed in the array b, where b(1,d ) is the lower and
        # b(2,d ) the upper bound in dimension d . On entry, n specifies the maximum number
        # of points that may be written to x. On exit, n must contain the actual number of
        # points in x.


        # Parameter relevant for Cuhre integration method
        # ------------------------------------------------

        default_opts['key'] = 0
        # *key* chooses the basic integration rule:
        # key = 7, 9, 11, 13 selects the cubature rule of degree key. Note that the degree-11
        # rule is available only in 3 dimensions, the degree-13 rule only in 2 dimensions.
        # For other values, the default rule is taken, which is the degree-13 rule in 2 dimensions,
        # the degree-11 rule in 3 dimensions, and the degree-9 rule otherwise.

        # Set instance attributes options
        for opt in default_opts:
            if opt in opts:
                setattr(self,opt,opts.pop(opt))
            else:
                setattr(self,opt,default_opts[opt])

        super(pyCubaIntegrator,self).__init__(integrands, **opts)

        supported_algorithms = ['Vegas', 'Suave', 'Divonne', 'Cuhre']
        if self.algorithm not in supported_algorithms:
            raise IntegratorError("The pyCubaIntegrator only supports the following algorithms: [ %s ]"%\
                                                                      (', '.join(supported_algorithms)))

        # After each integrations, results will be available from here
        self.last_integration_results = None

        # Set a stop flag which can be used to tell this python integrator to abort the run
        # if an abort occured.
        self.stopped = Value(c_bool,False)

        # Keep track of the number of iterations performed, it must be threadsafe so that
        # we provide a method for safely setting its value
        self.n_iterations_performed = Value('i',0)
        self.n_integrand_calls = Value('i',0)
        # The threadsafe value above will be forward to the variable below at the end of "integrate"
        self.tot_func_evals = 0

        # Also add a lock for when calling the observables, in order to avoid concurrency issues.
        self.observables_lock = threading.Lock()
        
        # Specify a path for the check_point file used by Cuba when running with the VEGAS
        # algorithm
        self.vegas_check_point_file = pjoin(self.state_file_folder,'cuba_vegas_check_point.dat')
        self.vegas_check_point_file_last_run = pjoin(self.state_file_folder,'cuba_vegas_check_point_last_run.dat')

    def set_stopped(self, value):
        with self.stopped.get_lock():
            self.stopped.value = value            

    def set_n_iterations_performed(self, n_iterations):
        """ Assign the value of the number of iterations performed in a thread-safe way."""
        with self.n_iterations_performed.get_lock():
            self.n_iterations_performed.value = n_iterations

    def increment_n_integrand_calls(self):
        """ Assign the value of the number of iterations performed in a thread-safe way."""
        with self.n_integrand_calls.get_lock():
            self.n_integrand_calls.value += 1

    def observables_normalization(self, n_integrand_calls):
        """ Given the number of integrand calls, return the appropriate overall normalization
        to apply to the observables. In the case of CUBA, only the numbber of iterations matter,
        which is kept track of internally."""
        
        if self.n_iterations_performed.value > 0:
            return 1.0/float(self.n_iterations_performed.value)
        else:
            return 1.0

    @staticmethod
    def wrapped_integrand_vegas(ndim, xx, ncomp, ff, userdata, nvec, core, jacobian_weight, iteration_number):
        """ Function to wrap the integrand to the pyCuba standards when using Vegas which
        provides the jacobian_weight and iteration number"""
        
        # Update the iteration count in a threadsafe manner
        if iteration_number.contents.value > IntegratorInstance.n_iterations_performed.value:
            IntegratorInstance.set_n_iterations_performed(iteration_number.contents.value)

        try:
            return pyCubaIntegrator.wrapped_integrand(
                ndim, xx, ncomp, ff, userdata, jacobian_weight[0])
        except (KeyboardInterrupt, SystemExit):
            IntegratorInstance.set_stopped(True)
            logger.critical('Integration with pyCuba aborted by user.')
            #sys.exit(1)
            return -999

    @staticmethod
    def wrapped_integrand_generic(ndim, xx, ncomp, ff, userdata):
        """ Generic dispatch of the integrand call"""
        try:        
            return pyCubaIntegrator.wrapped_integrand(ndim, xx, ncomp, ff, userdata)
        except (KeyboardInterrupt, SystemExit):
            IntegratorInstance.set_stopped(True)
            logger.critical('Integration with pyCuba aborted by user.')
            #sys.exit(1)
            return -999
    
    @staticmethod
    def wrapped_integrand(ndim, xx, ncomp, ff, userdata, jacobian_weight=1.0):
        """ Function to wrap the integrand to the pyCuba standards.
        IntegratorInstance is taken from a global which is set only during the
        time of the integration with the low-level language integrator implementation."""
        if IntegratorInstance.stopped.value:
            # Signal Cuba to stop
            return -999
        try:
            n_integrands_to_consider = ncomp.contents.value
            IntegratorInstance.increment_n_integrand_calls()
            for integrand_index, integrand in enumerate(IntegratorInstance.integrands[:n_integrands_to_consider]):
                dimensions = integrand.continuous_dimensions
                input_point = np.array([dimensions[i].lower_bound + xx[i]*(dimensions[i].upper_bound - dimensions[i].lower_bound) 
                                                                                         for i in range(ndim.contents.value)])
                try:
                    ff[integrand_index] = IntegratorInstance.integrands[integrand_index](
                        input_point, np.array([], dtype=int), 
                        integrator_jacobian=jacobian_weight, observables_lock=IntegratorInstance.observables_lock)
                except Exception as e:
                    IntegratorInstance.set_stopped(True)
                    logger.critical('The following exception was raised during the evaluation of the integrand:\n%s'%
                                                                    traceback.format_exc())
                    return -999
                # All integration methods here assume unit phase-space volume, so we need to rescale our evaluations here.
                ff[integrand_index] *= IntegratorInstance.integrands[integrand_index].get_dimensions().volume()

            return 0
        except (KeyboardInterrupt, SystemExit):
            IntegratorInstance.set_stopped(True)
            logger.critical('Integration with pyCuba aborted by user.')
            #sys.exit(1)
            return -999

    @staticmethod
    def wrapped_integrand_divonne(ndim, xx, ncomp, ff, phase):
        """ Wrapper for divonne which can specify the argument "phase".
            0, sampling of the points in xgiven,
            1, partitioning phase,
            2, final integration phase,
            3, refinement phase.
        This information might be useful if the integrand takes long to compute and a sufficiently
        accurate approximation of the integrand is available. The actual value of the integral is only
        of minor importance in the partitioning phase, which is instead much more dependent on
        the peak structure of the integrand to find an appropriate tessellation. An approximation
        which reproduces the peak structure while leaving out the fine details might hence be a
        perfectly viable and much faster substitute when phase < 2."""
        try:        
            return pyCubaIntegrator.wrapped_integrand(ndim, xx, ncomp, ff, None)
        except (KeyboardInterrupt, SystemExit):
            IntegratorInstance.set_stopped(True)
            logger.critical('Integration with pyCuba aborted by user.')
            #sys.exit(1)
            return -999

    def get_name(self):
        """ Returns the integrator name."""
        return self.__class__.__name__+'@%s'%self.algorithm

    def aggregate_results(self):
        """ Aggregate results from a MC run into a global summed contribution and error."""
        
        integral = 0.0
        error    = 0.0
        if self.last_integration_results is None:
            logger.critical('Abnormal termination of pyCuba, returning zero for the integral result.')
            return (integral, error)

        for res in self.last_integration_results['results']:
            integral += res['integral']
            error    += res['error']
        
        logger.debug('Cuba completed integration after a total of %d integrand calls (%d refine iterations).'%
                                           (self.n_integrand_calls.value, self.n_iterations_performed.value))
        logger.debug('Final results of the integrals are:\n')
        columns = '%-15s | %-20s | %-10s | %-10s | %-20s'
        logger.debug(columns%('integrand #', 'integral', 'error', 'rel. error', 'prob'))
        logger.debug('-'*84)
        columns = '%-15d | %-20.5e | %-10.2e | %-10.2e | %-20.5g'
        for i_integrand, integrand_result in enumerate(self.last_integration_results['results']):
            logger.debug(columns%(i_integrand+1,res['integral'],res['error'],
                abs(res['error']/res['integral']) if abs(res['integral']) != 0. else 0. ,res['prob']))
        logger.debug('\n')

        return (integral, error)

    def write_vegas_grid(self, source = None, destination = None, n_integrand_calls = None):
        """ Write the Vegas grid while keeping track on the total number of integrand calls
        already performed."""
        
        if any(el is None for el in [source, destination, n_integrand_calls]):
            raise IntegratorError("Function write_vegas_grid of pyCubaIntegrator requires all of its options to be set.")
        
        vegas_grid = { 'binary_grid'         : open(source,'rb').read(),
                       'n_integrand_calls'   : n_integrand_calls,
                       'last_integration_results' : self.last_integration_results }
        
        open(destination,'w').write(str(vegas_grid))

    def load_vegas_grid(self, source = None, destination = None):
        """ Load Vegas grid to destination file and retun the number of integrand_calls that
        lead to the grid. This is useful in particular because CUBA crashes if your start it
        from a grid that already contains more integrand calls than the value 'max_n_eval' set
        by the user."""

        if any(el is None for el in [source, destination]):
            raise IntegratorError("Function load_vegas_grid of pyCubaIntegrator requires all of its options to be set.")

        try:
            vegas_grid = eval(open(source,'r').read())
        except Exception as e:
            raise IntegratorError('pyCuba@Vegas could not read the grid from %s.'%source+
                ' Make sure it was produced by the pyCubaIntegrator and not Cuba directly.')            

        open(destination,'wb').write(vegas_grid['binary_grid'])
        self.last_integration_results = vegas_grid['last_integration_results']

        return vegas_grid['n_integrand_calls']


    def integrate(self):
        """ Return the final integral and error estimators. """

        for integrand in self.integrands:
            if len(integrand.discrete_dimensions)>0:
                raise IntegratorError("pyCuba integrator does not support discrete dimensions for now.")

        if len(set([len(integrand.continuous_dimensions) for integrand in self.integrands]))!=1:
            raise IntegratorError("pyCuba only supports multiple integrands with all the same number of dimensions.")

        # Make sure to reset the number of iterations performed back to 0.
        self.n_iterations_performed.value = 0
        # Likewise for the number of integrand calls
        self.n_integrand_calls.value = 0

        # Set the stopped flag to False
        self.set_stopped(False)
        
        # Make sure we can catch all interrupt signals
        signal.signal(signal.SIGINT, SignalHandler(self.stopped))

        if self.nb_core > 1:
            logger.warning(
"""\n\n-------
Parallelization with pyCuba is handled directly by the Cuba C++ library.
This is potentially dangerous because Python does not know that the wrapped integrands are called
in parallel and it may overwrite memory. In practice, this seems however not to be a problem for a wide
range of applications, especially since MadNkLO placed threadlocks where concurrent writing could occur.
To be safe however, consider making sure that your results are consistent with what is obtained with 
less statistics using one core only.
-------\n""")
            os.environ['CUBACORES'] = ""'%d'""%self.nb_core
            pycuba.set_cuba_nb_core(self.nb_core)
        else:
            os.environ['CUBACORES'] = '0'
            pycuba.set_cuba_nb_core(0)

        # Disable the observables if the integration algorithm is different than Vegas
        if self.algorithm != 'Vegas' and any(integrand.apply_observables for integrand in
                                                                          self.integrands):
            logger.warning('Plotting observables has been disabled when using a CUBA '+
                                              'integration algorithm different than Vegas')

        apply_observables_for_integrands_back_up = []
        if self.algorithm != 'Vegas':
            for i_integrand, integrand in enumerate(self.integrands):
                apply_observables_for_integrands_back_up.append(integrand.apply_observables)
                integrand.apply_observables = False

        n_dimensions = len(self.integrands[0].continuous_dimensions)
        if self.algorithm == 'Vegas':
            with tmpGlobal('IntegratorInstance',self):
                # Explicitly use the 'load_grids' option if you want to load from the check_point grid
                # Add a short_cut for the option 'last_run'
                grid_to_load = self.load_grids
                # If it was not done already make sure that any grid from previous runs is moved
                # to the 'last_run' grid name, hence making sure that it will be reused *only* if
                # it was specified by the user with the option '--last_run'.
                if os.path.exists(self.vegas_check_point_file):
                    shutil.move(self.vegas_check_point_file, self.vegas_check_point_file_last_run)
                if self.load_grids == 'last_run':
                    if os.path.exists(self.vegas_check_point_file_last_run):
                        grid_to_load = self.vegas_check_point_file_last_run
                    else:
                        logger.warning('No pyCuba@Vegas grid found from last run. Starting from scratch.')
                        grid_to_load = None
                else:
                    grid_to_load = self.load_grids
                if grid_to_load:
                    starting_n_integrand_calls = self.load_vegas_grid(
                        source      = grid_to_load,
                        destination = self.vegas_check_point_file )
                    logger.info("Starting pyCuba@Vegas from grid file '%s' (obtained from %d integrand calls)."%
                                                (grid_to_load,starting_n_integrand_calls))
                else:
                    starting_n_integrand_calls = 0
                survey_done = False
                refine_done = False
                apply_observables_for_integrands_back_up = []
                vegas_flags = (
                    2**0 * min(self.verbosity,3) + # Bits 0 and 1 encode the desired verbosity
                    2**4 * 1              + # We want to retain the state file after the run
                    2**5 * 0                # Use the state file if present
                )
                # Check if the survey was asked for
                if self.n_start_survey > 0 and self.max_eval_survey > 0:
                    logger.info('--------------------------------')
                    logger.info('%sNow running cuba@Vegas survey...%s'%(utils.bcolors.GREEN,utils.bcolors.ENDC))
                    logger.info('--------------------------------')
                    # During the survey, it is sub-optimal to fill in histograms, and we therefore disable here
                    # until the refine stage is reached.
                    # We must of course keep a back-up copy of the 'apply_observable' attribute of the
                    # integrand instances so as to be able to reinstate it before the refine.
                    for i_integrand, integrand in enumerate(self.integrands):
                        apply_observables_for_integrands_back_up.append(integrand.apply_observables)
                        integrand.apply_observables = False
                    try:
                        self.last_integration_results = pycuba.Vegas(
                            pyCubaIntegrator.wrapped_integrand_vegas, 
                            n_dimensions,
                            ncomp       = len(self.integrands),
                            flags       = vegas_flags,
                            userdata    = 0,
                            epsrel      = self.target_accuracy_survey,                  
                            nstart      = self.n_start_survey,
                            nincrease   = self.n_increase_survey,
                            nbatch      = self.n_batch,
                            seed        = self.seed,
                            mineval     = starting_n_integrand_calls + self.min_eval,
                            maxeval     = starting_n_integrand_calls + self.max_eval_survey,
                            gridno      = 0,
                            statefile   = self.vegas_check_point_file,
                            nvec        = 1
                        )
                        survey_done = True
                    except (KeyboardInterrupt, SystemExit):
                        logger.critical('Integration with pyCuba aborted by user.')
                        self.set_stopped(True)
                        survey_done = True
                
                # Now save the grids if asked for
                if survey_done and self.save_grids:
                    if os.path.exists(self.vegas_check_point_file):
                        logger.info("Cuba@Vegas grid after survey saved at '%s'"%self.save_grids)
                        self.write_vegas_grid(
                            source            = self.vegas_check_point_file, 
                            destination       = self.save_grids,
                            n_integrand_calls = starting_n_integrand_calls + self.n_integrand_calls.value)
                    else:
                        logger.warning("Cuba@Vegas grids could not be saved because "+
                                  "file '%s' cannot be found."%self.vegas_check_point_file)

                # Revert the apply_observable flag if the survey was performed
                if survey_done:
                    for i_integrand, integrand in enumerate(self.integrands):
                        apply_observables_for_integrands_back_up.append(integrand.apply_observables)
                        integrand.apply_observables = apply_observables_for_integrands_back_up[i_integrand]

                # Make sure to reset the iteration count to 0 and increment the number of
                # integrand calls we are starting from
                starting_n_integrand_calls += self.n_integrand_calls.value
                self.n_iterations_performed.value = 0
                
                # Check if the refine was asked for and if the survey didn't abort
                if self.n_start > 0 and self.max_eval > 0 and not self.stopped.value:
                    # Now run the refine
                    logger.info('--------------------------------')
                    logger.info('%sNow running cuba@Vegas refine...%s'%(utils.bcolors.GREEN,utils.bcolors.ENDC))
                    logger.info('--------------------------------')
                    try:
                        self.last_integration_results = pycuba.Vegas(
                            pyCubaIntegrator.wrapped_integrand_vegas, 
                            n_dimensions,
                            ncomp       = len(self.integrands),
                            flags       = vegas_flags,
                            userdata    = 0,
                            epsrel      = self.target_accuracy,                    
                            nstart      = self.n_start,
                            nincrease   = self.n_increase,
                            nbatch      = self.n_batch,
                            seed        = self.seed,
                            mineval     = starting_n_integrand_calls + self.min_eval,
                            maxeval     = starting_n_integrand_calls + self.max_eval,
                            gridno      = 0,
                            statefile   = self.vegas_check_point_file,
                            nvec        = 1
                        )
                        refine_done = True
                    except (KeyboardInterrupt, SystemExit):
                        logger.critical('Integration with pyCuba aborted by user.')
                        self.set_stopped(True)
                        refine_done = True
    
                # Now save the grids if asked for
                if refine_done and self.save_grids:
                    if os.path.exists(self.vegas_check_point_file):
                        logger.info("Cuba@Vegas grid after refine saved at '%s'"%self.save_grids)
                        self.write_vegas_grid(
                            source            = self.vegas_check_point_file, 
                            destination       = self.save_grids,
                            n_integrand_calls = starting_n_integrand_calls + self.n_integrand_calls.value)
                    else:
                        logger.warning("Cuba@Vegas grids could not be saved because "+
                                  "file '%s' cannot be found."%self.vegas_check_point_file)

                # Move the Vegas check_point file in case the run was interrupted 
                # and the user wants to restart from it
                if os.path.exists(self.vegas_check_point_file):
                    if self.keep_last_run_state_file:
                        logger.info('')
                        logger.info("If you want to run Cuba@Vegas from the last point where it stopped,"+
                            " launch again with the option:\n   %s--load_grids %s%s"%(
                                utils.bcolors.BLUE, self.vegas_check_point_file_last_run, utils.bcolors.ENDC))
                        logger.info('')
                        self.write_vegas_grid(
                                source            = self.vegas_check_point_file, 
                                destination       = self.vegas_check_point_file_last_run,
                                n_integrand_calls = starting_n_integrand_calls + self.n_integrand_calls.value)
                    os.remove(self.vegas_check_point_file)

            result = self.aggregate_results()

        elif self.algorithm == 'Suave':
            with tmpGlobal('IntegratorInstance',self):
                try:
                    self.last_integration_results = pycuba.Suave(
                        pyCubaIntegrator.wrapped_integrand_generic, 
                        n_dimensions,
                        ncomp       = len(self.integrands),
                        verbose     = self.verbosity,
                        userdata    = 0,
                        epsrel      = self.target_accuracy,                    
                        nnew        = self.n_new,
                        nmin        = self.n_min,
                        flatness    = self.flatness,
                        mineval     = self.min_eval,
                        maxeval     = self.max_eval
                    )
                except (KeyboardInterrupt, SystemExit):
                    self.set_stopped(True)
                    logger.critical('Integration with pyCuba aborted by user.')
            result = self.aggregate_results()

        elif self.algorithm == 'Divonne':
            with tmpGlobal('IntegratorInstance',self):
                try:
                    self.last_integration_results = pycuba.Divonne(
                        pyCubaIntegrator.wrapped_integrand_divonne, 
                        n_dimensions,
                        self.key1,
                        self.key2,
                        self.key3,
                        self.max_pass,
                        self.border,
                        self.maxchisq,
                        self.min_deviation,
                        epsrel          = self.target_accuracy,
                        ldxgiven        = self.ldxgiven,
                        xgiven          = self.x_given,
                        nextra          = self.n_extra,
                        peakfinder      = self.peak_finder,
                        ncomp           = len(self.integrands),
                        verbose         = self.verbosity,
                        userdata        = 0,
                        mineval     = self.min_eval,
                        maxeval     = self.max_eval)
                except (KeyboardInterrupt, SystemExit):
                    self.set_stopped(True)
                    logger.critical('Integration with pyCuba aborted by user.')
            result = self.aggregate_results()

        elif self.algorithm == 'Cuhre':
            with tmpGlobal('IntegratorInstance',self):
                try:
                    self.last_integration_results = pycuba.Cuhre(
                        pyCubaIntegrator.wrapped_integrand_generic, 
                        n_dimensions,
                        key         = self.key, 
                        mineval     = self.min_eval,
                        maxeval     = self.max_eval,
                        epsrel      = self.target_accuracy,                
                        ncomp       = len(self.integrands),
                        verbose     = self.verbosity,
                        userdata    = 0)
                except (KeyboardInterrupt, SystemExit):
                    self.set_stopped(True)
                    logger.critical('Integration with pyCuba aborted by user.')

            result = self.aggregate_results()

        # Restore the apply_observables attributes back to their original values
        if self.algorithm != 'Vegas':
            for i_integrand, integrand in enumerate(self.integrands):
                integrand.apply_observables = apply_observables_for_integrands_back_up[i_integrand]
        
        # Update the total of evaluations
        self.tot_func_evals = self.n_integrand_calls.value

        return result
        
if __name__ == '__main__':

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

    # Then the integrators
    # For a typical verbosity, chose 2
    verbosity = 0
    my_integrators = {}

    my_integrators['Naive MC'] = integrators.SimpleMonteCarloIntegrator(my_integrand, 
        n_iterations=10, n_points_per_iterations=1000, accuracy_target=None, verbosity = verbosity)

    my_integrators['VEGAS'] = pyCubaIntegrator(my_integrand, algorithm='Vegas', verbosity = verbosity)

    my_integrators['SUAVE'] = pyCubaIntegrator(my_integrand, algorithm='Suave', verbosity = verbosity)

    my_integrators['DIVONNE'] = pyCubaIntegrator(my_integrand, algorithm='Divonne', verbosity = verbosity)

    my_integrators['CUHRE'] = pyCubaIntegrator(my_integrand, algorithm='Cuhre', verbosity = verbosity)

    chosen_integrators = ['Naive MC', 'VEGAS', 'SUAVE', 'DIVONNE', 'CUHRE']

    for name in chosen_integrators:
        res = my_integrators[name].integrate()
        # Finally integrate
        print('\n'+'='*70)
        print('%-30s = %.4e +/- %.2e'%('%s :: Final result'%name, res[0], res[1]))
        print('='*70)

