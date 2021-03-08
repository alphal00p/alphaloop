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
import sys
import time
import os
import logging
import random
from multiprocessing import Process
from datetime import timedelta
#from aetypes import Enum

try:
    import vegas 
except ImportError:
    print("WARNING :: Vegas 3 could not be found, install it if needed.")

try:
    import numpy as np
except ImportError:
    print("ERROR :: Vegas 3 requires numpy, install it first.")
    sys.exit(1)

try:
    import gvar as gv
except ImportError:
    print("WARNING :: Vegas 3 requires gvar (gaussian variable handling package), install it if needed.")

import madgraph.various.misc as misc
import madgraph.iolibs.files as files
import alpha_loop.integrator.integrators as integrators
from alpha_loop.integrator.integrators import IntegratorError
import alpha_loop.integrator.integrands as integrands
import alpha_loop.integrator.functions as functions
import madgraph.various.cluster as cluster
import madgraph.iolibs.save_load_object as save_load_object

from madgraph import InvalidCmd, MadGraph5Error, MG5DIR
from multiprocessing import Value, Array, Manager

MPI_ACTIVE = False
MPI_SIZE = 1
MPI_RANK = 0 

if MPI_ACTIVE:
    from mpi4py import MPI
    mpi_comm = MPI.COMM_WORLD

numpy = np

class FakeLogger(object):

    def __init__(self):
        pass

    def debug(self,msg):
        print("DEBUG: %s"%msg)

    def info(self,msg):
        print("INFO: %s"%msg)
    
    def warning(self,msg):
        print("WARNING: %s"%msg)

    def critical(self,msg):
        print("CRITICAL: %s"%msg)

logger = logging.getLogger('madgraph.vegas3')
#logger = FakeLogger()

class ParallelWrappedIntegrand(vegas.BatchIntegrand):
    """ Class to wrap the integrand to the Vegas3 standards with parallelization."""
    
    def __init__(self, integrator, integrands, cluster_instance, start_time=None):
        " Save integrand; registers a cluster instance."
        self.start_time = time.time()
        self.integrator = integrator
        self.integrands = integrands
        self.cluster = cluster_instance
        self.start_time = start_time
        if not isinstance(cluster_instance,cluster.MultiCore):
            raise IntegratorError("Vegas3 integrator does not have cluster support "+
                                                                     "yet, only multicore")
        
    def __del__(self):
        " Standard cleanup. "
        pass

    @staticmethod
    def batch_integrand_process(integrands, inputs_list, result):
        """ Function to call the wrapped integrand for a batch of inputs"""
        
        for input_position, inputs in enumerate(inputs_list):
            summed_res = {}
            dims = integrands[0].continuous_dimensions
            fct_inputs = np.array(
                [dims[i].lower_bound + x*(dims[i].upper_bound - dims[i].lower_bound) 
                                                          for i, x in enumerate(inputs)])
            for i_integrand, integrand in enumerate(integrands):
                integrand_volume = integrand.get_dimensions().volume()
                res = integrand(fct_inputs, np.array([], dtype=int))
                if isinstance(res, dict):
                    res = dict( (k,v*integrand_volume) for k,v in res.items())
                else:
                    res = {'I': res*integrand_volume, 'I%d'%i_integrand: res*integrand_volume}
                
                summed_res = dict( (k,res.get(k,0.)+summed_res.get(k,0.)) for k in set(res.keys()).union(set(summed_res.keys())) )

            result[input_position] = summed_res

        return 0

    def batch_integrand(self, inputs_list, result):
        """ Integrand wrapper spawning an independent process for the computation of
        the batch of points listed in inputs_list."""
        p = Process(target=ParallelWrappedIntegrand.batch_integrand_process, args=(
                                                    self.integrands, inputs_list, result))
        p.start()
        p.join()
        return 0
        
    def wait_monitoring(self, Idle, Running, Done):
        if Idle+Running+Done == 0:
            return
        logger.debug('alphaLoop integration: %d Idle, %d Running, %d Done'%(Idle, Running, Done))

    def __call__(self, x):
        " Divide x into self.nproc chunks, feeding one to each process. "

        batch_size = x.shape[0]
        nx = batch_size // self.cluster.nb_core + 1
        inputs  = [x[i*nx : (i+1)*nx] for i in range(self.cluster.nb_core)]
#        print("T+%s: Sending a batch of %dx%d points. Eval so far: %d"%(
#            timedelta(seconds=time.time()-self.start_time),
#            self.cluster.nb_core, nx, self.integrator.n_function_evals ))
        # The results enttried will be accessed from independent processes and must
        # therefore be defined as shared memory
        
        #results = [Array('d', [0.]*len(input)) for input in inputs]

        local_manager = Manager()
        all_results = [ local_manager.list([{}]*len(input)) for input in inputs ]

        # launch evaluation of self.integrands for each chunk, in parallel
        if not self.start_time is None:
            logger.debug('Dispatching batch of %d points on %d cores: [%s]'%(
            batch_size,self.cluster.nb_core, misc.format_time(time.time()-self.start_time)))
        else:
            logger.debug('Dispatching batch of %d points on %d cores:'%(
                                                          batch_size,self.cluster.nb_core))            
        for position, input in enumerate(inputs):
            if len(input)==0: continue
            time.sleep(0.1)
            #self.cluster.submit(self.batch_integrand, [input, results[position]])            
            self.cluster.submit(self.batch_integrand, [input, all_results[position]])

        # Wait for all jobs to finish.
#        self.cluster.wait('', self.wait_monitoring,update_first=self.wait_monitoring)
        self.cluster.wait('', self.wait_monitoring)

        self.integrator.n_function_evals += batch_size
        #return np.concatenate([ np.array(result[:]) for result in results ])
        final_result = { k: 
            np.concatenate([ np.array([res[k] for res in results]) for results in all_results ])
            for k in all_results[0][0].keys()
        }
        return final_result

class Vegas3Integrator(integrators.VirtualIntegrator):

    def __init__(self, integrands, **opts):
        """ Initialize the Vegas3 integrator."""
       
        # General steering parameters
        default_opts= {
           'verbosity'          :   0,
           'target_accuracy'    :   1.0e-3,
           'seed'               :   None,
           'save_grids'         :   None,
           'load_grids'         :   None,
        }

        # Parameter for controlling VEGAS3
        # ------------------------------------------------
        
        # the number of iterations for the training session. 
        default_opts['n_iterations_survey'] = 10
        # the number of points per iteration in the training session. 
        default_opts['n_points_survey'] = 2000

        # the number of iterations for the production session. 
        default_opts['n_iterations_refine'] = 10
        # the number of points per iteration in the production session. 
        default_opts['n_points_refine'] = 10000

        # Steering parallelization
        if not MPI_ACTIVE:
            default_opts['cluster'] = cluster.onecore
        else:
            default_opts['cluster'] = cluster.MPICluster(MPI_RANK,MPI_SIZE)

        # Default batch size per node
        default_opts['batch_size'] = 1000

        # Set instance attributes options
        for opt in default_opts:
            if opt in opts:
                setattr(self,opt,opts.pop(opt))
            else:
                setattr(self,opt,default_opts[opt])
                
        super(Vegas3Integrator,self).__init__(integrands, **opts)

        if self.seed is not None:
            random.seed(self.seed)
            np.random.seed(self.seed)

        # Adjust batch size proportionally to the number of nodes
        self.batch_size *= self.cluster.nb_core

        # Save here the main instance for vegas3
        self.vegas3_integrator = None 
        
        if self.verbosity > 0:
            logger.level = logging.DEBUG
        else:
            logger.level = logging.INFO

        self.full_result = None

    def observables_normalization(self, n_integrand_calls):
        """ Given the number of integrand calls, return the appropriate overall normalization
        to apply to the observables."""
        return 1.0/float(self.n_iterations_refine)

    def wrapped_integrand(self, x_inputs):
        """ Function to wrap the integrand to the Vegas3 standards."""

        dims = self.integrands[0].continuous_dimensions
        if len(list(x_inputs))>len(dims):
            jacobian = x_inputs[-1]*dims.volume()
            x_inputs = list(x_inputs)[:-1]
        else:
            x_inputs = list(x_inputs)
            jacobian = 1.0*dims.volume()

        fct_inputs = np.array([dims[i].lower_bound + x*(dims[i].upper_bound - dims[i].lower_bound) 
                                                                         for i, x in enumerate(x_inputs)])
        
        all_results = {}
        found_non_dict = False
        summed_res = 0.
        for i_integrand, integrand in enumerate(self.integrands):
            integrand_volume = integrand.get_dimensions().volume()
            try:
                res = integrand(fct_inputs, np.array([], dtype=int), integrator_jacobian=jacobian)
            except AssertionError as err:
                traceback.print_tb(sys.exc_info()[-1])
                logger.warning('Assertion error encountered.')
                res = 0.
            if isinstance(res, dict):
                for k, v in res.items():
                    if k in all_results:
                        all_results[k] += v*integrand_volume
                    else:
                        all_results[k] = v*integrand_volume
            else:
                found_non_dict = True
                summed_res += res*integrand_volume
                all_results['I%d'%i_integrand] = res*integrand_volume
        if found_non_dict:
            all_results['I'] = summed_res
            
        self.n_function_evals +=1
        if self.n_function_evals%(max(self.curr_n_evals_per_iterations/2,1))==0 and MPI_RANK == 0:
            logger.debug('Evaluation #%d / %d*%d  (%.3g%%) [%s]'%(
                self.n_function_evals*MPI_SIZE, self.curr_n_iterations, self.curr_n_evals_per_iterations, 
                (100.0*self.n_function_evals*MPI_SIZE / (self.curr_n_iterations*self.curr_n_evals_per_iterations)),
                                            misc.format_time(time.time()-self.start_time)))
        
        return all_results #dict( ('I%d'%i, res) for i, res in enumerate(all_results) )

    def get_name(self):
        """ Returns the integrator name."""
        return 'Vegas3' 

    def integrate(self):
        """ Return the final integral and error estimators. """

        start_time = time.time()

        for integrand in self.integrands:
            if len(integrand.discrete_dimensions)>0:
                raise IntegratorError("Vegas3 integrator does not support discrete dimensions for now.")

        if len(set([len(integrand.continuous_dimensions) for integrand in self.integrands]))!=1:
            raise IntegratorError("Vegas3 only supports multiple integrands with all the same number of dimensions.")

        # During the survey, it is sub-optimal to fill in histograms, and we therefore disable here
        # until the refine stage is reached.
        # We must of course keep a back-up copy of the 'apply_observable' attribute of the
        # integrand instances so as to be able to reinstate it before the refine.
        if MPI_RANK == 0:
            apply_observables_for_integrands_back_up = []
            for i_integrand, integrand in enumerate(self.integrands):
                apply_observables_for_integrands_back_up.append(integrand.apply_observables)
                integrand.apply_observables = False

        n_dimensions = len(self.integrands[0].continuous_dimensions)
        # sync_ran is to decide if VEGAS3 random number generator should produce the same
        # number on different processors.
        if any(apply_observables_for_integrands_back_up) and self.cluster.nb_core==1:
#            if self.n_iterations_refine > 1:
#                logger.warning("Vegas3 can only run a single refine iteration when a fixed-order analysis is active.\n"+
#                               "The parameter 'n_iterations_refine' will consequently be forced to take the value 1.")
#                self.n_iterations_refine = 1
            self.vegas3_integrator = VegasWithJacobianInFunctionInput(n_dimensions * [[0., 1.]],
                analyzer        = vegas.reporter() if self.verbosity>1 else None, 
                nhcube_batch    = self.batch_size,
                # We must disable the hypercube optimisation for now when plotting observables
                max_nhcube      = 1,
                beta            = 0.,
                alpha           = 0.5,
                sync_ran        = True)
        else:
            self.vegas3_integrator = vegas.Integrator(n_dimensions * [[0., 1.]],
                analyzer        = vegas.reporter() if self.verbosity>1 else None, 
                nhcube_batch    = self.batch_size,
                # We must disable the hypercube optimisation for now when plotting observables
                max_nhcube      = 1e9,
                beta            = 0.75,
                alpha           = 0.5,
                sync_ran        = True)
        
        self.start_time = time.time()
        # Access/generate the wrapped integrand
        if False or self.cluster.nb_core==1:
            wrapped_integrand = self.wrapped_integrand
        else:
            wrapped_integrand = ParallelWrappedIntegrand(
                                      self, self.integrands, self.cluster, self.start_time)          

        self.tot_func_evals = 0
        result = None
        
        # Load pre-exisitng grids if specified. Note that no information on Vegas' hypercubes is supplied
        # here, only the standard integration grids along each dimension.
        if self.load_grids is not None:
            try:
                vegas_map = save_load_object.load_from_file(self.load_grids)
            except Exception as e:
                raise IntegratorError("Vegas3 could not load its integration grid from '%s'. Error: %s"%(
                                                                    self.load_grids,str(e)))
            logger.info('Vegas3 reading integration grids from file: %s'%self.load_grids)
            self.vegas3_integrator.set(map=vegas_map)
        
        if self.n_iterations_survey > 0 and self.n_points_survey > 0:
            # Train grid
            if MPI_RANK == 0:
                logger.debug("=================================")
                logger.debug("Vegas3 starting the survey stage.")
                logger.debug("=================================")
            self.n_function_evals = 0
            self.curr_n_iterations = self.n_iterations_survey
            self.curr_n_evals_per_iterations = self.n_points_survey
            result = self.vegas3_integrator(wrapped_integrand, 
                        nitn=self.n_iterations_survey, neval=self.n_points_survey, adapt=True)
    
            if MPI_RANK == 0:
                logger.debug('\n'+result.summary())
                
            self.tot_func_evals += self.n_function_evals
        else:
            logger.info("The survey step of Vegas3 integration is being skipped as per user's request.")
           
        # Save generated grids if asked for. Notice that no information on Vegas' hypercubes is saved here,
        # but just the standard integration grids.
        if self.save_grids is not None:
            try:
                vegas_map = save_load_object.save_to_file(self.save_grids, self.vegas3_integrator.map) 
            except Exception as e:
                raise IntegratorError("Vegas3 could not save its integration grid to '%s'. Error: %s"%(
                                                                    self.load_grids,str(e)))
            logger.info('Vegas3 saved its integration grids to file: %s'%self.save_grids)
           
        if self.n_iterations_refine > 0 and self.n_points_refine > 0:
            # Final integration
            if MPI_RANK == 0:
                logger.debug("=============================================")
                logger.debug("Vegas3 starting the refined integration stage")
                logger.debug("=============================================")
    
                # Restore the filling of the histograms for the refine, only if not running in parallel
                if self.cluster.nb_core==1:
                    for i_integrand, integrand in enumerate(self.integrands):
                        integrand.apply_observables=apply_observables_for_integrands_back_up[i_integrand]
                else:
                    if any(apply_observables_for_integrands_back_up):
                        logger.warning('Filling of the histograms of the fixed-order analysis will'+
                            ' now be disabled as this functionality is not yet available for parallel integration.')
    
            self.n_function_evals = 0
            self.curr_n_iterations = self.n_iterations_refine
            self.curr_n_evals_per_iterations = self.n_points_refine
            result = self.vegas3_integrator(wrapped_integrand, 
                        nitn=self.n_iterations_refine, neval=self.n_points_refine, adapt=False) 
            if MPI_RANK == 0:
                logger.debug('\n'+result.summary())
        else:
            logger.info("The refine step of Vegas3 integration is being skipped as per user's request.")

        if MPI_RANK == 0:
            # Make sure to restore original settings for the 'apply_observable' attribute
            for i_integrand, integrand in enumerate(self.integrands):
                integrand.apply_observables=apply_observables_for_integrands_back_up[i_integrand]
       
        self.tot_func_evals += self.n_function_evals
        if MPI_ACTIVE:
            # Aggregate the number of function evaluation from each MPI rank
            if MPI_RANK == 0:
                for rank in range(1,MPI_SIZE):
                    n_points_in_slave = mpi_comm.recv(source=rank,tag=11)
                    self.tot_func_evals += n_points_in_slave
            else:
                mpi_comm.send(self.n_function_evals, dest=0, tag=11)

        # Result is a list of instances of vegas.RAvg variablesr, with the following attributes: 
        # RAvg.mean : mean 
        # RAvg.sdev : standard dev.
        # RAvg.chi2 : Chi^2 of the weighted average
        # RAvg.dof  : number of degreeas of freedom
        # RAvg.Q    : p-value of the weighted average 
        # RAvg.itn_results : list of the integral estimates for each iteration
        # RAvg.summary() : summary of the integration
#        if False or self.cluster.nb_core==1:
#            summed_result = sum(result.values())
#        else:
#            summed_result = result
        
        self.full_result = result
        if 'I' not in result.keys():
            raise IntegratorError("The integrand did not meet the requirement of specifying the mandatory dictionary entry 'I'.")
        summed_result = result['I']

        integration_time = time.time() - start_time

        if MPI_RANK == 0:
            logger.debug("===============================================================")
            logger.debug('Total integration time: %s'%misc.format_time(integration_time))
            logger.debug('Vegas3 used a total of %d function evaluations, on %d core(s).'%
                                                  (self.tot_func_evals,self.cluster.nb_core))
            logger.debug('Vegas3 returned final results : %s'%summed_result)
            logger.debug("===============================================================")
        
        #out_data = open('0result.dat','a')
        #out_data.write('%.16e, '%(summed_result.mean)+'%.16e, '%(summed_result.sdev)+'%i\n'%(self.integrands[0].counter))
        #out_data.close()
        ##misc.sprint('\n%.16e, '%(summed_result.mean)+'%.16e, '%(summed_result.sdev)+'%i\n'%(self.integrands[0].counter))
        self.integrands[0].counter = 0

        return summed_result.mean, summed_result.sdev

    def show_grid(self, n_grid=40, shrink=False, axes=None):
        """ Shows the integration grid of Vegas3 in a neat matplotlib window.
        n_grid is the number of grid nodes in each direction to plot, 
        shrink is whether to show the entire range or just the n_grid nodes and
        axes can be a list of tuples specifying the list of directions to show."""

        if not self.vegas3_integrator:
            raise InvalidCmd("Vegas3 can only show grids after one integration has been performed.") 
        if not axes and not shrink:
            self.vegas3_integrator.map.show_grid(n_grid)
        else:
            raise IntegratorError("Options for show_grid not correctly handled yet.")
            self.vegas3_integrator.map.show_grid(n_grid, shrink, axes=axes)


class Channel():
    """ class that specifies a channel and its parameters when using the multi-channel integration technique."""
    def __init__(self,i,integrand,adaptive_wgt,selector_wgt):
        self.nr = i
        self.full_integrand = integrand
        self.adaptive_wgt = adaptive_wgt
        self.selector_wgt = selector_wgt
        self.dimension = len(integrand.continuous_dimensions)
        # TODO: Hypercubes of vegas not used here (see add_point)
        self.integrator = vegas.Integrator(self.dimension*[[0.,1.]])
        self.reset()
        return

    def integrand(self,random_variable,adaptive_wgts):
        """ Current integrand of channel. Adaptive_wgts change over course of integration. """
        res = self.full_integrand(random_variable, np.array([], dtype=int),
                                  adaptive_wgts=adaptive_wgts,
                                  channel_nr=self.nr)
        return res
    
    def update(self):
        """ Compute estimates of integral, integral^2, sdev, rel variance. """
        self.results_squared = [res**2 for res in self.results]
        assert(self.n == len(self.results))
        self.n = len(self.results)
        self.result = np.mean(self.results)
        self.result_squared = np.mean(self.results_squared)
        self.sdev = np.sqrt((self.result_squared-self.result**2)/self.n)
        self.rel_variance = (self.result_squared-self.result**2)/self.result**2
        return
    
    def reset(self):
        """ (Re)sets all data gathered while integrating to zero/empty. """
        self.results = []          # each evaluation of integrand*wgt
        self.random_variables = [] # = y in vegas manual
        self.grid_variables = []   # = x in vegas manual
        self.results_squared = []  # each evaluation of (integrand*wgt)^2
        self.wgts = []             # = wgt in vegas manual
        self.result = 0.           # current estimate of integrand (calculated in update)
        self.result_squared = 0.   # current estimate of integrand^2 (calculated in update)
        self.sdev = 0.             # current estimate of sdev (calculated in update)
        self.n = 0                 # number of total evaluation (==len(results))
        self.rel_variance = 0.     # current estimate of variance/result^2 (calculated in update)
        self.n_passed_cuts = 0     # number of evaluations that passed cuts
        return
        
    def add_point(self,adaptive_wgts):
        """ Adds an integration point sampled according to vegas3 grid and evaluates the integrand. """
        y = np.random.uniform(0., 1., (1, self.dimension))
        x = np.empty(y.shape, float)
        wgt = np.empty(y.shape[0], float)
        self.integrator.map.map(y,x,wgt)
        y = y[0]
        x= x[0]
        wgt = wgt[0]
        #misc.sprint('\n%s \n%s \n%e'%(x,y,wgt))
        self.random_variables += [y]
        self.grid_variables += [x]
        self.wgts += [wgt]
        try:
            res = self.integrand(x,adaptive_wgts)*wgt
        except AssertionError as err:
            traceback.print_tb(sys.exc_info()[-1])
            misc.sprint('Assertion error encountered during integration. Ignoring this point.')
            res = 0.
        self.results += [res]
        if res != 0.:
            self.n_passed_cuts += 1
        self.n += 1
        return
        
    def adjust_to_new_grid(self):
        # currently not used!
        """Integrand values can be reused during vegas3 grid updates.
        For a given x, there's an old y & old wgt from the old grid that is known
        and a new y & new wgt from the new grid that is unknown, yet.
        The xs are invariant by construction.
        Here the new y & new wgt are found and the integrand updated as integrand*=new_wgt/old_wgt."""
        grid = self.integrator.map.grid
        # collect new ys from old xs
        for j,x in enumerate(self.grid_variables):
            y_new = np.empty_like(x)
            for d in range(len(grid[:,0])):
                for i in range(len(grid[0,:])-1):
                    if grid[d,i]<=x[d]<grid[d,i+1]:
                        delta = (x[d]-grid[d,i])/(grid[d,i+1]-grid[d,i])
                        y_new[d] = (delta + i)*1./self.integrator.map.ninc#[d]
            # overwrite old ys with new ones
            #misc.sprint('\n',self.random_variables[j],'\n',y_new)
            self.random_variables[j]=y_new
        # check with built in function if I inverted correctly:
        check_x = np.empty_like(self.random_variables)
        new_wgt = np.empty(check_x.shape[0],np.float_)
        self.integrator.map.map(np.array(self.random_variables),check_x,new_wgt)
        #misc.sprint(check_x[0],self.grid_variables[0],self.integrator.map(np.array(self.random_variables[0])))
        assert(np.all(x==x_alt for x,x_alt in zip(self.grid_variables,check_x)))
        # adjust all saved results with the new wgts, so that they are reusable
        for i in range(len(self.results)):
            self.results[i] *= new_wgt[i]/self.wgts[i]
            #misc.sprint(new_wgt[i]/self.wgts[i])
            self.wgts[i] = new_wgt[i]
        return
    
    def header_string(self):
        msg = '{:<7}  {:>14}  {:>12}  {:>7}  {:<22}'.format('channel','#passed /  tot',
                                                              'adaptive wgt','rel var','estimate')
        return msg
    
    def summary_string(self):
        """ Returns a string with the current state of the integration. Formatted according to header_string. """
        msg = '{:<7}  {:>14}  {:>12}  {:>7}  {:>28}'.format(self.nr,'%4d / %4d'%(self.n_passed_cuts,self.n),
                                                               '%.10f'%self.adaptive_wgt,
                                                               '%3.2f'%self.rel_variance,
                                                               '%.6e +- %.6e'%(self.result,self.sdev))
        return msg

    def get_n_points(self,abs_acc):
        """ estimates the number of Monte Carlo points needed to arrive at an accuracy abs_acc of THIS CHANNEL. """
        if self.sdev < abs_acc:
            n_points = 0
        else:
            #var = self.sdev**2*self.n
            #n_improve = var/abs_target_accuracy**2 - self.n
            n_points = int((self.sdev**2/abs_acc**2-1.)*self.n)
        return n_points
    
class MultiChannelVegas3Integrator(integrators.VirtualIntegrator):
    """ Mutli-Channel integrator: Class that handles channels specified as instances of class Channel."""
    
    def __init__(self, integrands, adaptive_wgts=None, selector_wgts=None, **opts):
        """ specifies:
                adaptive_wgts: channel weights (usually denoted by \alpha), sum up to one
                selector_wgts: channel more often sampled if its selector wgt is large, sum up to one """

        super(MultiChannelVegas3Integrator,self).__init__(integrands, **opts)
        
        # TODO: #flattening_strategy = {
        #       channel_id: {'flattening_technique': 'diagram_fuction to call'/'jacobians',
        #       'use_alphas': True/False}}
        
        #np.random.seed(0)
        
        self.start_time = time.time()
        self.time = time.time()
        
        self.neval = opts.pop('n_points') #number of evaluations of first iteration, doubled after every itn
        self.rel_ch_acc = 0.1 #relative target accuracy of each channel
        self.v_acc = 0.01 #relative change of vegas grid
        self.a_acc = 0.01 #relative change of channel weights (alpha)
        self.min_nitn = 3
        self.max_nitn = opts.pop('n_iterations')
        assert(self.min_nitn<=self.max_nitn)
        self.rel_acc = 0.01 #relative target accuracy of final integral
        
        # channel properties
        if len(integrands)>1:
            misc.sprint('This multi channel integrator works with only one process/integrand for now.')
            raise
        self.n_channels = len(self.integrands[0].topologies_to_processes)
        if adaptive_wgts is None: #uniform
            adaptive_wgts = [1./self.n_channels]*self.n_channels
        if selector_wgts is None: #uniform
            selector_wgts = [1./self.n_channels]*self.n_channels
        self.adaptive_wgts = adaptive_wgts
        self.selector_wgts = selector_wgts
        
        
        # initialize channels
        integrand = self.integrands[0]
        self.channels = []
        for i in range(self.n_channels):
            self.channels += [Channel(i,integrand,self.adaptive_wgts[i],self.selector_wgts[i])]
        
        # variables for the full integrand   
        self.result = 0.
        self.sdev = 0.
        self.n = 0
        self.n_passed_cuts = 0
        self.rel_variance = 0.
        return
    
    def reset(self):
        """ Deletes all integration data, except for grid and channel weights (alpha) """
        for channel in self.channels:
            channel.reset()
        return
    
    def update(self):
        """ Combine the estimates from all the channels to give an estimate of the full integral. """
        for channel in self.channels:
            channel.update()
        self.result = sum([channel.result*channel.adaptive_wgt for channel in self.channels])
        self.n = sum([channel.n for channel in self.channels])
        self.n_passed_cuts = sum([channel.n_passed_cuts for channel in self.channels])
        W = sum([channel.result_squared*channel.adaptive_wgt for channel in self.channels])
        self.sdev = np.sqrt((W-self.result**2)/self.n)
        self.rel_variance = (W-self.result**2)/self.result**2 
        return
    
    def adapt_vegas(self,alpha):
        """ Adapts vegas grid of each channel with saved function evaluations
        and returns a measure for the change of the grid. """
        ch_wgt_change = []
        for channel in self.channels:
            old_wgts = np.array(channel.integrator.map.inc)
            channel.integrator.map.add_training_data(np.array(channel.random_variables),np.array(channel.results))
            channel.integrator.map.adapt(alpha=alpha)
            new_wgts = np.array(channel.integrator.map.inc)
            wgt_change = []
            wgt_change = [abs(new-old)/old for new,old in zip(new_wgts,old_wgts)]
            wgt_change = np.mean(wgt_change,axis=1) # sum within dimension (already normalized wgts)
            wgt_change = np.mean(wgt_change) #mean along each PS dimension
            ch_wgt_change += [wgt_change]
        return np.array(ch_wgt_change)
    
    
    def get_channels_with_points(self,n_points=None,rel_target_accuracy=None,training=-1):
        if n_points == None and rel_target_accuracy != None:
            self.update()
            if training == True:
                abs_target_accuracies = [rel_target_accuracy*channel.result for channel in self.channels]
                n_points = [channel.get_n_points(abs_target_accuracies[channel.nr]) for channel in self.channels]
                #if n_points > self.max_refine_neval:
                #    n_points = self.max_refine_neval
            elif training == False: 
                abs_target_accuracy = rel_target_accuracy*self.result
                n_points = self.get_n_points(abs_target_accuracy)
                n_points = [int(n_points*channel.selector_wgt) for channel in self.channels]
                abs_target_accuracy_alt = rel_target_accuracy*self.result/np.sqrt(self.n_channels)
                n_points_alt = [channel.get_n_points(abs_target_accuracy_alt) for channel in self.channels]
                if np.mean(n_points)>np.mean(n_points_alt):
                    n_points = n_points_alt
            else:
                logger.info('Set training either True or False.')
                raise
        elif n_points != None and rel_target_accuracy == None:
            n_points = [int(n_points*channel.selector_wgt) for channel in self.channels]
        else:
            self.selector_wgts = self.adaptive_wgts
            for channel in self.channels:
                channel.selector_wgt = channel.adaptive_wgt
            logger.info('Either choose n_points or rel_target_accuracy.')
            raise
        channels_with_points = {channel.nr: n_points[channel.nr] for channel in self.channels}
        return channels_with_points
    
        
    def get_n_points(self,abs_target_accuracy):
        """ estimates the number of Monte Carlo points needed to arrive at an accuracy abs_target_accuracy of FULL INTEGRAL. """
        if self.sdev < abs_target_accuracy:
            n_points = 0
        else:
            #var = self.sdev**2*self.n
            #n_improve = var/abs_target_accuracy**2 - self.n
            n_points = int((self.sdev**2/abs_target_accuracy**2-1.)*self.n)
        return n_points
    
    
    def adapt_wgts(self):
        """ One adaptation of adaptive_wgts. Per iteration there are n_vegas_adapt adaptations of the vegas grid.
            There are at least adapt_min_n_points points evaluations. If rel_target_accuracy is not reached,
            increase max n_max_improve_itns of n_improve_steps."""
        #now update the adaptive_wgts and normalize
        old_wgts = self.adaptive_wgts
        #misc.sprint(old_adaptive_wgts)
        all_adaptive_wgts = []
        for channel in self.channels:
            channel.adaptive_wgt *= np.sqrt(channel.result_squared)#**alpha
            all_adaptive_wgts += [channel.adaptive_wgt]
        adaptive_norm = sum(all_adaptive_wgts)
        for channel in self.channels:
            channel.adaptive_wgt *= 1./adaptive_norm
            all_adaptive_wgts[channel.nr] = channel.adaptive_wgt
        self.adaptive_wgts = all_adaptive_wgts
        new_wgts = self.adaptive_wgts
        wgt_change = [abs(new-old)/old for new,old in zip(new_wgts,old_wgts)]
        wgt_change = np.mean(wgt_change) #(already normalized wgts)
        return (wgt_change)
    
    def summary_string(self):
        """ Returns a string with the current state of the integration of all channels and their sum. """
        msg = '{:<7}  {:>14}  {:>12}  {:>7}  {:>28}'.format('full','%4d / %4d'%(self.n_passed_cuts,self.n),'-',
                                                             '%3.2f'%self.rel_variance,
                                                             '%.6e +- %.6e'%(self.result,self.sdev))
        return msg
    
    def print_summary(self):
        for channel in self.channels:  
            channel.update()    
            logger.info(channel.summary_string())
        logger.info('-'*70)
        self.update()
        logger.info(self.summary_string())
        return
        
    def get_channel_points(self,n_points):
        """ distributes a total of n_points on channels according to selector_wgts.
        returns a dict: {channel_nr: n_pts} """
        n_points = [int(n_points*channel.selector_wgt) for channel in self.channels]
        pts_per_channels = {channel.nr: n_points[channel.nr] for channel in self.channels}
        return pts_per_channels
    
    def get_refine_channel_points(self,rel_acc):
        """ estimates number of points needed to reach rel_acc. returns a dict: {channel_nr: n_pts} """
        abs_acc = rel_acc*self.result/np.sqrt(self.n_channels)
        abs_ch_acc = [abs_acc/channel.adaptive_wgt for channel in self.channels]
        n_points = [channel.get_n_points(abs_ch_acc[channel.nr]) for channel in self.channels]
        refine_pts_per_channels = {channel.nr: n_points[channel.nr] for channel in self.channels}
        return refine_pts_per_channels
    
    def survey(self, neval, max_nitn, min_nitn=5, a_acc=0.01, v_acc=0.01):
        """ survey stage of integration: girds and channel weights (alpha) are trained.
        while both the vegas grid and the channel weights have not stabilized (reach a_acc,v_acc),
        do an iteration with double the number of points neval.
        When adding points in parallel, print a full estimate after every "update_time" seconds.
        run a minimum and a maximum of iterations. """
        msg = '='*27 + ' training stage '+'='*27
        logger.info(msg,'$MG:color:GREEN')
        logger.info(self.channels[0].header_string())
        
        update_time = 10
        tot_it = 0
        a_it = 0
        v_it = 0
        while not (a_acc == -1 and v_acc == -1):
            self.reset()
            channel_points = self.get_channel_points(neval*2**(tot_it)) #double every time
            t0 = time.time()
            while sum(channel_points.values()) > 0:
                for nr in list(channel_points.keys()):
                    if channel_points[nr] == 0:
                        channel_points.pop(nr)
                    else:
                        self.channels[nr].add_point(self.adaptive_wgts)
                        channel_points[nr] -= 1
                if (time.time()-t0) > update_time:
                    total_time = time.time()-self.start_time
                    minutes = int(total_time/60)
                    time_string = '%d min %.2f sec'%(minutes,total_time-minutes*60)
                    self.update()
                    logger.info(self.summary_string()+'    '+'\033[94m'+time_string+'\033[0m')
                    t0 = time.time()
            self.print_summary()
            # if max is reached stop
            if tot_it == max_nitn:
                v_acc = -1
                a_acc = -1
            elif tot_it > 0 and tot_it >= min_nitn:
                if v_acc == -1 and a_acc != -1:
                    if a_change < a_acc:
                        a_acc = -1
                elif v_acc != -1 and a_acc == -1:
                    if mean_v_change < v_acc:
                        v_acc = -1
                else:
                    if mean_v_change < v_acc and a_change < a_acc:
                        v_acc = -1
                        a_acc = -1
            if v_acc != -1:
                v_change = self.adapt_vegas(alpha=.5)
                v_it += 1
                mean_v_change = np.mean(v_change)
                logger.info('#%d vegas grid changed by %.2f%%'%(v_it,mean_v_change*1e2),'$MG:color:GREEN')
                #logger.info('#%d vegas adapted: %s'%(v_it,' '.join(['(Ch. %d: %.4e)'%(nr,change) for nr,change in zip(xrange(self.n_channels),v_change)])),'$MG:color:GREEN')
            if a_acc != -1:
                a_change = self.adapt_wgts()
                a_it += 1
                logger.info('#%d: adaptive weights changed by %.2f%%'%(a_it,a_change*1e2),'$MG:color:GREEN') 
            tot_it += 1
        return
    
    def refine(self,rel_acc):
        """ refine stage of integration.
        estimates number of points needed to reach rel_acc.
        after every "update_time" seconds, prints out current estimate of full integrand
        and restimates number of points needed to reach rel_acc. """
        
        msg = '='*26 + ' refinement stage '+'='*26
        logger.info(msg,'$MG:color:GREEN')
        logger.info(self.channels[0].header_string())
        
        update_time = 10
        refine_channel_points = self.get_refine_channel_points(rel_acc=rel_acc)
        for nr,pts in list(refine_channel_points.items()):
            if pts != 0:
                logger.info('ESTIMATE: Refining channel %d with %d additional points.'%(nr,pts))
        t0 = time.time()
        while sum(refine_channel_points.values()) > 0:
            for nr in list(refine_channel_points.keys()):
                if refine_channel_points[nr] == 0:
                    refine_channel_points.pop(nr)
                else:
                    self.channels[nr].add_point(self.adaptive_wgts)
                    refine_channel_points[nr] -= 1
            if (time.time()-t0) > update_time:
                total_time = time.time()-self.start_time
                minutes = int(total_time/60)
                time_string = '%d min %.2f sec'%(int(total_time/60),total_time-minutes*60)
                self.update()
                logger.info(self.summary_string()+'    '+'\033[94m'+time_string+'\033[0m')
                t0 = time.time()
                refine_channel_points = self.get_refine_channel_points(rel_acc=self.rel_acc)
                for nr,pts in list(refine_channel_points.items()):
                    if pts != 0:
                        logger.info('ESTIMATE: Refining channel %d with %d additional points.'%(nr,pts))
        self.print_summary()
        return
    
    def integrate(self):
        """ combines survey and refine stage. and returns final result and Monte Carlo error. """
        
        self.survey(neval=self.neval,
                   max_nitn=self.max_nitn, min_nitn=self.min_nitn,
                   a_acc=self.a_acc, v_acc=self.v_acc)#,rel_ch_acc=self.rel_ch_acc)
        
        if self.result == 0:
            raise MadGraph5Error('Do training first.','$MG:color:BLUE')

        self.refine(rel_acc=self.rel_acc)
        
        misc.sprint('\n%.16e, '%(self.result)+'%.16e, '%(self.sdev)+'%i\n'%(self.n_passed_cuts))
        return self.result, self.sdev
        
        
# Wrapper around the original Vegas.Integrator which adds the jacobian of the point as
# the last input of the 1-dimensional array passed to the integrand function in the 
# __call__ function of the integrator.
# WARNING: does not work with vegas.Integrator.beta > 0 because Cython somehow makes the
# float attribute 'sum_sigf' of the integrand read-only, making the line:
#     self.sum_sigf = sum_sigf
# crash.
class VegasWithJacobianInFunctionInput(vegas.Integrator):

    def __init__(self, *args, **opts):
        super(VegasWithJacobianInFunctionInput, self).__init__(*args, **opts)

    def __call__(self, fcn, **kargs):
        """ Integrate integrand ``fcn``.

        A typical integrand has the form, for example::

            def f(x):
                return x[0] ** 2 + x[1] ** 4

        The argument ``x[d]`` is an integration point, where
        index ``d=0...`` represents direction within the
        integration volume.

        Integrands can be array-valued, representing multiple
        integrands: e.g., ::

            def f(x):
                return [x[0] ** 2, x[0] / x[1]]

        The return arrays can have any shape. Dictionary-valued
        integrands are also supported: e.g., ::

            def f(x):
                return {'a':x[0] ** 2, 'b':[x[0] / x[1], x[1] / x[0]]}


        Integrand functions that return arrays or dictionaries
        are useful for multiple integrands that are closely related,
        and can lead to substantial reductions in the errors for
        ratios or differences of the results.

        It is usually much faster to use |vegas| in batch
        mode, where integration points are presented to the
        integrand in batches. A simple batch integrand might
        be, for example::

            @vegas.batchintegrand
            def f(x):
                return x[:, 0] ** 2 + x[:, 1] ** 4

        where decorator ``@vegas.batchintegrand`` tells
        |vegas| that the integrand processes integration
        points in batches. The array ``x[i, d]``
        represents a collection of different integration
        points labeled by ``i=0...``. (The number of points is controlled
        |Integrator| parameter ``nhcube_batch``.) The batch index
        is always first.

        Batch integrands can also be constructed from classes
        derived from :class:`vegas.BatchIntegrand`.

        Batch mode is particularly useful (and fast) when the class
        derived from :class:`vegas.BatchIntegrand` is coded
        in Cython. Then loops over the integration points
        can be coded explicitly, avoiding the need to use
        :mod:`numpy`'s whole-array operators if they are not
        well suited to the integrand.

        Any |vegas| parameter can also be reset: e.g.,
        ``self(fcn, nitn=20, neval=1e6)``.
        """

        wgt      = numpy.empty(1, numpy.float_)
        mean     = numpy.empty(1, numpy.float_)
        var      = numpy.empty((1, 1), numpy.float_)

        neval_batch = self.nhcube_batch * self.min_neval_hcube
        self.fdv2 = numpy.empty(neval_batch, numpy.float_)

        firsteval = True

        self.neval_hcube = (
            numpy.zeros(self.nhcube_batch, numpy.intp) + self.min_neval_hcube
            )

        if kargs:
            self.set(kargs)

        # synchronize random numbers across all processes (mpi)
        if self.sync_ran:
            self.synchronize_random()

        # Put integrand into standard form
        fcn = vegas._vegas.VegasIntegrand(fcn)

        sigf = self.sigf
        for itn in range(self.nitn):
            # if self.minimize_mem:
            #     self.set()
            if self.analyzer is not None:
                self.analyzer.begin(itn, self)

            # initalize arrays that accumulate results for a single iteration
            mean[:] = 0.0
            var[:, :] = 0.0
            sum_sigf = 0.0
                        
            # iterate batch-slices of integration points
            for x, y, wgt, hcube in self.random_batch(
                yield_hcube=True, yield_y=True, fcn=fcn
                ):
                fdv2 = self.fdv2        # must be inside loop
               
                ##############################################################################################
                # START OF MODIFICATION FROM THE ORIGINAL IMPLEMENTATION OF VEGAS3
                ##############################################################################################
                # Append the jacobian of the integrator as the last entry of the inputs passed to the function
                x = numpy.asarray( [ numpy.append(one_x,[wgt[i]]) for i, one_x in enumerate(x)], numpy.float_)
                ##############################################################################################
                # EMD OF MODIFICATION FROM THE ORIGINAL IMPLEMENTATION OF VEGAS3
                ##############################################################################################

                # evaluate integrand at all points in x                
                fx = fcn.eval(x)

                if firsteval:
                    # allocate work arrays on first pass through;
                    # (needed a sample fcn evaluation in order to do this)
                    firsteval = False
                    wf = numpy.empty(fcn.size, numpy.float_)
                    sum_wf = numpy.empty(fcn.size, numpy.float_)
                    sum_wf2 = numpy.empty((fcn.size, fcn.size), numpy.float_)
                    mean = numpy.empty(fcn.size, numpy.float_)
                    var = numpy.empty((fcn.size, fcn.size), numpy.float_)
                    mean[:] = 0.0
                    var[:, :] = 0.0
                    result = vegas._vegas.VegasResult(fcn, weighted=self.adapt)

                # compute integral and variance for each h-cube
                # j is index of hcube within batch, i is absolute index
                j = 0
                for i in range(hcube[0], hcube[-1] + 1):
                    # iterate over h-cubes
                    sum_wf[:] = 0.0
                    sum_wf2[:, :] = 0.0
                    neval = 0
                    while j < len(hcube) and hcube[j] == i:
                        for s in range(fcn.size):
                            wf[s] = wgt[j] * fx[j, s]
                            sum_wf[s] += wf[s]
                            for t in range(s + 1):
                                sum_wf2[s, t] += wf[s] * wf[t]
                        fdv2[j] = (wf[0] * self.neval_hcube[i - hcube[0]]) ** 2
                        j += 1
                        neval += 1
                    for s in range(fcn.size):
                        mean[s] += sum_wf[s]
                        for t in range(s + 1):
                            var[s, t] += (
                                sum_wf2[s, t] * neval - sum_wf[s] * sum_wf[t]
                                ) / (neval - 1.)
                        if var[s, s] <= 0:
                            var[s, s] = mean[s] ** 2 * 1e-15
                    sigf2 = abs(sum_wf2[0, 0] * neval - sum_wf[0] * sum_wf[0])
                    if self.beta > 0 and self.adapt:
                        if not self.minimize_mem:
                            sigf[i] = sigf2 ** (self.beta / 2.)
                            sum_sigf += sigf[i]
                        else:
                            sum_sigf += sigf2 ** (self.beta / 2.)
                    if self.adapt_to_errors and self.adapt:
                        # replace fdv2 with variance
                        fdv2[j - 1] = sigf2
                        self.map.add_training_data(
                            y[j - 1:, :], fdv2[j - 1:], 1
                            )
                if (not self.adapt_to_errors) and self.adapt and self.alpha > 0:
                    self.map.add_training_data(y, fdv2, y.shape[0])

            for s in range(var.shape[0]):
                for t in range(s):
                    var[t, s] = var[s, t]

            # accumulate result from this iteration
            result.update(mean, var)

            if self.beta > 0 and self.adapt:
                self.sum_sigf = sum_sigf
            if self.alpha > 0 and self.adapt:
                self.map.adapt(alpha=self.alpha)
            if self.analyzer is not None:
                result.update_analyzer(self.analyzer)

            if result.converged(self.rtol, self.atol):
                break
        return result.result

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
    verbosity = 1
    my_integrators = {}

    my_integrators['Naive MC'] = integrators.SimpleMonteCarloIntegrator(my_integrand, 
        n_iterations=10, n_points_per_iterations=1000, accuracy_target=None, verbosity = verbosity)

    my_integrators['VEGAS3'] = Vegas3Integrator(my_integrand, verbosity = verbosity)

    chosen_integrators = ['Naive MC', 'VEGAS3']

    for name in chosen_integrators:
        res = my_integrators[name].integrate()
        # Finally integrate
        print('\n'+'='*70)
        print('%-30s = %.4e +/- %.2e'%('%s :: Final result'%name, res[0], res[1]))
        print('='*70)

