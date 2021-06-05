import traceback
import os
import logging
import math
import shutil
import numpy as np
import random
import sys
import multiprocessing
import socket
import numpy
import time
import pickle

if __name__ == '__main__':
    sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.path.pardir, os.path.pardir))

import madgraph.various.misc as misc
import alpha_loop.integrator.integrators as integrators
import alpha_loop.integrator.integrands as integrands
import alpha_loop.integrator.functions as functions
import alpha_loop.utils as utils

# Dask dependencies
import dask
from dask_jobqueue import HTCondorCluster
from distributed import Client, LocalCluster
from dask import delayed
from dask.distributed import progress
from dask.distributed import as_completed

import yaml
from yaml import Loader

class NoAliasDumper(yaml.SafeDumper):
    def ignore_aliases(self, data):
        return True

import datetime as dt
class MyFormatter(logging.Formatter):
    converter=dt.datetime.fromtimestamp
    def formatTime(self, record, datefmt=None):
        ct = self.converter(record.created)
        if datefmt:
            s = ct.strftime(datefmt)
        else:
            t = ct.strftime("%Y-%m-%d %H:%M:%S")
            s = "%s,%03d" % (t, record.msecs)
        return s
formatter = MyFormatter(fmt='%(name)-15s {}%(asctime)s{} %(message)s'.format(utils.bcolors.BLUE,utils.bcolors.END),datefmt='%Y-%m-%d,%H:%M:%S.%f')
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.INFO)
console_handler.setFormatter(formatter)

logger = logging.getLogger('madgraph.havana')
logger.handlers.clear()
logger.addHandler(console_handler)
logger.propagate = False

pjoin = os.path.join

class HavanaIntegratorError(integrators.IntegratorError):
    """Exception raised if an exception is triggered in Havana.""" 

class HavanaMockUp(object):

    def __init__(self, 
            n_dimensions, SG_ids=None, n_channels_per_SG=None, target_result=None, seed=None, n_integrands=1, 
            reference_integrand_index=0, phase='real', grid_file='havana_grid.yaml', 
            optimize_on_variance=True, max_prob_ratio=1000., fresh_integration=False
        ):

        self.n_dimensions = n_dimensions
        self.target_result = target_result
        self.n_integrands = n_integrands
        self.reference_integrand_index = reference_integrand_index
        self.phase = phase
        self.optimize_on_variance = optimize_on_variance
        self.max_prob_ratio = max_prob_ratio

        if self.phase=='real':
            self.reference_result_index = reference_integrand_index*2
        else:
            self.reference_result_index = reference_integrand_index*2+1

        self.n_discrete_dimension = 0
        if SG_ids is not None:
            self.n_discrete_dimension += 1
        if n_channels_per_SG is not None:
            self.n_discrete_dimension += 1

        self.discrete_grid = None
        self.SG_ids = SG_ids
        self.n_channels_per_SG = n_channels_per_SG

        self.grid_file = grid_file
        self.fresh_integration = fresh_integration

        if seed is not None:
            random.seed(seed)
            numpy.random.seed(seed)

        self.reset( clean = self.fresh_integration )

    def reset(self, clean=False):

        if not clean and os.path.exists(self.grid_file):
            
            logger.info("%sLoading Havana state and grids from %s.%s"%(
                utils.bcolors.RED, self.grid_file, utils.bcolors.END
            ))
            with open(self.grid_file,'r') as f:
                flat_record = yaml.load(f, Loader=Loader)
            self.current_integral_estimate = flat_record['current_integral_estimate']
            self.current_error_estimate = flat_record['current_error_estimate']
            self.wgt_sum = flat_record['wgt_sum']
            self.wgt_sum_squared = flat_record['wgt_sum_squared']
            self.n_points = flat_record['n_points']
            self.discrete_grid = flat_record['discrete_grid']
            self.n_iterations = flat_record['n_iterations']

        else:

            self.current_integral_estimate = [0.,]*(self.n_integrands*2)
            self.current_error_estimate = [0.,]*(self.n_integrands*2)
            self.wgt_sum = [0.,]*(self.n_integrands*2)
            self.wgt_sum_squared = [0.,]*(self.n_integrands*2)
            self.n_points = 0
            self.n_iterations = 0

            self.discrete_grid = None
            if self.SG_ids is not None:
                self.discrete_grid = []
                for i_SG, (SG_id, SG_name) in enumerate(self.SG_ids):
                    self.discrete_grid.append(
                        {
                            'SG_name'     : SG_name,
                            'SG_id'       : SG_id,
                            'n_points'    : 0,
                            'cum_wgt'     : 0.,
                            'cum_sqr_wgt' : 0.,
                            'integral_estimate': 0.,
                            'error_estimate'   : 0.,
                            'p'           : 1./len(self.SG_ids),
                            'p_cum'       : (1./len(self.SG_ids))*(i_SG+1),
                            'channel_grid'     :
                                None if self.n_channels_per_SG is None else [
                                    {
                                        'n_points'    : 0,
                                        'cum_wgt'     : 0.,
                                        'cum_sqr_wgt' : 0.,
                                        'integral_estimate': 0.,
                                        'error_estimate'   : 0.,
                                        'p'           : 1./self.n_channels_per_SG[i_SG],
                                        'p_cum'       : (1./self.n_channels_per_SG[i_SG])*(i_channel+1)
                                    } for i_channel in range(self.n_channels_per_SG[i_SG])
                                ]
                        }
                    )

    def dump_grid(self, file_name=None):

        flat_record = {
            'current_integral_estimate' : self.current_integral_estimate,
            'current_error_estimate' : self.current_error_estimate,
            'wgt_sum' : self.wgt_sum,
            'wgt_sum_squared' : self.wgt_sum_squared,
            'n_points' : self.n_points,
            'discrete_grid' : self.discrete_grid,
            'n_iterations' :self.n_iterations
        }
        with open(self.grid_file if file_name is None else file_name,'w') as f:
            f.write(yaml.dump(flat_record, Dumper=NoAliasDumper, default_flow_style=False))

    def sample(self, batch_size):

        sample = numpy.random.rand(batch_size,1+self.n_discrete_dimension+self.n_dimensions)
        for i_sample in range(batch_size):
            havana_jacobian = 1.0
            if self.discrete_grid is not None:
                selected_SG = None
                for i_SG, SG_bin in enumerate(self.discrete_grid):
                    if sample[i_sample][1]<=SG_bin['p_cum']:
                        selected_SG = i_SG
                        break
                if selected_SG is None:
                    raise HavanaIntegratorError("Could not select SG.")
                sample[i_sample][1] = float(self.discrete_grid[selected_SG]['SG_id'])
                havana_jacobian *= (1./self.discrete_grid[selected_SG]['p'])
                selected_channel = None
                if self.discrete_grid[selected_SG]['channel_grid'] is not None:
                    for i_channel, channel_bin in enumerate(self.discrete_grid[selected_SG]['channel_grid']):
                        if sample[i_sample][2]<=channel_bin['p_cum']:
                            selected_channel = i_channel
                            break
                    if selected_channel is None:
                        raise HavanaIntegratorError("Could not select channel.")
                    sample[i_sample][2] = float(selected_channel)
                    havana_jacobian *= (1./self.discrete_grid[selected_SG]['channel_grid'][selected_channel]['p'])

            sample[i_sample][0] = havana_jacobian
        
        return sample

    def accumulate_results(self, results):
        
        self.n_points += len(results)
        if self.n_points == 0:
            return

        for index in range(0,self.n_integrands*2):
            self.wgt_sum[index] += float(sum(res[index]*res[self.n_integrands*2] for res in results))
            self.wgt_sum_squared[index] += float(sum((res[index]*res[self.n_integrands*2])**2 for res in results))

            self.current_integral_estimate[index] = self.wgt_sum[index] / self.n_points
            self.current_error_estimate[index] = math.sqrt( ((self.wgt_sum_squared[index] / self.n_points) - self.current_integral_estimate[index]**2)/self.n_points )

        if self.discrete_grid is not None:
            for res in results:
                wgt = float(res[self.reference_result_index])
                wgt_squared = wgt**2
                SG_index = int(res[1+self.n_integrands*2])
                self.discrete_grid[SG_index]['n_points'] += 1
                self.discrete_grid[SG_index]['cum_wgt'] += wgt
                self.discrete_grid[SG_index]['cum_sqr_wgt'] += wgt_squared
                if self.discrete_grid[SG_index]['channel_grid'] is not None:
                    channel_index = int(res[2+self.n_integrands*2])
                    self.discrete_grid[SG_index]['channel_grid'][channel_index]['n_points'] += 1
                    self.discrete_grid[SG_index]['channel_grid'][channel_index]['cum_wgt'] += wgt
                    self.discrete_grid[SG_index]['channel_grid'][channel_index]['cum_sqr_wgt'] += wgt_squared

    def sync_grids(self, max_prob_ratio=None):
        
        self.sync_nested_grid(self.discrete_grid, max_prob_ratio=max_prob_ratio)

        if self.discrete_grid is not None:
            for SG_bin in self.discrete_grid:
                self.sync_nested_grid(SG_bin['channel_grid'], max_prob_ratio=max_prob_ratio)

    def update_iteration_count(self, dump_grids=True):

        self.n_iterations += 1

        if dump_grids:
            # Dump current grid and a copy to keep track of its evolution at each iteration.
            self.dump_grid()
            shutil.copy(self.grid_file, self.grid_file.replace('.yaml','_iteration_%d.yaml'%self.n_iterations) )

    def sync_nested_grid(self, nested_grid, max_prob_ratio=None):

        if nested_grid is None:
            return

        if max_prob_ratio is None:
            max_prob_ratio = self.max_prob_ratio

        for bin in nested_grid:
            if bin['n_points']==0:
                bin['integral_estimate'] = 0.
                bin['error_estimate'] = 0.
            else:
                bin['integral_estimate'] = ( bin['cum_wgt'] / bin['n_points'] )
                bin['error_estimate'] = math.sqrt( ((bin['cum_sqr_wgt'] / bin['n_points']) - bin['integral_estimate']**2)/bin['n_points'] )

            if self.optimize_on_variance:
                bin['p'] = bin['error_estimate']
            else:
                bin['p'] = abs(bin['integral_estimate'])

        max_p = max(bin['p'] for bin in nested_grid)
        if max_p == 0.:
            for bin in nested_grid:
                bin['p'] = 1./len(nested_grid)
        else:
            running_p_sum = 0.
            for bin in nested_grid:
                bin['p'] = max(bin['p'],max_p/max_prob_ratio) if max_prob_ratio is not None else bin['p']
                running_p_sum += bin['p']
                bin['p_cum'] = running_p_sum
            for bin in nested_grid:
                bin['p'] /= running_p_sum
                bin['p_cum'] /= running_p_sum

    def get_grid_summary(self, sort_by_variance=False, show_channel_grid=True):

        if self.discrete_grid is None:
            return 'N/A'

        sorted_bins = sorted([b for b in self.discrete_grid], key=lambda b: b['error_estimate'] if sort_by_variance else abs(b['integral_estimate']), reverse=True)
        res = ['Distribution of %d SGs:'%len(sorted_bins)]
        for bins_chunk in [sorted_bins[i_chunk:i_chunk+5] for i_chunk in range(0,len(sorted_bins),5)]:
            res.append('  '+' | '.join( '%-40s'%(
                '%s(I=%.5g(%.2g%%),p=%.2f%%)'%(
                    bin['SG_name'], bin['integral_estimate'], 
                    abs(bin['error_estimate']/bin['integral_estimate'])*100. if bin['integral_estimate']!=0. else 0.,
                    bin['p']*100.
                )
            ) for bin in bins_chunk ))

        if show_channel_grid:
            for SG_bin in sorted_bins:
                if SG_bin['channel_grid'] is not None:
                    res.append('Distribution of %d channels for SG %s:'%(len(SG_bin['channel_grid']), SG_bin['SG_name']))
                    sorted_channel_bins = sorted([(i_c, b) for i_c, b in enumerate(SG_bin['channel_grid'])], key=lambda b: b[1]['error_estimate'] if sort_by_variance else abs(b[1]['integral_estimate']), reverse=True)
                    for bins_chunk in [sorted_channel_bins[i_chunk:i_chunk+5] for i_chunk in range(0,len(sorted_channel_bins),5)]:
                        res.append('  '+' | '.join( '%-30s'%(
                            'C%d(I=%.5g(%.2g%%),p=%.2f%%)'%(
                                i_c, bin['integral_estimate'], 
                                abs(bin['error_estimate']/bin['integral_estimate'])*100. if bin['integral_estimate']!=0. else 0.,
                                bin['p']*100.
                            )
                        ) for i_c, bin in bins_chunk ))

        return '\n'.join(res)

    def get_summary(self):

        res = ["Result after %.1fM evaluations and %d iterations: %.5g +/- %.5g (%.2g%%)"%(
            self.n_points/1.0e6, self.n_iterations, self.current_integral_estimate[self.reference_result_index], self.current_error_estimate[self.reference_result_index], 
            0. if self.current_integral_estimate[self.reference_result_index]==0. else 
            (abs(self.current_error_estimate[self.reference_result_index]/self.current_integral_estimate[self.reference_result_index]))*100.
        )]
        if self.target_result is not None:
            res.append( ("%-{}s".format(len('Result after %d evaluations and %d iterations'%(self.n_points, self.n_iterations))))%("vs target")+
                ": %.5g del %.5g (%.2g%%)"%(
                    self.target_result, 
                    self.current_integral_estimate[self.reference_result_index]-self.target_result, 
                    0. if self.target_result==0. else ((self.current_integral_estimate[self.reference_result_index]-self.target_result)/self.target_result)*100.0
                ))
        for res_index in range(self.n_integrands*2):
            res.append('  %s %s[I%d] = %.5g +/- %.5g (%.2g%%)'%(
                '  ' if res_index!=self.reference_result_index else '->',
                'Re' if res_index%2==0 else 'Im',
                res_index//2,
                self.current_integral_estimate[res_index], self.current_error_estimate[res_index],
                0. if self.current_integral_estimate[res_index]==0. else 
                (abs(self.current_error_estimate[res_index]/self.current_integral_estimate[res_index]))*100.
            ))

        return '\n'.join(res)

    def get_current_estimate(self):

        return (self.current_integral_estimate[self.reference_result_index], self.current_error_estimate[self.reference_result_index])

def HavanaIntegrandWrapper(
        job_id,
        samples_batch,
        integrands_constructor_args,
        SG_ids_specified=True,
        channels_specified=True
    ):

    start_time = time.time()
    
    IO_pickle_out_name = None
    if isinstance(samples_batch, str):
        IO_pickle_out_name = os.path.join(os.path.dirname(samples_batch),'dask_output_job_%d.pkl'%job_id )
        sample_pickle_path = samples_batch
        samples_batch = pickle.load( open( sample_pickle_path, 'rb' ) )
        os.remove(sample_pickle_path)

    from alpha_loop.integrator.sampler import DaskHavanaALIntegrand

    integrands = [
        DaskHavanaALIntegrand(**integrand_constructor_args) for integrand_constructor_args in integrands_constructor_args
    ]
    all_res = numpy.ndarray(shape=(len(samples_batch),2*len(integrands)+len(samples_batch[0])), dtype=float)
    for i_sample, sample in enumerate(samples_batch):
        all_res[i_sample][2*len(integrands):] = sample

    if SG_ids_specified and channels_specified:
        SG_id_per_sample = [int(s[1]) for s in samples_batch]
        channel_id_per_sample = [int(s[2]) for s in samples_batch]
        xs_batch = [list(s[3:]) for s in samples_batch]
    elif SG_ids_specified and (not channels_specified):
        SG_id_per_sample = [int(s[1]) for s in samples_batch]
        channel_id_per_sample = None
        xs_batch = [list(s[2:]) for s in samples_batch]
    elif (not SG_ids_specified) and (not channels_specified):
        SG_id_per_sample = None
        channel_id_per_sample = None
        xs_batch = [list(s[1:]) for s in samples_batch]
    else:
        raise HavanaIntegratorError("Havana cannot integrate specific integration channels while summing over all supergraphs.")

    for i_integrand, integrand in enumerate(integrands):

        results = integrand(xs_batch, SG_id_per_sample = SG_id_per_sample, channel_id_per_sample = channel_id_per_sample)
        for i_res, res in enumerate(results):
            all_res[i_res][0+i_integrand*2] = res.real
            all_res[i_res][1+i_integrand*2] = res.imag

    if IO_pickle_out_name is not None:
        with open(IO_pickle_out_name, 'wb') as f:
            pickle.dump( all_res, f )
        all_res = IO_pickle_out_name

    return (job_id, time.time()-start_time, all_res)

class HavanaIntegrator(integrators.VirtualIntegrator):
    """ Steering of the havana integrator """
    
    _SUPPORTED_CLUSTER_ARCHITECTURES = [ 'dask_local', 'dask_condor' ]
    _DEBUG = False

    def __init__(self, integrands,
                 cross_section_set=None,
                 all_supergraphs=None,
                 run_workspace=None,
                 accuracy_target=None,
                 n_iterations=None,
                 n_start = 10000,
                 n_increase = 10000,
                 verbosity = logging.INFO,
                 n_max = None,
                 seed = None,
                 n_workers = None,
                 batch_size = 1e6,
                 cluster_type = 'dask_local',
                 dask_local_options = None,
                 dask_condor_options = None,
                 target_result = None,
                 MC_over_SGs = True,
                 MC_over_channels = True,
                 selected_SGs = None, # Means ALL
                 stream_monitor = True,
                 havana_optimize_on_variance=True, 
                 havana_max_prob_ratio=1000.,
                 phase = 'real',
                 show_SG_grid = True,
                 show_channel_grid = False,
                 show_grids_sorted_by_variance = False,
                 dump_havana_grids = True,
                 fresh_integration = False,
                 pickle_IO = False,
                 **opts):

        """ Initialize the simplest MC integrator."""
        
        self.accuracy_target = accuracy_target
        self.n_iterations = n_iterations
        self.cross_section_set = cross_section_set
        self.all_supergraphs = all_supergraphs
        if not self._DEBUG:
            self.stream_monitor = stream_monitor
        else:
            self.stream_monitor = False
        self.phase = phase
        self.run_workspace = run_workspace
        self.havana_optimize_on_variance = havana_optimize_on_variance
        self.havana_max_prob_ratio = havana_max_prob_ratio
        self.show_SG_grid = show_SG_grid
        self.show_channel_grid = show_channel_grid
        self.show_grids_sorted_by_variance = show_grids_sorted_by_variance
        self.dump_havana_grids = dump_havana_grids
        self.fresh_integration = fresh_integration
        self.pickle_IO = pickle_IO

        self.n_start = n_start
        self.n_max = n_max
        self.n_increase = n_increase

        logger.setLevel(verbosity)

        if n_workers is None:
            self.n_workers = multiprocessing.cpu_count()
        else:
            self.n_workers = n_workers

        self.dask_local_options = {'n_workers' : self.n_workers}
        if sys.platform == "darwin":
            self.dask_local_options['interface'] = 'lo0'
        if dask_local_options is not None:
            self.dask_local_options.update(dask_local_options)

        self.dask_condor_options = {
            'memory'  : '2000MB',
            'disk'    : '1000MB',
            'IO_port' : 8786,
            'n_cores_per_worker' : 1,
            'job_flavour' : 'tomorrow'
        }
        if dask_condor_options is not None:
            self.dask_condor_options.update(dask_condor_options)

        if cluster_type not in self._SUPPORTED_CLUSTER_ARCHITECTURES:
            raise HavanaIntegratorError("Integrator %s only support the following cluster_types: %s, not '%s'."%(
                self.get_name(), self._SUPPORTED_CLUSTER_ARCHITECTURES, cluster_type
            ))
        
        self.batch_size = batch_size
        self.cluster_type = cluster_type

        self.target_result = target_result
        self.MC_over_SGs = MC_over_SGs
        self.MC_over_channels = MC_over_channels
        self.selected_SGs = selected_SGs

        if (self.selected_SGs is not None) and (not self.MC_over_SGs):
            raise HavanaIntegratorError("A selection of supergraphs can only be specified when Monte-Carlo'ing over them.")

        self.seed = seed
        self.havana = None

        super(HavanaIntegrator, self).__init__(integrands, **opts)
        #misc.sprint(self.integrands)
        
    def dask_integrate(self, client):

        n_dimensions = len(self.integrands[0].dimensions.get_continuous_dimensions())
        if not all(len(integrand.dimensions.get_continuous_dimensions())==n_dimensions for integrand in self.integrands):
            raise HavanaIntegratorError("In Havana implementations, all integrands must have the same dimension.")

        if not self.MC_over_SGs:
            SG_ids = None
        else:
            if self.selected_SGs is None:
                SG_ids = [ (i_SG, SG['name']) for i_SG, SG in enumerate(self.cross_section_set['topologies']) ]
            else:
                SG_ids = []
                for i_SG, SG in enumerate(self.cross_section_set['topologies']):
                    if SG['name'] in self.selected_SGs:
                        SG_ids.append( (i_SG, SG['name']) )                    

                if len(SG_ids)!=len(self.selected_SGs) or len(SG_ids)==0:
                    raise HavanaIntegratorError("Not all specified SG names were found in the cross_section set.")

        if not self.MC_over_channels:
            n_channels_per_SG = None
        else:
            n_channels_per_SG = []
            for (i_SG, SG_name) in SG_ids:
                n_channels_per_SG.append(len(self.all_supergraphs[SG_name]['multi_channeling_bases']))

        self.havana = HavanaMockUp(
            n_dimensions, 
            SG_ids=SG_ids, 
            n_channels_per_SG=n_channels_per_SG, 
            target_result=self.target_result, 
            seed=self.seed,
            n_integrands = len(self.integrands),
            reference_integrand_index = 0,
            phase = self.phase,
            grid_file = pjoin(self.run_workspace,'havana_grid.yaml'),
            optimize_on_variance = self.havana_optimize_on_variance,
            max_prob_ratio = self.havana_max_prob_ratio,
            fresh_integration = self.fresh_integration
        )

        # Now perform the integration
        
        logger.info("Staring alphaLoop integration with Havana, lean back and enjoy...")
        logger.info("Visit http://localhost:8787/status to follow Dask jobs from a dashboard.")
        start_time = time.time()

        current_iteration = 1
        current_n_points = self.n_start
        current_step = self.n_increase
        n_tot_points = 0

        integrands_constructor_args = [
            integrand.get_constructor_arguments() for integrand in self.integrands
        ]
        scattered_integrands_constructor_args = client.scatter(integrands_constructor_args)

        job_ID = 0

        self.canvas = None

        n_jobs_total_completed = 0
        cumulative_job_time = 0.
        cumulative_processing_time = 0.
        cumulative_IO_time = 0.
        
        while True:
            
            n_remaining_points = current_n_points

            n_points_for_this_iteration = 0

            futures = []
            logger.info("Preparing samples and submission for iteration #%d"%current_iteration)
            if self._DEBUG: logger.info("Starting submission now...")
            while n_remaining_points > 0:
                this_n_points = min(n_remaining_points, self.batch_size)
                n_remaining_points -= this_n_points   
                t2 = time.time()           
                job_ID += 1
                if self._DEBUG: logger.info("Generating sample for job #%d..."%job_ID)
                this_sample = self.havana.sample(this_n_points)
                if self._DEBUG: logger.info("Done with sample generation.")
                cumulative_processing_time = time.time() - t2
    
                t0 = time.time()
                if self.pickle_IO:
                    if self._DEBUG: logger.info("Dumping input of job %d to pickle..."%job_ID)
                    IO_path = pjoin(self.run_workspace,'dask_input_job_%d.pkl'%job_ID)
                    with open(IO_path, 'wb') as f:
                        pickle.dump( this_sample, f )
                    if self._DEBUG: logger.info("Done with pickle dump.")
                    this_scattered_sample = str(IO_path)
                else:
                    # Scattering data is not really necessary as I believe this is the right thing to do only for data common to multiple jobs. 
                    # But I don't like the warning when not doing so, so I'll put it.
                    if self._DEBUG: logger.info("Scattering sampling input for job #%d..."%job_ID)
                    this_scattered_sample = client.scatter(this_sample)
                    if self._DEBUG: logger.info("Done with scattering sampling input.")

                if self._DEBUG: logger.info("Submitting job %d."%job_ID)
                futures.append(client.submit(HavanaIntegrandWrapper, 
                    job_ID,
                    this_scattered_sample,
                    scattered_integrands_constructor_args,
                    self.MC_over_SGs,
                    self.MC_over_channels
                ))
                if self._DEBUG: logger.info("Done with job submission.")
                cumulative_IO_time += time.time()-t0

            logger.info("Now running iteration #%d"%current_iteration)

            #print("DONE SUBMITTING")
            n_submitted = len(futures)
            n_done = 0

            self.update_status(
                    start_time, n_jobs_total_completed, n_submitted, 0,
                    n_tot_points, n_points_for_this_iteration, current_n_points, cumulative_IO_time,
                    cumulative_processing_time, cumulative_job_time, current_iteration)

            try:

                if self._DEBUG: logger.info("Now waiting for a job to complete...")
                for completed_future in as_completed(futures):
                    
                    t1 = time.time()
                    if self._DEBUG: logger.info("Now receiving result from a job...")
                    job_ID, run_time, results = completed_future.result()
                    if self.pickle_IO:
                        if self._DEBUG: logger.info("Loading pickled results for job %d..."%job_ID)
                        results_pickle_path = results
                        results = pickle.load( open( results_pickle_path, 'rb' ) )
                        os.remove(results_pickle_path)
                        if self._DEBUG: logger.info("Finished loading pickled results of %d samples."%len(results))
                    else:
                        if self._DEBUG: logger.info("Finished receiving %d sample result from job %d..."%(len(results),job_ID))

                    cumulative_IO_time += time.time()-t1
                    cumulative_job_time += run_time

                    t2 = time.time()
                    #logger.debug("Job #%d completed in %.0f s for %d points."%(job_ID, run_time, len(results) ))
                    if self._DEBUG: logger.info("Accumulating new results from job %d into havana..."%job_ID)
                    self.havana.accumulate_results(results)
                    max_prob_ratio_for_this_iteration = max( self.havana_max_prob_ratio / max( 100. - 10.*(current_iteration-1) , 1.), 3.)
                    if self._DEBUG: logger.info("Syncing Havana grids with the new results...")
                    self.havana.sync_grids( max_prob_ratio = max_prob_ratio_for_this_iteration)
                    if self._DEBUG: logger.info("Done updating Havana grids...")

                    for res in results:
                        xs = res[1+len(self.integrands)*2:]
                        for index in range(0,len(self.integrands)):
                            if self.phase=='real':
                                wgt = res[index*2]
                            else:
                                wgt = res[1+index*2]
                            havana_jacobian = res[len(self.integrands)*2]
                            self.integrands[index].update_evaluation_statistics(xs, wgt*havana_jacobian)
                    
                    n_done += 1
                    n_jobs_total_completed += 1
                    n_tot_points += len(results)
                    n_points_for_this_iteration += len(results)
                    
                    self.update_status(
                        start_time, n_jobs_total_completed, n_submitted, n_done,
                        n_tot_points, n_points_for_this_iteration, current_n_points, cumulative_IO_time,
                        cumulative_processing_time, cumulative_job_time, current_iteration)

                    if self._DEBUG: logger.info("Now waiting for a job to complete...")

                self.havana.update_iteration_count(dump_grids=self.dump_havana_grids)
                cumulative_processing_time = time.time() - t2

            except KeyboardInterrupt:
                logger.warning("Aborting Havana integration now.")
                break

            if self.n_iterations is not None and current_iteration >= self.n_iterations:
                logger.info("Max number of iterations %d reached."%self.n_iterations) 
                break               

            if self.n_max is not None and n_tot_points >= self.n_max:
                logger.info("Max number of sample points %d reached."%self.n_max)
                break

            res_int, res_error = self.havana.get_current_estimate()
            if self.accuracy_target is not None and (res_error/abs(res_int) if res_int!=0. else 0.)<self.accuracy_target:
                logger.info("Target accuracy of %.2g reached."%self.accuracy_target)
                break

            current_n_points += current_step
            current_step += self.n_increase
            current_iteration += 1

            # Keep printout of last step of this iteration.
            self.canvas = None

        return self.havana.get_current_estimate()

    def update_status( self,
        start_time, n_jobs_total_completed, n_submitted, n_done,
        n_tot_points, n_points_for_this_iteration, current_n_points, cumulative_IO_time,
        cumulative_processing_time, cumulative_job_time, current_iteration_number
    ):

        curr_integration_time = time.time()-start_time

        # Update the run stat line and monitoring canvas
        monitoring_report = []
        monitoring_report.append( '| Jobs: completed = %d (avg %s) and running = %d/%d on %d workers.\n| Total n_pts: %.1fM ( %s ms / pt on one core, %s pts / s overall ). n_pts for this iteration #%d: %.2fM/%.2fM (%.1f%%).'%(
            n_jobs_total_completed, '%.3g min/job'%((cumulative_job_time/n_jobs_total_completed)/60.) if n_jobs_total_completed>0 else 'N/A', n_submitted-n_done, n_submitted, self.n_workers,
            n_tot_points/(1.e6), '%.3g'%((curr_integration_time/n_tot_points)*self.n_workers*1000.) if n_tot_points>0 else 'N/A', 
            ('%.1fK'%(n_tot_points/curr_integration_time/1000.) if n_tot_points>0 else 'N/A'),
            current_iteration_number,
            n_points_for_this_iteration/(1.e6), current_n_points/(1.e6), (float(n_points_for_this_iteration)/current_n_points)*100.
        ))
        wall_time = curr_integration_time
        IO_time = (cumulative_IO_time / wall_time)*100.
        processing_time = (cumulative_processing_time / wall_time)*100.
        cpu_hours = (cumulative_job_time / 3600.)
        parallelisation_efficiency = (cumulative_job_time / (wall_time*self.n_workers) )*100.
        monitoring_report.append('| Wall time: %.1f h, IO: %.2f%%, processing: %.2f%%, jobs: %.1f CPU-hours, efficiency: %.2f%%'%(
            wall_time/3600., IO_time, processing_time, cpu_hours, parallelisation_efficiency
        ))
        monitoring_report.append('')
        monitoring_report.append( self.havana.get_summary() )

        if self.show_SG_grid:
            monitoring_report.append( self.havana.get_grid_summary( 
                sort_by_variance=self.show_grids_sorted_by_variance, 
                show_channel_grid= (self.MC_over_channels and self.show_channel_grid)
            ) )
        monitoring_report.append('')

        if self.canvas is None:
            self.canvas = utils.Canvas( '\n'.join(monitoring_report), stream=self.stream_monitor, overwrite=True )
        else:
            self.canvas.print('\n'.join(monitoring_report))


    def integrate(self):
        """ Return the final integral and error estimates."""

        logger.info("Now deploying the Dask cluster...")
        if self.cluster_type == 'dask_local':
            
            with LocalCluster(**self.dask_local_options) as cluster:                
                client = Client(cluster)
                return self.dask_integrate(client)

        elif self.cluster_type == 'dask_condor':

            with HTCondorCluster(
                    cores=self.dask_condor_options['n_cores_per_worker'],
                    memory=self.dask_condor_options['memory'],
                    disk=self.dask_condor_options['disk'],
                    death_timeout = '60',
                    nanny = False,
                    scheduler_options={
                        'port': self.dask_condor_options['IO_port'],
                        'host': socket.gethostname()
                    },
                    job_extra={
                        'log'    : 'dask_job_output.log',
                        'output' : 'dask_job_output.out',
                        'error'  : 'dask_job_output.err',
                        'should_transfer_files'   : 'Yes',
                        'when_to_transfer_output' : 'ON_EXIT',
                        '+JobFlavour' : self.dask_condor_options['job_flavour'],
                    },
                    extra = [ '--worker-port {}'.format(self.dask_condor_options['IO_port']) ]
            ) as cluster:

                logger.info("Scaling condor cluster to %d workers..."%self.n_workers)
                cluster.scale(self.n_workers)
                logger.info("Initialising condor cluster...")
                logger.debug("Condor submission script:\n%s"%cluster.job_script())
                client = Client(cluster)
                return self.dask_integrate(client)

        else:
            raise HavanaIntegratorError("Cluster architecture %s not supported."%self.cluster_type)