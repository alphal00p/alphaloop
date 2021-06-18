from asyncio.tasks import create_task
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
import subprocess
import glob
import re

alphaloop_basedir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.path.pardir, os.path.pardir))
if alphaloop_basedir not in sys.path:
    sys.path.append(alphaloop_basedir)
#if __name__ == '__main__':
#    sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.path.pardir, os.path.pardir))

import madgraph.various.misc as misc
import alpha_loop.integrator.integrators as integrators
import alpha_loop.integrator.integrands as integrands
import alpha_loop.integrator.functions as functions
import alpha_loop.utils as utils
from alpha_loop.utils import bcolors
from alpha_loop.integrator.worker import HavanaIntegrandWrapper, ALStandaloneIntegrand, HavanaMockUp, Havana

import asyncio
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


def fast_remove(path, trashcan):
    if trashcan is None:
        if os.path.isdir(path):
            shutil.rmtree(path)
        else:
            os.remove(path)
    else:
        shutil.move(path, pjoin(trashcan, os.path.basename(path)))

class HavanaIntegratorError(Exception):
    pass

class AL_cluster(object):

    _SUPPORTED_CLUSTER_ARCHITECTURES = ['local', 'condor']
    _JOB_CHECK_FREQUENCY = 0.3
    _FORWARD_WORKER_OUTPUT = True

    def __init__(self, architecture, n_workers, n_cores_per_worker, run_workspace, process_results_callback, run_id, monitor_callback=None, cluster_options=None, keep=False, debug=False, trashcan=None):

        self.n_workers = n_workers
        self.n_cores_per_worker = n_cores_per_worker
        self.run_workspace = run_workspace
        self.available_worker_ids = asyncio.Queue()
        self.n_jobs_in_queue = 0
        self.jobs_queue = asyncio.Queue()
        self.workers = { }
        self.all_worker_hooks = []
        self.active = True

        self.process_results_callback = process_results_callback
        self.monitor_callback = monitor_callback
        self.architecture = architecture
        self.cluster_options = {
            'RAM_required' : 1024,
            'job_flavour'  : 'tomorrow'
        }
        if cluster_options is not None:
            self.cluster_options.update(cluster_options)
        
        self.debug  = debug
        self.keep   = keep

        self.run_id = run_id
        self.current_status = None

        if self.architecture == 'local':
            self.work_monitoring_path = pjoin(self.run_workspace,'run_%d_follow_workers'%self.run_id)
            self.worker_monitoring_file = open(self.work_monitoring_path,'a')
        else:
            self.work_monitoring_path = None
            self.worker_monitoring_file = None

        self.trashcan = trashcan

        if self.cluster_options is None:
            self.cluster_options = {}

        if self.architecture not in self._SUPPORTED_CLUSTER_ARCHITECTURES:
            raise HavanaIntegratorError("Cluster architecture not supported by ALcluster: %s"%architecture)
    
    async def update_status(self):

        all_statuses = [ v['status'] for v in self.workers.values() ]
        new_status = {
            'FREE' : all_statuses.count('FREE'),
            'PENDING' : all_statuses.count('PENDING'),
            'RUNNING' : all_statuses.count('RUNNING')
        }
        if self.current_status is None or new_status!=self.current_status:
            self.current_status = new_status
            if self.monitor_callback is not None:
                await self.monitor_callback(self.current_status)

    async def send_job(self, worker_id, payload):

        input_job_path = pjoin(self.run_workspace, 'run_%d_job_input_worker_%d.pkl'%(self.run_id, worker_id))
        input_job_path_done = pjoin(self.run_workspace, 'run_%d_job_input_worker_%d.done'%(self.run_id, worker_id))
        if os.path.exists(input_job_path):
            raise HavanaIntegratorError("Attempted to send a job to a busy worker, this should never happen.")

        with open(input_job_path,'wb') as f:
            pickle.dump( payload, f )
        with open(input_job_path_done,'w') as f:
            f.write("Marker file for job input dump done.")

        self.workers[worker_id]['status'] = 'RUNNING'
        await self.update_status()

    async def follow_workers(self):
        while self.active:
            try:
                if self.debug: logger.info("Scanning for new job results...")
                io_files = [os.path.basename(f) for f in glob.glob(pjoin(self.run_workspace, 'run_%d_job_output_worker_*.done'%self.run_id))]
                output_worker_ids = []
                for io_file in io_files:
                    output_match = re.match(r'^run_(\d+)_job_output_worker_(?P<worker_id>\d+).done$',io_file)
                    if not output_match:
                        raise HavanaIntegratorError("Unexpected IO file found: %s"%io_file)
                    worker_id = int(output_match.group('worker_id'))
                    if os.path.exists(pjoin(self.run_workspace, 'run_%d_job_output_worker_%d.pkl'%(self.run_id,worker_id))):
                        output_worker_ids.append(worker_id)

                if self.debug and len(output_worker_ids)>0: logger.info("Found new available results from workers: %s"%str(output_worker_ids))
                for output_worker_id in output_worker_ids:
                    output_path = pjoin(self.run_workspace, 'run_%d_job_output_worker_%d.pkl'%(self.run_id, output_worker_id))
                    output_path_done = pjoin(self.run_workspace, 'run_%d_job_output_worker_%d.done'%(self.run_id, output_worker_id))
                    if self.debug: logger.info("Deserialising the following job result: %s"%str(output_path))
                    job_result = pickle.load( open( output_path, 'rb' ) )
                    if self.debug: logger.info("Done.")
                    input_path = pjoin(self.run_workspace, 'run_%d_job_input_worker_%d.pkl'%(self.run_id, output_worker_id))
                    input_path_done = pjoin(self.run_workspace, 'run_%d_job_input_worker_%d.done'%(self.run_id, output_worker_id))
                    if not os.path.exists(input_path):
                        raise HavanaIntegratorError("Matching input file for worker %d not found."%output_worker_id)

                    fast_remove(input_path, self.trashcan)
                    fast_remove(input_path_done, self.trashcan)
                    fast_remove(output_path, self.trashcan)
                    fast_remove(output_path_done, self.trashcan)
                    
                    await self.available_worker_ids.put(output_worker_id)
                    self.workers[output_worker_id]['status'] = 'FREE'
                    await asyncio.sleep(0.01)

                    await self.update_status()

                    if job_result != 'pong':
                        if self.debug: logger.info("Processing results from worker: %d"%output_worker_id)
                        await self.process_results_callback(job_result)
                        if self.debug: logger.info("Done.")
                        #asyncio.create_task(self.process_results_callback(job_result))

                if self.architecture == 'local':
                    await asyncio.sleep(self._JOB_CHECK_FREQUENCY)
                else:
                    await asyncio.sleep(self._JOB_CHECK_FREQUENCY*15.)
            except KeyboardInterrupt as e:
                break

    async def deploy(self, wait_for_one_worker=True):

        if self.architecture == 'local':

            for i_worker in range(self.n_workers):

                worker_id_min = i_worker*self.n_cores_per_worker
                worker_id_max = (i_worker+1)*self.n_cores_per_worker-1

                cmd = [ sys.executable, pjoin(alphaloop_basedir,'alpha_loop','integrator','worker.py'), 
                        '--run_id', str(self.run_id), '--worker_id_min', str(worker_id_min), '--worker_id_max', str(worker_id_max), 
                        '--workspace_path', str(self.run_workspace), '--timeout', '%.2f'%self._JOB_CHECK_FREQUENCY]
                self.all_worker_hooks.append(subprocess.Popen(
                        cmd,
                        cwd=self.run_workspace,
                        stdout=(self.worker_monitoring_file if self._FORWARD_WORKER_OUTPUT else subprocess.DEVNULL),
                        stderr=(self.worker_monitoring_file if self._FORWARD_WORKER_OUTPUT else subprocess.DEVNULL)
                    ))
                    
        elif self.architecture == 'condor':
            
            self.condor_deploy()

        for i_worker in range(self.n_workers):

            worker_id_min = i_worker*self.n_cores_per_worker
            worker_id_max = (i_worker+1)*self.n_cores_per_worker-1

            for worker_id in range(worker_id_min, worker_id_max+1):
                self.workers[worker_id] = {
                    'hook' : self.all_worker_hooks[-1],
                    'status' : 'PENDING'
                }
                await self.send_job( worker_id, 'ping')

        await self.update_status()
    
        asyncio.create_task(self.follow_workers())

        # Wait for at least one worker to become active
        if wait_for_one_worker:
            available_worker_id = await self.available_worker_ids.get()
            await self.available_worker_ids.put(available_worker_id)

        asyncio.create_task(self.job_submitter())

    def condor_deploy(self):
        
        if not os.path.exists(pjoin(self.run_workspace,'run_%d_condor_logs'%self.run_id)):
            os.mkdir(pjoin(self.run_workspace,'run_%d_condor_logs'%self.run_id))

        with open(pjoin(self.run_workspace,'run_%d_condor_worker_arguments.txt'%self.run_id),'w') as f:
            f.write('\n'.join(
                '%d, %d'%(i_worker*self.n_cores_per_worker, (i_worker+1)*self.n_cores_per_worker-1) for i_worker in range(self.n_workers)
            ))
        
        libscsdir_path = pjoin(alphaloop_basedir,'libraries','scs','out','libscsdir.so')
        if not os.path.exists(libscsdir_path):
            raise HavanaIntegratorError("Could not find libscsdir.so at '%s'. Make sure it is present and compiled."%libscsdir_path)

        with open(pjoin(self.run_workspace,'run_%d_condor_submission.sub'%self.run_id),'w') as f:
            f.write(
"""executable         = %(python)s
arguments             = %(worker_script)s --run_id %(run_id)d --worker_id_min $(worker_id_min) --worker_id_max $(worker_id_max) --workspace_path %(workspace)s --timeout %(timeout)s
environment           = "LD_PRELOAD=%(libscsdir_path)s"
output                = %(workspace)s/run_%(run_id)d_condor_logs/worker_$(worker_id_min)_$(worker_id_max).out
error                 = %(workspace)s/run_%(run_id)d_condor_logs/worker_$(worker_id_min)_$(worker_id_max).err
log                   = %(workspace)s/run_%(run_id)d_condor_logs/worker_$(worker_id_min)_$(worker_id_max).log
RequestCpus           = %(n_cpus_per_worker)d
RequestMemory         = %(requested_memory_in_MB)d
RequestDisk           = DiskUsage
should_transfer_files = Yes
when_to_transfer_output = ON_EXIT
+JobFlavour           = "%(job_flavour)s"
queue worker_id_min,worker_id_max from %(workspace)s/run_%(run_id)d_condor_worker_arguments.txt
"""%{
        'python' : sys.executable,
        'libscsdir_path' : libscsdir_path,
        'run_id' : self.run_id,
        'worker_script' : pjoin(alphaloop_basedir,'alpha_loop','integrator','worker.py'),
        'workspace' : self.run_workspace,
        'timeout' : '%.2f'%(self._JOB_CHECK_FREQUENCY*15.),
        'job_flavour' : self.cluster_options['job_flavour'],
        'n_cpus_per_worker' : self.n_cores_per_worker,
        'requested_memory_in_MB' : self.cluster_options['RAM_required']
    }
            )

        submit_proc = subprocess.Popen(['condor_submit',pjoin(self.run_workspace,'run_%d_condor_submission.sub'%self.run_id)], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        submit_stdout, submit_stderr = submit_proc.communicate()
        if submit_proc.returncode != 0:
            raise HavanaIntegratorError('Could not successfully submit condor jobs. Error:\n%s'%(submit_stderr.decode('utf-8')))
        else:
            submit_stdout_decoded = submit_stdout.decode('utf-8')
            matched_stdout = re.match(r'(.|\n)*cluster\s(?P<cluster_id>(\d+))(.|\n)*',submit_stdout_decoded)
            if matched_stdout:
                cluster_id = int(matched_stdout.group('cluster_id'))
                self.all_worker_hooks.append(cluster_id)
                logger.info("A total of %d cluster workers were submitted to condor cluster id %d"%(self.n_workers,cluster_id))
            else:
                raise HavanaIntegratorError('Could not interpret response from condor cluster:\n%s'%submit_stdout_decoded)

    def condor_kill(self):
        if len(self.all_worker_hooks)==0:
            return
        try:
            kill_proc = subprocess.Popen(['condor_rm',]+[str(hook) for hook in self.all_worker_hooks], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            kill_stdout, kill_stderr = kill_proc.communicate()
            if kill_proc.returncode != 0:
                logger.warning('Could not successfully clean condor jobs. Error:\n%s'%(kill_stderr.decode('utf-8')))
            else:
                logger.info('Successfully clean condor jobs:\n%s'%(kill_stdout.decode('utf-8')))
        except Exception as e:
            logger.warning("Failed to clean condor run. Error: %s"%str(e))

    async def submit_job(self, payload, blocking=False):

        if blocking:
            available_worker_id = await self.available_worker_ids.get()
            await self.send_job( available_worker_id, payload)
        else:
            self.n_jobs_in_queue += 1
            if self.debug: logger.info("Adding a new job to the queue (now pending=%d)..."%self.n_jobs_in_queue)
            await self.jobs_queue.put(payload)
            await self.update_status()

    async def job_submitter(self):
        while self.active:
            try:
                payload = await self.jobs_queue.get()
                if self.debug: logger.info("Job submitter picked up a new job from the queue and it is now waiting for an available worker...")
                available_worker_id = await self.available_worker_ids.get()
                self.n_jobs_in_queue -=1
                if self.debug: logger.info("Job submitter found a free worker (#%d) to dispatch the job (now pending=%d)..."%(available_worker_id,self.n_jobs_in_queue))
                await self.send_job( available_worker_id, payload)
                await self.update_status()
            except KeyboardInterrupt as e:
                break

    def terminate(self):
        
        logger.info("Cleaning up cluster run #%d..."%self.run_id)
        self.active = False

        if self.architecture == 'local':
            for hook in self.all_worker_hooks:
                try:
                    hook.kill()
                except Exception as e:
                    pass
            try:
                self.worker_monitoring_file.close()
                if not self.keep:
                    if os.path.exists(self.work_monitoring_path):
                        os.remove(self.work_monitoring_path)
            except Exception as e:
                pass

        elif self.architecture == 'condor':
            if not self.keep:
                if os.path.exists(pjoin(self.run_workspace,'run_%d_condor_logs'%self.run_id)):
                    shutil.rmtree(pjoin(self.run_workspace,'run_%d_condor_logs'%self.run_id))
                
                if os.path.exists(pjoin(self.run_workspace,'run_%d_condor_submission.sub'%self.run_id)):
                    os.remove(pjoin(self.run_workspace,'run_%d_condor_submission.sub'%self.run_id))
                if os.path.exists(pjoin(self.run_workspace,'run_%d_condor_worker_arguments.txt'%self.run_id)):
                    os.remove(pjoin(self.run_workspace,'run_%d_condor_worker_arguments.txt'%self.run_id))

            self.condor_kill()

        io_files = \
            [f for f in glob.glob(pjoin(self.run_workspace, 'run_%d_job_*.pkl'%self.run_id))]+\
            [f for f in glob.glob(pjoin(self.run_workspace, 'run_%d_job_*.done'%self.run_id))]+\
            [f for f in glob.glob(pjoin(self.run_workspace, 'run_%d_job_*.bin'%self.run_id))]+\
            [f for f in glob.glob(pjoin(self.run_workspace, 'run_%d_job_*.yaml'%self.run_id))]
        for io_file in io_files:
            os.remove(io_file)

class HavanaIntegrator(integrators.VirtualIntegrator):
    """ Steering of the havana integrator """
    
    _SUPPORTED_CLUSTER_ARCHITECTURES = [ 'local', 'condor' ]
    _DEBUG = False
    _USE_HAVANA_MOCKUP = False

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
                 n_cores_per_worker = None,
                 batch_size = 1e6,
                 cluster_type = 'local',
                 local_options = None,
                 condor_options = None,
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
                 keep=False,
                 run_id = None,
                 havana_starting_n_bins = 128,
                 havana_n_points_min = 1000,
                 havana_learning_rate = 1.5,
                 havana_bin_increase_factor_schedule = None,
                 show_selected_phase_only = False,
                 show_all_information_for_all_integrands = False,
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

        self.run_id = run_id
        self.run_workspace = pjoin(run_workspace,'run_%d'%self.run_id)
        
        # Clean up of possibly previously crashed run
        io_files = \
            [f for f in glob.glob(pjoin(self.run_workspace, 'run_%d_job_*.pkl'%self.run_id))]+\
            [f for f in glob.glob(pjoin(self.run_workspace, 'run_%d_job_*.done'%self.run_id))]+\
            [f for f in glob.glob(pjoin(self.run_workspace, 'run_%d_job_*.bin'%self.run_id))]+\
            [f for f in glob.glob(pjoin(self.run_workspace, 'run_%d_job_*.yaml'%self.run_id))]
        for io_file in io_files:
            os.remove(io_file)

        # Setup a trashcan
        self.trashcan = pjoin(self.run_workspace,'run_%s_trashcan'%self.run_id)
        if os.path.exists(self.trashcan):
            shutil.rmtree(self.trashcan)
        if os.path.exists(self.trashcan+'_being_removed'):
            shutil.rmtree(self.trashcan+'_being_removed')
        os.mkdir(self.trashcan)

        self.havana_optimize_on_variance = havana_optimize_on_variance
        self.havana_max_prob_ratio = havana_max_prob_ratio
        self.show_SG_grid = show_SG_grid
        self.show_channel_grid = show_channel_grid
        self.show_grids_sorted_by_variance = show_grids_sorted_by_variance
        self.show_selected_phase_only = show_selected_phase_only
        self.show_all_information_for_all_integrands = show_all_information_for_all_integrands
        self.dump_havana_grids = dump_havana_grids
        self.fresh_integration = fresh_integration

        self.n_start = n_start
        self.n_max = n_max
        self.n_increase = n_increase

        self.havana_starting_n_bins = havana_starting_n_bins
        self.havana_n_points_min = havana_n_points_min
        self.havana_learning_rate = havana_learning_rate
        if havana_bin_increase_factor_schedule is None:
            # This is a reasonable default: double the number of bins every 5 iterations for the first 20
            self.havana_bin_increase_factor_schedule = tuple(sum([[1,1,1,1,2],]*4,[]))
        else:
            self.havana_bin_increase_factor_schedule = havana_bin_increase_factor_schedule

        self.keep = keep

        self.cluster_status_summary = None

        logger.setLevel(verbosity)

        # Only useful later in run_interface for final report
        self.tot_func_evals = 0

        if n_workers is None:
            self.n_workers = multiprocessing.cpu_count()
        else:
            self.n_workers = n_workers
        if n_cores_per_worker is None:
            self.n_cores_per_worker = 1
        else:
            self.n_cores_per_worker = n_cores_per_worker
        self.n_cpus = self.n_workers*self.n_cores_per_worker

        self.local_options = {'n_workers' : self.n_workers}
        if local_options is not None:
            self.local_options.update(local_options)

        self.condor_options = {
            'memory'  : '2000MB',
            'disk'    : '1000MB',
            'IO_port' : 8786,
            'n_cores_per_worker' : 1,
            'job_flavour' : 'tomorrow',
            'RAM_required' : 1024
        }

        if condor_options is not None:
            self.condor_options.update(condor_options)

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

        self.seed = seed
        if not self._USE_HAVANA_MOCKUP and self.seed is None:
            # We must use a definite seed when not using havana mockup.
            self.seed = random.randint(1,1001)

        self.havana = None

        super(HavanaIntegrator, self).__init__(integrands, **opts)

    def update_status( self,
        start_time, n_jobs_total_completed, n_submitted, n_done,
        n_tot_points, n_points_for_this_iteration, current_n_points, cumulative_IO_time,
        cumulative_processing_time, cumulative_job_time, current_iteration_number,
        n_jobs_for_this_iteration=None
    ):

        curr_integration_time = time.time()-start_time

        # Update the run stat line and monitoring canvas
        monitoring_report = []
        if self.cluster_status_summary is not None:
            monitoring_report.append( '\n'.join('| %s'%line for line in self.cluster_status_summary.split('\n')) )

        monitoring_report.append( '| Jobs completed: overall = %d (avg %s), and for this iteration = %d/%d/%d on %d cpus.\n| Total n_pts: %s%.1fM%s ( %s%s ms / pt%s on one core, %s %s pts / s%s overall ). n_pts for this iteration %s#%d: %.2fM/%.2fM (%.1f%%)%s.'%(
            n_jobs_total_completed, '%.3g min/job'%((cumulative_job_time/n_jobs_total_completed)/60.) if n_jobs_total_completed>0 else 'N/A', n_done, n_submitted, n_jobs_for_this_iteration if n_jobs_for_this_iteration is not None else n_submitted,
            self.n_cpus, bcolors.GREEN, n_tot_points/(1.e6),bcolors.END, bcolors.GREEN, '%.3g'%((curr_integration_time/n_tot_points)*self.n_cpus*1000.) if n_tot_points>0 else 'N/A', bcolors.END,
            bcolors.GREEN, ('%.1fK'%(n_tot_points/curr_integration_time/1000.) if n_tot_points>0 else 'N/A'), bcolors.END, bcolors.BLUE,
            current_iteration_number,
            n_points_for_this_iteration/(1.e6), current_n_points/(1.e6), (float(n_points_for_this_iteration)/current_n_points)*100., bcolors.END
        ))
        wall_time = curr_integration_time
        IO_time = (cumulative_IO_time / wall_time)*100.
        processing_time = (cumulative_processing_time / wall_time)*100.
        cpu_hours = (cumulative_job_time / 3600.)
        parallelisation_efficiency = (cumulative_job_time / (wall_time*self.n_cpus) )*100.
        monitoring_report.append('| Wall time: %.3f h, IO: %.2f%%, processing: %.2f%%, jobs: %.3f CPU-hours, efficiency: %s%.2f%%%s'%(
            wall_time/3600., IO_time, processing_time, cpu_hours, bcolors.GREEN if parallelisation_efficiency > 90. else bcolors.RED, parallelisation_efficiency, bcolors.END
        ))
        monitoring_report.append('')
        monitoring_report.append( self.havana.get_summary() )

        if self.show_SG_grid:
            grid_summary = self.havana.get_grid_summary( 
                sort_by_variance=self.show_grids_sorted_by_variance, 
                show_channel_grid= (self.MC_over_channels and self.show_channel_grid),
                show_selected_phase_only = self.show_selected_phase_only,
                show_all_information_for_all_integrands = self.show_all_information_for_all_integrands
            )
            if grid_summary != '':
                monitoring_report.append( grid_summary )
        monitoring_report.append('')

        if self.canvas is None:
            self.canvas = utils.Canvas( '\n'.join(monitoring_report), stream=self.stream_monitor, overwrite=True )
        else:
            self.canvas.print('\n'.join(monitoring_report))

    async def process_job_result(self, job_result):
        if self.exit_now:
            return
        try:
            job_ID, run_time, results = job_result
            if self._DEBUG: logger.info("Processing result of %d samples from job #%d."%(len(results), job_ID))

            self.cumulative_job_time += run_time

            t2 = time.time()
            if self._DEBUG: logger.info("Accumulating new results from job %d into havana..."%job_ID)

            havana_constructor_args = results['havana_updater']
            havana_constructor_args['integrands'] = self.integrands

            if self._USE_HAVANA_MOCKUP:
                havana_updater = HavanaMockUp(**havana_constructor_args)
            else:
                havana_updater = Havana(**havana_constructor_args)
                # Clean-up job output
                for grid_file_path in havana_constructor_args['flat_record']['grid_files']:
                    fast_remove(grid_file_path, self.trashcan)
            
            self.havana.accumulate_from_havana_grid(havana_updater)
            
            if self._USE_HAVANA_MOCKUP:
                max_prob_ratio_for_this_iteration = max( self.havana_max_prob_ratio / max( 100. - 10.*(self.current_iteration-1) , 1.), 3.)
                if self._DEBUG: logger.info("Syncing Havana grids with the new results...")
                self.havana.sync_grids( max_prob_ratio = max_prob_ratio_for_this_iteration )
                if self._DEBUG: logger.info("Done updating Havana grids...")

                for index in range(0,len(self.integrands)):
                    self.integrands[index].n_evals.value += results['sample_statistics'][index]['n_evals']
                    self.integrands[index].n_evals_failed.value += results['sample_statistics'][index]['n_evals_failed']
                    self.integrands[index].n_zero_evals.value += results['sample_statistics'][index]['n_zero_evals']
                    if (results['sample_statistics'][index]['max_eval_positive_xs'] is not None) and results['sample_statistics'][index]['max_eval_positive'] > self.integrands[index].max_eval_positive.value:
                        self.integrands[index].max_eval_positive.value = results['sample_statistics'][index]['max_eval_positive']
                        self.integrands[index].max_eval_positive_xs[:] = results['sample_statistics'][index]['max_eval_positive_xs']
                    if (results['sample_statistics'][index]['max_eval_negative_xs'] is not None) and results['sample_statistics'][index]['max_eval_negative'] < self.integrands[index].max_eval_negative.value:
                        self.integrands[index].max_eval_negative.value = results['sample_statistics'][index]['max_eval_negative']
                        self.integrands[index].max_eval_negative_xs[:] = results['sample_statistics'][index]['max_eval_negative_xs']

            self.n_points_for_this_iteration += havana_updater.n_points
            self.n_tot_points += havana_updater.n_points


            self.cumulative_processing_time = time.time() - t2

            self.n_jobs_awaiting_completion -= 1

        except KeyboardInterrupt:
            self.exit_now = True

    async def process_al_cluster_status_update(self, new_status):
        
        self.cluster_status_summary = 'Status of %d workers: %d pending, %d available, %s%d%s running. A total of %d jobs are pending.'%(
            new_status['PENDING']+new_status['FREE']+new_status['RUNNING'],
            new_status['PENDING'], new_status['FREE'], bcolors.GREEN, new_status['RUNNING'], bcolors.END,
            self.al_cluster.n_jobs_in_queue
        )
        if self.canvas is not None and self.stream_monitor:
            self.canvas.print('| %s'%self.cluster_status_summary)

    def handle_exception(self, loop, context):
        self.exit_now = True
        logger.critical("The following exception happened in a coroutine:\n%s\nAborting integration now."%str(context))

    @staticmethod
    async def async_rm(rm_target, is_file=False):
        if os.path.isdir(rm_target):
            shutil.rmtree(rm_target)
        else:
            os.remove(rm_target)

    async def async_integrate(self):

        n_dimensions = len(self.integrands[0].dimensions.get_continuous_dimensions())
        if not all(len(integrand.dimensions.get_continuous_dimensions())==n_dimensions for integrand in self.integrands):
            raise HavanaIntegratorError("In Havana implementations, all integrands must have the same dimension.")

        if not self.MC_over_SGs and self.MC_over_channels:
            raise HavanaIntegratorError("Havana cannot perform a Monte-Carlo over channels but not over supergraphs.")

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
            total_number_of_integration_channels = sum(
                len(self.all_supergraphs[SG['name']]['multi_channeling_bases']) if len(self.all_supergraphs[SG['name']]['optimal_channel_ids'])==0 
                    else len(self.all_supergraphs[SG['name']]['optimal_channel_ids']) for SG in self.cross_section_set['topologies']
            )
        else:
            n_channels_per_SG = []
            for (i_SG, SG_name) in SG_ids:
                n_channels_per_SG.append(
                    len(self.all_supergraphs[SG_name]['multi_channeling_bases']) if len(self.all_supergraphs[SG_name]['optimal_channel_ids'])==0 
                    else len(self.all_supergraphs[SG_name]['optimal_channel_ids'])
                )
        
        if SG_ids is None:
            total_number_of_integration_channels = sum(
                len(self.all_supergraphs[SG['name']]['multi_channeling_bases']) if len(self.all_supergraphs[SG['name']]['optimal_channel_ids'])==0 
                    else len(self.all_supergraphs[SG['name']]['optimal_channel_ids']) for SG in self.cross_section_set['topologies']
            )
        else:
            total_number_of_integration_channels = 0
            for (i_SG, SG_name) in SG_ids:
                total_number_of_integration_channels += (
                    len(self.all_supergraphs[SG_name]['multi_channeling_bases']) if len(self.all_supergraphs[SG_name]['optimal_channel_ids'])==0 
                    else len(self.all_supergraphs[SG_name]['optimal_channel_ids'])
                )

        if self._USE_HAVANA_MOCKUP:
            self.havana = HavanaMockUp(
                n_dimensions, 
                self.integrands,
                SG_ids=SG_ids, 
                n_channels_per_SG=n_channels_per_SG, 
                target_result=self.target_result, 
                seed=self.seed,
                reference_integrand_index = 0,
                phase = self.phase,
                grid_file = pjoin(self.run_workspace,'havana_grid.yaml'),
                optimize_on_variance = self.havana_optimize_on_variance,
                max_prob_ratio = self.havana_max_prob_ratio,
                fresh_integration = self.fresh_integration
            )
        else:
            self.havana = Havana(
                n_dimensions, 
                self.integrands,
                SG_ids=SG_ids, 
                n_channels_per_SG=n_channels_per_SG, 
                target_result=self.target_result, 
                seed=self.seed,
                reference_integrand_index = 0,
                phase = self.phase,
                grid_file = pjoin(self.run_workspace,'havana_grid.yaml'),
                optimize_on_variance = self.havana_optimize_on_variance,
                max_prob_ratio = self.havana_max_prob_ratio,
                fresh_integration = self.fresh_integration,
                n_bins = self.havana_starting_n_bins,
                n_points_min = self.havana_n_points_min,
                learning_rate = self.havana_learning_rate,
                bin_increase_factor_schedule = self.havana_bin_increase_factor_schedule,
                alpha_loop_path = alphaloop_basedir,
            )

        # Now perform the integration
        
        logger.info("Staring alphaLoop integration with Havana and cluster '%s' as run #%d"%(self.cluster_type, self.run_id))
        logger.info("Number of SG considered: %s (%s)"%(
            len(self.cross_section_set['topologies']) if SG_ids is None else len(SG_ids), "Monte-Carlo'ed over" if self.MC_over_SGs else 'explicitly summed per sample'
        ))
        logger.info("Total number of integration channels considered: %s (%s)"%(
            total_number_of_integration_channels, "Monte-Carlo'ed over" if self.MC_over_channels else 'explicitly summed per sample'
        ))
        logger.info("%s%s%s"%(bcolors.GREEN, "Lean back and enjoy...", bcolors.END))

        self.current_iteration = 1
        current_n_points = self.n_start
        current_step = self.n_increase
        self.n_tot_points = 0

        integrands_constructor_args = [
            integrand.get_constructor_arguments() for integrand in self.integrands
        ]

        job_ID = 0

        self.canvas = None

        self.n_jobs_awaiting_completion = 0

        n_jobs_total_completed = 0
        self.cumulative_job_time = 0.
        self.cumulative_processing_time = 0.
        self.cumulative_IO_time = 0.
        self.exit_now = False

        self.n_points_for_this_iteration = 0

        cluster_options = None
        if self.cluster_type == 'condor':
            cluster_options = self.condor_options
        self.al_cluster = AL_cluster(
            self.cluster_type, self.n_workers, self.n_cores_per_worker, self.run_workspace, self.process_job_result, self.run_id,
            monitor_callback=self.process_al_cluster_status_update, cluster_options=cluster_options, keep=self.keep, debug=self._DEBUG, trashcan=self.trashcan
        )

        try:
            asyncio.get_event_loop().set_exception_handler(self.handle_exception)

            logger.info("Now deploying cluster '%s' and waiting for at least one job to become active..."%self.cluster_type)
            await self.al_cluster.deploy()
            logger.info("Cluster '%s' now active."%self.cluster_type)
            start_time = time.time()

            while True:

                n_remaining_points = current_n_points

                n_jobs_for_this_iteration = math.ceil(current_n_points/float(self.batch_size))
                self.n_points_for_this_iteration = 0
                
                sampling_grid_constructor_arguments = self.havana.get_constructor_arguments()

                n_submitted = 0
                n_done = 0
                last_n_done = 0

                logger.info("Now running iteration %s#%d%s"%(bcolors.GREEN, self.current_iteration, bcolors.END))
                self.update_status(
                    start_time, n_jobs_total_completed, n_submitted, n_done,
                    self.n_tot_points, self.n_points_for_this_iteration, current_n_points, self.cumulative_IO_time,
                    self.cumulative_processing_time, self.cumulative_job_time, self.current_iteration, n_jobs_for_this_iteration=n_jobs_for_this_iteration
                )

                show_waiting = False
                while True:
                    if self.exit_now:
                        break

                    timeout = 0.3
                    # Add a threshold of 5% of workers to insure smooth rollover
                    if n_remaining_points > 0 and self.n_jobs_awaiting_completion < self.n_cpus+max(int(self.n_cpus/20.),2):

                        this_n_points = min(n_remaining_points, self.batch_size)
                        n_remaining_points -= this_n_points   
                        t2 = time.time()           
                        job_ID += 1
                        job_payload = {
                                'job_id' : job_ID,
                                'integrands_constructor_args' : integrands_constructor_args,
                                'SG_ids_specified' : self.MC_over_SGs,
                                'channels_specified' : self.MC_over_channels
                            }

                        havana_grid_constructor_arguments = dict(sampling_grid_constructor_arguments)
                        if self.seed:
                            havana_grid_constructor_arguments['seed'] = int(self.seed+100000*job_ID)
                        else:
                            havana_grid_constructor_arguments['seed'] = None
                        havana_grid_constructor_arguments['fresh_integration'] = False
                        job_payload['phase'] = self.phase
                        job_payload['havana_grid'] = havana_grid_constructor_arguments
                        job_payload['havana_mockup'] = self._USE_HAVANA_MOCKUP
                        job_payload['n_points_to_sample'] = this_n_points

                        self.cumulative_processing_time = time.time() - t2
                        t0 = time.time()
                        if self._DEBUG: logger.info("Submitting job %d."%job_ID)
                        await self.al_cluster.submit_job(
                            job_payload,
                            blocking=False
                        )
                        if self._DEBUG: logger.info("Done with submission of job %d."%job_ID)
                        self.cumulative_IO_time += time.time()-t0

                        self.n_jobs_awaiting_completion += 1
                        n_submitted += 1
                        timeout = 0.01
                        show_waiting = True
                        self.update_status(
                            start_time, n_jobs_total_completed, n_submitted, n_done,
                            self.n_tot_points, self.n_points_for_this_iteration, current_n_points, self.cumulative_IO_time,
                            self.cumulative_processing_time, self.cumulative_job_time, self.current_iteration, n_jobs_for_this_iteration=n_jobs_for_this_iteration)

                    else:
                        if show_waiting:
                            show_waiting = False
                            if self._DEBUG: logger.info("Now waiting for jobs to complete and freeing workers...")


                    await asyncio.sleep(timeout)
                    n_done = n_submitted - self.n_jobs_awaiting_completion

                    if last_n_done != n_done:
                        n_jobs_total_completed += n_done-last_n_done
                        last_n_done = n_done

                        self.update_status(
                            start_time, n_jobs_total_completed, n_submitted, n_done,
                            self.n_tot_points, self.n_points_for_this_iteration, current_n_points, self.cumulative_IO_time,
                            self.cumulative_processing_time, self.cumulative_job_time, self.current_iteration, n_jobs_for_this_iteration=n_jobs_for_this_iteration)

                    if n_remaining_points == 0 and self.n_jobs_awaiting_completion == 0:
                        break
                
                shutil.move(self.trashcan, self.trashcan+'_being_removed')
                os.mkdir(self.trashcan)
                asyncio.create_task(self.async_rm(self.trashcan+'_being_removed'))

                if self.exit_now:
                    break

                t2 = time.time()
                if self._USE_HAVANA_MOCKUP and self._DEBUG: logger.info("Updating Havana grids for iteration #%d..."%self.current_iteration)
                self.havana.update_iteration_count(dump_grids=self.dump_havana_grids)
                if not self._USE_HAVANA_MOCKUP:
                    for i_itg, integrand in enumerate(self.integrands):
                        grid_index = ( i_itg*2 if self.phase=='real' else ( i_itg*2 + 1 ) )
                        avg, err, chi_sq, max_eval_negative, max_eval_positive, n_evals, n_zero_evals = self.havana.havana_grids[grid_index].get_current_estimate()
                        integrand.n_evals.value = n_evals
                        integrand.n_evals_failed.value = 0
                        integrand.n_zero_evals.value = n_zero_evals
                        integrand.max_eval_positive.value = max_eval_positive
                        #TODO propagate this info too
                        integrand.max_eval_positive_xs[:] = [0.,]*len(integrand.max_eval_positive_xs)
                        integrand.max_eval_negative.value = max_eval_negative
                        #TODO propagate this info too
                        integrand.max_eval_negative_xs[:] = [0.,]*len(integrand.max_eval_negative_xs)

                if self._USE_HAVANA_MOCKUP and self._DEBUG: logger.info("Done updating Havana grids in %.3g ms."%((time.time()-t2)*1000.))
                self.cumulative_processing_time += time.time() - t2

                if self.n_iterations is not None and self.current_iteration >= self.n_iterations:
                    logger.info("Max number of iterations %d reached."%self.n_iterations) 
                    break               

                if self.n_max is not None and self.n_tot_points >= self.n_max:
                    logger.info("Max number of sample points %d reached."%self.n_max)
                    break

                res_int, res_error = self.havana.get_current_estimate()
                if self.accuracy_target is not None and (res_error/abs(res_int) if res_int!=0. else 0.)<self.accuracy_target:
                    logger.info("Target accuracy of %.2g reached."%self.accuracy_target)
                    break

                current_n_points += current_step
                current_step += self.n_increase
                self.current_iteration += 1

                # Keep printout of last step of this iteration.
                self.canvas = None

        except KeyboardInterrupt as e:
            logger.warning("Aborting Havana integration now.")

        except Exception as e:
            logger.critical("Exception '%s' occurred during run:\n%s\nAborting integration now."%(str(e),str(traceback.format_exc())))

        if self.exit_now:
            logger.warning("Aborting Havana integration now.")

        if not self._USE_HAVANA_MOCKUP:
            for i_itg, integrand in enumerate(self.integrands):
                grid_index = ( i_itg*2 if self.phase=='real' else ( i_itg*2 + 1 ) )
                avg, err, chi_sq, max_eval_negative, max_eval_positive, n_evals, n_zero_evals = self.havana.havana_grids[grid_index].get_current_estimate()
                integrand.n_evals.value = n_evals
                integrand.n_evals_failed.value = 0
                integrand.n_zero_evals.value = n_zero_evals
                integrand.max_eval_positive.value = max_eval_positive
                #TODO propagate this info too
                integrand.max_eval_positive_xs[:] = [0.,]*len(integrand.max_eval_positive_xs)
                integrand.max_eval_negative.value = max_eval_negative
                #TODO propagate this info too
                integrand.max_eval_negative_xs[:] = [0.,]*len(integrand.max_eval_negative_xs)

    def integrate(self):
        """ Return the final integral and error estimates."""

        final_res = None

        logger.info("Now deploying the running cluster...")
        
        if self.cluster_type in ['local','condor']:

            self.al_cluster = None
            try:
                asyncio.get_event_loop().run_until_complete(self.async_integrate())
            except KeyboardInterrupt:
                logger.info("Havana integration aborted by user.")

            if self.al_cluster is not None:
                self.al_cluster.terminate()

            self.tot_func_evals = self.havana.get_n_evals()
            final_res = self.havana.get_current_estimate()

        else:
            raise HavanaIntegratorError("Cluster architecture %s not supported."%self.cluster_type)

        if not self.keep:
            try:
                shutil.rmtree(self.run_workspace)
            except:
                pass

        return final_res
