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
import socket
import uuid
import tempfile
from pprint import pprint, pformat

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
from alpha_loop.integrator.worker import HavanaIntegrandWrapper, ALStandaloneIntegrand, HavanaMockUp, Havana, run_job

import asyncio
import subprocess
import yaml
from yaml import Loader


try:
    import redis
    from redis import Redis
    import rq
except ImportError as e:
    # This is fine if not using the redis-based parallelisation implementation
    pass

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


trash_counter=0
def fast_remove(path, trashcan):
    global trash_counter
    if trashcan is None:
        if os.path.isdir(path):
            shutil.rmtree(path)
        else:
            os.remove(path)
    else:
        trash_counter += 1
        shutil.move(path, pjoin(trashcan, '%d_%s'%(trash_counter,os.path.basename(path))))

class HavanaIntegratorError(Exception):
    pass

class AL_cluster(object):

    _SUPPORTED_CLUSTER_ARCHITECTURES = ['local', 'condor']
    _JOB_CHECK_FREQUENCY = 0.3
    _FORWARD_WORKER_OUTPUT = True

    def __init__(self, architecture, n_workers, n_cores_per_worker, run_workspace, process_results_callback, run_id, 
                monitor_callback=None, cluster_options=None, keep=False, debug=False, trashcan=None, 
                use_redis=False, redis_port=8786, redis_hostname=None, redis_max_job_time=None, external_redis=False,
                redis_queue_name = None, worker_resources_path = None, standalone_workers = None
                ):

        self.n_workers = n_workers
        self.n_cores_per_worker = n_cores_per_worker
        self.run_workspace = run_workspace
        self.available_worker_ids = asyncio.Queue()
        self.n_jobs_in_queue = 0
        self.jobs_queue = asyncio.Queue()
        self.workers = { }
        self.all_worker_hooks = []
        self.active = True
        self.worker_resources_path = worker_resources_path
        self.standalone_workers = standalone_workers

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
        
        self.use_redis = use_redis
        self.redis_server_process = None
        self.redis_port = redis_port
        self.external_redis = external_redis
        self.rq_path = None
        self.redis_initialization = None
        self.redis_queue = None
        self.redis_connection = None
        self.redis_submitter_hostname = redis_hostname
        self.redis_queue_name = redis_queue_name
        self.jobs_for_current_iteration = {}
        self.redis_max_job_time = redis_max_job_time
        self.n_redis_job_failed_for_this_iteration = 0
        if self.use_redis:
            if self.redis_submitter_hostname is None:
                self.redis_submitter_hostname = socket.gethostname()
            self.redis_initialization = asyncio.Event()
            self.rq_path = shutil.which('rq')
            if self.rq_path is None:
                raise HavanaIntegratorError("The executable 'rq' could not be found within the current PATH environment. Make sure to install it with pip install rq.")

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
        
        if self.use_redis:
            #workers = rq.Worker.all(connection=self.redis_connection)
            workers = rq.worker.Worker.all(queue=self.redis_queue)
            # Possible worker states are suspended, started, busy and idle
            worker_statuses = [worker.state for worker in workers]
            n_idle_workers = worker_statuses.count('idle')
            n_running_workers = worker_statuses.count('busy')
            new_status = {
                'FREE' : n_idle_workers,
                'PENDING' : len(worker_statuses)-n_idle_workers-n_running_workers,
                'RUNNING' : n_running_workers
            }
            new_status['jobs'] = {}
            new_status['jobs']['started'] = self.redis_queue.started_job_registry.count
            new_status['jobs']['deferred'] = self.redis_queue.deferred_job_registry.count
            new_status['jobs']['finished_and_not_processed'] = self.redis_queue.finished_job_registry.count
            new_status['jobs']['failed'] = self.n_redis_job_failed_for_this_iteration
            #new_status['jobs']['failed'] = self.redis_queue.failed_job_registry.count
            new_status['jobs']['scheduled'] = self.redis_queue.scheduled_job_registry.count
        else:
            all_statuses = [ v['status'] for v in self.workers.values() ]
            new_status = {
                'FREE' : all_statuses.count('FREE'),
                'PENDING' : all_statuses.count('PENDING'),
                'RUNNING' : all_statuses.count('RUNNING')
            }

        if self.current_status is None or new_status!=self.current_status:
            self.current_status = new_status
            if self.monitor_callback is not None:
                cluster_status_str = 'Status of %d workers: %d pending, %d available, %s%d%s running. '%(
                            new_status['PENDING']+new_status['FREE']+new_status['RUNNING'],
                            new_status['PENDING'], new_status['FREE'], bcolors.GREEN, new_status['RUNNING'], bcolors.END   
                        )
                if 'jobs' in new_status:
                    n_jobs_tot = sum(new_status['jobs'].values())
                    cluster_status_str += '\nStatus of %d jobs: %s%d%s scheduled, %s%d%s submitted for this iteration, %s%d%s started, %s%d%s finished and not processed'%(
                        n_jobs_tot, 
                        bcolors.BLUE, new_status['jobs']['scheduled'], bcolors.END,
                        bcolors.BLUE, len(self.jobs_for_current_iteration), bcolors.END,
                        bcolors.GREEN, new_status['jobs']['started'], bcolors.END,
                        bcolors.RED, new_status['jobs']['finished_and_not_processed'], bcolors.END
                    )
                    if new_status['jobs']['deferred'] != 0:
                        cluster_status_str += ', %s%d%s deferred'%(bcolors.BLUE, new_status['jobs']['deferred'], bcolors.END)
                    if new_status['jobs']['failed'] != 0:
                        cluster_status_str += ', %s%d%s failed for this iteration'%(bcolors.RED, new_status['jobs']['failed'], bcolors.END)
                else:
                    cluster_status_str += 'A total of %d jobs are pending.'%(self.n_jobs_in_queue)
                await self.monitor_callback(cluster_status_str)

    async def send_redis_job(self, payload):

        # result_ttl is how long (in seconds) to keep the job (if successful) and its results
        # ttl can also be supplied, and it if it is it provides a maximum queued time (in seconds) of the job before it's discarded.
        # job = rq.job.Job.create( func=run_job, result_ttl=3600, connection=self.redis_connection, retry=rq.Retry(3),
        #   args=(payload,),
        #   kwargs={}
        # )
        #self.redis_queue.enqueue_job( job )

        # One could add retry=rq.Retry(3) below to add an automated retry, but that is not necessary here

        sent_job = self.redis_queue.enqueue( 
            run_job, result_ttl=3600, connection=self.redis_connection, job_timeout=self.redis_max_job_time,
            args=(payload,self.run_id),
            kwargs={}
        )
        if 'job_id' in payload:
            self.jobs_for_current_iteration[sent_job.id] = {
                'alpha_loop_job_id' : payload['job_id'],
                'sent_time' : time.time()
            }
        # if payload=='ping':
        #     print("Just sent REDIS job with payload: %s"%pformat(payload))
        #     sent_job_record = rq.job.Job.fetch(sent_job.id, connection=self.redis_connection)
        #     print("Sent job record immediately after: %s, %s"%(
        #         sent_job_record.get_status(), str(sent_job_record.result)
        #     ))
        #     await asyncio.sleep(1.)
        #     sent_job_record = rq.job.Job.fetch(sent_job.id, connection=self.redis_connection)
        #     print("Sent job record after one second: %s, %s"%(
        #         sent_job_record.get_status(), str(sent_job_record.result)
        #     ))

        await self.update_status()

    async def send_many_redis_jobs(self, payloads):

        with self.redis_queue.connection.pipeline() as pipe:
            jobs = self.redis_queue.enqueue_many(
                [
                    rq.Queue.prepare_data(
                        run_job, result_ttl=3600, timeout=self.redis_max_job_time,
                        args=(payload,self.run_id),
                        kwargs={}
                    ) for payload in payloads
                ],
                pipeline=pipe
            )
            pipe.execute()

            sent_time = time.time()
            for sent_job, payload in zip(jobs, payloads):
                if 'job_id' in payload:
                    self.jobs_for_current_iteration[sent_job.id] = {
                        'alpha_loop_job_id' : payload['job_id'],
                        'sent_time' : sent_time
                    }

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

        if 'job_id' in payload:
            self.jobs_for_current_iteration[payload['job_id']] = {'alpha_loop_job_id': payload['job_id'], 'sent_time': time.time()}

        await self.update_status()

    def start_new_iteration(self):
        # Make sure to reset the internal redis jobs records when starting a new iteration, so as to be sure
        # to not send the job result from an older iteration for processing.
        if self.use_redis:
            for job_id in self.jobs_for_current_iteration:
                try:
                    job_to_delete = rq.job.Job.fetch(job_id, connection=self.redis_connection)
                    job_to_delete.delete()
                except Exception as e:
                    logger.warning("Failed to delete redis job %s. Exception: %s"%(job_id, str(e)))
            self.n_redis_job_failed_for_this_iteration = 0
        else:
            for _ in range(self.n_jobs_in_queue):
                try:
                    self.jobs_queue.get_nowait()
                    self.jobs_queue.task_done()
                except Exception as e:
                    pass
            self.n_jobs_in_queue = 0

        self.jobs_for_current_iteration = {}
                
    async def follow_redis_workers(self):

        # Possible job statuses are: queued, started, deferred, finished, stopped, scheduled and failed
        last_status = None
        job_ids_already_handled = set([])
        while self.active:

            try:

                if self.architecture == 'local':
                    await asyncio.sleep(self._JOB_CHECK_FREQUENCY)
                else:
                    await asyncio.sleep(self._JOB_CHECK_FREQUENCY*5.)

                await self.update_status()

                # This no longer necessary really
                # if self.redis_max_job_time is not None:
                #     curr_time = time.time()
                #     jobs_to_slow = [ (job_id, job_info) for job_id, job_info in self.jobs_for_current_iteration.items() if (curr_time-job_info['sent_time'])>self.redis_max_job_time ]
                #     for job_id, job_info in jobs_to_slow:
                #         del self.jobs_for_current_iteration[job_id]
                #         # Notify that this job failed
                #         await self.process_results_callback(job_info['alpha_loop_job_id'])
                #         job_ids_already_handled.add(job_id)
                #         try:
                #             job_to_delete = rq.job.Job.fetch(job_id, connection=self.redis_connection)
                #             job_to_delete.delete()
                #         except Exception as e:
                #             logger.warning("Failed to delete redis job %s. Exception: %s"%(job_id, str(e)))

                new_status = {
                    'started' : self.redis_queue.started_job_registry.count,
                    'deferred' : self.redis_queue.deferred_job_registry.count,
                    'finished_and_not_processed' : self.redis_queue.finished_job_registry.count,
                    'failed' : self.redis_queue.failed_job_registry.count,
                    'scheduled' : self.redis_queue.scheduled_job_registry.count
                }

                #pprint(new_status)
                if last_status is not None and new_status==last_status and new_status['finished_and_not_processed']==0 and new_status['failed']==0:
                    continue

                # For now we should only have to act on the finished statuses, the failed ones will be resubmitted automatically
                jobs_to_remove = []
                if new_status['finished_and_not_processed']>0:
                    finished_jobs = rq.job.Job.fetch_many(self.redis_queue.finished_job_registry.get_job_ids(), connection=self.redis_connection)
                    for finished_job in finished_jobs:
                        do_process_result = True
                        job_result = None
                        job_id = finished_job.id
                        if job_id in job_ids_already_handled:
                            logger.warning("The following job ID '%s' appeared on multiple occasions in the list of finished jobs, this should not happen."%str(job_id))
                            do_process_result = False
                        else:
                            finished_job = rq.job.Job.fetch(job_id, connection=self.redis_connection)
                            if finished_job.get_status()!='finished':
                                logger.warning("The following redis job '%s' was in the register of finished jobs but its status is '%s'."%(str(job_id), finished_job.get_status()))
                                do_process_result = False
                            else:
                                job_result = finished_job.result
                                if job_result is None:
                                    logger.warning("The following redis job '%s' was in the register of finished jobs but its result is still None.")
                                    do_process_result = False

                        if do_process_result:
                            if job_result == 'pong':
                                self.redis_initialization.set()
                            else:
                                if job_id in self.jobs_for_current_iteration:
                                    await self.process_results_callback(job_result)
                                job_ids_already_handled.add(job_id)
                        
                        if job_id in self.jobs_for_current_iteration:
                            del self.jobs_for_current_iteration[job_id]
                        jobs_to_remove.append(finished_job)
                
                if new_status['failed']>0:
                    failed_jobs = rq.job.Job.fetch_many(self.redis_queue.failed_job_registry.get_job_ids(), connection=self.redis_connection)
                    for failed_job in failed_jobs:
                        job_id = failed_job.id
                        logger.warning("%sRedis job '%s' failed with the following execution info:\n%s%s"%('\n'*50,str(job_id), str(failed_job.exc_info),'\n'*5))
                        if job_id in job_ids_already_handled:
                            logger.warning("The following job ID '%s' appeared on multiple occasions in the list of failed jobs, this should not happen."%str(job_id))
                        else:
                            # Notify that this job failed
                            if job_id in self.jobs_for_current_iteration:
                                self.n_redis_job_failed_for_this_iteration += 1
                                await self.process_results_callback(self.jobs_for_current_iteration[job_id]['alpha_loop_job_id'])
                            job_ids_already_handled.add(job_id)

                        # Remove the failed job so that it should no longer show up upon future queries
                        try:
                            self.redis_queue.failed_job_registry.remove(job_id, delete_job=True)
                        except Exception as e:
                            logger.warning("Failed to remove rq job %s, even though it should have been successful: %s."%(str(job_id), str(e)))
                            pass
                        if job_id in self.jobs_for_current_iteration:
                            del self.jobs_for_current_iteration[job_id]
                        jobs_to_remove.append(failed_job)

                # Remove the completed job so that it should no longer show up upon future queries
                for job in jobs_to_remove:
                    try:
                        #self.redis_queue.finished_job_registry.remove(job_id, delete_job=True)
                        job.delete()
                    except Exception as e:
                        logger.warning("Failed to remove rq job %s, even though it should have been successful: %s."%(str(job.id), str(e)))
                        pass

                new_status = {
                    'started' : self.redis_queue.started_job_registry.count,
                    'deferred' : self.redis_queue.deferred_job_registry.count,
                    'finished_and_not_processed' : self.redis_queue.finished_job_registry.count,
                    'failed' : self.redis_queue.failed_job_registry.count,
                    'scheduled' : self.redis_queue.scheduled_job_registry.count
                }
                last_status = new_status

            except KeyboardInterrupt as e:
                break

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
                        if job_result[0] in self.jobs_for_current_iteration:
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

        if self.use_redis:
            redis_server_cmd = ['redis-server','--port','%d'%self.redis_port,'--protected-mode','no','--save','""','--appendonly','no']
            if self.external_redis:
                warned = False
                while True:
                    try:
                        self.redis_connection = Redis(host=self.redis_submitter_hostname, port=self.redis_port)
                        # Test if the connection is established
                        self.redis_connection.ping()
                    except redis.exceptions.RedisError as e:
                        if not warned:
                            warned = True
                            logger.info("You specified using an external instance of a redis server, but alphaLoop could not connect to any because of: %s."%str(e))
                            logger.info("You can now manually start a redis instance with the following command for example:\n\n%s\n\nand alphaLoop will attempt to connect to it again every 3 seconds."%(
                                ' '.join(redis_server_cmd)
                            ))
                        time.sleep(3.)
                        continue
                    break
            else:
                redis_server_directory = pjoin(self.run_workspace,'redis_server')
                if os.path.exists(redis_server_directory):
                    shutil.rmtree(redis_server_directory)
                os.mkdir(redis_server_directory)
                # Start the redis server
                logger.info("Starting up a redis server in %s"%redis_server_directory)
                redis_monitor_file = open(pjoin(redis_server_directory,'redis_server_live_log.txt'),'a')
                self.redis_server_process = subprocess.Popen(redis_server_cmd, cwd=redis_server_directory, stdout=redis_monitor_file, stderr=redis_monitor_file)
                # Let the server boot (there may be more elegant ways of doing this, but w/e)
                time.sleep(1.)
                self.redis_connection = Redis(host=self.redis_submitter_hostname, port=self.redis_port)
            
            # Now use a unique queue
            if self.redis_queue_name is None:
                existing_queue_names = [q.name for q in rq.Queue.all(self.redis_connection)]
                suffix = 0
                this_queue_name = 'run_%d'%self.run_id
                while this_queue_name+'_%d'%suffix in existing_queue_names:
                    suffix += 1
                this_queue_name = this_queue_name+'_%d'%suffix
            else:
                this_queue_name = self.redis_queue_name

            self.redis_queue = rq.Queue(this_queue_name, connection=self.redis_connection)

        if self.architecture == 'local':

            for i_worker in range(self.n_workers):

                worker_id_min = i_worker*self.n_cores_per_worker
                worker_id_max = (i_worker+1)*self.n_cores_per_worker-1

                if self.use_redis:
                    cmd = [ sys.executable, pjoin(alphaloop_basedir,'alpha_loop','integrator','redis_worker.py'), 
                            '--run_id', str(self.run_id), '--worker_id_min', str(worker_id_min), '--worker_id_max', str(worker_id_max), 
                            '--workspace_path', str(self.run_workspace), '--redis_server', 'redis://%s:%d'%(self.redis_submitter_hostname, self.redis_port),
                            '--rq_path', self.rq_path, '--work_monitoring_path', ('none' if not self.work_monitoring_path else self.work_monitoring_path),
                            '--redis_queue_name', self.redis_queue.name
                        ]
                else:
                    cmd = [ sys.executable, pjoin(alphaloop_basedir,'alpha_loop','integrator','worker.py'), 
                            '--run_id', str(self.run_id), '--worker_id_min', str(worker_id_min), '--worker_id_max', str(worker_id_max), 
                            '--workspace_path', str(self.run_workspace), '--timeout', '%.2f'%self._JOB_CHECK_FREQUENCY]

                worker_env = os.environ.copy()
                if self.worker_resources_path:
                    worker_env['LD_PRELOAD'] = pjoin(self.worker_resources_path,'libs','libscsdir.so')
                else:
                    worker_env['LD_PRELOAD'] = pjoin(alphaloop_basedir,'libraries','scs','out','libscsdir.so')
                #worker_env['LD_PRELOAD'] = ' '.join([
                #    pjoin(alphaloop_basedir,'libraries','scs','out','libscsdir.so'),
                #    '/usr/lib/gcc/x86_64-linux-gnu/9/libstdc++.so',
                #    '/usr/lib/gcc/x86_64-linux-gnu/9/libquadmath.so',
                #    '/scratch/hirschva/mp++/mppp-0.26_build/lib/libmp++.so'
                #])
                #print("worker_env=%s"%worker_env)
                self.all_worker_hooks.append(subprocess.Popen(
                        cmd,
                        cwd=self.run_workspace,
                        stdout=(self.worker_monitoring_file if self._FORWARD_WORKER_OUTPUT else subprocess.DEVNULL),
                        stderr=(self.worker_monitoring_file if self._FORWARD_WORKER_OUTPUT else subprocess.DEVNULL),
                        env = worker_env
                    ))
                    
        elif self.architecture == 'condor':
            
            self.condor_deploy()

        if not self.use_redis:

            for i_worker in range(self.n_workers):

                worker_id_min = i_worker*self.n_cores_per_worker
                worker_id_max = (i_worker+1)*self.n_cores_per_worker-1

                for worker_id in range(worker_id_min, worker_id_max+1):
                    self.workers[worker_id] = {
                        'hook' : self.all_worker_hooks[-1],
                        'status' : 'PENDING'
                    }
                    await self.send_job( worker_id, 'ping')
        else:

            await self.send_redis_job('ping')

        await self.update_status()

        if self.use_redis:
            asyncio.create_task(self.follow_redis_workers())
        else:
            asyncio.create_task(self.follow_workers())

        # Wait for at least one worker to become active
        if wait_for_one_worker:
            if self.use_redis:
                # Find a way to lock on initialisation
                logger.info("Waiting for the redis parallelisation scheme to be properly setup...")
                await self.redis_initialization.wait()
            else:
                available_worker_id = await self.available_worker_ids.get()
                await self.available_worker_ids.put(available_worker_id)

        await self.update_status()

        if not self.use_redis:
            asyncio.create_task(self.job_submitter())

    def condor_deploy(self):
        
        if not os.path.exists(pjoin(self.run_workspace,'run_%d_condor_logs'%self.run_id)):
            os.mkdir(pjoin(self.run_workspace,'run_%d_condor_logs'%self.run_id))

        with open(pjoin(self.run_workspace,'run_%d_condor_worker_arguments.txt'%self.run_id),'w') as f:
            f.write('\n'.join(
                '%d, %d'%(i_worker*self.n_cores_per_worker, (i_worker+1)*self.n_cores_per_worker-1) for i_worker in range(self.n_workers)
            ))
        
        if self.worker_resources_path:
            libscsdir_path = pjoin(self.worker_resources_path,'libs','libscsdir.so')
        else:
            libscsdir_path = pjoin(alphaloop_basedir,'libraries','scs','out','libscsdir.so')
            if not os.path.exists(libscsdir_path):
                raise HavanaIntegratorError("Could not find libscsdir.so at '%s'. Make sure it is present and compiled."%libscsdir_path)

        if self.worker_resources_path:
            worker_run_workspace_dir = pjoin(self.worker_resources_path,'ProcessOutput','run_workspace')
        else:
            worker_run_workspace_dir = self.run_workspace

        with open(pjoin(self.run_workspace,'run_%d_condor_submission.sub'%self.run_id),'w') as f:
            
            format_dict = {
                'python' : sys.executable,
                'libscsdir_path' : libscsdir_path,
                'run_id' : self.run_id,
                'workspace' : self.run_workspace,
                'worker_run_workspace' : worker_run_workspace_dir,
                'timeout' : '%.2f'%(self._JOB_CHECK_FREQUENCY*15.),
                'job_flavour' : self.cluster_options['job_flavour'],
                'n_cpus_per_worker' : self.n_cores_per_worker,
                'requested_memory_in_MB' : self.cluster_options['RAM_required'],
                'redis_server' : 'redis://%s:%d'%(self.redis_submitter_hostname, self.redis_port),
                'rq_path' : self.rq_path,
                'transfer_input_files' : ( "" if self.worker_resources_path is None else pjoin(self.standalone_workers,os.path.basename(self.worker_resources_path)) )
            }
            if self.use_redis:
                format_dict['worker_script'] = pjoin(alphaloop_basedir,'alpha_loop','integrator','redis_worker.py')
                format_dict['redis_queue_name'] = self.redis_queue.name
                f.write(
    """executable         = %(python)s
    arguments             = %(worker_script)s --run_id %(run_id)d --worker_id_min $(worker_id_min) --worker_id_max $(worker_id_max) --workspace_path %(worker_run_workspace_dir)s --redis_server %(redis_server)s --rq_path %(rq_path)s --work_monitoring_path none --redis_queue_name %(redis_queue_name)s
    environment           = "LD_PRELOAD=%(libscsdir_path)s"
    output                = %(workspace)s/run_%(run_id)d_condor_logs/worker_$(worker_id_min)_$(worker_id_max).out
    error                 = %(workspace)s/run_%(run_id)d_condor_logs/worker_$(worker_id_min)_$(worker_id_max).err
    log                   = %(workspace)s/run_%(run_id)d_condor_logs/worker_$(worker_id_min)_$(worker_id_max).log
    RequestCpus           = %(n_cpus_per_worker)d
    RequestMemory         = %(requested_memory_in_MB)d
    RequestDisk           = DiskUsage
    should_transfer_files = Yes
    transfer_input_files  = %(transfer_input_files)s
    transfer_output_files = ""
    when_to_transfer_output = ON_EXIT
    +JobFlavour           = "%(job_flavour)s"
    queue worker_id_min,worker_id_max from %(workspace)s/run_%(run_id)d_condor_worker_arguments.txt
    """%format_dict)
            else:
                format_dict['worker_script'] = pjoin(alphaloop_basedir,'alpha_loop','integrator','worker.py')
                f.write(
    """executable         = %(python)s
    arguments             = %(worker_script)s --run_id %(run_id)d --worker_id_min $(worker_id_min) --worker_id_max $(worker_id_max) --workspace_path %(worker_run_workspace_dir)s --timeout %(timeout)s
    environment           = "LD_PRELOAD=%(libscsdir_path)s"
    output                = %(workspace)s/run_%(run_id)d_condor_logs/worker_$(worker_id_min)_$(worker_id_max).out
    error                 = %(workspace)s/run_%(run_id)d_condor_logs/worker_$(worker_id_min)_$(worker_id_max).err
    log                   = %(workspace)s/run_%(run_id)d_condor_logs/worker_$(worker_id_min)_$(worker_id_max).log
    RequestCpus           = %(n_cpus_per_worker)d
    RequestMemory         = %(requested_memory_in_MB)d
    RequestDisk           = DiskUsage
    should_transfer_files = Yes
    transfer_input_files  = %(transfer_input_files)s
    transfer_output_files = ""
    when_to_transfer_output = ON_EXIT
    +JobFlavour           = "%(job_flavour)s"
    queue worker_id_min,worker_id_max from %(workspace)s/run_%(run_id)d_condor_worker_arguments.txt
    """%format_dict)

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

    async def submit_many_jobs(self, payloads):
        if not self.use_redis:
            raise HavanaIntegratorError("Submitting many jobs at once is not supported nor useful when not using redis.")
        await self.send_many_redis_jobs(payloads)

    async def submit_job(self, payload, blocking=False):

        if self.use_redis:
            await self.send_redis_job(payload)
            return

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

        if self.use_redis:
            # Shutdown workers
            try:
                workers = rq.worker.Worker.all(queue=self.redis_queue)
                for worker in workers:
                    try:
                        rq.command.send_shutdown_command(self.redis_connection, worker.name)
                    except Exception as e:
                        pass
            except Exception as e:
                pass
            time.sleep(0.5)
    
        if self.redis_server_process:
            logger.info("Terminating redis server...")
            try:
                self.redis_server_process.terminate()
                time.sleep(0.5)
                self.redis_server_process.kill()

                redis_server_directory = pjoin(self.run_workspace,'redis_server')
                if os.path.exists(redis_server_directory):
                    shutil.rmtree(redis_server_directory)

            except Exception as e:
                pass

        io_files = \
            [f for f in glob.glob(pjoin(self.run_workspace, 'run_%d_job_*.pkl'%self.run_id))]+\
            [f for f in glob.glob(pjoin(self.run_workspace, 'run_%d_job_*.done'%self.run_id))]+\
            [f for f in glob.glob(pjoin(self.run_workspace, 'run_%d_job_*.bin'%self.run_id))]+\
            [f for f in glob.glob(pjoin(self.run_workspace, 'run_%d_job_*.yaml'%self.run_id))]+\
            [f for f in glob.glob(pjoin(self.run_workspace, '*trashcan*'))]
        for io_file in io_files:
            if os.path.isdir(io_file):
                shutil.rmtree(io_file)
            else:
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
                 use_optimal_integration_channels = True,
                 use_redis = False,
                 redis_max_job_time = None,
                 max_iteration_time = None,
                 external_redis = False,
                 redis_port=8786,
                 redis_hostname=None,
                 bulk_redis_enqueuing=0,
                 redis_queue = None,
                 write_common_grid_inputs_to_disk=True,
                 standalone_workers = False,
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
            [f for f in glob.glob(pjoin(self.run_workspace, 'run_%d_job_*.yaml'%self.run_id))]+\
            [f for f in glob.glob(pjoin(self.run_workspace, '*trashcan*'))]
        for io_file in io_files:
            if os.path.isdir(io_file):
                shutil.rmtree(io_file)
            else:
                os.remove(io_file)

        # Setup a trashcan
        self.trashcan = pjoin(self.run_workspace,'run_%s_trashcan'%self.run_id)
        if os.path.exists(self.trashcan):
            shutil.rmtree(self.trashcan)
        if os.path.exists(self.trashcan+'_being_removed'):
            shutil.rmtree(self.trashcan+'_being_removed')
        os.mkdir(self.trashcan)

        self.max_iteration_time = max_iteration_time

        self.havana_optimize_on_variance = havana_optimize_on_variance
        self.havana_max_prob_ratio = havana_max_prob_ratio
        self.show_SG_grid = show_SG_grid
        self.show_channel_grid = show_channel_grid
        self.show_grids_sorted_by_variance = show_grids_sorted_by_variance
        self.show_selected_phase_only = show_selected_phase_only
        self.show_all_information_for_all_integrands = show_all_information_for_all_integrands
        self.dump_havana_grids = dump_havana_grids
        self.fresh_integration = fresh_integration
        self.use_optimal_integration_channels = use_optimal_integration_channels
        self.use_redis = use_redis
        self.external_redis = external_redis
        self.redis_max_job_time = redis_max_job_time

        self.n_start = n_start
        self.n_max = n_max
        self.n_increase = n_increase
        self.redis_port = redis_port
        self.redis_hostname = redis_hostname
        self.redis_queue = redis_queue
        self.bulk_redis_enqueuing = bulk_redis_enqueuing
        self.write_common_grid_inputs_to_disk = write_common_grid_inputs_to_disk
        self.standalone_workers = standalone_workers
        self.worker_resources_path = None

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
            if isinstance(job_result,int):
                if self._DEBUG: logger.info("Job #%d failed and will thus be ignored."%job_result)
                self.n_jobs_awaiting_completion -= 1
                return

            job_ID, run_time, results = job_result
            if self._DEBUG: logger.info("Processing result from job #%d."%job_ID)

            self.cumulative_job_time += run_time

            t2 = time.time()
            if self._DEBUG: logger.info("Accumulating new results from job %d into havana..."%job_ID)

            havana_constructor_args = results['havana_updater']
            havana_constructor_args['integrands'] = self.integrands

            if self._USE_HAVANA_MOCKUP:
                havana_updater = HavanaMockUp(**havana_constructor_args)
            else:
                havana_updater = Havana(**havana_constructor_args)
                if havana_constructor_args['flat_record']['on_disk']:
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

            if self._DEBUG: logger.info("Done processing new results from job %d into havana..."%job_ID)

        except KeyboardInterrupt:
            self.exit_now = True

    async def process_al_cluster_status_update(self, new_status):
        
        self.cluster_status_summary = new_status
        if self.canvas is not None and self.stream_monitor:
            self.canvas.print('\n'.join('| %s'%line for line in self.cluster_status_summary.split('\n')))

    def handle_exception(self, loop, context):
        self.exit_now = True
        logger.critical("The following exception happened in a coroutine:\n%s\nAborting integration now."%str(context))

    @staticmethod
    def custom_rm(rm_target):
        if os.path.isdir(rm_target):
            shutil.rmtree(rm_target)
        else:
            os.remove(rm_target)

    def prepare_standalone_directory(self):
        """ Prepare a standalone directory with all necessary resources for running workers so that he can be submitted and copy onto the worker nodes."""
        
        itg_args = self.integrands[0].get_constructor_arguments()
        
        if os.path.isdir(self.standalone_workers) and os.path.isdir(pjoin(self.standalone_workers,'ProcessOutput','run_workspace')):
            logger.info("Recycling existing resources directory '%s'"%self.standalone_workers)
            run_workspace = pjoin(self.standalone_workers,'ProcessOutput','run_workspace')
            if not os.path.isdir(pjoin(run_workspace,itg_args['run_dir'])):
                os.mkdir(pjoin(run_workspace,itg_args['run_dir']))
            for itg in self.integrands:
                yaml_hyperparams = itg.get_constructor_arguments()['run_hyperparameters_filename']
                shutil.copy(
                    pjoin(itg_args['run_workspace'], itg_args['run_dir'], yaml_hyperparams),
                    pjoin(run_workspace, itg_args['run_dir'], yaml_hyperparams)
                )
            return self.standalone_workers

        elif (not os.path.isdir(self.standalone_workers)) and os.path.isdir(os.path.abspath(pjoin(self.standalone_workers, os.path.pardir))):
            os.mkdir(self.standalone_workers)
            standalone_dir = self.standalone_workers
        elif os.path.isdir(self.standalone_workers) and not any(os.scandir(self.standalone_workers)):
            standalone_dir = self.standalone_workers
        else:
            standalone_dir = pjoin(self.standalone_workers,'alphaLoop_%s'%str(uuid.uuid4()))
            os.mkdir(standalone_dir)

        # Copy all Python resources
        #python_path = os.path.abspath(pjoin(sys.executable,os.path.pardir,os.path.pardir))
        #logger.info("Copying Python resources...")
        #shutil.copytree(python_path,pjoin(standalone_dir,'Python'))

        # Copy relevant process output resources
        logger.info("Copying relevant process output resources into '%s'..."%standalone_dir)
        process_output_path = pjoin(standalone_dir,'ProcessOutput')
        os.mkdir(process_output_path)
        shutil.copytree(itg_args['rust_input_folder'], pjoin(process_output_path,'Rust_inputs') )
        shutil.copy(
            itg_args['cross_section_set_file_path'], 
            pjoin(process_output_path,'Rust_inputs',os.path.basename(itg_args['cross_section_set_file_path'])) 
        )
        shutil.copytree(os.path.join(itg_args['rust_input_folder'], os.path.pardir,'lib'), pjoin(process_output_path,'lib') )

        # Copy relevant alphaLoop resources
        logger.info("Copying relevant alphaLoop resources...")
        alphaLoop_path = pjoin(standalone_dir,'alpha_loop')
        os.mkdir(alphaLoop_path)
        shutil.copy(pjoin(itg_args['alpha_loop_path'],'ltd.so'),pjoin(alphaLoop_path,'ltd.so'))

        # Copy relevant run_workspace resources
        logger.info("Copying relevant run workspace resources...")
        run_workspace = pjoin(process_output_path,'run_workspace')
        os.mkdir(run_workspace)
        os.mkdir(pjoin(run_workspace,itg_args['run_dir']))
        for itg in self.integrands:
            yaml_hyperparams = itg.get_constructor_arguments()['run_hyperparameters_filename']
            shutil.copy(
                pjoin(itg_args['run_workspace'], itg_args['run_dir'], yaml_hyperparams),
                pjoin(run_workspace, itg_args['run_dir'], yaml_hyperparams)
            )

        # Copy relevant libraries
        logger.info("Copying relevant library resources...")
        libs_path = pjoin(standalone_dir,'libs')
        os.mkdir(libs_path)
        found_it = False
        for ext in ['so', 'dylib']:
            libscsdir_path = pjoin(itg_args['alpha_loop_path'],'libraries','scs','out','libscsdir.%s'%ext)
            if os.path.exists(libscsdir_path):
                shutil.copy(libscsdir_path,pjoin(libs_path,'libscsdir.%s'%ext))
                found_it = True
                break
        if not found_it:
            raise HavanaIntegratorError("Could not find libscsdir.so at '%s'. Make sure it is present and compiled."%(pjoin(itg_args['alpha_loop_path'],'libraries','scs','out','libscsdir.[so|dylib]')))

        return standalone_dir

    async def async_integrate(self):

        global trash_counter

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
                len(self.all_supergraphs[SG['name']]['multi_channeling_bases']) if (len(self.all_supergraphs[SG['name']]['optimal_channel_ids'])==0 or not self.use_optimal_integration_channels)
                    else len(self.all_supergraphs[SG['name']]['optimal_channel_ids']) for SG in self.cross_section_set['topologies']
            )
        else:
            n_channels_per_SG = []
            for (i_SG, SG_name) in SG_ids:
                n_channels_per_SG.append(
                    len(self.all_supergraphs[SG_name]['multi_channeling_bases']) if (len(self.all_supergraphs[SG_name]['optimal_channel_ids'])==0 or not self.use_optimal_integration_channels)
                    else len(self.all_supergraphs[SG_name]['optimal_channel_ids'])
                )
        
        if SG_ids is None:
            total_number_of_integration_channels = sum(
                len(self.all_supergraphs[SG['name']]['multi_channeling_bases']) if (len(self.all_supergraphs[SG['name']]['optimal_channel_ids'])==0 or not self.use_optimal_integration_channels)
                    else len(self.all_supergraphs[SG['name']]['optimal_channel_ids']) for SG in self.cross_section_set['topologies']
            )
        else:
            total_number_of_integration_channels = 0
            for (i_SG, SG_name) in SG_ids:
                total_number_of_integration_channels += (
                    len(self.all_supergraphs[SG_name]['multi_channeling_bases']) if (len(self.all_supergraphs[SG_name]['optimal_channel_ids'])==0 or not self.use_optimal_integration_channels)
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
            if os.path.isfile(pjoin(self.run_workspace,'run_description.txt')):
                with open(pjoin(self.run_workspace,'run_description.txt'), 'r') as f:
                    run_description = f.read()
            else:
                run_description = None
            if os.path.isfile(pjoin(self.run_workspace,'integrand_descriptions.txt')):
                with open(pjoin(self.run_workspace,'integrand_descriptions.txt'), 'r') as f:
                    integrand_descriptions = f.read().split('\n')
            else:
                integrand_descriptions = None
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
                run_description = run_description,
                integrand_descriptions = integrand_descriptions,
                run_id = self.run_id
            )

        if self.n_max < 0:
            logger.info("Latest result:\n\n%s\n%s"%(
                self.havana.get_summary(),
                self.havana.get_grid_summary(
			sort_by_variance=self.show_grids_sorted_by_variance,
                	show_channel_grid= (self.MC_over_channels and self.show_channel_grid),
                	show_selected_phase_only = self.show_selected_phase_only,
                	show_all_information_for_all_integrands = self.show_all_information_for_all_integrands
                )
            ))
            return

        integrands_constructor_args = [
            integrand.get_constructor_arguments() for integrand in self.integrands
        ]

        # Now perform the integration
        if self.standalone_workers:
            resources_dir_name = self.prepare_standalone_directory()
            if self.cluster_type == 'condor':
                # On condor, resources will be copied to the local directory
                self.worker_resources_path = pjoin('.',os.path.basename(resources_dir_name))
            else:
                # Locally, we keep resources at the original location
                self.worker_resources_path = os.path.abspath(resources_dir_name)
            # Now overwrite the integrands_constructor_args so as to specify local resources instead.
            for a in integrands_constructor_args:
                a['alpha_loop_path'] = pjoin(self.worker_resources_path,'alpha_loop')
                a['run_workspace'] = pjoin(self.worker_resources_path,'ProcessOutput','run_workspace')
                a['cross_section_set_file_path'] = pjoin(self.worker_resources_path,'ProcessOutput','Rust_inputs',os.path.basename(a['cross_section_set_file_path']))
                a['rust_input_folder'] = pjoin(self.worker_resources_path,'ProcessOutput','Rust_inputs')

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
            monitor_callback=self.process_al_cluster_status_update, cluster_options=cluster_options, keep=self.keep, debug=self._DEBUG, 
            trashcan=self.trashcan, use_redis=self.use_redis, redis_max_job_time = self.redis_max_job_time, external_redis = self.external_redis,
            redis_hostname=self.redis_hostname, redis_port=self.redis_port, redis_queue_name = self.redis_queue, worker_resources_path = self.worker_resources_path,
            standalone_workers = self.standalone_workers
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
                
                # Because these input grids will be common to all jobs, it is often more efficient to write it once and for all on disk and have the jobs read it from there
                write_common_grid_inputs_to_disk = True
                # except if explicitly asked not to do this
                if self.use_redis and (not self.write_common_grid_inputs_to_disk):
                    write_common_grid_inputs_to_disk = False
                sampling_grid_constructor_arguments = self.havana.get_constructor_arguments( on_disk = write_common_grid_inputs_to_disk )

                n_submitted = 0
                n_done = 0
                last_n_done = 0

                logger.info("Now running iteration %s#%d%s"%(bcolors.GREEN, self.current_iteration, bcolors.END))
                iteration_start_time = time.time()
                # Notify the cluster that a new iteration starts now, so that all potentially still running/pending jobs will be ignore when returning,
                # as they would pertain to an older iteration
                self.al_cluster.start_new_iteration()
                self.update_status(
                    start_time, n_jobs_total_completed, n_submitted, n_done,
                    self.n_tot_points, self.n_points_for_this_iteration, current_n_points, self.cumulative_IO_time,
                    self.cumulative_processing_time, self.cumulative_job_time, self.current_iteration, n_jobs_for_this_iteration=n_jobs_for_this_iteration
                )

                # Collect jobs here in case we want to do bulk submission
                job_payloads = []

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
                        if self.use_redis and self.bulk_redis_enqueuing>0:
                            job_payloads.append(job_payload)
                            if len(job_payloads)>=self.bulk_redis_enqueuing:
                                if self._DEBUG: logger.info("Submitting batch of %d jobs."%(len(job_payloads)))
                                await self.al_cluster.submit_many_jobs(job_payloads)
                                job_payloads.clear()
                                if self._DEBUG: logger.info("Done with submission of batch of %s jobs."%(len(job_payloads)))
                        else:
                            if self._DEBUG: logger.info("Submitting job %d."%job_ID)
                            await self.al_cluster.submit_job(
                                job_payload,
                                blocking=False
                            )
                            if self._DEBUG: logger.info("Done with submission of job %d."%job_ID)
                        self.cumulative_IO_time += time.time()-t0

                        self.n_jobs_awaiting_completion += 1
                        n_submitted += 1
                        timeout = 0.001
                        show_waiting = True
                        self.update_status(
                            start_time, n_jobs_total_completed, n_submitted, n_done,
                            self.n_tot_points, self.n_points_for_this_iteration, current_n_points, self.cumulative_IO_time,
                            self.cumulative_processing_time, self.cumulative_job_time, self.current_iteration, n_jobs_for_this_iteration=n_jobs_for_this_iteration)

                    else:
                        if self.use_redis and len(job_payloads)>0:
                            if self._DEBUG: logger.info("Submitting batch of %d jobs."%(len(job_payloads)))
                            await self.al_cluster.submit_many_jobs(job_payloads)
                            job_payloads.clear()
                            if self._DEBUG: logger.info("Done with submission of batch of %s jobs."%(len(job_payloads)))
                            self.update_status(
                                start_time, n_jobs_total_completed, n_submitted, n_done,
                                self.n_tot_points, self.n_points_for_this_iteration, current_n_points, self.cumulative_IO_time,
                                self.cumulative_processing_time, self.cumulative_job_time, self.current_iteration, n_jobs_for_this_iteration=n_jobs_for_this_iteration)

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

                    if self.n_points_for_this_iteration>0:
                        if n_remaining_points == 0 and self.n_jobs_awaiting_completion == 0:
                            break
                        if self.max_iteration_time is not None:
                            curr_iteration_time = time.time()-iteration_start_time
                            if curr_iteration_time >  self.max_iteration_time:
                                logger.warning("Max iteration time reached: %.5g > %.5g. Forcing the termination of this iteration #%d now."%(
                                    curr_iteration_time, self.max_iteration_time, self.current_iteration
                                ))
                                break
                
                # Force the number of jobs awaiting completion to be 0
                self.n_jobs_awaiting_completion = 0

                trash_counter += 1
                shutil.move(self.trashcan, '%s_%d_%s'%(self.trashcan,trash_counter,'being_removed'))
                os.mkdir(self.trashcan)
                subprocess.Popen('rm -rf %s_%d_%s'%(self.trashcan,trash_counter,'being_removed'), shell=True, stdout=subprocess.PIPE)
                if self.exit_now:
                    break

                t2 = time.time()
                if self._USE_HAVANA_MOCKUP and self._DEBUG: logger.info("Updating Havana grids for iteration #%d..."%self.current_iteration)
                self.havana.update_iteration_count(dump_grids=self.dump_havana_grids)

                if not self._USE_HAVANA_MOCKUP:
                    self.havana.save_state(pjoin(self.run_workspace,'latest_results.yaml'))
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
        
        if self.n_max > 0:
            logger.info("Now deploying the running cluster...")
        
        if self.cluster_type in ['local','condor']:

            self.al_cluster = None
            try:
                asyncio.get_event_loop().run_until_complete(self.async_integrate())
            except KeyboardInterrupt:
                logger.info("Havana integration aborted by user.")

            if self.al_cluster is not None:
                self.al_cluster.terminate()

            if not self._USE_HAVANA_MOCKUP:
                self.havana.save_state(pjoin(self.run_workspace,'latest_results.yaml'))

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
