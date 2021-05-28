#!/usr/bin/env python3
#
# sudo -H python3 -m pip install bokeh dask dask-jobqueue
#
import dask
from dask_jobqueue import HTCondorCluster
from distributed import Client, LocalCluster
from dask import delayed
from dask.distributed import progress
from dask.distributed import as_completed
from dask import delayed
import numpy
import time
import random
import math
import logging
from dask.distributed import progress

class HavanaMockUp(object):

    def __init__(self, dimensions, target_result=None, seed=None):
        self.dimensions = dimensions
        self.wgt_sum = 0.
        self.wgt_sum_squared = 0.
        self.n_points = 0
        self.target_result = target_result

        self.current_integral_estimate = 0.
        self.current_error_estimate = 0.

        if seed is not None:
            random.seed(seed)
            numpy.random.seed(seed)

    def get_copy(self, seed=None):
        return HavanaMockUp(self.dimensions, target_result=self.target_result, seed=seed)
    
    def sample(self, batch_size=None):
        # if batch_size is None:
        #     return (tuple([random.random() for _ in range(self.dimensions)]), 1.0)
        # else:
        #     return [ (tuple([random.random() for _ in range(self.dimensions)]), 1.0) for _ in range(batch_size) ]
        if batch_size is None:
            sample = numpy.ndarray(shape=(self.dimensions+1), dtype=float)
            sample[1:] = numpy.random.rand(self.dimensions)
            sample[0] = 1.0
            return sample
        else:
            sample = numpy.random.rand(batch_size,self.dimensions+1)
            for i_batch in range(batch_size):
                sample[i_batch][0] = 1.0
            return sample

    def accumulate_results(self, results):
        
        self.n_points += len(results)
        if self.n_points == 0:
            return

        self.wgt_sum += sum(res[0]*res[1] for res in results)
        self.wgt_sum_squared += sum((res[0]*res[1])**2 for res in results)

        self.current_integral_estimate = self.wgt_sum / self.n_points
        self.current_error_estimate = math.sqrt( ((self.wgt_sum_squared / self.n_points) - self.current_integral_estimate**2)/self.n_points )
    
    def get_summary(self):
        res = ["Result after %d evaluations: %.5g +/- %.5g (%.2g%%)"%(
            self.n_points, self.current_integral_estimate, self.current_error_estimate, 
            0. if self.current_integral_estimate==0. else (abs(self.current_error_estimate)/self.current_integral_estimate)*100.
        )]
        if self.target_result is not None:
            res.append( "vs target                      : %.5g del %.5g (%.2g%%)"%(
                    self.target_result, 
                    self.current_integral_estimate-self.target_result, 
                    0. if self.target_result==0. else ((self.current_integral_estimate-self.target_result)/self.target_result)*100.0
                ))
        return '\n'.join(res)

def integrand(job_id, sample_points):
    start_time = time.time()
    all_res = numpy.ndarray(shape=(len(sample_points),len(sample_points[0])+1), dtype=float)
    for i_sample, sample_point in enumerate(sample_points):
        sample_jacobian_wgt = sample_point[0]
        sample_point_coordinates = sample_point[1:]
        integrand_wgt = math.exp(-sum((pi/(1.-pi))**2 for pi in sample_point_coordinates))
        for pi in sample_point_coordinates:
            integrand_wgt *= 1./(1-pi)**2
        all_res[i_sample][0] = integrand_wgt
        all_res[i_sample][1] = sample_jacobian_wgt
        all_res[i_sample][2:] = sample_point_coordinates

    return (job_id, time.time()-start_time, all_res)

def run(client, n_dims, n_jobs, n_iterations, batch_size):

    havana = HavanaMockUp(
        n_dims,
        target_result = (math.sqrt(math.pi)/2.)**n_dims,
        seed = 1
    )

    global_job_ID = 0
    for iteration_number in range(n_iterations):

        # print("Preparing %d jobs for iteration #%d ..."%(n_jobs,iteration_number))
        # jobs = [ delayed(integrand)(global_job_ID+job_id, havana.sample(batch_size)) for job_id in range(n_jobs) ]
        # global_job_ID += n_jobs
    
        # print("Submitting %d jobs for iteration #%d ..."%(n_jobs,iteration_number))
        # futures = client.compute(jobs)
        # print("Done. Now waiting for job completion...")
        
        print("")
        futures = []
        for job_id in range(n_jobs):
            print("Generating input data for job %d/%d for iteration #%d ..."%(job_id+1,n_jobs,iteration_number))
            a_sample = havana.sample(batch_size)
            # Scattering data is not really necessary as I believe this is the right thing to do only for data common to multiple jobs. But I don't like the warning so I'll put it.
            print("Scattering input data for job %d/%d for iteration #%d ..."%(job_id+1,n_jobs,iteration_number))
            scattered_sample = client.scatter(a_sample)
            print("Submitting job %d/%d for iteration #%d ..."%(job_id+1,n_jobs,iteration_number))
            futures.append(client.submit(integrand, global_job_ID+job_id, scattered_sample))
        print("")
        print("Now waiting for completion of %d jobs at iteration #%d ... "%(n_jobs, iteration_number))
        print("")

        # Passive monitoring of jobs can be done with the widget below.
        #progress(futures)

        try:
            for completed_future in as_completed(futures):
                job_ID, run_time, results = completed_future.result()
                print("Job #%d completed in %.0f s for %d points."%(
                    job_ID, run_time, len(results)
                ))
                havana.accumulate_results(results)
                print("Integral thus far at iteration %d:\n%s"%(iteration_number, havana.get_summary()))
        except KeyboardInterrupt:
            print("Aborting computation now.")

def MC_example(
        n_dims = 3,
        n_jobs = 8,
        n_iterations = 3,
        batch_size = int(1.0e3),
        n_workers = 16
    ):

    with LocalCluster(n_workers=n_workers) as cluster:
        
        client = Client(cluster)
        run(client, n_dims, n_jobs, n_iterations, batch_size)

def CONDOR_MC_example(
        n_dims = 3,
        n_jobs = 8,
        n_iterations = 3,
        batch_size = int(1.0e3),
        n_workers = 16,
        io_port = 8786
    ):

    with HTCondorCluster(
            cores=1,
            memory='2000MB',
            disk='1000MB',
            death_timeout = '60',
            nanny = False,
            scheduler_options={
                'port': io_port,
                'host': socket.gethostname()
            },
            job_extra={
                'log'    : 'dask_job_output.log',
                'output' : 'dask_job_output.out',
                'error'  : 'dask_job_output.err',
                'should_transfer_files'   : 'Yes',
                'when_to_transfer_output' : 'ON_EXIT',
                '+JobFlavour' : "espresso",
            },
            extra = [ '--worker-port {}'.format(io_port) ]
    ) as cluster:

        print("Scaling cluster ...")
        cluster.scale(n_workers)

        print("Initialising cluster ...")
        print(cluster.job_script())
        client = Client(cluster)

        run(client, n_dims, n_jobs, n_iterations, batch_size)

if __name__ == "__main__":

    # Local multicore example
    MC_example(
        n_dims = 3,
        n_jobs = 4,
        n_iterations = 3,
        batch_size = int(1.0e5),
        n_workers = 2
    )

    # Cluster (lxplus) example
    # CONDOR_MC_example(
    #     n_dims = 3,
    #     n_jobs = 4,
    #     n_iterations = 3,
    #     batch_size = int(1.0e5),
    #     n_workers = 2
    # )