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

# classdistributed.deploy.local.LocalCluster(
#     name=None, n_workers=None, threads_per_worker=None, processes=True, loop=None, 
#     start=None, host=None, ip=None, scheduler_port=0, silence_logs=30, 
#     dashboard_address=':8787', worker_dashboard_address=None, diagnostics_port=None, 
#     services=None, worker_services=None, service_kwargs=None, asynchronous=False, 
#     security=None, protocol=None, blocked_handlers=None, interface=None, 
#     worker_class=None, scheduler_kwargs=None, **worker_kwargs)

# export DASK_DISTRIBUTED__SCHEDULER__ALLOWED_FAILURES=20
# export DASK_DISTRIBUTED__COMM__TIMEOUTS__CONNECT=5
# export DASK_DISTRIBUTED__COMM__RETRY__COUNT=10

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
            return (numpy.random.rand(self.dimensions), 1.0)
        else:
            return [ (numpy.random.rand(self.dimensions), 1.0) for _ in range(batch_size) ]

    def accumulate_results(self, results):
        
        self.n_points += len(results)
        if self.n_points == 0:
            return

        self.wgt_sum += sum(res[1]*res[2] for res in results)
        self.wgt_sum_squared += sum((res[1]*res[2])**2 for res in results)

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

@dask.delayed
def integrand(job_id, sample_points):
    #print("Start job %d"%job_id)
    start_time = time.time()
    all_res = []
    for sample_point, sample_jacobian_wgt in sample_points:
        # Slow down evaluation
        for _ in range(10000):
            integrand_wgt = math.exp(-sum((pi/(1.-pi))**2 for pi in sample_point))
            for pi in sample_point:
                integrand_wgt *= 1./(1-pi)**2
        all_res.append((sample_point, sample_jacobian_wgt, integrand_wgt))

    return (job_id, time.time()-start_time, all_res)

def MC_example(
        n_dims = 3,
        n_jobs = 8,
        n_iterations = 3,
        batch_size = int(1.0e3),
        n_workers = 16
    ):

    havana = HavanaMockUp(
        n_dims,
        target_result = (math.sqrt(math.pi)/2.)**n_dims,
        seed = 1
    )

    cluster = LocalCluster(
        n_workers=n_workers,
        interface='lo0' # This is the magic line that's supposed to help with MacOS. Try removing it on unix too.
    )
    client = Client(cluster)

    global_job_ID = 0
    for iteration_number in range(n_iterations):

        print("Preparing %d jobs for iteration #%d ..."%(n_jobs,iteration_number))
        jobs = [ integrand(global_job_ID+job_id, havana.sample(batch_size)) for job_id in range(n_jobs) ]
        global_job_ID += n_jobs
    
        print("Submitting %d jobs for iteration #%d ..."%(n_jobs,iteration_number))
        futures = client.compute(jobs)
        print("Done. Now waiting for job completion...")

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

@dask.delayed
def integrand_in_place(job_id, batch_size, havana):
    #print("Start job %d"%job_id)
    start_time = time.time()
    all_res = []
    if isinstance(havana, int):
        local_havana = HavanaMockUp(
            havana,
            target_result = (math.sqrt(math.pi)/2.)**havana,
            seed = job_id
        )
    else:
        local_havana = havana
    for _ in range(batch_size):
        sample_point, sample_jacobian_wgt = local_havana.sample()
        # No slow-down here
        for _ in range(1):
            integrand_wgt = math.exp(-sum((pi/(1.-pi))**2 for pi in sample_point))
            for pi in sample_point:
                integrand_wgt *= 1./(1-pi)**2
        all_res.append((sample_point, sample_jacobian_wgt, integrand_wgt))

    return (job_id, time.time()-start_time, all_res)

def MC_example_in_place_generation(
        n_dims = 3,
        n_jobs = 8,
        n_iterations = 3,
        batch_size = int(1.0e3),
        n_workers = 16
    ):

    havana = HavanaMockUp(
        n_dims,
        target_result = (math.sqrt(math.pi)/2.)**n_dims,
        seed = 1
    )

    cluster = LocalCluster(
        n_workers=n_workers,
        interface='lo0' # This is the magic line that's supposed to help with MacOS. Try removing it on unix too.
    )
    client = Client(cluster)

    global_job_ID = 0
    for iteration_number in range(n_iterations):

        print("Preparing %d jobs for iteration #%d ..."%(n_jobs,iteration_number))
        jobs = []
        for job_ID in range(n_jobs):
            local_havana = havana.get_copy(seed= (global_job_ID+job_ID) )
            # jobs.append(
            #     integrand_in_place(global_job_ID+job_ID, batch_size, local_havana )
            # )
            jobs.append(
                integrand_in_place(global_job_ID+job_ID, batch_size, n_dims )
            )
        global_job_ID += n_jobs
    
        print("Submitting %d jobs for iteration #%d ..."%(n_jobs,iteration_number))
        futures = client.compute(jobs)
        print("Done. Now waiting for job completion...")

        try:
            for completed_future in as_completed(futures):
                print("one job completed!")
                job_ID, run_time, results = completed_future.result()
                print("Job #%d completed in %.0f s for %d points."%(
                    job_ID, run_time, len(results)
                ))
                havana.accumulate_results(results)
                print("Integral thus far at iteration %d:\n%s"%(iteration_number, havana.get_summary()))
        except KeyboardInterrupt:
            print("Aborting computation now.")

@dask.delayed
def worker(job_id, data):
    print("Start job %d"%job_id)
    time.sleep(len(data))
    return (job_id, "I received: %s"%data)

def simple_example():
    #cluster = HTCondorCluster(cores=1, memory='100MB', disk='100MB')
    #cluster.scale(2)
    #print(cluster.job_script())
    cluster = LocalCluster()
    client = Client(cluster)

    #two_jobs = [ delayed(worker)(arg) for arg in ['D','CC'] ]
    two_jobs = [ worker(job_id, arg) for job_id, arg in enumerate(['D','CC']) ]

    print("Submitting jobs...")
    all_results = client.compute(two_jobs)
    print("Done. Now waiting for job completion...")
    for completed_future in as_completed(all_results):
        res = completed_future.result()
        print("Job #%d finished with result: %s"%(res[0],res[1]))
    # for i, res in enumerate(all_results):
    #     print("Now waiting on result #%d"%i)
    #     print(res.result())

if __name__ == "__main__":
    #simple_example()
    MC_example_in_place_generation(
        n_dims = 3,
        n_jobs = 12,
        n_iterations = 5,
        n_workers = 4,
        batch_size = int(1.0e3)
    )
    # MC_example(
    #     n_dims = 3,
    #     n_jobs = 16,
    #     n_iterations = 5,
    #     batch_size = int(1.0e3),
    #     n_workers = 16
    # )