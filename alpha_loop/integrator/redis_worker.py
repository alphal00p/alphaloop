#!/usr/bin/env python3
import subprocess
import sys
import os
from argparse import ArgumentParser

# python3 run_rq_worker.py worker --with-scheduler --url redis://lxplus733.cern.ch:8786 --path /afs/cern.ch/work/v/vjhirsch/private/MG_v3_0_2_py3/PLUGIN/alphaloop/alpha_loop/integrator

def run_redis_worker(args):
    run_id, worker_id, work_monitoring_path, worker_path, redis_server, rq_path = args
    rq_cmd = [rq_path, 'worker', '--with-scheduler', '--url', redis_server, '--path', worker_path]
    print("Starting redis worker with command: %s"%(' '.join(rq_cmd)))
    if work_monitoring_path != 'none':
        with open(work_monitoring_path, 'a') as f:
            rq_job = subprocess.Popen(rq_cmd, stdout=f, stderr=f)
            res_stdout, res_stderr = rq_job.communicate()
    else:
        rq_job = subprocess.Popen(rq_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        res_stdout, res_stderr = rq_job.communicate()
    print("Redis worker #%d terminated with stdout:"%worker_id)
    print(res_stdout.decode('utf-8'))
    print("and stderr:")
    print(res_stderr.decode('utf-8'))
    sys.stdout.flush()

if __name__ == '__main__':

    print("Redis worker was called with: %s"%str(' '.join(sys.argv)))
    sys.stdout.flush()
    parser = ArgumentParser(prog='al_worker')

    parser.add_argument("--run_id", dest="run_id", type=int)
    parser.add_argument("--worker_id_min", dest="worker_id_min", type=int)
    parser.add_argument("--worker_id_max", dest="worker_id_max", type=int)
    parser.add_argument("--workspace_path", dest="workspace_path", type=str)
    parser.add_argument("--redis_server", dest="redis_server", type=str)
    parser.add_argument("--rq_path", dest="rq_path", type=str)
    parser.add_argument("--work_monitoring_path", dest="work_monitoring_path", type=str)
    
    args = parser.parse_args()

    worker_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),os.path.pardir,os.path.pardir))

    if args.worker_id_min==args.worker_id_max:
        print("Starting the following worker %d for run #%d for process '%s' and redis server: '%s'."%(
           args.worker_id_min, args.run_id, worker_path, args.redis_server))
        sys.stdout.flush()
        try:
            run_redis_worker(tuple([args.run_id,args.worker_id_min, args.work_monitoring_path, worker_path, args.redis_server, args.rq_path]))
        except Exception as e:
            print("Worker %d finished (%s)."%(args.worker_id_min, str(e)))
            sys.stdout.flush()
    else:
        print("Starting the following range of workers %d->%d for run #%d for process '%s' and redis server: '%s'."%(
           args.worker_id_min, args.worker_id_max, args.run_id, worker_path, args.redis_server))
        sys.stdout.flush()
        from multiprocessing import Pool
        try:
            with Pool(args.worker_id_max-args.worker_id_min+1) as p:
                p.map( run_redis_worker, [ (args.run_id, worker_id, args.work_monitoring_path, worker_path, args.redis_server, args.rq_path) for worker_id in range(args.worker_id_min, args.worker_id_max+1)] )
        except Exception as e:
            print("Worker %d->%d finished (%s)."%(args.worker_id_min, args.worker_id_max, str(e)))
            sys.stdout.flush()