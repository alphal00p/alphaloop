class WorkerException(Exception):
    pass

class ALStandaloneIntegrand(object):

    #TODO make this path specific to the current run (using the unique run id of the HavanaIntegrator session)
    _dask_run_hyperparameters_filename = 'cluster_run_hyperparameters.yaml'

    def __init__(self, n_integration_dimensions, alpha_loop_path, run_workspace, rust_input_folder, 
                        cross_section_set_file_name, E_cm, n_dimensions_per_SG_id=None, frozen_momenta=None):

        import sys
        import os
        pjoin = os.path.join

        self.constructor_arguments = {
            'alpha_loop_path' : alpha_loop_path,
            'run_workspace' : run_workspace,
            'cross_section_set_file_name' : cross_section_set_file_name,
            'rust_input_folder' : rust_input_folder,
            'n_integration_dimensions' : n_integration_dimensions,
            'n_dimensions_per_SG_id' : n_dimensions_per_SG_id,
            'frozen_momenta' : frozen_momenta,
            'E_cm' : E_cm
        }
        self.alpha_loop_path = alpha_loop_path
        self.run_workspace = run_workspace
        self.rust_input_folder = rust_input_folder
        self.cross_section_set_file_name = cross_section_set_file_name
        self.frozen_momenta = frozen_momenta
        self.E_cm = E_cm

        # This dimensions specification is not really used for now, only the length of the set of 
        # continuous and discrete dimensions matters for now.
        self.n_integration_dimensions = n_integration_dimensions
        self.n_dimensions_per_SG_id = n_dimensions_per_SG_id

        self.frozen_momenta = frozen_momenta

        # Now start a rust worker
        try:
            if self.alpha_loop_path not in sys.path:
                sys.path.insert(0, self.alpha_loop_path)
            # Import the rust bindings
            from ltd import CrossSection
        except ImportError as e:
            raise WorkerException("Could not import the rust back-end 'ltd' module in '%s'. Compile it first with:\n"%self.alpha_loop_path+
                " ./make_lib\nfrom within the alphaLoop directory.")

        #os.environ['MG_NUMERATOR_PATH'] = proc_path if proc_path.endswith('/') else '%s/'%proc_path
        if not os.path.isfile(pjoin(self.run_workspace, self._dask_run_hyperparameters_filename)):
            raise WorkerException("Could not find hyperparameter file at %s."%pjoin(self.run_workspace, self._dask_run_hyperparameters_filename))

        try:
            self.rust_worker = CrossSection(
                pjoin(self.rust_input_folder, cross_section_set_file_name),
                pjoin(self.run_workspace, self._dask_run_hyperparameters_filename),
                cross_section_set = True
            )
        except:
            raise WorkerException("Failed to load the rust backend API in the Dask integrand.")

    def __call__(self, samples_batch, SG_id_per_sample = None, channel_id_per_sample = None):

        # Pad xs with frozen momenta if necessary
        # TODO either do this internally in Rust or at least compute the quantity below as a single batch to avoid Python call overhead
        if self.frozen_momenta is not None:
            overall_inv_jac = [1.,]*len(samples_batch)
            for i_sample, sample in enumerate(samples_batch):
                n_loop_vs = len(sample)//3
                for i_v, v in enumerate([v[1:] for v in self.frozen_momenta['out']]):
                    xx, xy, xz, inv_jac = self.rust_worker.inv_parameterize(list(v), n_loop_vs+i_v, self.E_cm**2)
                    sample.extend([xx, xy, xz])
                    overall_inv_jac[i_sample] *= inv_jac*(2.0*math.pi)**4

        if self.n_dimensions_per_SG_id is not None and SG_id_per_sample is not None:
            for i_sample in range(len(samples_batch)):
                samples_batch[i_sample] = samples_batch[i_sample][:self.n_dimensions_per_SG_id[SG_id_per_sample[i_sample]]]

        if SG_id_per_sample is None and channel_id_per_sample is None:
            all_res = self.rust_worker.evaluate_integrand_batch(samples_batch)
        elif (SG_id_per_sample is not None) and channel_id_per_sample is None:
            all_res = self.rust_worker.evaluate_integrand_batch(samples_batch, sg_ids=SG_id_per_sample)
        elif (SG_id_per_sample is not None) and (channel_id_per_sample is not None):
            all_res = self.rust_worker.evaluate_integrand_batch(samples_batch, sg_ids=SG_id_per_sample, channel_ids=channel_id_per_sample )
        else:
            raise WorkerException("AlphaLoop integrand cannot sum over supergraphs but evaluate a single specific channel.")

        if self.frozen_momenta is not None:
            return [complex(*r)*overall_inv_jac[i_r] for i_r, r in enumerate(all_res)]
        else:
            return [complex(*r) for r in all_res]

    def get_constructor_arguments(self):
        import copy
        return copy.deepcopy(self.constructor_arguments)

def HavanaIntegrandWrapper(
        job_id,
        samples_batch,
        integrands_constructor_args,
        SG_ids_specified=True,
        channels_specified=True,
        preconstructed_integrands=None
    ):

    import time
    start_time = time.time()

    #import random
    #time.sleep(20+10*random.random())

    import os
    import numpy
    
    IO_pickle_out_name = None
    if isinstance(samples_batch, str):
        import pickle
        IO_pickle_out_name = os.path.join(os.path.dirname(samples_batch),'dask_output_job_%d.pkl'%job_id )
        sample_pickle_path = samples_batch
        samples_batch = pickle.load( open( sample_pickle_path, 'rb' ) )
        os.remove(sample_pickle_path)
    
    if preconstructed_integrands is None:
        integrands = [
            ALStandaloneIntegrand(**integrand_constructor_args) for integrand_constructor_args in integrands_constructor_args
        ]
    else:
        integrands =preconstructed_integrands

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
        raise Exception("Havana cannot integrate specific integration channels while summing over all supergraphs.")

    for i_integrand, integrand in enumerate(integrands):

        #for _ in range(10):
        results = integrand(xs_batch, SG_id_per_sample = SG_id_per_sample, channel_id_per_sample = channel_id_per_sample)
        for i_res, res in enumerate(results):
            all_res[i_res][0+i_integrand*2] = res.real
            all_res[i_res][1+i_integrand*2] = res.imag

    if IO_pickle_out_name is not None:
        import pickle
        with open(IO_pickle_out_name, 'wb') as f:
            pickle.dump( all_res, f )
        all_res = IO_pickle_out_name

    #print('run_time=%s'%str(time.time()-start_time))
    return (job_id, time.time()-start_time, all_res)

def main(arg_tuple):

    run_id, worker_id, run_workspace = arg_tuple
    import pickle
    import time
    import os

    standalone_integrands = None
    last_integrands_contructor_args = None

    input_path = os.path.join(run_workspace,'run_%d_job_input_worker_%d.pkl'%(run_id, worker_id))
    input_path_done = os.path.join(run_workspace,'run_%d_job_input_worker_%d.done'%(run_id, worker_id))
    output_path = os.path.join(run_workspace,'run_%d_job_output_worker_%d.pkl'%(run_id, worker_id))
    output_path_done = os.path.join(run_workspace,'run_%d_job_output_worker_%d.done'%(run_id, worker_id))
    while True:

        try:
            if (not os.path.exists(input_path_done) or not os.path.exists(input_path)) or os.path.exists(output_path_done):
                time.sleep(0.2)
                continue
        
            input_payload = pickle.load( open( input_path, 'rb' ) )
            
            output = None
            if input_payload == 'ping':
                output = 'pong'
            else:
                if last_integrands_contructor_args is None or input_payload['integrands_constructor_args'] != last_integrands_contructor_args:
                    last_integrands_contructor_args = input_payload['integrands_constructor_args']
                    standalone_integrands = [
                        ALStandaloneIntegrand(**integrand_constructor_args) for integrand_constructor_args in input_payload['integrands_constructor_args']
                    ]
                call_options = dict(input_payload)
                call_options['preconstructed_integrands'] = standalone_integrands
                output = HavanaIntegrandWrapper(**call_options)

            with open(output_path,'wb') as f:
                pickle.dump( output, f )
            with open(output_path_done,'w') as f:
                f.write("Marker file for job output done.")

        except Exception as e:
            print("The following exception was encountered in run #%d and worker #%d: %s"%(run_id, worker_id, str(e)))
            break
        except KeyboardInterrupt as e:
            #print("Keyboard interrupts in run #%d and worker #%d."%(run_id, worker_id))
            break

if __name__ == '__main__':
    import sys

    run_id = int(sys.argv[1])
    worker_id_min = int(sys.argv[2])
    worker_id_max = int(sys.argv[3])
    run_workspace = str(sys.argv[4])

    if worker_id_min==worker_id_max:
        print("Starting the following worker %d for run #%d for process '%s'."%(
           worker_id_min, run_id, run_workspace))
        sys.stdout.flush()
        try:
            main(tuple([run_id,worker_id_min,run_workspace]))
        except Exception as e:
            print("Worker %d finished (%s)."%(worker_id_min, str(e)))
    else:
        print("Starting the following range of workers %d->%d for run #%d for process '%s'."%(
           worker_id_min, worker_id_max, run_id, run_workspace))
        sys.stdout.flush()
        from multiprocessing import Pool
        try:
            with Pool(worker_id_max-worker_id_min+1) as p:
                    p.map( main, [ (run_id, worker_id, run_workspace) for worker_id in range(worker_id_min, worker_id_max+1)] )
        except Exception as e:
            print("Worker %d->%d finished (%s)."%(worker_id_min, worker_id_max, str(e)))