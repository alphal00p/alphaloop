import math 

class HavanaIntegratorError(Exception):
    """Exception raised if an exception is triggered in Havana.""" 

class HavanaMockUp(object):

    def __init__(self, 
            n_dimensions, SG_ids=None, n_channels_per_SG=None, target_result=None, seed=None, n_integrands=1, 
            reference_integrand_index=0, phase='real', grid_file='havana_grid.yaml', 
            optimize_on_variance=True, max_prob_ratio=1000., fresh_integration=False,
            flat_record = None
        ):

        import numpy
        import random

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

        self.seed = seed
        if seed is not None:
            random.seed(seed)
            numpy.random.seed(seed)

        self.reset( clean = self.fresh_integration, flat_record=flat_record )

    def get_constructor_arguments(self):
        
        import copy
        return {
            'n_dimensions' : self.n_dimensions,
            'SG_ids' : self.SG_ids,
            'n_channels_per_SG' : self.n_channels_per_SG,
            'target_result' : self.target_result,
            'seed' : self.seed,
            'n_integrands' : self.n_integrands,
            'reference_integrand_index' : self.reference_integrand_index,
            'phase' : self.phase,
            'grid_file' : self.grid_file,
            'optimize_on_variance' : self.optimize_on_variance,
            'max_prob_ratio' :self.max_prob_ratio,
            'fresh_integration' : self.fresh_integration,
            'flat_record' : copy.deepcopy(self.get_flat_record())
        }

    def reset(self, clean=False, flat_record = None):

        import os

        if not clean and ( (flat_record is not None) or os.path.exists(self.grid_file) ):
            
            if flat_record is None:
                import yaml
                with open(self.grid_file,'r') as f:
                    flat_record = yaml.load(f, Loader=yaml.Loader)

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

    def get_flat_record(self):

        return {
            'current_integral_estimate' : self.current_integral_estimate,
            'current_error_estimate' : self.current_error_estimate,
            'wgt_sum' : self.wgt_sum,
            'wgt_sum_squared' : self.wgt_sum_squared,
            'n_points' : self.n_points,
            'discrete_grid' : self.discrete_grid,
            'n_iterations' :self.n_iterations
        }

    def dump_grid(self, file_name=None):

        flat_record = self.get_flat_record()

        with open(self.grid_file if file_name is None else file_name,'w') as f:
            import yaml
            class NoAliasDumper(yaml.SafeDumper):
                def ignore_aliases(self, data):
                    return True
            f.write(yaml.dump(flat_record, Dumper=NoAliasDumper, default_flow_style=False))

    def sample(self, batch_size):

        import numpy

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

    def accumulate_from_havana_grid(self, grid):

        if grid.n_points == 0:
            return

        self.n_points += grid.n_points
        for index in range(0,self.n_integrands*2):
            self.wgt_sum[index] += grid.wgt_sum[index]
            self.wgt_sum_squared[index] += grid.wgt_sum_squared[index]

            self.current_integral_estimate[index] = self.wgt_sum[index] / self.n_points
            self.current_error_estimate[index] = math.sqrt( ((self.wgt_sum_squared[index] / self.n_points) - self.current_integral_estimate[index]**2)/self.n_points )

        if self.discrete_grid is not None:
            for SG_index in range(len(self.discrete_grid)):
                self.discrete_grid[SG_index]['n_points'] += grid.discrete_grid[SG_index]['n_points']
                self.discrete_grid[SG_index]['cum_wgt'] += grid.discrete_grid[SG_index]['cum_wgt']
                self.discrete_grid[SG_index]['cum_sqr_wgt'] += grid.discrete_grid[SG_index]['cum_sqr_wgt']
                if self.discrete_grid[SG_index]['channel_grid'] is not None:
                    for channel_index in range(len(self.discrete_grid[SG_index]['channel_grid'])):
                        self.discrete_grid[SG_index]['channel_grid'][channel_index]['n_points'] += grid.discrete_grid[SG_index]['channel_grid'][channel_index]['n_points']
                        self.discrete_grid[SG_index]['channel_grid'][channel_index]['cum_wgt'] += grid.discrete_grid[SG_index]['channel_grid'][channel_index]['cum_wgt']
                        self.discrete_grid[SG_index]['channel_grid'][channel_index]['cum_sqr_wgt'] += grid.discrete_grid[SG_index]['channel_grid'][channel_index]['cum_sqr_wgt']

    def accumulate_results(self, results):
        
        if isinstance(results, HavanaMockUp):
            return self.accumulate_from_havana_grid(results)

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
            import shutil
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
            return ''

        sorted_bins = sorted([b for b in self.discrete_grid], key=lambda b: b['error_estimate'] if sort_by_variance else abs(b['integral_estimate']), reverse=True)
        res = ['Distribution of %d SGs:'%len(sorted_bins)]
        for bins_chunk in [sorted_bins[i_chunk:i_chunk+5] for i_chunk in range(0,len(sorted_bins),5)]:
            res.append('  '+' | '.join( '%-45s'%(
                '%s(I=%.5g(%.2g%%),p=%.2f%%,n=%.2gM)'%(
                    bin['SG_name'], bin['integral_estimate'], 
                    abs(bin['error_estimate']/bin['integral_estimate'])*100. if bin['integral_estimate']!=0. else 0.,
                    bin['p']*100., bin['n_points']/1000000.
                )
            ) for bin in bins_chunk ))

        if show_channel_grid:
            for SG_bin in sorted_bins:
                if SG_bin['channel_grid'] is not None:
                    res.append('Distribution of %d channels for SG %s:'%(len(SG_bin['channel_grid']), SG_bin['SG_name']))
                    sorted_channel_bins = sorted([(i_c, b) for i_c, b in enumerate(SG_bin['channel_grid'])], key=lambda b: b[1]['error_estimate'] if sort_by_variance else abs(b[1]['integral_estimate']), reverse=True)
                    for bins_chunk in [sorted_channel_bins[i_chunk:i_chunk+5] for i_chunk in range(0,len(sorted_channel_bins),5)]:
                        res.append('  '+' | '.join( '%-35s'%(
                            'C%d(I=%.5g(%.2g%%),p=%.2f%%,n=%.2gM)'%(
                                i_c, bin['integral_estimate'], 
                                abs(bin['error_estimate']/bin['integral_estimate'])*100. if bin['integral_estimate']!=0. else 0.,
                                bin['p']*100., bin['n_points']/1000000.
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

    def get_n_evals(self):
        return self.n_points

class WorkerException(Exception):
    pass

class ALStandaloneIntegrand(object):

    #TODO make this path specific to the current run (using the unique run id of the HavanaIntegrator session)
    _run_hyperparameters_filename = 'cluster_run_hyperparameters.yaml'

    def __init__(self, n_integration_dimensions, alpha_loop_path, run_workspace, rust_input_folder, 
                        cross_section_set_file_name, E_cm, run_dir='', n_dimensions_per_SG_id=None, frozen_momenta=None):

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
        hyperparameters_path = pjoin(self.run_workspace, run_dir, self._run_hyperparameters_filename)
        if not os.path.isfile(hyperparameters_path):
            raise WorkerException("Could not find hyperparameter file at %s."%hyperparameters_path)

        try:
            self.rust_worker = CrossSection(
                pjoin(self.rust_input_folder, cross_section_set_file_name),
                hyperparameters_path,
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

    run_id, worker_id, run_workspace, timeout = arg_tuple
    import pickle
    import time
    import os
    import sys

    standalone_integrands = None
    last_integrands_contructor_args = None

    input_path = os.path.join(run_workspace,'run_%d_job_input_worker_%d.pkl'%(run_id, worker_id))
    input_path_done = os.path.join(run_workspace,'run_%d_job_input_worker_%d.done'%(run_id, worker_id))
    output_path = os.path.join(run_workspace,'run_%d_job_output_worker_%d.pkl'%(run_id, worker_id))
    output_path_done = os.path.join(run_workspace,'run_%d_job_output_worker_%d.done'%(run_id, worker_id))
    while True:

        try:
            if (not os.path.exists(input_path_done) or not os.path.exists(input_path)) or os.path.exists(output_path_done):
                time.sleep(timeout)
                continue
            
            print("Worker #%d received a new job."%(worker_id))
            t0=time.time()
            input_payload = pickle.load( open( input_path, 'rb' ) )

            output = None
            if input_payload == 'ping':
                output = 'pong'
            else:
                if last_integrands_contructor_args is None or input_payload['integrands_constructor_args'] != last_integrands_contructor_args:
                    last_integrands_contructor_args = input_payload['integrands_constructor_args']
                    standalone_integrands = [
                        ALStandaloneIntegrand(run_dir='run_%d'%run_id, **integrand_constructor_args) for integrand_constructor_args in input_payload['integrands_constructor_args']
                    ]

                call_options = {
                    'job_id' : input_payload['job_id'],
                    'integrands_constructor_args' : input_payload['integrands_constructor_args'],
                    'SG_ids_specified' : input_payload['SG_ids_specified'],
                    'channels_specified' : input_payload['channels_specified'],
                    'preconstructed_integrands' : standalone_integrands
                }
                if not input_payload['worker_generates_samples']:

                    call_options['samples_batch'] = input_payload['samples_batch']
                    print("Worker #%d deserialized input for job #%d in %.5gs."%(worker_id, call_options['job_id'], time.time()-t0 ))
                    t0 = time.time()
                    output = HavanaIntegrandWrapper(**call_options)
                    print("Worker #%d completed computation of job #%d in %.5gs."%(worker_id, call_options['job_id'], time.time()-t0 ))

                else:

                    print("Worker #%d generates a sample of %d points for job #%d..."%(worker_id, input_payload['n_points_to_sample'], call_options['job_id']))
                    t0 = time.time()
                    havana_sampler = HavanaMockUp(**input_payload['havana_grid'])
                    this_sample = havana_sampler.sample(input_payload['n_points_to_sample'])
                    sample_generation_time = time.time()-t0
                    print("Worker #%d completed sample generation in %.5s."%(worker_id, sample_generation_time))
                    call_options['samples_batch'] = this_sample
                    job_id, run_time, samples_evaluated = HavanaIntegrandWrapper(**call_options)
                    havana_updater_constructor_args = dict(input_payload['havana_grid'])
                    havana_updater_constructor_args['fresh_integration'] = True
                    havana_updater = HavanaMockUp(**havana_updater_constructor_args)
                    havana_updater.accumulate_results(samples_evaluated)

                    is_phase_real = input_payload['phase']=='real'
                    n_integrands = len(standalone_integrands)
                    sample_statistics = [ {
                        'n_evals' : 0,
                        'n_evals_failed' : 0,
                        'n_zero_evals' : 0,
                        'max_eval_positive' : 0.,
                        'max_eval_positive_xs' : None,
                        'max_eval_negative' : 0.,
                        'max_eval_negative_xs' : None,
                    } for _ in range(n_integrands) ]

                    t0 = time.time()
                    for res in samples_evaluated:
                        for index in range(0,n_integrands):
                            if is_phase_real:
                                wgt = res[index*2]
                            else:
                                wgt = res[1+index*2]
                            sample_statistics[index]['n_evals'] += 1
                            if wgt is None:
                                sample_statistics[index]['n_evals_failed'] += 1
                                continue
                            if wgt == 0.:
                                sample_statistics[index]['n_zero_evals'] += 1
                                continue
                            
                            # Include havana jacobian 
                            wgt *= res[n_integrands*2]

                            if wgt>0. and wgt > sample_statistics[index]['max_eval_positive']:
                                sample_statistics[index]['max_eval_positive'] = wgt
                                sample_statistics[index]['max_eval_positive_xs'] = res[1+n_integrands*2:]

                            if wgt<0. and wgt < sample_statistics[index]['max_eval_negative']:
                                sample_statistics[index]['max_eval_negative'] = wgt
                                sample_statistics[index]['max_eval_negative_xs'] = res[1+n_integrands*2:]

                    havana_updater_constructor_arguments = havana_updater.get_constructor_arguments()
                    havana_updater_constructor_arguments['fresh_integration'] = False
                    complete_result = {
                        'havana_updater' : havana_updater_constructor_arguments,
                        'sample_statistics' : sample_statistics
                    }
                    post_processing_time = time.time()-t0
                    output = ( job_id, run_time+sample_generation_time+post_processing_time, complete_result)

            print("Worker #%d now writing job output."%(worker_id))
            t0 = time.time()
            with open(output_path,'wb') as f:
                pickle.dump( output, f )
            with open(output_path_done,'w') as f:
                f.write("Marker file for job output done.")
            print("Worker #%d finished serialisation and dump of job output in %.5gs."%( worker_id, time.time()-t0 ))
            sys.stdout.flush()

        except Exception as e:
            print("The following exception was encountered in run #%d and worker #%d: %s"%(run_id, worker_id, str(e)))
            break
        except KeyboardInterrupt as e:
            #print("Keyboard interrupts in run #%d and worker #%d."%(run_id, worker_id))
            break

if __name__ == '__main__':
    import sys

    from argparse import ArgumentParser

    parser = ArgumentParser(prog='al_worker')

    parser.add_argument("--run_id", dest="run_id", type=int)
    parser.add_argument("--worker_id_min", dest="worker_id_min", type=int)
    parser.add_argument("--worker_id_max", dest="worker_id_max", type=int)
    parser.add_argument("--workspace_path", dest="workspace_path", type=str)
    parser.add_argument("--timeout", dest="timeout", default=0.2, type=float)
    args = parser.parse_args()

    if args.worker_id_min==args.worker_id_max:
        print("Starting the following worker %d for run #%d for process '%s'."%(
           args.worker_id_min, args.run_id, args.workspace_path))
        sys.stdout.flush()
        try:
            main(tuple([args.run_id,args.worker_id_min,args.workspace_path, args.timeout]))
        except Exception as e:
            print("Worker %d finished (%s)."%(args.worker_id_min, str(e)))
    else:
        print("Starting the following range of workers %d->%d for run #%d for process '%s'."%(
           args.worker_id_min, args.worker_id_max, args.run_id, args.workspace_path))
        sys.stdout.flush()
        from multiprocessing import Pool
        try:
            with Pool(args.worker_id_max-args.worker_id_min+1) as p:
                    p.map( main, [ (args.run_id, worker_id, args.workspace_path, args.timeout) for worker_id in range(args.worker_id_min, args.worker_id_max+1)] )
        except Exception as e:
            print("Worker %d->%d finished (%s)."%(args.worker_id_min, args.worker_id_max, str(e)))