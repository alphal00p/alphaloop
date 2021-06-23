import math 
import tempfile

class bcolors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    WARNING = YELLOW
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    RED = '\033[91m'
    END = ENDC

class HavanaIntegratorError(Exception):
    """Exception raised if an exception is triggered in Havana.""" 

class Havana(object):
    """ A wrapper class around the Havana rust bindings."""

    def __init__(self, 
            n_dimensions, integrands, n_bins=128, n_points_min=1000, SG_ids=None, n_channels_per_SG=None, target_result=None, seed=None, 
            reference_integrand_index=0, phase='real', grid_file='havana_grid.yaml', 
            optimize_on_variance=True, max_prob_ratio=1000., fresh_integration=False,
            flat_record = None, alpha_loop_path=None, learning_rate=1.5, bin_increase_factor_schedule=None,
            **opts
        ):

        self.alpha_loop_path = alpha_loop_path

        self.n_dimensions = n_dimensions
        self.target_result = target_result
        self.integrands = integrands
        self.n_integrands = len(self.integrands)
        self.n_bins = n_bins
        self.n_points_min = n_points_min
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

        self.SG_ids = SG_ids
        self.n_channels_per_SG = n_channels_per_SG
        self.seed = seed

        self.havana_grids = None

        self.learning_rate = learning_rate
        if bin_increase_factor_schedule is None:
            self.bin_increase_factor_schedule = []
        else:
            self.bin_increase_factor_schedule = bin_increase_factor_schedule

        self.grid_file = grid_file
        self.fresh_integration = fresh_integration

        self.reset( clean = self.fresh_integration, flat_record=flat_record )
    

    def get_constructor_arguments(self, dump_format='bin', tmp_folder=None):
        """ Is a tmp folder is specified, it will be used to temporary write out grid from rust, 
        but it will be transmitted over network."""
        
        import copy
        return {
            'n_dimensions' : self.n_dimensions,
            'SG_ids' : self.SG_ids,
            'n_channels_per_SG' : self.n_channels_per_SG,
            'target_result' : self.target_result,
            'seed' : self.seed,
            'integrands' : None, # will be filled in at runtime
            'reference_integrand_index' : self.reference_integrand_index,
            'phase' : self.phase,
            'n_bins' : self.n_bins,
            'n_points_min' : self.n_points_min,
            'grid_file' : self.grid_file,
            'optimize_on_variance' : self.optimize_on_variance,
            'max_prob_ratio' :self.max_prob_ratio,
            'learning_rate' : self.learning_rate,
            'bin_increase_factor_schedule' : self.bin_increase_factor_schedule,
            'flat_record' : self.get_flat_record(dump_format=dump_format,tmp_folder=tmp_folder)
        }


    def reset(self, clean=False, flat_record = None):

        import os
        import sys

        # Access havana bindings
        try:
            if self.alpha_loop_path not in sys.path:
                sys.path.insert(0, self.alpha_loop_path)
            # Import the havana bindings, in the form of the ltd shared object library that contains the following submodule:
            #import ltd.havana as havana
            from ltd.havana import Havana, Sample, GridConstructor, ContinuousGridConstructor, DiscreteGridConstructor
        except ImportError as e:
            raise HavanaIntegratorError("Could not import the rust back-end 'ltd' module in '%s'. Compile it first with:\n"%self.alpha_loop_path+
                " ./make_lib\nfrom within the alphaLoop directory.")

        if not clean and ( (flat_record is not None) or os.path.exists(self.grid_file) ):
            
            if flat_record is None:
                import yaml
                with open(self.grid_file,'r') as f:
                    flat_record = yaml.load(f, Loader=yaml.Loader)
            self.n_points = flat_record['n_points']
            self.n_iterations = flat_record['n_iterations']
            if flat_record['tmp_folder'] is None:
                self.havana_grids = []
                for grid_file_path in flat_record['grid_files']:
                    with open(grid_file_path, 'rb') as f:
                        self.havana_grids.append(
                            Havana.load_grid(f.read(), seed=self.seed, format=grid_file_path.split('.')[-1]) for grid_file_path in flat_record['grid_files']
                        )
            else:
                self.havana_grids = []
                for serialized_grid in flat_record['grid_files']:
                    self.havana_grids.append( Havana.load_grid(serialized_grid, seed=self.seed, format='bin') )
        else:

            self.n_points = 0
            self.n_iterations = 0
            self.havana_grids = []
            for i_itg, integrand in enumerate(self.integrands):
                # Add one havana grid for the real and imaginary part
                for phase in ['real','imag']:
                    if self.SG_ids is None:
                        a_grid = GridConstructor(cgc=ContinuousGridConstructor(self.n_dimensions, self.n_bins, self.n_points_min))
                    elif self.n_channels_per_SG is None:                    
                        a_grid=GridConstructor(dgc=DiscreteGridConstructor(
                            [len(self.SG_ids),],
                            [
                                GridConstructor(cgc=ContinuousGridConstructor(integrand.n_dimensions_per_SG_id[SG_id], self.n_bins, self.n_points_min))
                                for (SG_id, SG_name) in self.SG_ids 
                            ],
                            self.max_prob_ratio 
                        ))
                    else:
                        a_grid=GridConstructor(dgc=DiscreteGridConstructor(
                            [len(self.SG_ids),],
                            [
                                GridConstructor(dgc=DiscreteGridConstructor(
                                    [n_channels,],
                                    [
                                        GridConstructor(cgc=ContinuousGridConstructor(integrand.n_dimensions_per_SG_id[SG_id], self.n_bins, self.n_points_min))
                                        for i_channel in range(n_channels) 
                                    ],
                                    self.max_prob_ratio 
                                ))
                                for n_channels, (SG_id, SG_name) in zip(self.n_channels_per_SG, self.SG_ids) 
                            ],
                            self.max_prob_ratio 
                        ))
                    self.havana_grids.append(Havana(a_grid, seed=self.seed))


    def get_havana_grid_filenames(self, top_level_grid_filename, dump_format='bin'):

        import os
        subgrid_file_paths = []
        for i_itg in range(0,len(self.havana_grids)//2):
            for i_phase, phase in enumerate(['real','imag']):
                grid_base_path = os.path.dirname(top_level_grid_filename)
                grid_name = '.'.join(s for s in os.path.basename(top_level_grid_filename).split('.')[:-1])
                grid_suffix = os.path.basename(top_level_grid_filename).split('.')[-1]
                subgrid_file_paths.append(
                    os.path.join(grid_base_path,'%s_itg_%d_%s.%s'%(
                        grid_name, i_itg+1, phase, dump_format
                    ))
                )
        return subgrid_file_paths

    def get_flat_record(self, file_name=None, dump_format='bin', tmp_folder=None):
        import time
        import sys
        import os

        havana_grid_file_paths = self.get_havana_grid_filenames(self.grid_file if file_name is None else file_name, dump_format=dump_format)
        if tmp_folder is None:
            for havana_grid, grid_file_path in zip(self.havana_grids, havana_grid_file_paths):
                bytesvec = havana_grid.dump_grid(format=dump_format)
                with open(grid_file_path, 'wb') as f:
                    f.write(bytesvec)
            grid_files = havana_grid_file_paths
        else:
            grid_files = []
            for havana_grid, grid_file_path in zip(self.havana_grids, havana_grid_file_paths):
                bytesvec = havana_grid.dump_grid(format=dump_format)
                grid_files.append( bytesvec )
        return {
            'tmp_folder' : tmp_folder,
            'n_points' : self.n_points,
            'grid_files' : grid_files,
            'n_iterations' :self.n_iterations
        }

    def dump_grid(self, file_name=None, dump_format='yaml'):

        flat_record = self.get_flat_record(file_name=file_name, dump_format=dump_format)

        with open(self.grid_file if file_name is None else file_name,'w') as f:
            import yaml
            class NoAliasDumper(yaml.SafeDumper):
                def ignore_aliases(self, data):
                    return True
            f.write(yaml.dump(flat_record, Dumper=NoAliasDumper, default_flow_style=False))

    def sample(self, batch_size, in_place=True):

        self.havana_grids[self.reference_result_index].sample(batch_size)
        self.n_points += batch_size

        # Return None if we use an Havana in-place sampling approach
        if in_place:
            return self
        else:
            return self.havana_grids[self.reference_result_index].get_samples()

    def accumulate_from_havana_grid(self, grid):

        if grid.n_points == 0:
            return

        self.n_points += grid.n_points

        for index in range(0,self.n_integrands*2):
            self.havana_grids[index].merge(grid.havana_grids[index])

    def accumulate_results(self, results):
        
        if isinstance(results, Havana):
            return self.accumulate_from_havana_grid(results)

        raise NotImplementedError("Currently the Havana wrapper only supports accumulating results by merging other havana grids.")

    def sync_grids(self, max_prob_ratio=None):
        
        # Havana grids should only be updated at the end of an iteration round.
        # The current estimate gets automatically synced every time a new grid is merged into this.
        pass

    def update_iteration_count(self, dump_grids=True):
        
        import os

        self.n_iterations += 1

        # Sync the grid at this stage
        if len(self.bin_increase_factor_schedule)<self.n_iterations:
            bin_increase_factor = 1
        else:
            bin_increase_factor = int(self.bin_increase_factor_schedule[self.n_iterations-1])
        self.n_bins = bin_increase_factor*self.n_bins
        for havana_grid in self.havana_grids:
            havana_grid.update(alpha = self.learning_rate, new_bin_length=self.n_bins, train_on_avg=(not self.optimize_on_variance) )
        
        if dump_grids:
            # Dump current grid and a copy to keep track of its evolution at each iteration.
            self.dump_grid()
            file_path_for_this_iteration = os.path.join(
                os.path.dirname(self.grid_file),
                '%s_iteration_%d.%s'%(
                    '.'.join(os.path.basename(self.grid_file).split('.')[:-1]),
                    self.n_iterations,
                    os.path.basename(self.grid_file).split('.')[-1]
                )
            )
            self.dump_grid(file_name=file_path_for_this_iteration)

    def get_grid_summary(self, sort_by_variance=True, show_channel_grid=True, 
                            show_all_information_for_all_integrands=False, show_selected_phase_only = False):

        if (self.havana_grids is None):
            return ''

        from prettytable import PrettyTable
        pt = PrettyTable()
        all_results = []
        reference_result_index_here=None
        for i_havana_grid, havana_grid in enumerate(self.havana_grids):
            if show_selected_phase_only and ((self.phase=='real' and i_havana_grid%2==1) or (self.phase=='imag' and i_havana_grid%2==0)):
                continue
            avg, err, chi_sq, max_eval_negative, max_eval_positive, n_evals, n_zero_evals = havana_grid.get_current_estimate()
            integrand_result = [
                    {
                        'i_SG'             : -1,
                        'SG_name'          : 'Total',
                        'avg'              : avg,
                        'err'              : err,
                        'chi_sq_per_dof'   : chi_sq/self.n_iterations if self.n_iterations > 0 else 0.,
                        'max_wgt_infl'     : 0. if n_evals==0 else max(abs(max_eval_negative),max_eval_positive)/(err*n_evals),
                        'n_evals'          : n_evals,
                        'n_zero_evals'     : n_zero_evals,
                        'p'                : None
                    }                
            ]
            if self.SG_ids is not None:
                integrand_result.extend([ 
                        {
                            'i_SG'             : i_SG,
                            'SG_name'          : self.SG_ids[i_SG][1],
                            'avg'              : avg,
                            'err'              : err,
                            'chi_sq_per_dof'   : chi_sq/self.n_iterations if self.n_iterations > 0 else 0.,
                            'max_wgt_infl'     : 0. if n_evals==0 else max(abs(max_eval_negative),max_eval_positive)/(err*n_evals),
                            'n_evals'          : n_evals,
                            'n_zero_evals'     : n_zero_evals
                        }
                        for i_SG, (avg, err, chi_sq, max_eval_negative, max_eval_positive, n_evals, n_zero_evals) in enumerate(havana_grid.get_top_level_accumulators()) 
                    ]
                )
                cdfs = havana_grid.get_top_level_cdfs()
                pdfs = [ cdfs[i]-(0. if i==0 else cdfs[i-1]) for i in range(len(cdfs)) ]
                for i_SG, pdf in enumerate(pdfs):
                    integrand_result[i_SG+1]['p'] = pdf

            if i_havana_grid==self.reference_result_index:
                reference_result_index_here = len(all_results)
            all_results.append(('%s[I%d]'%('Re' if i_havana_grid%2==0 else 'Im', i_havana_grid//2),integrand_result) )

        # Place the reference result first
        all_results = [all_results[reference_result_index_here],]+[r for i_r, r in enumerate(all_results) if i_r!=reference_result_index_here]

        # Now sort SGs
        if self.SG_ids is not None:
            sorted_SG_ids = sorted( 
                [ 
                    ( (sg['err'] if sort_by_variance else abs(sg['avg'])), sg['i_SG']) for sg in all_results[0][1][1:]
                ], key=lambda el: el[0], reverse=True 
            )
            sorted_all_results = []
            for integrand_name, integrand_results in all_results:
                sorted_all_results.append(
                    (integrand_name, 
                        [integrand_results[0],]+
                        [
                            integrand_results[i_SG+1] for metric, i_SG in sorted_SG_ids
                        ]
                    )
                )
            all_results = sorted_all_results

        def format_entry(r, short=False, show_p=False):
            str_res = '%s%-12.6g +- %-9.3g%-8s%s'%(
                bcolors.RED if (r['avg']==0. or abs(r['err']/r['avg'])>0.01) else bcolors.GREEN,
                r['avg'],r['err'],
                '(%.2g%%)'%(0. if r['avg']==0. else abs(r['err']/r['avg'])*100.),
                bcolors.END
            )
            if short:
                return str_res
            str_res += ' | chi2=%s%.3g%s'%(
                bcolors.RED if r['chi_sq_per_dof'] > 5. else bcolors.GREEN,
                r['chi_sq_per_dof'],
                bcolors.END
            )
            str_res += '|mwi=%s%.3g%s'%(
                bcolors.RED if r['max_wgt_infl'] > 5. else bcolors.GREEN,
                r['max_wgt_infl'],
                bcolors.END
            )
            str_res += '|n=%.1fM(%.2g%%)'%(
                r['n_evals']/1.0e6,
                0. if self.n_points == 0 else (r['n_evals']/self.n_points)*100.
            )
            str_res += '|n0=%.1fM(%.2g%%)'%(
                r['n_zero_evals']/1.0e6,
                0. if r['n_evals'] == 0 else (r['n_zero_evals']/r['n_evals'])*100.
            )
            if show_p and r['p'] is not None:
                str_res += '|p=%.2g%%'%(r['p']*100.0)
            return str_res

        pt.add_column("Result",[ entry['SG_name'] for entry in all_results[0][1] ])

        for i_res, (integrand_name, integrand_results) in enumerate(all_results):
            pt.add_column(integrand_name, [ format_entry(integrand_result, short=((i_res!=0) and not show_all_information_for_all_integrands), show_p=(i_res==0) ) for integrand_result in integrand_results] )
            pt.align[integrand_name] = 'l'

        return pt.get_string()

    def get_summary(self):

        avg, err, chi_sq, max_eval_negative, max_eval_positive, n_evals, n_zero_evals = self.havana_grids[self.reference_result_index].get_current_estimate()

        # One can do a consistency check of comparing self.n_points and n_evals
        #print(avg, err, chi_sq, max_eval_negative, max_eval_positive, n_evals, n_zero_evals)
        #print(err,self.n_points)
        max_wgt_infl = 0. if self.n_points==0 else max(abs(max_eval_negative),max_eval_positive)/(err*self.n_points)

        res = ["Result after %.1fM evaluations and %d iterations: %s%.6g +/- %.4g (%.2g%%)%s, chi2=%s%.3g%s, max_wgt_infl=%s%.3g%s, zero_evals=%.3g%%"%(
            self.n_points/1.0e6, self.n_iterations, bcolors.RED if avg==0. or abs(err/avg)>0.01 else bcolors.GREEN, avg, err, 
            0. if avg==0. else (abs(err/avg))*100., bcolors.END,
            bcolors.RED if self.n_iterations==0. or chi_sq/self.n_iterations > 5. else bcolors.GREEN, chi_sq/self.n_iterations if self.n_iterations>0 else 0., bcolors.END, 
            bcolors.RED if max_wgt_infl>5. else bcolors.GREEN, max_wgt_infl, bcolors.END, 
            (n_zero_evals/float(self.n_points))*100. if self.n_points>0 else 0.
        )]
        if self.target_result is not None:
            res.append( ("%-{}s".format(len('Result after %d evaluations and %d iterations'%(self.n_points, self.n_iterations))))%("vs target")+
                ": %.5g del %.5g (%.2g%%)"%(
                    self.target_result, 
                    avg-self.target_result, 
                    0. if self.target_result==0. else ((avg-self.target_result)/self.target_result)*100.0
                ))
        for res_index, havana_grid in enumerate(self.havana_grids):
            avg, err, chi_sq, max_eval_negative, max_eval_positive, n_evals, n_zero_evals = havana_grid.get_current_estimate()
            max_wgt_infl = 0. if self.n_points==0 else max(abs(max_eval_negative),max_eval_positive)/(err*self.n_points)
            res.append('  %s %s[I%d] = %.5g +/- %.5g (%.2g%%), chi2=%.3g, max_wgt_infl=%.3g, zero_evals=%.3g%%'%(
                '  ' if res_index!=self.reference_result_index else '->',
                'Re' if res_index%2==0 else 'Im',
                res_index//2,
                avg, err, 0. if avg==0. else (abs(err/avg))*100.,
                chi_sq, max_wgt_infl, (n_zero_evals/float(self.n_points))*100. if self.n_points>0 else 0.
            ))

        return '\n'.join(res)

    def get_current_estimate(self):
        avg, err, chi_sq, max_eval_negative, max_eval_positive, n_evals, n_zero_evals = self.havana_grids[self.reference_result_index].get_current_estimate()
        return (avg, err)

    def get_n_evals(self):
        avg, err, chi_sq, max_eval_negative, max_eval_positive, n_evals, n_zero_evals = self.havana_grids[self.reference_result_index].get_current_estimate()
        return n_evals

class HavanaMockUp(object):

    def __init__(self, 
            n_dimensions, integrands, SG_ids=None, n_channels_per_SG=None, target_result=None, seed=None, 
            reference_integrand_index=0, phase='real', grid_file='havana_grid.yaml', 
            optimize_on_variance=True, max_prob_ratio=1000., fresh_integration=False,
            flat_record = None, **opts
        ):

        import numpy
        import random

        self.n_dimensions = n_dimensions
        self.target_result = target_result
        self.integrands = integrands
        self.n_integrands = len(self.integrands)
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

    def get_constructor_arguments(self, **opts):
        
        import copy
        return {
            'n_dimensions' : self.n_dimensions,
            'SG_ids' : self.SG_ids,
            'n_channels_per_SG' : self.n_channels_per_SG,
            'target_result' : self.target_result,
            'seed' : self.seed,
            'integrands' : None, # Will be filled at runtime
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

    def get_grid_summary(self, sort_by_variance=True, show_channel_grid=True, **opts):

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

    def __init__(self, n_integration_dimensions, alpha_loop_path, run_workspace, rust_input_folder, 
                        cross_section_set_file_path, E_cm, run_dir='', n_dimensions_per_SG_id=None, frozen_momenta=None,
                        run_hyperparameters_filename=None):

        import sys
        import os

        self.constructor_arguments = {
            'alpha_loop_path' : alpha_loop_path,
            'run_workspace' : run_workspace,
            'cross_section_set_file_path' : cross_section_set_file_path,
            'rust_input_folder' : rust_input_folder,
            'n_integration_dimensions' : n_integration_dimensions,
            'n_dimensions_per_SG_id' : n_dimensions_per_SG_id,
            'frozen_momenta' : frozen_momenta,
            'run_dir' : run_dir,
            'E_cm' : E_cm,
            'run_hyperparameters_filename' : run_hyperparameters_filename
        }
        self.alpha_loop_path = alpha_loop_path
        self.run_workspace = run_workspace
        self.rust_input_folder = rust_input_folder
        self.cross_section_set_file_path = cross_section_set_file_path
        self.frozen_momenta = frozen_momenta
        self.E_cm = E_cm

        # This dimensions specification is not really used for now, only the length of the set of 
        # continuous and discrete dimensions matters for now.
        self.n_integration_dimensions = n_integration_dimensions
        self.n_dimensions_per_SG_id = n_dimensions_per_SG_id

        self.frozen_momenta = frozen_momenta

        self.run_hyperparameters_filename = run_hyperparameters_filename

        # Now start a rust worker
        try:
            if self.alpha_loop_path not in sys.path:
                sys.path.insert(0, self.alpha_loop_path)
            # Import the rust bindings
            from ltd import CrossSection
        except ImportError as e:
            raise WorkerException("Could not import the rust back-end 'ltd' module in '%s'. Compile it first with:\n"%self.alpha_loop_path+
                " ./make_lib\nfrom within the alphaLoop directory.")

        os.environ['MG_NUMERATOR_PATH'] = os.path.abspath(os.path.join(self.rust_input_folder, os.path.pardir))
        hyperparameters_path = os.path.join(self.run_workspace, run_dir, self.run_hyperparameters_filename)
        if not os.path.isfile(hyperparameters_path):
            raise WorkerException("Could not find hyperparameter file at %s."%hyperparameters_path)

        try:
            self.rust_worker = CrossSection(
                cross_section_set_file_path,
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
        preconstructed_integrands=None,
        havana_mockup=False
    ):

    import time
    start_time = time.time()

    if preconstructed_integrands is None:
        integrands = [
            ALStandaloneIntegrand(**integrand_constructor_args) for integrand_constructor_args in integrands_constructor_args
        ]
    else:
        integrands =preconstructed_integrands

    #import random
    #time.sleep(20+10*random.random())
    if not havana_mockup:
        # In this case samples_batch is a complete HavanaInstance with all samples generated in-place
        havana_sampler = samples_batch
        for i_integrand, integrand in enumerate(integrands):
            if i_integrand != havana_sampler.reference_integrand_index:
                integrand.rust_worker.evaluate_integrand_havana(
                    havana_sampler.havana_grids[havana_sampler.reference_result_index],
                    havana_updater_re = havana_sampler.havana_grids[i_integrand*2],
                    havana_updater_im = havana_sampler.havana_grids[i_integrand*2 + 1]
                )
            else:
                if havana_sampler.phase == 'real':
                    integrand.rust_worker.evaluate_integrand_havana(
                        havana_sampler.havana_grids[havana_sampler.reference_result_index],
                        havana_updater_im = havana_sampler.havana_grids[i_integrand*2 + 1],
                        real_phase = True
                    )
                else:
                    integrand.rust_worker.evaluate_integrand_havana(
                        havana_sampler.havana_grids[havana_sampler.reference_result_index],
                        havana_updater_re = havana_sampler.havana_grids[i_integrand*2],
                        real_phase = False
                    )

        return (job_id, time.time()-start_time, havana_sampler)
    
    else:

        import os
        import numpy

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

        #print('run_time=%s'%str(time.time()-start_time))
        return (job_id, time.time()-start_time, all_res)

def run_job(input_payload, run_id, worker_id=-1, cache=None, output_path=None):

    if input_payload == 'ping':
        return 'pong'

    import time
    import os
    import sys

    if cache is not None:
        if cache['last_integrands_contructor_args'] is None or input_payload['integrands_constructor_args'] != cache['last_integrands_contructor_args']:
            standalone_integrands = [
                ALStandaloneIntegrand(**integrand_constructor_args) for integrand_constructor_args in input_payload['integrands_constructor_args']
            ]
            cache['standalone_integrands'] = standalone_integrands
            cache['last_integrands_contructor_args'] = input_payload['integrands_constructor_args']
        else:
            standalone_integrands = cache['standalone_integrands']
    else:
        standalone_integrands = [
            ALStandaloneIntegrand(**integrand_constructor_args) for integrand_constructor_args in input_payload['integrands_constructor_args']
        ]

    call_options = {
        'job_id' : input_payload['job_id'],
        'integrands_constructor_args' : input_payload['integrands_constructor_args'],
        'SG_ids_specified' : input_payload['SG_ids_specified'],
        'channels_specified' : input_payload['channels_specified'],
        'preconstructed_integrands' : standalone_integrands
    }

    print("Worker #%d processes a sample of %d points for job #%d..."%(worker_id, input_payload['n_points_to_sample'], call_options['job_id']))
    t0 = time.time()
    havana_constructor_arguments = input_payload['havana_grid']
    havana_constructor_arguments['integrands'] = standalone_integrands
    if input_payload['havana_mockup']:
        havana_sampler = HavanaMockUp(**havana_constructor_arguments)
    else:
        if havana_constructor_arguments['flat_record']['tmp_folder'] is None:
            print("Worker #%d loads the havana sampler grid from files %s"%(
                worker_id, ', '.join(os.path.basename(fp) for fp in havana_constructor_arguments['flat_record']['grid_files'])
            ))
        else:
            print("Worker #%d loads the havana sampler grid from serialized records using temporary file."%(worker_id))
        t0 = time.time()
        havana_constructor_arguments['flat_record']['n_points'] = 0
        if output_path is None:
            assert(havana_constructor_arguments['flat_record']['tmp_folder'] is not None)
            havana_constructor_arguments['grid_file'] =  os.path.join(
                havana_constructor_arguments['flat_record']['tmp_folder'],
                'run_%d_job_output_%d_grids.yaml'%(run_id, input_payload['job_id'] )
            )
        else:
            havana_constructor_arguments['grid_file'] =  os.path.join(
                os.path.dirname(output_path),
                'run_%d_job_output_%d_grids.yaml'%(run_id, input_payload['job_id'] )
            )
        havana_sampler = Havana(**havana_constructor_arguments)
        havana_loading_time = time.time()-t0
        print("Worker #%d loaded the havana grids in %.5fs"%(worker_id, havana_loading_time))

        
    this_sample = havana_sampler.sample(input_payload['n_points_to_sample'])
    sample_generation_time = time.time()-t0
    print("Worker #%d completed sample generation in %.5fs."%(worker_id, sample_generation_time))
    call_options['samples_batch'] = this_sample
    call_options['havana_mockup'] = input_payload['havana_mockup']
    t0 = time.time()
    job_id, run_time, samples_evaluated = HavanaIntegrandWrapper(**call_options)
    print("Worker #%d evaluated all samples in %.5s."%(worker_id, time.time()-t0))

    if not input_payload['havana_mockup']:
        # In this case, samples_evaluated is a complete Havana instance with all evaluations accumulated in-place
        havana_updater = samples_evaluated
        print("Worker #%d now dumps havana updater grid ..."%(worker_id))
        t0 = time.time()
        havana_updater_constructor_arguments = havana_updater.get_constructor_arguments(tmp_folder=havana_constructor_arguments['flat_record']['tmp_folder'])
        havana_dumping_time = time.time()-t0
        if havana_constructor_arguments['flat_record']['tmp_folder'] is None:
            print("Worker #%d dumped Havana updater grids in %.5fs to files %s"%(worker_id, 
                havana_dumping_time,
                ', '.join(os.path.basename(fp) for fp in havana_updater_constructor_arguments['flat_record']['grid_files'])
            ))
        else:
            print("Worker #%d serialised results as a havana grid using a temporary local file."%(worker_id,))
        complete_result = {
            'havana_updater' : havana_updater_constructor_arguments,
        }
        return ( job_id, havana_loading_time+run_time+sample_generation_time+havana_dumping_time, complete_result )

    else:
        t0 = time.time()
        havana_updater_constructor_args = dict(input_payload['havana_grid'])
        havana_updater_constructor_args['fresh_integration'] = True
        havana_updater_constructor_args['integrands'] = standalone_integrands
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
        print("Worker #%d has run_time=%.2f s, sample_generation_time=%.2f s, post_processing_time=%.2f s"%(
            worker_id, run_time, sample_generation_time, post_processing_time
        ))
        
        return ( job_id, run_time+sample_generation_time+post_processing_time, complete_result )

def main(arg_tuple):

    run_id, worker_id, run_workspace, timeout = arg_tuple
    import pickle
    import time
    import os
    import sys

    cache = {
        'standalone_integrands' : None,
        'last_integrands_contructor_args' : None
    }

    input_path = os.path.join(run_workspace,'run_%d_job_input_worker_%d.pkl'%(run_id, worker_id))
    input_path_done = os.path.join(run_workspace,'run_%d_job_input_worker_%d.done'%(run_id, worker_id))
    output_path = os.path.join(run_workspace,'run_%d_job_output_worker_%d.pkl'%(run_id, worker_id))
    output_path_done = os.path.join(run_workspace,'run_%d_job_output_worker_%d.done'%(run_id, worker_id))
    while True:

        try:
            if (not os.path.exists(input_path_done) or not os.path.exists(input_path)) or os.path.exists(output_path_done):
                time.sleep(timeout)
                continue

            t_overall = time.time()
            print("Worker #%d received a new job at time t=%s"%(worker_id, time.time()))
            sys.stdout.flush()
            t0=time.time()
            input_payload = pickle.load( open( input_path, 'rb' ) )
            pickling_time = time.time() - t0
            print("Worker #%d unpickle job payload in %.5f s."%(worker_id, pickling_time))

            output = run_job(input_payload, run_id, worker_id=worker_id, cache=cache, output_path = output_path)
            if output=='pong':
                with open(output_path,'wb') as f:
                    pickle.dump( output, f )
                with open(output_path_done,'w') as f:
                    f.write("Marker file for job output done.")
                continue

            print("Worker #%d now writing job output."%(worker_id))
            cumulative_measured_timing = output[1]
            # Overwrite with actual timing
            output = (output[0], time.time()-t_overall, output[2])
            t0 = time.time()
            with open(output_path,'wb') as f:
                pickle.dump( output, f )
            with open(output_path_done,'w') as f:
                f.write("Marker file for job output done.")
            print("Worker #%d finished serialisation and dump of job output in %.5gs."%( worker_id, time.time()-t0 ))
            print("Worker #%d exact timing vs cumulative one: %.5fs vs %.5fs."%(worker_id, time.time()-t_overall, cumulative_measured_timing) )
            print("Worker #%d completed job at time t=%s"%(worker_id, time.time()))
            sys.stdout.flush()

        except Exception as e:
            print("The following exception was encountered in run #%d and worker #%d: %s"%(run_id, worker_id, str(e)))
            break
        except KeyboardInterrupt as e:
            #print("Keyboard interrupts in run #%d and worker #%d."%(run_id, worker_id))
            break


def example_job(arg):
    import time
    time.sleep(3.)
    return 'DONE '+str(arg)

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
        except ImportError as e:
            #VHTOFIX
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
