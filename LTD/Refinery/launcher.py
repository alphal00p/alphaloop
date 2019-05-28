#!/usr/bin/env python
import os
import sys
import copy
import subprocess
import time
import shutil
import math

pjoin = os.path.join
root_path = os.path.dirname(os.path.realpath( __file__ ))
sys.path.insert(0, pjoin(root_path,os.pardir))
sys.path.insert(0, pjoin(root_path,os.pardir,os.pardir))

rust_executable_path = os.path.abspath(pjoin(root_path,os.pardir,os.pardir,'rust_backend','target','release','ltd'))
import yaml
from yaml import Loader, Dumper

import ltd_commons
import ltd_utils
from ltd_utils import Colour

_RUN_LOCALLY = True
_CLEAN = False
_SILENCE = False
_N_CORES = 8
_CLEAN = False 
_NAME = None
_IDS = None

_N_REFINES = 10
_REFINE_N_POINTS = int(1e7)

_N_ITERATIONS = 10
_N_START = int(1e6)
_N_INCREASE = int(1e6)

_TOPOLOGY = None
_SEED_START = 1
_PHASE = 'real'

_WALL_TIME = 24 # in hours 

general_hyperparams = copy.deepcopy(ltd_commons.hyperparameters)
general_hyperparams['General']['statistics_interval'] = 100000000
general_hyperparams['General']['screen_log_core'] = 1
general_hyperparams['General']['absolute_precision'] = 1.0e+5
general_hyperparams['General']['relative_precision'] = 5.
general_hyperparams['General']['integration_statistics'] = True 
general_hyperparams['General']['numerical_instability_check'] = True 
general_hyperparams['General']['unstable_point_warning_percentage'] = 10.
general_hyperparams['General']['num_digits_different_for_inconsistency'] = 10.
general_hyperparams['General']['minimal_precision_for_returning_result'] = 2.

general_hyperparams['Parameterization']['mode'] = 'spherical'
general_hyperparams['Parameterization']['mapping'] = 'log'
general_hyperparams['Parameterization']['b'] = 0.1

general_hyperparams['Integrator']['integrator'] = 'vegas'
general_hyperparams['Integrator']['n_start'] = int(1e5)
general_hyperparams['Integrator']['n_max'] = int(1e6)
general_hyperparams['Integrator']['n_increase'] = int(1e5)
general_hyperparams['Integrator']['seed'] = 0

general_hyperparams['Integrator']['eps_rel'] = 1.0e-99
general_hyperparams['Integrator']['eps_abs'] = 0.
general_hyperparams['Integrator']['border'] = 1.0e-3
general_hyperparams['Integrator']['maxpass'] = 5
general_hyperparams['Integrator']['maxchisq'] = 10.
general_hyperparams['Integrator']['mindeviation'] = 0.25

general_hyperparams['Integrator']['reset_vegas_integrator'] = True
general_hyperparams['Integrator']['use_only_last_sample'] = False

general_hyperparams['input_rescaling'] = [
            [[0., 1.], [0., 1.], [0., 1.]],
            [[0., 1.], [0., 1.], [0., 1.]],
            [[0., 1.], [0., 1.], [0., 1.]],
            [[0., 1.], [0., 1.], [0., 1.]],
        ]

# shift the loop momenta. the first component is a rescaling of the radius
general_hyperparams['shifts'] = [
            [1., 0., 0., 0.],
            [1.1, 0., 0., 0.],
            [1.2, 0., 0., 0.],
            [1.3, 0., 0., 0.],
        ]

def load_results_from_yaml(log_file_path):
    """Load a full-fledged scan from a yaml dump"""

    raw_data = yaml.load(open(log_file_path,'r'), Loader=Loader)
    processed_data = {'run_results': []}
    for entry_name, entry_value in raw_data:
        if entry_name == 'run_result':
            processed_data['run_results'].append(entry_value)
        else:
            processed_data[entry_name] = entry_value
    return processed_data

def run_topology(topo, dir_name, run_options, result_path, job_name_suffix='', local=True):
    """ Run topology of specified index and directory locally or on a SLURM scheduled cluster."""

    run_options = dict(run_options)
    output_log_path = run_options.pop('output_log_path', '%s/LTD_runs/logs/%s%s.out'%(os.environ.get('SCRATCH','/tmp'), dir_name, job_name_suffix))
    error_log_path = run_options.pop('error_log_path', '%s/LTD_runs/logs/%s%s.err'%(os.environ.get('SCRATCH','/tmp'), dir_name, job_name_suffix))

    cmd = [ rust_executable_path, 
            '-t','%s'%topo.name, 
            '-f','%s'%pjoin(root_path,dir_name,'hyperparameters.yaml'), 
            '-l','%s'%pjoin(root_path,dir_name,'topologies.yaml'),
            '-c','%d'%_N_CORES
    ]
    cmd.extend(sum([['--%s'%opt,'%s'%str(value)] for opt, value in sorted(run_options.items())],[]))
  
    if os.path.exists(result_path):
        os.remove(result_path)

    #print(' '.join(cmd))
    if _RUN_LOCALLY:
    	with ltd_utils.Silence(active=_SILENCE):
        	subprocess.call(cmd, cwd=pjoin(root_path,dir_name))
        	# Sleep for a very short while to allow output file flushing
        	time.sleep(0.3)
    else:
        submission_script = open(pjoin(root_path,'submission_template.run'),'r').read()
        open(pjoin(root_path,'submitter.run'),'w').write(submission_script%{
		'job_name' : '%s%s'%(topo.name, job_name_suffix),
                'n_hours' : _WALL_TIME,
                'n_cpus_per_task' :_N_CORES,
                'output' : output_log_path,
                'error' : error_log_path,
                'executable_line' : ' '.join(cmd)
	})
        subprocess.call(['sbatch','submitter.run'], cwd=root_path)
        return    

    if not os.path.isfile(result_path):
        print("Error: Run did not successfully complete as the results yaml dump '%s' could not be found."%result_path )
    else:
        result = yaml.load(open(result_path,'r'), Loader=Loader)
        analytic_result = topo.analytic_result.real if _PHASE=='real' else topo.analytic_result.imag
        print ">> Analytic result : %.16e"%analytic_result
        print ">> LTD result      : %.16e"%(result['result'][0])
        print ">> LTD error       : %.16e"%(result['error'][0])
        n_sigmas = abs((analytic_result-result['result'][0])/result['error'][0])
        print ">> LTD discrepancy : %s%.2f sigmas%s"%(Colour.GREEN if n_sigmas <=3. else Colour.RED,  n_sigmas, Colour.END)
        print ">> LTD rel. discr. : %.2g%%"%(100.0*abs((analytic_result-result['result'][0])/analytic_result))
        print ">> LTD n_points    : %dM"%(int(result['neval']/1.e6))        
        return result

def gather_result(topo, dir_name, clean=False):
    """ Combine all results into an 'final_result' file."""
   
    analytic_result = topo.analytic_result.real if \
                    abs(topo.analytic_result.real)>abs(topo.analytic_result.imag) else topo.analytic_result.imag
    
    all_res_files = []
    all_results = {}
    
    all_res_lines = []    
    n_eval_tot = 0
    
    for filename in os.listdir(pjoin(root_path,dir_name)):
        if not filename.startswith('refine_'):
            continue
        if filename.startswith('refine_grid') or filename.startswith('refine_runs'):
            continue
        refine_id = int(filename.split('_')[1])
        result_path = pjoin(root_path,dir_name,filename)
        all_res_files.append((refine_id,result_path)) 
    all_res_files = sorted(all_res_files, key=lambda el: el[0])
   
    for refine_id, result_path in all_res_files:
        result = yaml.load(open(result_path,'r'), Loader=Loader)
        n_sigmas = abs((analytic_result-result['result'][0])/result['error'][0])
        n_eval_tot += result['neval']
        these_res_lines = [
            "Results from refine #%d of topology %s"%(refine_id, topo.name), 
            ">> Analytic result : %.16e"%analytic_result,
            ">> LTD result      : %.16e"%(result['result'][0]),
            ">> LTD error       : %.16e"%(result['error'][0]),
            ">> LTD discrepancy : %s%.2f sigmas%s"%(Colour.GREEN if n_sigmas <=3. else Colour.RED,  n_sigmas, Colour.END),
            ">> LTD rel. discr. : %.2g%%"%(100.0*abs((analytic_result-result['result'][0])/analytic_result)),
            ">> LTD n_points    : %dM"%(int(result['neval']/1.e6)),
        ]
        if not _SILENCE:
            print '\n'.join(these_res_lines)
        all_res_lines.extend(these_res_lines)
        all_results[refine_id] = (result['result'][0], result['error'][0])

    # Now aggregate the results nicely and dump them in final_result.dat
    final_central_value = 0.0
    final_error = 0.0
    for (central, error) in all_results.values():
        final_central_value += central/(error**2)
        final_error += 1./(error**2)
    final_central_value /= final_error
    final_error = 1./ math.sqrt(final_error)

    n_sigmas = abs((analytic_result-final_central_value)/final_error)
    final_res_lines = [
            "%.16e %.16e"%(final_central_value,final_error),
            "Final result for topology %s"%topo.name, 
            ">> Analytic result : %.16e"%analytic_result,
            ">> LTD result      : %.16e"%final_central_value,
            ">> LTD error       : %.16e"%final_error,
            ">> LTD discrepancy : %s%.2f sigmas%s"%(Colour.GREEN if n_sigmas <=3. else Colour.RED,  n_sigmas, Colour.END),
            ">> LTD rel. discr. : %.2g%%"%(100.0*abs((analytic_result-final_central_value)/analytic_result)),
            ">> LTD n_points    : %dM"%(n_eval_tot/1.e6),
    ]
    if not _SILENCE:
        print '='*50
        print '\n'.join(final_res_lines[1:])
        print '='*50

    final_result = open(pjoin(root_path,dir_name,'final_result.dat'),'w')
    final_result.write('\n'.join(final_res_lines+these_res_lines))
    final_result.close()

    if _CLEAN:
        for _, filepath in all_res_files:
            os.remove(filepath)

if __name__ == '__main__':

    n_cores_user_set = False
    # Parse options
    processed_args = []
    for arg in sys.argv[1:]:
        if arg.startswith('--'):
            try:
                key, value = arg.split('=')
            except ValueError:
                key = args
                value = None
            key = key[2:]        
        else:
            processed_args.append(arg)
            continue

        if key=='dir':
            _NAME = value
        elif key=='cores':
            _N_CORES = int(value)
            n_cores_user_set = True
        elif key=='cores':
            _N_CORES = int(value)
            n_cores_user_set = True
        elif key=='quiet':
            _SILENCE = True
        elif key=='topology':
            _TOPOLOGY = value
        elif key=='seed':
            _SEED_START = int(value)
        elif key=='phase':
            _PHASE = value
        elif key=='clean':
            _CLEAN = True
        elif key=='n_iterations':
            _N_ITERATIONS = int(value)
        elif key=='n_refines':
            _N_REFINES = int(value)
        elif key=='refine_n_points':
            _REFINE_N_POINTS = int(value)
        elif key=='n_start':
            _N_START = int(value)
        elif key=='n_increase':
            _N_INCREASE = int(value)
        elif key=='cluster':
	    _RUN_LOCALLY = False
            if not n_cores_user_set:
                _N_CORES = 36
        else:
            print('Error: options %s not reckognised.'%arg)
            sys.exit(1)

    if _TOPOLOGY is None:
        print('Error: you must specify a topology with the --topology option.')
        sys.exit(1)
    if _NAME is None:
        print('Error: you must specify a run directory with the --dir option.')
        sys.exit(1)

    topology = ltd_utils.TopologyCollection.import_from(os.path.join(root_path, _NAME,'topologies.yaml'))[_TOPOLOGY]
    
    this_params = copy.deepcopy(general_hyperparams)
    this_params['Integrator']['integrated_phase'] = _PHASE

    rust_run_options = {}

    if _RUN_LOCALLY:
        running = "running"
    else:
        running = "launching job for"

    subcommands = [arg for arg in sys.argv[1:] if not arg.startswith('--')]
    if len(subcommands)<1:
        subcommand = 'survey'
    else:
        subcommand = subcommands[0]

    def compute_n_max(n_start, n_increase, n_iterations):
        n_max = 0
        for i_iteration in range(n_iterations):
            n_max += n_start + i_iteration*n_increase
        return n_max

    if subcommand == 'survey':

        this_params['Integrator']['n_start'] = _N_START
        this_params['Integrator']['n_increase'] = _N_INCREASE
        n_max = compute_n_max(_N_START, _N_INCREASE, _N_ITERATIONS)
        this_params['Integrator']['n_max'] = n_max
        this_params['Integrator']['keep_state_file'] = True
        this_params.export_to(os.path.join(root_path, _NAME, 'hyperparameters.yaml'))

        rust_run_options['log_file_prefix'] = pjoin(root_path, _NAME, 'integration_statistics', 'survey_')
        rust_run_options['res_file_prefix'] = pjoin(root_path, _NAME, 'survey_')
        if os.path.exists(pjoin(root_path, _NAME, 'survey_grid_%s_state.dat'%topology.name)):
            os.remove(pjoin(root_path, _NAME, 'survey_grid_%s_state.dat'%topology.name))
        rust_run_options['state_filename_prefix'] = pjoin(root_path, _NAME, 'survey_grid_')
        rust_run_options['seed'] = _SEED_START

        print "Now %s survey stage for %s (n_start=%dM, n_increase=%dM, n_iteration=%dM => n_max=%dM)"%(
            running, _TOPOLOGY, int(_N_START/1.e6), int(_N_INCREASE/1.e6), int(_N_ITERATIONS/1.e6), int(n_max/1.e6))
        run_topology(topology, _NAME, rust_run_options, pjoin(root_path,_NAME,'survey_%s_res.dat'%topology.name), 
                     job_name_suffix='_survey', local=True)

    elif subcommand == 'refine':
        
        this_params['Integrator']['n_start'] = _REFINE_N_POINTS
        this_params['Integrator']['n_increase'] = 0
        this_params['Integrator']['n_max'] = _REFINE_N_POINTS
        this_params['Integrator']['keep_state_file'] = False

        this_params.export_to(os.path.join(root_path, _NAME, 'hyperparameters.yaml'))

        for i_refine in range(1, _N_REFINES+1):
            # Copy the grid for the corresponding refine run
            shutil.copy(pjoin(root_path,_NAME,'survey_grid_%s_state.dat'%topology.name),
                        pjoin(root_path,_NAME,'refine_grid_%d_%s_state.dat'%(i_refine, topology.name)))
            rust_run_options['log_file_prefix'] = pjoin(root_path,_NAME,'integration_statistics', 'refine_%d_'%i_refine)
            rust_run_options['res_file_prefix'] = pjoin(root_path, _NAME, 'refine_%d_'%i_refine)
            rust_run_options['state_filename_prefix'] = pjoin(root_path, _NAME, 'refine_grid_%d_'%i_refine)            
            rust_run_options['seed'] = _SEED_START + i_refine 
            print "Now %s #%d refine for %s with %dM points."%(running, i_refine, _TOPOLOGY, int(_REFINE_N_POINTS/1.e6))
            run_topology(topology, _NAME, rust_run_options, pjoin(root_path,_NAME,'refine_%d_%s_res.dat'%(i_refine, topology.name)), 
                     job_name_suffix='_refine_%d'%(i_refine), local=True)
            if os.path.exists(pjoin(root_path,_NAME,'refine_grid_%d_%s_state.dat'%(i_refine, topology.name))):
                os.remove(pjoin(root_path,_NAME,'refine_grid_%d_%s_state.dat'%(i_refine, topology.name)))

    elif subcommand == 'gather':
        gather_result(topology, _NAME, clean=_CLEAN)

    else:
        print "Unknown subcommand '%s'"%subcommand
        sys.exit(1)
