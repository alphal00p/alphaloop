#!/usr/bin/env python
import os
import sys
import copy
import subprocess
import time
import math
import glob

pjoin = os.path.join
root_path = os.path.dirname(os.path.realpath( __file__ ))
sys.path.insert(0, pjoin(root_path,os.pardir))
sys.path.insert(0, pjoin(root_path,os.pardir,os.pardir))

rust_executable_path = os.path.abspath(pjoin(root_path,os.pardir,'rust_backend','target','release','ltd'))
pyvegas_executable_path = os.path.abspath(pjoin(root_path,os.pardir,'pyvegas.py'))
import yaml
from yaml import Loader, Dumper

import ltd_commons
import ltd_utils
from ltd_utils import Colour

from pprint import pprint, pformat

_PYVEGAS = False
_INTEGRATOR = 'vegas'
_RUN_LOCALLY = True
_CLEAN = False
_COLLECT = False
_SILENCE = False
_N_CORES = 8
_CLEAN = False 
_PREFIX = ''
_IDS = None
_TOPOLOGY = 'Box'
_RESULTS = None
_TARGET_ACCURACY = 1.0e-2
_INCREMENT = int(1e6)
_N_START = int(1e5)
_N_INCREASE = int(1e4)
_WALL_TIME = 24 # in hours 
_ACCOUNT = 'eth5e'

_PHASES = ['real','imag']

pyvegas_hyperparms = {'survey_iter': 5, 'survey_neval': int(1e5), 'refine_iter': 10, 'refine_neval': int(1e6), 'phase': 'imag'}

general_hyperparams = copy.deepcopy(ltd_commons.hyperparameters)
general_hyperparams['General']['multi_channeling'] = True
general_hyperparams['Integrator']['reset_vegas_integrator'] = True
#general_hyperparams['General']['absolute_precision'] = 1.0e+5
#general_hyperparams['Integrator']['integrator'] = 'vegas'
#general_hyperparams['Integrator']['n_start'] = int(1e5)
#general_hyperparams['Integrator']['n_max'] = int(1e6)
#general_hyperparams['Integrator']['n_increase'] = int(1e5)
#general_hyperparams['Integrator']['seed'] = 0
#general_hyperparams['Integrator']['state_filename'] = None
#general_hyperparams['Parameterization']['mode'] = 'spherical'
#general_hyperparams['Parameterization']['mapping'] = 'linear'
#general_hyperparams['Parameterization']['b'] = 1.0
general_hyperparams['General']['integration_statistics'] = False
#general_hyperparams['General']['statistics_interval'] = 100000000
#general_hyperparams['General']['screen_log_core'] = 1
#general_hyperparams['General']['numerical_instability_check'] = True 
#general_hyperparams['General']['unstable_point_warning_percentage'] = 200.
#general_hyperparams['General']['num_digits_different_for_inconsistency'] = 200.
#general_hyperparams['General']['minimal_precision_for_returning_result'] = 2.
#general_hyperparams['Integrator']['eps_rel'] = 1.0e-3
#general_hyperparams['Integrator']['eps_abs'] = 0.
#general_hyperparams['Integrator']['border'] = 1.0e-3
#general_hyperparams['Integrator']['maxpass'] = 5
#general_hyperparams['Integrator']['maxchisq'] = 10.
#general_hyperparams['Integrator']['mindeviation'] = 0.25


#general_hyperparams['input_rescaling'] = [
#            [[0., 1.], [0., 1.], [0., 1.]],
#            [[0., 1.], [0., 1.], [0., 1.]],
#            [[0., 1.], [0., 1.], [0., 1.]],
#            [[0., 1.], [0., 1.], [0., 1.]],
#        ]

# shift the loop momenta. the first component is a rescaling of the radius
#general_hyperparams['shifts'] = [
#            [1., 0., 0., 0.],
#            [1., 0., 0., 0.],
#            [1., 0., 0., 0.],
#            [1., 0., 0., 0.],
#        ]

def load_results_from_yaml(log_file_path):
    """Load a full-fledged scan from a yaml dump"""

    raw_data = yaml.load(open(log_file_path,'r'), Loader=Loader)
    processed_data = {'channel_data': []}
    for entry_name, entry_value in raw_data:
        if entry_name == 'channel_data':
            processed_data['channel_data'].append(entry_value)
        else:
            processed_data[entry_name] = entry_value
    return processed_data

def run_topology(topo,dir_name, index, n_hours, local=True):
    """ Run topology of specified index and directory locally or on a SLURM scheduled cluster."""
    if _RUN_LOCALLY:
       print("Now running %s topology with channel #%d"%(dir_name, index))
    else:
       print("Now launch job for %s topology with channel #%d"%(dir_name, index))
        
    if os.path.exists(pjoin(root_path,dir_name,'amplitudes.yaml')):
        amplitude_path = pjoin(root_path,dir_name,'amplitudes.yaml')
    else:
        if not os.path.exists(pjoin(root_path,'amplitudes.yaml')):
            print("ERROR: make sure the file 'amplitude.yaml' can be found at '%s'."%pjoin(root_path,'amplitudes.yaml'))
            sys.exit(1)
        amplitude_path = pjoin(root_path,'amplitudes.yaml')

    if _PYVEGAS:
        cmd = [ pyvegas_executable_path,
                '%d'%topo.n_loops,
                '-t','%s'%_TOPOLOGY,
                '-hf','%s'%pjoin(root_path,dir_name,'hyperparameters_channel_%d.yaml'%index),
                '-f','%s'%pjoin(root_path,dir_name,'topologies.yaml'),
                '-c','%d'%_N_CORES,
                '--si','%d'%pyvegas_hyperparms['survey_iter'],
                '--sn','%d'%pyvegas_hyperparms['survey_neval'],
                '--ri','%d'%pyvegas_hyperparms['refine_iter'],
                '--rn','%d'%pyvegas_hyperparms['refine_neval'],
                '--phase',pyvegas_hyperparms['phase'],
                '--out',pjoin(root_path, dir_name, topo.name + '_res.dat')
        ]
    else:
        cmd = [ rust_executable_path, 
                '-t','%s'%_TOPOLOGY,
                '-f','%s'%pjoin(root_path,dir_name,'hyperparameters_channel_%d.yaml'%index),
                '-l','%s'%pjoin(root_path,dir_name,'topologies.yaml'),
                '-c','%d'%_N_CORES,
                '-p','%s'%amplitude_path
        ]

    try:
        os.remove(pjoin(root_path,dir_name,'scan_%d_state.dat'%index))
    except:
        pass
    if _RUN_LOCALLY:
        # print(' '.join(cmd),pjoin(root_path,dir_name))        
    	with ltd_utils.Silence(active=_SILENCE):
            subprocess.call(cmd, cwd=pjoin(root_path,dir_name))
            # Sleep for a very short while to allow output file flushing
            time.sleep(1.0)
    else:
        submission_script = open(pjoin(root_path,os.pardir,'submission_template.run'),'r').read()
        open(pjoin(root_path,dir_name,'submitter_%d.run'%index),'w').write(submission_script%{
		'job_name' : 'channel_%d_%s'%(index, _TOPOLOGY),
                'n_hours' : n_hours,
                'n_cpus_per_task' :_N_CORES,
                'output' : '%s/LTD_runs/logs/channel_%d_%s.out'%(os.environ['SCRATCH'], index, _TOPOLOGY),
                'error' : '%s/LTD_runs/logs/channel_%d_%s.err'%(os.environ['SCRATCH'], index, _TOPOLOGY),
                'account' : _ACCOUNT,
                'executable_line' : ' '.join(cmd)
	})
        time.sleep(3.0)
        subprocess.call(['sbatch',pjoin(root_path,dir_name,'submitter_%d.run'%index)], cwd=root_path)
        return    

    result_path = pjoin(root_path,dir_name,'channel_%d_%s_res.dat'%(index,_TOPOLOGY))
    if not os.path.isfile(result_path):
        print(("Error: Run did not successfully complete as the results yaml dump '%s' could not be found."%result_path ))
    else:
        str_data = open(result_path,'r').read()
        # Some weird I/O issue of rust
        str_data = '...'.join(str_data.split('...')[:1])
        open(result_path,'w').write(str_data)
        result = yaml.load(open(result_path,'r'), Loader=Loader)
        return (result['result'][0],result['error'][0])

def get_n_channels(topology):
    n_propagators_per_loop_line = [len(ll.propagators) for ll in topology.loop_lines]
    n_channels = 0
    for cs in topology.ltd_cut_structure:
        n_channels_for_this_cut = 1
        for i_ll, cut_sign in enumerate(cs):
            if cut_sign != 0:
                n_channels_for_this_cut *= n_propagators_per_loop_line[i_ll]
        n_channels += n_channels_for_this_cut
    return n_channels

def combine_results(results, verbose=True, individual_channel_results=True):

    max_error = None
    max_central = None
    central = 0.
    error = 0.
    n_tot_points = 0
    for i_channel, channel_result in enumerate(results['channel_results']):
        if verbose and individual_channel_results:
            print(("Result for channel %d = %.8e +/- %.3e (n_points=%.1fM)"%(
                i_channel, channel_result[0][0], channel_result[0][1], channel_result[1]/1.0e6)))
        central += channel_result[0][0]
        error += channel_result[0][1]**2
        n_tot_points += channel_result[1]
        if max_error is None or channel_result[0][1] > max_error[0]:
            max_error = (channel_result[0][1], i_channel)
        if max_central is None or abs(channel_result[0][0]) > abs(max_central[0]):
            max_central = (channel_result[0][0], i_channel)
    error = math.sqrt(error)

    res = {
        'central_value' : central,
        'error'         : error,
        'max_error'     : max_error,
        'max_central'   : max_central,
        'n_tot_points'  : n_tot_points
    }
    if verbose:
        n_sigmas = abs(res['central_value']-results['analytic_result'])/res['error'] if res['error'] != 0. else 0.
        print('%sresult: %.8e +/- %.3e (%.2g%%) vs %.8e (n_points=%.1fM) n_sigmas=%.2g %s'%(
            Colour.RED if n_sigmas > 3.0 else Colour.GREEN,
            res['central_value'], res['error'], 
            (res['error']/abs(res['central_value']))*100.0 if res['central_value']!=0. else 0.,
            results['analytic_result'],
            res['n_tot_points']/1.0e6,
            n_sigmas,
            Colour.END
        ))
        if res['max_central'] is not None:
            print(('Max contribution (%.4e) from channel #%d'%(res['max_central'][0], res['max_central'][1])))
        if res['max_error'] is not None:
            print(('Max error (%.4e) from channel #%d'%(res['max_error'][0], res['max_error'][1])))

    return res

if __name__ == '__main__':

    n_cores_user_set = False
    # Parse options
    processed_args = []
    for arg in sys.argv[1:]:
        if arg.startswith('--'):
            try:
                key, value = arg.split('=')
            except ValueError:
                key = arg
                value = None
            key = key[2:]        
        else:
            processed_args.append(arg)
            continue

        if key=='prefix':
            _PREFIX = value
        elif key=='collect':
            _COLLECT = True
        elif key=='wall_time':
            _WALL_TIME = int(value)
        elif key=='target_accuracy':
            _TARGET_ACCURACY = float(value)
        elif key=='topology':
            _TOPOLOGY = value
        elif key=='increment':
            _INCREMENT = int(eval(value.replace('M','*1000000')))
        elif key=='n_start':
            _N_START = int(eval(value.replace('M','*1000000')))
        elif key=='integrator':
            _INTEGRATOR = value
        elif key=='n_increase':
            _N_INCREASE = int(eval(value.replace('M','*1000000')))
        elif key=='phases':
            _PHASES = eval(value)
        elif key=='results':
            _RESULTS = value
        elif key=='ids':
            _IDS = eval(value)

        elif key=='cores':
            _N_CORES = int(value)
            n_cores_user_set = True
        else:
            print(('Error: options %s not reckognised.'%arg))
            sys.exit(1)

    dir_name = processed_args[0]

    if _COLLECT:
        topology = ltd_utils.TopologyCollection.import_from(os.path.join('%s_%s'%(dir_name, _PHASES[0]),'topologies.yaml'))[_TOPOLOGY]
        channel_results = {}
        channel_results[_TOPOLOGY] = {
                'real' : {
                    'analytic_result' : topology.analytic_result.real,
                    'channel_results' : [],
                },
                'imag' : {
                    'analytic_result' : topology.analytic_result.imag,                    
                    'channel_results' : [],        
                },
            }
        n_channels = get_n_channels(topology)
        for phase in _PHASES:
            run_name = '%s_%s'%(dir_name, phase)
            for channel_ID in range(n_channels):
                result_path = pjoin(root_path,run_name,'channel_%d_%s_res.dat'%(channel_ID,_TOPOLOGY))
                if not os.path.isfile(result_path):
                    print(("Error: Run did not successfully complete as the results yaml dump '%s' could not be found."%result_path ))
                    sys.exit(1)
                str_data = open(result_path,'r').read()
                # Some weird I/O issue of rust
                str_data = '...'.join(str_data.split('...')[:1])
                open(result_path,'w').write(str_data)
                result = yaml.load(open(result_path,'r'), Loader=Loader)
                channel_results[_TOPOLOGY][phase]['channel_results'].append( [ (result['result'][0],result['error'][0]), result['neval']] )
        
        if not _RESULTS:
            _RESULTS = dir_name+'_results.dat'

        save_path = _RESULTS
        print("Saving results to %s"%save_path)
        open(save_path,'w').write(pformat(channel_results))
        
    elif _RESULTS and os.path.exists(_RESULTS):
        print("Loading results from '%s'."%(_RESULTS))
        channel_results = eval(open(_RESULTS,'r').read())
    else:
        if not _RESULTS:
            _RESULTS = dir_name+'_results.dat'
        channel_results = {}

    if not all(os.path.exists(dir_name+'_%s'%phase) for phase in _PHASES):
        print("The first argument must be an existing directory path when suffixed with phases.")
        sys.exit(1)

    if any(arg=='quiet' for arg in processed_args):
        _SILENCE = True
        
    if any(arg=='cluster' for arg in processed_args):
        _RUN_LOCALLY = False
        if not n_cores_user_set:
            _N_CORES = 36

    if any(arg=='clean' for arg in processed_args):
        _CLEAN = True

    if any(arg=='gather' for arg in processed_args):
        
        print("")
        for topology, results in list(channel_results.items()):
            for phase, phase_results in list(results.items()):
                print("")
                print(">>> Final results for topology %s and phase %s:"%(topology, phase))
                combined_results = combine_results(phase_results, verbose=True)

    if any(arg=='run' for arg in processed_args):
            
        # First refresh configuration files
        hyperparams = copy.deepcopy(general_hyperparams)
       
        topology = ltd_utils.TopologyCollection.import_from(os.path.join('%s_%s'%(dir_name, _PHASES[0]),'topologies.yaml'))[_TOPOLOGY]
        
        if _TOPOLOGY not in channel_results:
            channel_results[_TOPOLOGY] = {
                'real' : {
                    'analytic_result' : topology.analytic_result.real if topology.analytic_result else 0.,
                    'channel_results' : [],
                },
                'imag' : {
                    'analytic_result' : topology.analytic_result.imag if topology.analytic_result else 0., 
                    'channel_results' : [],        
                },
            }

        for phase in _PHASES:

            run_name = '%s_%s'%(dir_name, phase)
            if _CLEAN:
                for f in glob.glob(pjoin(run_name,'*_state.dat')):
                    os.remove(f)
                for f in glob.glob(pjoin(run_name,'*_res.dat')):
                    os.remove(f)
                for f in glob.glob(pjoin(run_name,'integration_statistics','*')):
                    os.remove(f)
            
            n_channels = get_n_channels(topology)
            
            hyperparams['Integrator']['integrated_phase'] = phase

            mode = 'REFINE'
            if len(channel_results[_TOPOLOGY][phase]['channel_results'])==0 or \
                any(res is None for res in channel_results[_TOPOLOGY][phase]['channel_results']):
                if not _CLEAN and os.path.isfile(_RESULTS):
                    print("Loading results from '%s'."%(_RESULTS))
                    channel_results = eval(open(_RESULTS,'r').read())
                if len(channel_results[_TOPOLOGY][phase]['channel_results'])==0 or \
                   any(res is None for res in channel_results[_TOPOLOGY][phase]['channel_results']):
                    mode = 'INITIALISE'
            
            try:

                if mode == 'INITIALISE' or not _RUN_LOCALLY:
                    print(("Now initialising first results for the channels of topology %s and phase '%s'."%(_TOPOLOGY, phase)))
                    if len(channel_results[_TOPOLOGY][phase]['channel_results'])==0:
                        channel_results[_TOPOLOGY][phase]['channel_results'] = [None,]*n_channels
                    hyperparams['Integrator']['integrator'] = _INTEGRATOR
                    hyperparams['Integrator']['eps_rel'] = 1.0e-99
                    hyperparams['Integrator']['n_start'] = _N_START
                    hyperparams['Integrator']['n_max'] = _INCREMENT 
                    hyperparams['Integrator']['n_increase'] = _N_INCREASE
                    hyperparams['Integrator']['reset_vegas_integrator'] = False
                    if not os.path.isdir(pjoin(run_name,'integration_statistics')):
                        os.makedirs(pjoin(run_name,'integration_statistics'))
                    hyperparams['General']['log_file_prefix'] = pjoin(root_path,run_name,'integration_statistics')+'/'
                    
                    for channel_ID in range(n_channels):
                        if channel_results[_TOPOLOGY][phase]['channel_results'][channel_ID] is not None:
                            print("Channel %d of topology %s already in the database."%(channel_ID, _TOPOLOGY))                
                            continue 
                        hyperparams['General']['res_file_prefix'] = pjoin(root_path,run_name,'channel_%d_'%channel_ID)
                        hyperparams['Integrator']['state_filename_prefix'] = pjoin(root_path,run_name,'channel_%d_'%channel_ID)
                        hyperparams['General']['multi_channeling_channel'] = channel_ID
                        hyperparams.export_to(os.path.join(root_path, run_name,'hyperparameters_channel_%d.yaml'%channel_ID))
                        print("Now running channel %d/%d for the first time (%.5gM points required)."%(channel_ID,n_channels,_INCREMENT/1.0e6))
                        result = run_topology(topology,run_name, channel_ID, _WALL_TIME, _RUN_LOCALLY)
                        if result:
                            channel_results[_TOPOLOGY][phase]['channel_results'][channel_ID] = [result,_INCREMENT]
                            print("Result for channel #%d = %.7e +/- %.3e (n_points=%.1fM)"%(channel_ID, result[0], result[1], _INCREMENT//1.0e6))
                    if _RUN_LOCALLY:
                        combined_results = combine_results(channel_results[_TOPOLOGY][phase], verbose=True)
                        mode = 'REFINE'
                    else:
                        print("Now wait for cluster results.")
                        sys.exit()

                if mode == 'REFINE' and _RUN_LOCALLY:
                    hyperparams['Integrator']['integrator'] = _INTEGRATOR
                    hyperparams['Integrator']['eps_rel'] = 1.0e-99
                    hyperparams['Integrator']['n_start'] = _N_START
                    hyperparams['Integrator']['n_max'] = _INCREMENT 
                    hyperparams['Integrator']['n_increase'] = _N_INCREASE
                    if not os.path.isdir(pjoin(run_name,'integration_statistics')):
                        os.makedirs(pjoin(run_name,'integration_statistics'))
                    hyperparams['General']['log_file_prefix'] = pjoin(root_path,run_name,'integration_statistics')+'/'

                    print(("Now refining channels of topology %s and phase '%s'."%(_TOPOLOGY, phase)))
                    hyperparams['Integrator']['reset_vegas_integrator'] = False 
                    combined_results = combine_results(channel_results[_TOPOLOGY][phase], verbose=False)
                    if phase == 'real':
                        analytic_result = topology.analytic_result.real
                    else:
                        analytic_result = topology.analytic_result.imag
                    while (combined_results['error'] / abs(combined_results['central_value'])) > _TARGET_ACCURACY:
                        n_sigmas = abs(combined_results['central_value']-analytic_result)/combined_results['error']
                        print("%s Current result: %.8e +/- %.8e (%.3g%%) vs %.8e with n_tot_points=%d (n_sigmas=%.2g) %s"%(
                            Colour.GREEN if n_sigmas < 3.0 else Colour.RED,
                            combined_results['central_value'], combined_results['error'], 
                            (combined_results['error'] / abs(combined_results['central_value']))*100.0,
                            analytic_result,
                            combined_results['n_tot_points'],
                            n_sigmas,
                            Colour.END
                        ))
                        channel_to_refine = combined_results['max_error'][1]  
                        n_points_required_for_this_run = channel_results[_TOPOLOGY][phase]['channel_results'][channel_to_refine][1]+_INCREMENT
                        print(">>> Now refining channel #%d (with current error %.3g%% out of %.3g%% (%.2g%%)), n_points_required=%.5gM"%(
                            channel_to_refine,
                            (combined_results['max_error'][0] / abs(combined_results['central_value']))*100.0,
                            (combined_results['error'] / abs(combined_results['central_value']))*100.0,
                            (   ((combined_results['max_error'][0] / abs(combined_results['central_value'])) /
                                (combined_results['error'] / abs(combined_results['central_value'])) )*100.0
                            ),
                            n_points_required_for_this_run/1.0e6
                        ))
                        old_result = [
                            channel_results[_TOPOLOGY][phase]['channel_results'][channel_to_refine][0][0],
                            channel_results[_TOPOLOGY][phase]['channel_results'][channel_to_refine][0][1]
                        ]
                        hyperparams['Integrator']['n_max'] = n_points_required_for_this_run
                                
                        hyperparams['General']['res_file_prefix'] = pjoin(root_path,run_name,'channel_%d_'%channel_to_refine)
                        hyperparams['Integrator']['state_filename_prefix'] = pjoin(root_path,run_name,'channel_%d_'%channel_to_refine)
                        hyperparams['General']['multi_channeling_channel'] = channel_to_refine
                        hyperparams.export_to(os.path.join(root_path, run_name,'hyperparameters_channel_%d.yaml'%channel_to_refine))
                        result = run_topology(topology,run_name, channel_to_refine, _WALL_TIME, _RUN_LOCALLY)
                        channel_results[_TOPOLOGY][phase]['channel_results'][channel_to_refine][0] = result
                        channel_results[_TOPOLOGY][phase]['channel_results'][channel_to_refine][1] += _INCREMENT
                        n_points = channel_results[_TOPOLOGY][phase]['channel_results'][channel_to_refine][1]
                        print("OLD result for channel #%d = %.7e +/- %.3e (n_points=%.1fM)"%(
                                        channel_to_refine, old_result[0], old_result[1], (n_points-_INCREMENT)//1.0e6))
                        print("NEW result for channel #%d = %.7e +/- %.3e (n_points=%.1fM)"%(
                                        channel_to_refine, result[0], result[1], n_points//1.0e6))
                        combined_results = combine_results(channel_results[_TOPOLOGY][phase], verbose=False)

                    combined_results = combine_results(channel_results[_TOPOLOGY][phase], verbose=True)
                    print("Final result: %.8e +/- %.8e (%.3f%% required: %.3f%%) vs %.8e"%(
                            combined_results['central_value'], combined_results['error'], 
                            (combined_results['error'] / abs(combined_results['central_value']))*100.0,
                            _TARGET_ACCURACY*100.0,
                            analytic_result
                        ))

            except KeyboardInterrupt:
                save_path = _RESULTS
                print("Saving results to %s"%save_path)
                open(save_path,'w').write(pformat(channel_results))
                sys.exit()
            
            save_path = _RESULTS
            print("Saving results to %s"%save_path)
            open(save_path,'w').write(pformat(channel_results))
        
            print("")
            print("+"*80)
            print("="*80)
            print("+"*80)
            for topology_name, results in list(channel_results.items()):
                for phase, phase_results in list(results.items()):
                    print("")                
                    print(">>> Final results for topology %s and phase %s:"%(topology_name, phase))
                    combined_results = combine_results(phase_results, verbose=True)

