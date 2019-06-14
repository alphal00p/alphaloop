#!/usr/bin/env python
import os
import sys
import copy
import subprocess
import time

pjoin = os.path.join
root_path = os.path.dirname(os.path.realpath( __file__ ))
sys.path.insert(0, pjoin(root_path,os.pardir))
sys.path.insert(0, pjoin(root_path,os.pardir,os.pardir))

rust_executable_path = os.path.abspath(pjoin(root_path,os.pardir,os.pardir,'rust_backend','target','release','ltd'))
pyvegas_executable_path = os.path.abspath(pjoin(root_path,os.pardir,os.pardir,'pyvegas.py'))
import yaml
from yaml import Loader, Dumper

import ltd_commons
import ltd_utils
from ltd_utils import Colour

_PYVEGAS = False
_RUN_LOCALLY = True
_CLEAN = False
_SILENCE = False
_N_CORES = 8
_CLEAN = False 
_PREFIX = ''
_IDS = None

_WALL_TIMES = {'box': 1, 'doublebox':10, 'triplebox':24}

pyvegas_hyperparms = {'survey_iter': 5, 'survey_neval': int(1e5), 'refine_iter': 10, 'refine_neval': int(1e6), 'phase': 'imag'}

general_hyperparams = copy.deepcopy(ltd_commons.hyperparameters)
general_hyperparams['General']['absolute_precision'] = 1.0e+5
general_hyperparams['Integrator']['integrator'] = 'vegas'
general_hyperparams['Integrator']['n_start'] = int(1e5)
general_hyperparams['Integrator']['n_max'] = int(1e6)
general_hyperparams['Integrator']['n_increase'] = int(1e5)
general_hyperparams['Integrator']['seed'] = 0
general_hyperparams['Integrator']['state_filename'] = None
general_hyperparams['Parameterization']['mode'] = 'spherical'
general_hyperparams['Parameterization']['mapping'] = 'log'
general_hyperparams['Parameterization']['b'] = 0.1
general_hyperparams['General']['integration_statistics'] = False
general_hyperparams['General']['statistics_interval'] = 100000000
general_hyperparams['General']['log_file_prefix'] = pjoin(root_path,'integration_statistics')+'/'
general_hyperparams['General']['screen_log_core'] = 1
general_hyperparams['General']['numerical_instability_check'] = True 
general_hyperparams['General']['unstable_point_warning_percentage'] = 200.
general_hyperparams['General']['num_digits_different_for_inconsistency'] = 200.
general_hyperparams['General']['minimal_precision_for_returning_result'] = 2.
general_hyperparams['Integrator']['eps_rel'] = 1.0e-3
general_hyperparams['Integrator']['eps_abs'] = 0.
general_hyperparams['Integrator']['border'] = 1.0e-3
general_hyperparams['Integrator']['maxpass'] = 5
general_hyperparams['Integrator']['maxchisq'] = 10.
general_hyperparams['Integrator']['mindeviation'] = 0.25


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

def run_topology(topo,dir_name, index, n_hours, local=True):
    """ Run topology of specified index and directory locally or on a SLURM scheduled cluster."""
    if _RUN_LOCALLY:
       print "Now running %s topology #%d"%(dir_name, index)
    else:
       print "Now launch job for %s topology #%d"%(dir_name, index)

    print(pjoin(root_path, dir_name, topo.name + '_res.dat'))

    if _PYVEGAS:
        cmd = [ pyvegas_executable_path,
                '%d'%topo.n_loops,
                '-t','scan_%d'%index,
                '-hf','%s'%pjoin(root_path,dir_name,'hyperparameters.yaml'),
                '-f','%s'%pjoin(root_path,dir_name,'topologies.yaml'),
                '-c','%d'%_N_CORES,
                '--si','%d'%pyvegas_hyperparms['survey_iter'],
                '--sn','%d'%pyvegas_hyperparms['survey_neval'],
                '--ri','%d'%pyvegas_hyperparms['refine_neval'],
                '--rn','%d'%pyvegas_hyperparms['refine_neval'],
                '--phase',pyvegas_hyperparms['phase'],
                '--out',pjoin(root_path, dir_name, topo.name + '_res.dat')
        ]
    else:
        cmd = [ rust_executable_path, 
                '-t','scan_%d'%index,
                '-f','%s'%pjoin(root_path,dir_name,'hyperparameters.yaml'),
                '-l','%s'%pjoin(root_path,dir_name,'topologies.yaml'),
                '-c','%d'%_N_CORES
        ]

    if _RUN_LOCALLY:
    	with ltd_utils.Silence(active=_SILENCE):
            subprocess.call(cmd, cwd=pjoin(root_path,dir_name))
            # Sleep for a very short while to allow output file flushing
            time.sleep(0.3)
    else:
        submission_script = open(pjoin(root_path,'submission_template.run'),'r').read()
        open(pjoin(root_path,'submitter.run'),'w').write(submission_script%{
		'job_name' : '%s_scan_%d'%(dir_name, index),
                'n_hours' : n_hours,
                'n_cpus_per_task' :_N_CORES,
                'output' : '%s/LTD_runs/logs/%s_scan_%d.out'%(os.environ['SCRATCH'], dir_name, index),
                'error' : '%s/LTD_runs/logs/%s_scan_%d.err'%(os.environ['SCRATCH'], dir_name, index),
                'executable_line' : ' '.join(cmd)
	})
        subprocess.call(['sbatch','submitter.run'], cwd=root_path)
        return    

    result_path = pjoin(root_path,dir_name,'scan_%d_res.dat'%index)
    if not os.path.isfile(result_path):
        print("Error: Run did not successfully complete as the results yaml dump '%s' could not be found."%result_path )
    else:
        result = yaml.load(open(result_path,'r'), Loader=Loader)
        analytic_result = topo.analytic_result.real if \
                    abs(topo.analytic_result.real)>abs(topo.analytic_result.imag) else topo.analytic_result.imag
        print ">> Analytic result : %.16e"%analytic_result
        print ">> LTD result      : %.16e"%(result['result'][0])
        print ">> LTD error       : %.16e"%(result['error'][0])
        n_sigmas = abs((analytic_result-result['result'][0])/result['error'][0])
        print ">> LTD discrepancy : %s%.2f sigmas%s"%(Colour.GREEN if n_sigmas <=3. else Colour.RED,  n_sigmas, Colour.END)
        print ">> LTD rel. discr. : %.2g%%"%(100.0*abs((analytic_result-result['result'][0])/analytic_result))
        print ">> LTD n_points    : %dM"%(int(result['neval']/1.e6))        
        return result

def gather_result(all_topologies, dir_name, clean=False):
    """ Combine all results in scan_# into a file ready for gnuplotting."""
    
    data_stream = open(pjoin(root_path,dir_name,'ltd_results.dat'),'w')
    data_stream_no_line = open(pjoin(root_path,dir_name,'ltd_results_no_line.dat'),'w')
 
    x_values = []
    for line in open(pjoin(root_path,dir_name,'analytic_result.dat'),'r'):
        x_values.append(float(line.split(' ')[0]))

    n_scans = len(all_topologies)
    for index in range(1, n_scans+1):
        if _IDS is not None and index not in _IDS: continue
        result_path = pjoin(root_path,dir_name,'scan_%d_res.dat'%index)
        if not os.path.isfile(result_path):
            print("Error: Run did not successfully complete as the results yaml dump '%s' could not be found."%result_path )
            continue
        
        result = yaml.load(open(result_path,'r'), Loader=Loader)
        topo = all_topologies['scan_%d'%index]
        analytic_result = topo.analytic_result.real if \
                    abs(topo.analytic_result.real)>abs(topo.analytic_result.imag) else topo.analytic_result.imag

        if not _SILENCE:
            print("Results for %s topology #%d"%(dir_name, index))            
            print ">> Analytic result : %.16e"%analytic_result
            print ">> LTD result      : %.16e"%(result['result'][0])
            print ">> LTD error       : %.16e"%(result['error'][0])
            n_sigmas = abs((analytic_result-result['result'][0])/result['error'][0])
            print ">> LTD discrepancy : %s%.2f sigmas%s"%(Colour.GREEN if n_sigmas <=3. else Colour.RED,  n_sigmas, Colour.END)
            print ">> LTD rel. discr. : %.2g%%"%(100.0*abs((analytic_result-result['result'][0])/analytic_result))
            print ">> LTD n_points    : %dM"%(int(result['neval']/1.e6))
        
        line = '%.16e %.16e %.16e %d %.16e'%(x_values[index-1], result['result'][0], result['error'][0], result['neval'], analytic_result)
        if index==1:
            data_stream.write(line)
            data_stream_no_line.write(line+'\n')
        else:
            data_stream.write('\n%s'%line)
            data_stream_no_line.write('\n%s'%line+'\n')
    data_stream.close()
    data_stream_no_line.close()

    if _CLEAN:
        for index in range(1, n_scans+1):
            if os.path.isfile(result_path):
                os.remove(result_path)

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

        if key=='prefix':
            _PREFIX = value
        elif key=='ids':
            _IDS = eval(value)
        elif key=='cores':
            _N_CORES = int(value)
            n_cores_user_set = True
        else:
            print('Error: options %s not reckognised.'%arg)
            sys.exit(1)


    if any('quiet' in arg for arg in processed_args):
        _SILENCE = True
        
    if any('cluster' in arg for arg in processed_args):
	_RUN_LOCALLY = False
        if not n_cores_user_set:
            _N_CORES = 36

    if any('clean' in arg for arg in processed_args):
        _CLEAN = True

    if any('gather' in arg for arg in processed_args):

        if any('singlebox' in arg for arg in processed_args):     
            # Gather box results
            topologies = ltd_utils.TopologyCollection.import_from(os.path.join(root_path, '%sbox'%_PREFIX,'topologies.yaml')) 
            gather_result(topologies, '%sbox'%_PREFIX, _CLEAN)

        if any('doublebox' in arg for arg in processed_args):     
            # Gather box results
            topologies = ltd_utils.TopologyCollection.import_from(os.path.join(root_path, '%sdoublebox'%_PREFIX,'topologies.yaml')) 
            gather_result(topologies, '%sdoublebox'%_PREFIX, _CLEAN)

        if any('triplebox' in arg for arg in processed_args):     
            # Gather box results
            topologies = ltd_utils.TopologyCollection.import_from(os.path.join(root_path, '%striplebox'%_PREFIX,'topologies.yaml')) 
            gather_result(topologies, '%striplebox'%_PREFIX, _CLEAN)

    if any('run' in arg for arg in processed_args):

        if any('singlebox' in arg for arg in processed_args):
            # Run box
            
            # First refresh configuration files
            box_hyperparams = copy.deepcopy(general_hyperparams)

#            box_hyperparams['General']['relative_precision'] = -10.
#            box_hyperparams['General']['absolute_precision'] = 1.0e+5
##            box_hyperparams['General']['relative_precision'] = 99.
##            box_hyperparams['General']['absolute_precision'] = 1.0e-99
##            box_hyperparams['Integrator']['integrator'] = 'cuhre'
##            box_hyperparams['Integrator']['n_max'] = int(1e7)

            box_hyperparams['General']['relative_precision'] = 3.
            box_hyperparams['General']['absolute_precision'] = 1.0e-4
            box_hyperparams['Integrator']['integrator'] = 'divonne'
            box_hyperparams['Integrator']['n_start'] = int(1e6)
            box_hyperparams['Integrator']['n_increase'] = int(1e6)
            box_hyperparams['Integrator']['n_max'] = int(1e8)

            box_hyperparams['Integrator']['integrated_phase'] = 'imag'
            box_hyperparams['General']['res_file_prefix'] = pjoin(root_path,'%sbox'%_PREFIX)+'/'            
            box_hyperparams.export_to(os.path.join(root_path, '%sbox'%_PREFIX,'hyperparameters.yaml'))

            # Get topologies
            topologies = ltd_utils.TopologyCollection.import_from(os.path.join(root_path, '%sbox'%_PREFIX,'topologies.yaml'))
            n_topologies = len(topologies)
            for index in range(1, n_topologies+1):
                if _IDS is not None and index not in _IDS: continue
                run_topology(topologies['scan_%d'%index],'%sbox'%_PREFIX, index, _WALL_TIMES['box'], _RUN_LOCALLY)


        if any('doublebox' in arg for arg in processed_args):
            # Run box
            
            # First refresh configuration files
            doublebox_hyperparams = copy.deepcopy(general_hyperparams)
            doublebox_hyperparams['General']['relative_precision'] = 3.        
            doublebox_hyperparams['General']['absolute_precision'] = 1.0e-6
#            doublebox_hyperparams['General']['relative_precision'] = 99.
#            doublebox_hyperparams['General']['absolute_precision'] = 1.0e-99
            doublebox_hyperparams['Integrator']['integrator'] = 'divonne'
            doublebox_hyperparams['Integrator']['n_start'] = int(1e6)
            doublebox_hyperparams['Integrator']['n_increase'] = int(1e6)            
            doublebox_hyperparams['Integrator']['n_max'] = int(1e8)
            doublebox_hyperparams['Integrator']['seed'] = 1
            doublebox_hyperparams['Integrator']['integrated_phase'] = 'real'
            doublebox_hyperparams['Integrator']['eps_rel'] = 1.0e-3
            doublebox_hyperparams['General']['res_file_prefix'] = pjoin(root_path,'%sdoublebox'%_PREFIX)+'/'

#            doublebox_hyperparams['General']['minimal_precision_for_returning_result'] = 2.
#            doublebox_hyperparams['Integrator']['state_filename'] = '/users/hirschva/MG5/git_madnklo/PLUGIN/pynloop/LTD/boxes_scan/experiment_doublebox/experiment_state.dat'

            doublebox_hyperparams.export_to(os.path.join(root_path, '%sdoublebox'%_PREFIX,'hyperparameters.yaml'))

            # Get topologies
            topologies = ltd_utils.TopologyCollection.import_from(os.path.join(root_path, '%sdoublebox'%_PREFIX,'topologies.yaml'))
            n_topologies = len(topologies)
            for index in range(1, n_topologies+1):
                if _IDS is not None and index not in _IDS: continue
                run_topology(topologies['scan_%d'%index],'%sdoublebox'%_PREFIX, index, _WALL_TIMES['doublebox'], _RUN_LOCALLY)

        if any('triplebox' in arg for arg in processed_args):
            # Run box
            
            # First refresh configuration files
            triplebox_hyperparams = copy.deepcopy(general_hyperparams)
            triplebox_hyperparams['General']['relative_precision'] = 3.      
            triplebox_hyperparams['General']['absolute_precision'] = 1.0e-8
#            triplebox_hyperparams['General']['relative_precision'] = 99.
#            triplebox_hyperparams['General']['absolute_precision'] = 1.0e-99
            triplebox_hyperparams['Integrator']['integrator'] = 'divonne'
            triplebox_hyperparams['Integrator']['n_start'] = int(1e7)
            triplebox_hyperparams['Integrator']['n_increase'] = int(1e6)            
            triplebox_hyperparams['Integrator']['n_max'] = int(1e8)
            triplebox_hyperparams['Integrator']['seed'] = 1
            triplebox_hyperparams['Integrator']['integrated_phase'] = 'imag'
            triplebox_hyperparams['Integrator']['eps_rel'] = 1.0e-3

            triplebox_hyperparams['General']['minimal_precision_for_returning_result'] = 2.

            triplebox_hyperparams['General']['res_file_prefix'] = pjoin(root_path,'%striplebox'%_PREFIX)+'/' 
            triplebox_hyperparams.export_to(os.path.join(root_path, '%striplebox'%_PREFIX,'hyperparameters.yaml'))

            # Get topologies
            topologies = ltd_utils.TopologyCollection.import_from(os.path.join(root_path, '%striplebox'%_PREFIX,'topologies.yaml'))
            n_topologies = len(topologies)
            for index in range(1, n_topologies+1):
                if _IDS is not None and index not in _IDS: continue
                run_topology(topologies['scan_%d'%index],'%striplebox'%_PREFIX, index, _WALL_TIMES['triplebox'],  _RUN_LOCALLY)
