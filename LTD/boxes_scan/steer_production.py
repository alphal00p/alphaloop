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

general_hyperparams = copy.deepcopy(ltd_commons.hyperparameters)
general_hyperparams['General']['absolute_precision'] = 1.0e+5
general_hyperparams['Integrator']['integrator'] = 'vegas'
general_hyperparams['Integrator']['n_start'] = int(1e7)
general_hyperparams['Integrator']['n_max'] = int(1e10)
general_hyperparams['Integrator']['n_increase'] = int(1e6)
general_hyperparams['Parameterization']['mode'] = 'spherical'
general_hyperparams['Parameterization']['mapping'] = 'log'
general_hyperparams['Parameterization']['b'] = 1.0
general_hyperparams['General']['integration_statistics'] = False
general_hyperparams['General']['log_file_prefix'] = pjoin(root_path,'integration_statistics')
general_hyperparams['General']['screen_log_core'] = 0

general_hyperparams.export_to(os.path.join(root_path, 'box','hyperparameters.yaml'))
 
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

def run_topology(topo,dir_name, index, local=True):
    """ Run topology of specified index and directory locally or on a SLURM scheduled cluster."""
    print("Now running %s topology #%d"%(dir_name, index))
    cmd = [ rust_executable_path, 
            '-t','scan_%d'%index, 
            '-f','%s'%pjoin(root_path,dir_name,'hyperparameters.yaml'), 
            '-l','%s'%pjoin(root_path,dir_name,'topologies.yaml'),
            '-c','%d'%_N_CORES
    ]
    
    #print(' '.join(cmd))
    with ltd_utils.Silence(active=_SILENCE):
        subprocess.call(cmd, cwd=pjoin(root_path,dir_name))
        # Sleep for a very short while to allow output file flushing
        time.sleep(0.3)

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
    
    x_values = []
    for line in open(pjoin(root_path,dir_name,'analytic_result.dat'),'r'):
        x_values.append(float(line.split(' ')[0]))

    n_scans = len(all_topologies)
    for index in range(1, n_scans+1):
        result_path = pjoin(root_path,dir_name,'scan_%d_res.dat'%index)
        if not os.path.isfile(result_path):
            print("Error: Run did not successfully complete as the results yaml dump '%s' could not be found."%result_path )
            data_stream.close()
            return
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
        else:
            data_stream.write('\n%s'%line)
    data_stream.close()

    if _CLEAN:
        for index in range(1, n_scans+1):
            if os.path.isfile(result_path):
                os.remove(result_path)

if __name__ == '__main__':


    if any('quiet' in arg for arg in sys.argv[1:]):
        _SILENCE = True

    if any('clean' in arg for arg in sys.argv[1:]):
        _CLEAN = True

    if any('gather' in arg for arg in sys.argv[1:]):

        if any('box' in arg for arg in sys.argv[1:]):     
            # Gather box results
            topologies = ltd_utils.TopologyCollection.import_from(os.path.join(root_path, 'box','topologies.yaml')) 
            gather_result(topologies, 'box', _CLEAN)

    if any('run' in arg for arg in sys.argv[1:]):

        if any('box' in arg for arg in sys.argv[1:]):
            # Run box
            
            # First refresh configuration files
            box_hyperparams = copy.deepcopy(general_hyperparams)
            box_hyperparams['General']['absolute_precision'] = 1.0e+5
            box_hyperparams['Integrator']['integrator'] = 'cuhre'
            box_hyperparams['Integrator']['n_max'] = int(1e6)
            box_hyperparams.export_to(os.path.join(root_path, 'box','hyperparameters.yaml'))

            # Get topologies
            topologies = ltd_utils.TopologyCollection.import_from(os.path.join(root_path, 'box','topologies.yaml'))
            n_topologies = len(topologies)
            for index in range(1, n_topologies+1):
                run_topology(topologies['scan_%d'%index],'box', index, _RUN_LOCALLY)

