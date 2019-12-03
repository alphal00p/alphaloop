#!/usr/bin/env python2

import sys
import os
import itertools
import subprocess
import time
import glob
import copy
import shutil

pjoin = os.path.join
root_path = os.path.dirname(os.path.realpath( __file__ ))
sys.path.insert(0, root_path)
sys.path.insert(0, pjoin(root_path,'LTD'))

rust_executable_path = os.path.abspath(pjoin(root_path,'rust_backend','target','release','ltd'))
pyvegas_executable_path = os.path.abspath(pjoin(root_path,'pyvegas.py'))
import yaml
from yaml import Loader, Dumper

import ltd_commons
import ltd_utils
from ltd_utils import Colour

from pprint import pprint, pformat

_WORK_DIR = pjoin(root_path, 'scan_dir')
_N_CORES = 8
_SILENCE = False

class HyperParamProber(object):
    
    _DUMMY_VEGAS_STATE_FILE_PATH = 'dummy_vegas_'
    _TOPOLOGIES_FILE = 'LTD/topologies.yaml'
    def __init__(self, 
                 rust_config_file_path = pjoin(_WORK_DIR,'default_hyperparameters.yaml'),
                 topology = 'Hexagon_3s',
                 log_stream = None,
                 additional_run_options = None
                ):
        self._rust_config_file_path = rust_config_file_path
        self._rust_executable = rust_executable_path
        self.topology = topology

        try:
            self._default_run_configuration = yaml.load(open(self._rust_config_file_path,'r'), Loader=Loader)
        except Exception as e:
            print("ERROR: Could not parse yaml configuration file at '%s'."%self._rust_config_file_path)
            raise e
        # Apply additional run options 
        # Save the VEGAS state to a dummy file that will need to be removed after each run
        self._default_run_configuration['state_filename'] = self._DUMMY_VEGAS_STATE_FILE_PATH
        self.log_stream = log_stream
        # Assign some default additional run optionis
        self._default_run_configuration['General']['res_file_prefix'] = _WORK_DIR+'/'
        self._default_run_configuration['Integrator']['state_filename_prefix'] = pjoin(_WORK_DIR,self._DUMMY_VEGAS_STATE_FILE_PATH)
        self._default_run_configuration['General']['log_file_prefix'] = pjoin(_WORK_DIR,'integration_statistics')+'/'
        if not os.path.isdir(pjoin(_WORK_DIR, 'integration_statistics')):
            os.makedirs(pjoin(_WORK_DIR, 'integration_statistics'))
        if not os.path.isfile(pjoin(_WORK_DIR, 'amplitudes.yaml')):
            shutil.copy(pjoin(root_path,'LTD','amplitudes.yaml'),pjoin(_WORK_DIR, 'amplitudes.yaml'))

        # Now save the default yaml configuration to the log stream if specified.
        if log_stream is not None:
            log_stream.write(yaml.dump([('default_configuration',self._default_run_configuration),], Dumper=Dumper))

    def clean(self):
        try:
            for fpath in glob.glob(pjoin(_WORK_DIR,self._DUMMY_VEGAS_STATE_FILE_PATH+'*.dat')):
                os.remove(fpath)
            for fpath in glob.glob(pjoin(_WORK_DIR,'integrration_statistics','*.log')):
                os.remove(fpath)             
        except Exception as e:
            pass

    def run_rust(self):
        """ Steer a run of rust integration with the current yaml config setup and return the result."""
        
        # Remove some temporary files if necessary
        self.clean() 

        cmd = [ self._rust_executable, 
                '-t','%s'%self.topology,
                '-f','%s'%pjoin(_WORK_DIR, 'this_run_hyperparameters.yaml'),
                '-l','%s'%pjoin(_WORK_DIR,'topologies.yaml'),
                '-c','%d'%_N_CORES,
                '-p','%s'%pjoin(_WORK_DIR, 'amplitudes.yaml')
        ]
        #print "Now running command: %s"%(' '.join(cmd))
    	with ltd_utils.Silence(active=_SILENCE):
            subprocess.call(cmd, cwd=_WORK_DIR)            
            # Sleep for a very short while to allow output file flushing
            time.sleep(1.0)

        # Now we should be able to return the result of the run
        try:
            result_path = pjoin(_WORK_DIR,'%s_res.dat'%self.topology)
            str_data = open(result_path,'r').read()
            # Some weird I/O issue of rust
            str_data = '...'.join(str_data.split('...')[:1])
            open(result_path,'w').write(str_data)
            result = yaml.load(open(result_path,'r'), Loader=Loader)
            os.remove(result_path)
            return result
        except Exception as e:
            raise BaseException("ERROR: Could not load rust integration result from file: %s. Error: %s"%(
                result_path, str(e)
            ))

    def probe(self, parameter_values, correlated_variations = False):

        # Compute all combination of hyper parameters to test, those not specified are always considered default
        # We always start by a default run
        all_combinations = [{},]
        if not correlated_variations:
            for param, values in parameter_values.items():
                all_combinations.extend(
                    {param: v} for v in values 
                )
        else:
            all_combinations.extend([ dict(param_values) for param_values in 
                            itertools.product( *[[(param, v) for v in values] for 
                                            param,values in parameter_values.items()]) ])

        #print all_combinations
        #stop

        # Store all results produced below
        all_results = []
        try:
            for param_combination in all_combinations:
                if len(param_combination)==0:
                    print("Running the default setup:")
                else:
                    print("Now considering the following combination of parameters: %s"%
                                        ', '.join('%s=%s'%(p,v) for p,v in param_combination.items()))
                this_run_config = copy.deepcopy(self._default_run_configuration)
                for p, v in param_combination.items():
                    elems = p.split('>')
                    sub_dict = this_run_config
                    try:
                        for elem in elems[:-1]:
                            sub_dict=sub_dict[elem]
                        sub_dict[elems[-1]] = v
                    except:
                        raise BaseException("ERROR: Could not find parameter named '%s' in the yaml config file."%p)
                # Now write this configuration file
                open(pjoin(_WORK_DIR,'this_run_hyperparameters.yaml'),'w').write(yaml.dump(this_run_config, Dumper=Dumper))
                try:
                    run_result = self.run_rust()
                except Exception as e:
                    print('An error occurred when running rust: %s'%str(e))
                    raise e
                # Save the result into the logstream if specified
                if log_stream is not None:
                    log_stream.write(yaml.dump([('run_result', (param_combination,run_result) ), ], Dumper=Dumper))
                print("Result of the run:")
                pprint(run_result)
                all_results.append((param_combination,run_result))

        except KeyboardInterrupt:
            # Remove some temporary files if necessary
            try:
                pass
            ##    os.remove(self._rust_config_file_path+'_this_run')
            except:
                pass
            try:
                os.remove(self._DUMMY_VEGAS_STATE_FILE_PATH)
            except:
                pass

            print("Run interrupted by user, partial results will be returned")

        return {
            'default_configuration' : self._default_run_configuration,
            'run_results'           : all_results
        }

class ResultsAnalyser(object):
    
    def __init__(self, results):
        self.default_configuration = results['default_configuration']
        self.results = results['run_results']
        self.run_time = results['run_time']

    def get_param(self, param):

        elems = param.split('>')
        subdict = self.default_configuration
        for elem in elems[:-1]:
            subdict = subdict[elem]
        return subdict[elems[-1]]

    def get_name(self, name):
        return name.split('>')[-1]
    
    def analyse(self, n_optimal_to_show=5):

        # List all parameters studied
        all_params = sorted(list(set(sum((r[0].keys() for r in self.results),[]))))
        print('='*100)
        print("Statistics considered for topology '%s':"%self.default_configuration['General']['topology'])
        stats_params = ['Integrator>n_start','Integrator>n_increase','Integrator>n_max','Integrator>integrator']
        for param in stats_params:
            print('  %-40s = %s' % (param, self.get_param(param)))
        print('-'*100)
        print("run time: %.3e [h]"%(self.run_time/3600.0))
        print('-'*100)
        print("Default parameter values:")
        for param in all_params:
            print('  %-40s = %.3e'%(self.get_name(param), self.get_param(param)))
        print('-'*100)
        print("Parameters studied: %s"%(', '.join(all_params)))
        print('-'*100)
        default_result = None
        for r in self.results:
            if r[0]=={}:
                default_result = r
                break
        if default_result is None:
            print('WARNING: Could not find default result.')
        else:
            print('%-40s%.3e +- %.3e'%('Default -> ',default_result[1]['result'][0],default_result[1]['error'][0]))
        print('-'*100)            
        for param in all_params:
            print param
            sorted_param_values = sorted(
                ( [r for r in self.results if param in r[0]]+
                  [({param : self.get_param(param) },default_result[1])]),
                key=lambda el: el[1]['error'][0])
            print("First %d optimal values for parameter '%s':"%(min(n_optimal_to_show,len(sorted_param_values)), param))
            for i_optimal, result in enumerate(sorted_param_values[:n_optimal_to_show]):
                print("%-30s%-40s -> %.3e +- %.3e"%(
                    '  | ranked #%d%s: '%(
                        i_optimal+1, ' @DEFAULT' if result[0][param]==self.get_param(param) else ''),
                    ', '.join('%s=%g'%(self.get_name(p),result[0][p]) for p in sorted(result[0].keys()) ),
                    result[1]['result'][0],
                    result[1]['error'][0]
                ))
        print('-'*100)
        sorted_all_param_values = sorted(self.results, key=lambda el: el[1]['error'][0])
        print("Overall %d best results:"%min(n_optimal_to_show,len(sorted_all_param_values)))
        for i_optimal, result in enumerate(sorted_all_param_values[:n_optimal_to_show]):
            print("%-30s%-40s -> %.3e +- %.3e"%(
                '  | ranked #%d%s: '%(
                    i_optimal+1, ' @DEFAULT' if len(result[0]) == 0 else ''),
                ', '.join('%s=%g'%(self.get_name(p),result[0][p]) for p in sorted(result[0].keys()) ),
                result[1]['result'][0],
                result[1]['error'][0]
            ))
        print('='*100)

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

if __name__ == '__main__':

    default_multiplication_factors = [0.0001,0.001,0.01,0.05,0.1,0.25,0.5,0.75,1.5,2.,4.,10.,50.,100.,1000.,10000.]
    default_multiplication_factors = [0.5,1.0]
    parameter_values_to_test = {
        'Deformation>scaling>expansion_threshold' : [-0.15*mult for mult in default_multiplication_factors],
        'Deformation>fixed>M_ij' : [1.0*mult for mult in default_multiplication_factors],
    }

    rust_config_file_path = None
    topology = 'Hexagon_3s'
    n_optimal_to_show = 5
    load_results_from = None
    save_results_to = None
    correlated_scan = False

    integrator = 'cuhre'
    n_max = int(1.0e6)
    n_start = int(1.0e5)
    n_increase = int(1.0e5)

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

        if key in ['rust_config_file_path', 'rcfp']:
            rust_config_file_path = value
        elif key in ['topology','t']:
            topology = value
        elif key in ['correlated']:
            correlated_scan = True
        elif key in ['n_max']:
            n_max = int(eval(value.replace('M','*1.0e6')))
        elif key in ['n_start']:
            n_start = int(eval(value.replace('M','*1.0e6')))
        elif key in ['n_increase']:
            n_increase = int(eval(value.replace('M','*1.0e6')))
        elif key in ['integrator']:
            integrator = value
        elif key in ['quiet']:
            _SILENCE = True
        elif key in ['work_dir','wd']:
            _WORK_DIR = os.path.abspath(value)
        elif key in ['cores','c']:
            _N_CORES = int(value)
        elif key in ['load_results_from','load']:
            load_results_from = value
        elif key in ['save_results_to','save']:
            save_results_to = value
        elif key in ['n_optimal_to_show','show']:
            n_optimal_to_show = eval(value)
        elif key.endswith('_values'):
            parameter_values_to_test[key[:-7]] = eval(value)
        else:
            print "Unrecognized option: %s"%value
            sys.exit(1)
   
    if not save_results_to:
        save_results_to = pjoin(_WORK_DIR,'hyperparameters_scan.yaml')

    if load_results_from is None:

        if not os.path.isfile(pjoin(_WORK_DIR,'topologies.yaml')):
            if os.path.isfile(pjoin(root_path,'LTD','topologies','%s.yaml'%topology)):
                shutil.copy(
                    pjoin(root_path,'LTD','topologies','%s.yaml'%topology),
                    pjoin(_WORK_DIR,'topologies.yaml')
                )
            else:
                shutil.copy(
                    pjoin(root_path,'LTD','topologies.yaml'),
                    pjoin(_WORK_DIR,'topologies.yaml')
                )

        if not os.path.isfile(pjoin(_WORK_DIR, 'topologies.yaml')):
            shutil.copy(pjoin(root_path,'LTD','amplitudes.yaml'),pjoin(_WORK_DIR, 'amplitudes.yaml'))


        if not rust_config_file_path:
            hyperparams = copy.deepcopy(ltd_commons.hyperparameters)
            hyperparams['General']['topology'] = topology            
            if n_max:
                hyperparams['Integrator']['n_max'] = n_max
            if n_start:
                hyperparams['Integrator']['n_start'] = n_start
            if n_increase:
                hyperparams['Integrator']['n_increase'] = n_increase
            if integrator:
                hyperparams['Integrator']['integrator'] = integrator

            hyperparams.export_to(pjoin(_WORK_DIR,'default_hyperparameters.yaml')) 
        else:
            shutil.copy(rust_config_file_path, pjoin(_WORK_DIR,'default_hyperparameters.yaml'))

        log_stream = open(save_results_to,'w')
        prober = HyperParamProber(
            topology=topology,
            log_stream=log_stream
        )
        start_time = time.time()
        test_results = prober.probe(parameter_values_to_test, correlated_variations=correlated_scan)
        test_results['run_time'] = time.time()-start_time
        log_stream.write(yaml.dump([('run_time', test_results['run_time']), ], Dumper=Dumper))
        log_stream.close()
    else:
        try:
            test_results = load_results_from_yaml(load_results_from)
        except Exception as e:
            print("ERROR: Could not hyperparameters scan results from file '%s'. Error: %s"%(
                load_results_from, str(e)
                ))
            raise e

    #pprint(test_results)
    analyser = ResultsAnalyser(test_results)
    analyser.analyse(n_optimal_to_show=n_optimal_to_show)
