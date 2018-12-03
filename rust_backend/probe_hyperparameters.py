#!/usr/bin/env python

import sys
import os
import shutil
import subprocess
import yaml
from yaml import Loader, Dumper
from pprint import pprint

pjoin = os.path.join 

class HyperParamProber(object):
    
    _DUMMY_VEGAS_STATE_FILE_PATH = 'dummy_vegas_state_file.dat'
    def __init__(self, 
                 rust_config_file_path = pjoin(os.getcwd(),'config.yaml'),
                 rust_executable_path = pjoin(os.getcwd(),'target','release','integrator'),
                 topology = 'box',
                 additional_run_options = {
                     'min_eval'   : 0,
                     'nstart'     : 1000000,
                     'nincrease'  : 5000,
                     'max_eval'   : 50000000,
                     'cores'      : 0,
                     'seed'       : 1,
                     'refine_n_runs' : 0
                 }
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
        for option, value in additional_run_options.items():
            if option not in self._default_run_configuration:
                raise BaseException("Could not find parameter '%s' in the the config yaml file."%option)
            self._default_run_configuration[option] = value

    def run_rust(self):
        """ Steer a run of rust integration with the current yaml config setup and return the result."""

        # Remove some temporary files if necessary
        try:
            os.remove(self._DUMMY_VEGAS_STATE_FILE_PATH)
        except:
            pass 
        subprocess.call(' '.join([
            self._rust_executable,'-t=%s'%self.topology,'--config=%s_this_run'%self._rust_config_file_path,
        ]),shell=True)
        # Remove some temporary files if necessary
        try:
            os.remove(self._rust_config_file_path+'_this_run')
        except:
            pass
        try:
            os.remove(self._DUMMY_VEGAS_STATE_FILE_PATH)
        except:
            pass
        # Now we should be able to return the result of the run
        try:
            result = yaml.load(open(pjoin(os.getcwd(),'%s_res.dat'%self.topology),'r'), Loader=Loader)
            os.remove(pjoin(os.getcwd(),'%s_res.dat'%self.topology))
            return result
        except Exception as e:
            raise BaseException("ERROR: Could not load rust integration result from file: %s. Error: %s"%(
                '%s_res.dat'%self.topology,
                str(e)
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
            all_combinations = [ dict(param_values) for param_values in 
                            itertools.product( [(param, v) for v in values] for 
                                            param,values in parameter_values.items()) ]

        # Store all results produced below
        all_results = []
        try:
            for param_combination in all_combinations:
                if len(param_combination)==0:
                    print("Running the default setup:")
                else:
                    print("Now considering the following combination of parameters: %s"%
                                        ', '.join('%s=%s'%(p,v) for p,v in param_combination.items()))
                this_run_config = dict(self._default_run_configuration)
                for p, v in param_combination.items():
                    if p not in this_run_config:
                        raise BaseException("ERROR: Could not find parameter named '%s' in the yaml config file."%p)
                    this_run_config[p] = v
                # Now write this configuration file
                open(self._rust_config_file_path+'_this_run','w').write(yaml.dump(this_run_config, Dumper=Dumper))
                try:
                    run_result = self.run_rust()
                except Exception as e:
                    print('An error occurred when running rust: %s'%str(e))
                    raise e
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
                os.remove(_DUMMY_VEGAS_STATE_FILE_PATH)
            except:
                pass

            print("Run interrupted by user, partial results will be returned")

        return all_results

class ResultsAnalyser(object):
    
    def __init__(self, results):
        self.results = results

    def analyse(self, n_optimal_to_show=5):

        # List all parameters studied
        all_params = sorted(list(set(sum((r[0].keys() for r in self.results),[]))))
        print('='*40)
        print("Parameters studied: %s"%(', '.join(all_params)))
        print('-'*40)
        default_result = None
        for r in self.results:
            if r[0]=={}:
                default_result = r
                break
        if default_result is None:
            print('WARNING: Could not find default result.')
        else:
            print('Default MC error: %.3e'%default_result[1]['error'][0])
        print('-'*40)            
        for param in all_params:
            sorted_param_values = sorted([r for r in self.results if param in r[0]], key=lambda el: el[1]['error'][0]) 
            print("First %d optimal values for parameter '%s':"%(n_optimal_to_show, param))
            for i_optimal, result in enumerate(sorted_param_values[:n_optimal_to_show]):
                print("  | optimal values #%d: %s -> MC_error: %.3e"%(
                    i_optimal, 
                    ', '.join('%s=%g'%(p,result[0][p]) for p in sorted(result[0].keys()) ),
                    result[1]['error'][0]
                ))
        print('-'*40)
        print("Overall %d best results:"%n_optimal_to_show)
        sorted_all_param_values = sorted(self.results, key=lambda el: el[1]['error'][0]) 
        for i_optimal, result in enumerate(sorted_all_param_values[:n_optimal_to_show]):
            print("  | optimal values #%d: %s -> MC_error: %.3e"%(
                i_optimal, 
                ', '.join('%s=%g'%(p,result[0][p]) for p in sorted(result[0].keys()) ),
                result[1]['error'][0]
            ))
        print('='*40)


if __name__ == '__main__':
    prog_options = ' '.join(sys.argv[1:]).split('--')

    default_multiplication_factors = [0.0001,0.001,0.01,0.05,0.1,0.25,0.5,0.75,1.5,2.,4.,10.,50.,100.,1000.,10000.]
    parameter_values_to_test = {
        'm1_fac' : [0.035*mult for mult in default_multiplication_factors],
        'm2_fac' : [0.7*mult for mult in default_multiplication_factors],
        'm3_fac' : [0.035*mult for mult in default_multiplication_factors],
        'gamma1' : [0.7*mult for mult in default_multiplication_factors],
    }
    additional_run_options = {
                     'min_eval'   : 0,
                     'nstart'     : 1000000,
                     'nincrease'  : 5000,
                     'max_eval'   : 50000000,
                     'cores'      : 0,
                     'seed'       : 1,
                     'refine_n_runs' : 0
                 }
    rust_config_file_path = pjoin(os.getcwd(),'config.yaml')
    topology = 'box'
    n_optimal_to_show = 5
    load_results_from = None
    save_results_to = pjoin(os.getcwd(),'hyperparameters_scan.dat')
    for prog_option in prog_options:
        if prog_option.strip()=='':
            continue
        try:
            option_name, option_value = prog_option.strip().split('=')
        except ValueError:
            option_name, option_value = prog_option, None

        if option_name in ['rust_config_file_path', 'rcfp']:
            rust_config_file_path = option_value
        elif option_name in ['topology','t']:
            topology = option_value
        elif option_name in ['load_results_from','load']:
            load_results_from = option_value
        elif option_name in ['save_results_to','save']:
            save_results_to = option_value
        elif option_name in ['n_optimal_to_show','show']:
            n_optimal_to_show = eval(option_value)
        elif option_name.endswith('_values'):
            parameter_values_to_test[option_name[:-7]] = eval(option_value)
        else:
            # Unknown options will be considered by default as additional value to give to rust
            additional_run_options[option_name] = eval(option_value)
    
    prober = HyperParamProber(
        rust_config_file_path = rust_config_file_path,
        topology = topology,
        additional_run_options = additional_run_options
    )
    if load_results_from is None:
        test_results = prober.probe(parameter_values_to_test, correlated_variations=False)
        # Dump results to a file for possible later use
        open(save_results_to,'w').write(yaml.dump(test_results, Dumper=Dumper))
    else:
        try:
            test_results = yaml.load(open(load_results_from,'r'), Loader=Loader)
        except Exception as e:
            print("ERROR: Could not hyperparameters scan results from file '%s'. Error: %s"%(
                load_results_from, str(e)
                ))
            raise e

    #pprint(test_results)
    analyser = ResultsAnalyser(test_results)
    analyser.analyse(n_optimal_to_show=n_optimal_to_show)
