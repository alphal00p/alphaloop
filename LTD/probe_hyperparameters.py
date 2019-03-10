#!/usr/bin/env python

import sys
import os
import itertools
import subprocess
import time
import yaml
from yaml import Loader, Dumper
from pprint import pprint
import ltd_commons
import copy
import multiprocessing

pjoin = os.path.join 

root_path = os.path.dirname(os.path.realpath( __file__ ))
rust_back_end_path = os.path.abspath(pjoin(root_path,os.path.pardir,'rust_backend'))
class HyperParamProber(object):
    
    _DUMMY_VEGAS_STATE_FILE_PATH = pjoin(root_path,'dummy_vegas_state_file.dat')
    _CURRENT_RUST_CONFIG_PATH = pjoin(root_path,'hyperparameter_scan_current_config.yaml')
    def __init__(self, 
                 rust_config_file_path = pjoin(root_path,'hyperparameters.yaml'),
                 rust_executable_path = pjoin(rust_back_end_path,'target','release','ltd'),
                 topology = 'box',
                 log_stream = None,
                 additional_run_options = None,
                 silence_rust = False
                ):
        self._rust_config_file_path = rust_config_file_path
        self._rust_executable = rust_executable_path
        self.topology = topology
        self.silence_rust = silence_rust

        try:
            self._default_run_configuration = ltd_commons.HyperParameters.import_from(self._rust_config_file_path)
        except Exception as e:
            print("ERROR: Could not parse yaml configuration file at '%s'."%self._rust_config_file_path)
            raise e

        # Apply additional run options 
        # Save the VEGAS state to a dummy file that will need to be removed after each run
        # self._default_run_configuration['Integrator.save_state_file'] = self._DUMMY_VEGAS_STATE_FILE_PATH
        self.log_stream = log_stream
        
        for option, value in additional_run_options.items():
            if not self.set_option_in_configuration(self._default_run_configuration,option,value):
                raise BaseException("Could not find parameter '%s' in the the config yaml file."%option)
        # And the topology as well
        self.set_option_in_configuration(self._default_run_configuration,'General.topology',self.topology)
        
        # Now save the default yaml configuration to the log stream if specified.
        if log_stream is not None:
            log_stream.write(yaml.dump([('default_configuration',self._default_run_configuration),], Dumper=Dumper))

    def set_option_in_configuration(self, config, option, value):
        """ Sets the specified option on the config file. Option names are organised in categories like this:
               CategoryA.SubCategoryC.option_name
        """

        option_path = option.split('.')
        category_dict = config
        for name in option_path[:-1]:
            try:
                category_dict = category_dict[name]
            except KeyError:
                return False
        try:
            category_dict[option_path[-1]] = value
        except KeyError:
            return False
        
        return True

    def run_rust(self):
        """ Steer a run of rust integration with the current yaml config setup and return the result."""

        # Remove some temporary files if necessary
        try:
            os.remove(self._DUMMY_VEGAS_STATE_FILE_PATH)
        except:
            pass 
        # Cores are currently not read from the config file, so specify them here for now.
        cmd = ' '.join([
            self._rust_executable,'--topology %s'%self.topology,
            '--config %s'%self._CURRENT_RUST_CONFIG_PATH,
            '--cores %d'%multiprocessing.cpu_count()
        ])
        if self.silence_rust:
            subprocess.call(cmd,shell=True,stdout=open(os.devnull, 'w'))
        else:
            subprocess.call(cmd,shell=True)

        # Remove some temporary files if necessary
        try:
            os.remove(self._CURRENT_RUST_CONFIG_PATH)
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

                this_run_config = copy.deepcopy(self._default_run_configuration)

                for p, v in param_combination.items():
                    if not self.set_option_in_configuration(this_run_config, p, v):
                        raise BaseException("ERROR: Could not find parameter named '%s' in the yaml config file."%p)
                
                # Now write this configuration file
                this_run_config.export_to(self._CURRENT_RUST_CONFIG_PATH) 
                
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
                os.remove(self._CURRENT_RUST_CONFIG_PATH)
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

    def get_entry_from_config(self, config, entry):
        """ Fetches an entry from a config file, using the notation
            CategoryA.SubCategoryC.option_name
        """

        entry_path = entry.split('.')
        category_dict = config
        for name in entry_path[:-1]:
            try:
                category_dict = category_dict[name]
            except KeyError:
                raise BaseException('Cannot find category named "%s" in config file.'%name) 
        try:
            return category_dict[entry_path[-1]]
        except KeyError:
            raise BaseException('Cannot find entry named "%s" in config file.'%entry_path[-1]) 

    def analyse(self, n_optimal_to_show=5):

        # List all parameters studied
        all_params = sorted(list(set(sum((r[0].keys() for r in self.results),[]))))
        print('='*80)
        print("Statistics considered for topology '%s':"%self.default_configuration['General']['topology'])
        stats_params = ['n_start','n_increase','n_max']
        for param in stats_params:
            print('  %-20s = %d' % (param, self.default_configuration['Integrator'][param]))
        print('-'*80)
        print("run time: %.3e [h]"%(self.run_time/3600.0))
        print('-'*80)
        print("Default parameter values:")
        for param in all_params:
            print param
            print('  %-20s = %.3e'%(param, self.get_entry_from_config(self.default_configuration,param)))
        print('-'*80)
        print("Parameters studied: %s"%(', '.join(all_params)))
        print('-'*80)
        default_result = None
        for r in self.results:
            if r[0]=={}:
                default_result = r
                break
        if default_result is None:
            print('WARNING: Could not find default result.')
        else:
            print('%-50s%.3e +- %.3e'%('Default -> ',default_result[1]['result'][0],default_result[1]['error'][0]))
        print('-'*40)            
        for param in all_params:
            sorted_param_values = sorted(
                ( [r for r in self.results if param in r[0]]+
                  [({param : self.get_entry_from_config(self.default_configuration,param) },default_result[1])]),
                key=lambda el: el[1]['error'][0])
            print("First %d optimal values for parameter '%s':"%(min(n_optimal_to_show,len(sorted_param_values)), param))
            for i_optimal, result in enumerate(sorted_param_values[:n_optimal_to_show]):
                print("%-30s%-40s -> %.3e +- %.3e"%(
                    '  | ranked #%d%s: '%(
                        i_optimal+1, ' @DEFAULT' if result[0][param]==self.get_entry_from_config(self.default_configuration,param) else ''),
                    ', '.join('%s=%g'%(p,result[0][p]) for p in sorted(result[0].keys()) ),
                    result[1]['result'][0],
                    result[1]['error'][0]
                ))
        print('-'*80)
        sorted_all_param_values = sorted(self.results, key=lambda el: el[1]['error'][0])
        print("Overall %d best results:"%min(n_optimal_to_show,len(sorted_all_param_values)))
        for i_optimal, result in enumerate(sorted_all_param_values[:n_optimal_to_show]):
            print("%-30s%-40s -> %.3e +- %.3e"%(
                '  | ranked #%d%s: '%(
                    i_optimal+1, ' @DEFAULT' if len(result[0]) == 0 else ''),
                ', '.join('%s=%g'%(p,result[0][p]) for p in sorted(result[0].keys()) ),
                result[1]['result'][0],
                result[1]['error'][0]
            ))
        print('='*80)


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
    prog_options = ' '.join(sys.argv[1:]).split('--')

    default_multiplication_factors = [0.0001,0.001,0.01,0.05,0.1,0.25,0.5,0.75,1.5,2.,4.,10.,50.,100.,1000.,10000.]
    parameter_values_to_test = {
        'Deformation.generic.M_ij' : [10., 100.],
        'Deformation.rodrigo.a_ij' : [2.*mult for mult in default_multiplication_factors],
        #TODO Add Deformation.rodrigo.a_ij here and implement the corresponding quantity in rust. 
        # < 0. would means auto-scaled and > 0. means set fixed.
    }
    additional_run_options = {
                     'General.deformation_strategy' : 'generic',
                     'General.topology'             : 'Box',
                     'Integrator.n_increase'        : 0,
                     'Integrator.n_max'             : 1000000,
                     'Integrator.n_start'           : 10000,
                     'Integrator.integrated_phase'  : 'real',
                     # Parameters below not yet available as hyperparameters
                     # 'Integrator.cores'             : 0,
                     # 'Integrator.seed'              : 1,
                     # 'Integrator.refine_n_runs'     : 0,
                     # 'Integrator.save_state_file'   : 'vegas_state.dat',
                 }
    rust_config_file_path = pjoin(root_path,'hyperparameters.yaml')
    topology = 'Box'
    n_optimal_to_show = 5
    load_results_from = None
    save_results_to = pjoin(os.getcwd(),'hyperparameters_scan.yaml')
    silence_rust = False
    correlated_variations = False
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
        elif option_name in ['correlated_variations', 'cv']:
            if option_value is None:
                correlated_variations = True
            else:
                correlated_variations = bool(eval(option_value))
        elif option_name in ['silence_rust','sr']:
            if option_value is None:
                silence_rust = True
            else:
                silence_rust = bool(eval(option_value))
        elif option_name.endswith('_values'):
            parameter_values_to_test[option_name[:-7]] = eval(option_value)
        else:
            # Unknown options will be considered by default as additional value to give to rust
            additional_run_options[option_name] = eval(option_value)

    if load_results_from is None:
        log_stream = open(save_results_to,'w')
        prober = HyperParamProber(
            rust_config_file_path=rust_config_file_path,
            topology=topology,
            additional_run_options=additional_run_options,
            log_stream=log_stream,
            silence_rust=silence_rust
        )
        start_time = time.time()
        test_results = prober.probe(parameter_values_to_test, correlated_variations=correlated_variations)
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
