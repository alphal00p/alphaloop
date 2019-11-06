#!/usr/bin/env python3
import subprocess
from tabulate import tabulate
from uncertainties import ufloat
import time
import yaml
import json
from tqdm import tqdm
import os
import argparse
from math import sqrt
import sys
from pprint import pprint, pformat

file_path = os.path.dirname(os.path.realpath( __file__ ))
pjoin = os.path.join

W  = '\033[0m'  # white (normal)
R  = '\033[31m' # red
G  = '\033[32m' # green
O  = '\033[33m' # orange
B  = '\033[34m' # blue
P  = '\033[35m' # purple

_VERBOSITY = 0

class Units:
    K = 1000
    M = 1000*K
    B = 1000*M

class BenchmarkRun(dict):

    def __init__(self,*args,**opts):
        if opts:
            self.update(opts)
        if len(args)==1:
            self['topology'] = args[0]

    def __call__(self,n_cores=1):
        return self.get_rust_result(self['topology'], self['phase'], self['samples'], n_cores, self['integrator'],
            [
                '--n_start', str(self['n_start']),
                '--n_increase', str(self['n_increase']),               
            ]
        )

    @classmethod
    def get_rust_result(cls,topology, phase, num_samples, cores, integrator, ltd_extra_args=None):

        # Check if there exists a yaml specification file for that topology
        topology_resource = pjoin(file_path,"LTD","topologies.yaml")
        if os.path.isfile(pjoin(file_path,"LTD","topologies","%s.yaml"%topology)):
            topology_resource = pjoin(file_path,"LTD","topologies","%s.yaml"%topology)
    
        # get analytical result
        analytical_result = (0.,0.)
        n_loops = 0
        with open(topology_resource, 'r') as f:
            topologies = yaml.safe_load(f)
            try:
                for t in topologies:
                    if t['name'] == topology:
                        analytical_result = (t['analytical_result_real'], t['analytical_result_imag'])
                        n_loops = t['n_loops'] 
            except:
                print("Could not extract analytic result and/or load information for topology: %s"%topology)
                sys.exit(1)
    
        if analytical_result[0]==analytical_result[1]==0.:
            print("Topology '%s' does not exist or does not specify an analytical result. The benchmark tool is meant to be used for topologies with a target result."%topology)
            sys.exit(1)
    
        # remove info from old runs
        try:
            os.remove(pjoin(file_path,"rust_backend","%s_res.dat" % topology))
        except:
            pass
        try:
            os.remove(pjoin(file_path,"rust_backend","%s_state.dat" % topology))
        except:
            pass
       
        cargo_options = ["run", "--release", "--bin", "ltd"]
        if _VERBOSITY<=1:
            cargo_options.append("--quiet")

        ltd_options = [ "-t", topology]
        ltd_options.extend(["-s", str(num_samples)])
        ltd_options.extend(["-c", str(cores)])
        ltd_options.extend(["-l", str(topology_resource)])
        if integrator.lower()=='auto':
            if n_loops>1:
                ltd_options.extend(["--integrator", 'vegas'])
            else:
                ltd_options.extend(["--integrator", 'cuhre'])
        else:
            ltd_options.extend(["--integrator", integrator])

        if ltd_extra_args:
            ltd_options.extend(ltd_extra_args)
    
        # TODO: check for errors
        subprocess_options = {
            'stderr': subprocess.STDOUT,
            'cwd' : pjoin(file_path,"rust_backend")
        }
        if _VERBOSITY>1: 
            print("Running command: %s"%(' '.join(str(opt) for opt in ["cargo",]+cargo_options+["--",]+ltd_options)))
            subprocess_options['stdout'] = subprocess.PIPE
       
        if _VERBOSITY>1:
            rust_process = subprocess.Popen(
                ["cargo",]+cargo_options+["--",]+ltd_options, 
                **subprocess_options
            )
            for line in iter(rust_process.stdout.readline, b''):
                sys.stdout.write(line.decode())
        else:
            output = subprocess.check_output(
                ["cargo",]+cargo_options+["--",]+ltd_options, 
                **subprocess_options
            ).strip().decode()

        # read the output file
        with open(pjoin(file_path,"rust_backend","%s_res.dat" % topology), 'r') as f:
            rust_result = yaml.safe_load(f)

        # get git revision
        git_revision = subprocess.check_output(["git", "describe", "--always"]).strip().decode()

        # get git diff
        git_diff = subprocess.check_output(["git", "diff"]).strip().decode()

        integral_result = [None,None]
        error = [None,None]
        if phase == 'both':
            integral_result = rust_result['result']
            error = rust_result['error']
        elif phase == 'imag':
            integral_result[1] = rust_result['result'][0]
            error[1] = rust_result['error'][0]
        elif phase == 'real':
            integral_result[0] = rust_result['result'][0]
            error[0] = rust_result['error'][0]

        result = {
            'revision': git_revision,
            'diff': git_diff,
            'num_samples': num_samples,
            'topology': topology,
            'result': tuple(integral_result),
            'error': tuple(error),
            'analytical_result': analytical_result,
            'time': time.time()
        }
        return result

class BenchmarkRun1loop(BenchmarkRun):

    def __init__(self,*args,**opts):
        super(BenchmarkRun1loop, self).__init__(*args,**opts)
        self['phase'] = self.get('phase','both')        
        self['samples'] = self.get('samples',10*Units.M)
        self['integrator'] = self.get('integrator','cuhre')
        self['n_start'] = self.get('n_start',100*Units.K)
        self['n_increase'] = self.get('n_increase',100*Units.K)

class BenchmarkRun2loop(BenchmarkRun):

    def __init__(self,*args,**opts):
        super(BenchmarkRun2loop, self).__init__(*args,**opts)
        self['phase'] = self.get('phase','both')        
        self['samples'] = self.get('samples',100*Units.M)
        self['integrator'] = self.get('integrator','vegas')
        self['n_start'] = self.get('n_start', Units.M)
        self['n_increase'] = self.get('n_increase',100*Units.K)

class BenchmarkRun3loop(BenchmarkRun):

    def __init__(self,*args,**opts):
        super(BenchmarkRun3loop, self).__init__(*args,**opts)
        self['phase'] = self.get('phase','both')        
        self['samples'] = self.get('samples',Units.B)
        self['integrator'] = self.get('integrator','vegas')
        self['n_start'] = self.get('n_start', Units.M)
        self['n_increase'] = self.get('n_increase',100*Units.K)

class BenchmarkRunNloop(BenchmarkRun):

    def __init__(self,*args,**opts):
        super(BenchmarkRunNloop, self).__init__(*args,**opts)
        self['phase'] = self.get('phase','both')        
        self['samples'] = self.get('samples',10*Units.B)
        self['integrator'] = self.get('integrator','vegas')
        self['n_start'] = self.get('n_start', Units.M)
        self['n_increase'] = self.get('n_increase',Units.M)

class Benchmark(list):

    _ALL_1LOOP_TOPOLOGIES = [
        
        # Dario customised ones
        "Pentagon_1s",
        "Pentagon_2s",
        "Pentagon_3s",
        "Hexagon_1s",
        "Hexagon_2s",
        "Hexagon_3s",
        "Hexagon_4s",

        # Zeno customised ones
        "Pentagon_10E_1s",
        "Pentagon_6E_4s",
        "Pentagon_8E_5s",
        "Hexagon_6E_4s",
        "Hexagon_10E_4s",
        "Hexagon_9E_4s",
        "Hexagon_10E_7s",
        "Hexagon_10E_5s",
        "Hexagon_6E_2s",

    ]

    _ALL_2LOOP_TOPOLOGIES = [
        # Our PRL topology
        "T3_DoubleBox_Weinzierl",

        # Weinzierl's topologies
        "T2_6P_2L_Weinzierl_A",
        "T2_6P_2L_Weinzierl_B",
        "T2_6P_2L_Weinzierl_C",
        "T2_6P_2L_Weinzierl_D",
        "T2_6P_2L_Weinzierl_E",
        "T2_6P_2L_Weinzierl_F",

    ]

    def __init__(self, name, manual_specifications=None):
        self.name = name
        if manual_specifications is not None:
            self.extend(manual_specifications)
        
        # Now fill in the above
        for benchmark_name in self.name.split('+'):
            if benchmark_name == 'manual':
                continue
            if not hasattr(self, "get_%s"%benchmark_name):
                print("Benchmark named '%s' not recognized."%benchmark_name)
                sys.exit(1)
            additional_runs = eval("self.get_%s()"%benchmark_name)
            self.extend([r for r in additional_runs if r['topology'] not in [_['topology'] for _ in self]])

    #TODO Below please all contribute to create nice benchmarks that balance speed and sensitivity!

    def get_1loop(self):
        
        res = []
        for topology in self._ALL_1LOOP_TOPOLOGIES:
            res.append(BenchmarkRun1loop(topology))

        return res

    def get_quick(self):
        
        res = []
        res.append(BenchmarkRun1loop("Box_3E"))
        res.append(BenchmarkRun1loop("Box_4E"))

        return res

    def get_Zeno(self):
        res = []
        res.append(BenchmarkRun1loop("Box_3E"))
        res.append(BenchmarkRun1loop("Box_4E"))
        return res

    def get_Dario(self):
        res = []
        res.append(BenchmarkRun1loop("Box_3E"))
        res.append(BenchmarkRun1loop("Box_4E"))
        return res

    def get_Ben(self):
        res = []
        res.append(BenchmarkRun1loop("Box_3E"))
        res.append(BenchmarkRun1loop("Box_4E"))
        return res

    def get_Valentin(self):
        res = []
        
        # 1-loop topologies
        res.extend(self.get_1loop())

        # 2-loop topologies
        res.append(BenchmarkRun2loop("T2_6P_2L_Weinzierl_A", n_start=100*Units.K, n_increase=10*Units.K, samples=100*Units.M))
        res.append(BenchmarkRun2loop("T2_6P_2L_Weinzierl_B", n_start=100*Units.K, n_increase=10*Units.K, samples=100*Units.M))
        res.append(BenchmarkRun2loop("T2_6P_2L_Weinzierl_C", n_start=100*Units.K, n_increase=10*Units.K, samples=100*Units.M))
        res.append(BenchmarkRun2loop("T2_6P_2L_Weinzierl_D", n_start=100*Units.K, n_increase=10*Units.K, samples=100*Units.M))        
        res.append(BenchmarkRun2loop("T2_6P_2L_Weinzierl_E", n_start=100*Units.K, n_increase=10*Units.K, samples=100*Units.M))
        res.append(BenchmarkRun2loop("T2_6P_2L_Weinzierl_F", n_start=100*Units.K, n_increase=10*Units.K, samples=100*Units.M))
        res.append(BenchmarkRun2loop("T3_DoubleBox_Weinzierl",n_start=100*Units.K, n_increase=10*Units.K, samples=100*Units.M))
        
        # 3-loop topologies
        res.append(BenchmarkRun3loop("T4_TripleBox_Weinzierl",n_start=100*Units.K, n_increase=10*Units.K, samples=100*Units.M))
 
        # 4-loop topologies
        res.append(BenchmarkRun3loop("T4_Quadruple_Box_Weinzierl",n_start=100*Units.K, n_increase=10*Units.K, samples=100*Units.M))

        return res


def get_history(history_path):
    historical_data = []

    try:
        with open(history_path, 'r') as f:
            historical_data = json.load(f)
    except:
        pass

    return historical_data

def save_to_history(samples, output_path):
    historical_data = []
    try:
        with open(output_path, 'r') as f:
            historical_data = json.load(f)
    except:
        pass

    historical_data.extend(samples)

    with open(output_path, 'w') as f:
        json.dump(historical_data, f)

def get_score_for_sample(sample, number_of_samples):
    scores = {'real': None, 'imag': None}
    for i_phase, phase in enumerate(['real', 'imag']):
        if sample['result'][i_phase] is not None:
            scores[phase] = {'accuracy': None, 'precision': None, 'percentage': None}
            scores[phase]['accuracy'] = abs(sample['analytical_result'][i_phase] - sample['result'][i_phase]) / sample['error'][i_phase]
            if abs(sample['analytical_result'][i_phase]) != 0.:
                scores[phase]['precision'] = sample['error'][i_phase] * sqrt(number_of_samples) / abs(sample['analytical_result'][i_phase])
                scores[phase]['percentage'] = 100.0*(abs(sample['analytical_result'][i_phase]-sample['result'][i_phase]) / abs(sample['analytical_result'][i_phase]))        
            else:
                scores[phase]['precision'] = None
                scores[phase]['percentage'] = None

    return scores

def render_data(samples, number_of_samples, sort=False):
    """ Render the data in a table. """
    data = []

    def sort_by_accuracy(sample):
        score = get_score_for_sample(sample, number_of_samples)
        if score['imag'] is None:
            return score['real']['accuracy']
        elif score['real'] is None:
            return score['imag']['accuracy']
        else:
            return min(score['real']['accuracy'],score['imag']['accuracy'])

    if sort:
        samples = sorted(samples, key=sort_by_accuracy)

    for sample in samples:
        score = get_score_for_sample(sample, number_of_samples)

        for i_phase, (phase,phase_name) in enumerate( [('real', 'Real'), ('imag', 'Imag')]):
            if sample['result'][i_phase] is None:
                continue

            (accuracy, precision, percentage) = score[phase]['accuracy'], score[phase]['precision'], score[phase]['percentage']
            
            data.append(
                [sample['topology'] + ' ' + phase_name, "{:,}".format(int(sample['num_samples'])),
                    ufloat(sample['result'][i_phase], sample['error'][i_phase]), 
                    sample['analytical_result'][i_phase],
                    R + str(accuracy) + W if accuracy > 2.0 else G + str(accuracy) + W,
                    precision if precision is not None else 'N/A',
                    (R + '%.2g'%percentage + W if percentage > 1.0 else G + '%.2g'%percentage + W) if percentage is not None else 'N/A',
                    sample['revision'], sample['diff'] == '']
            )
    
    print(tabulate(data, ['Topology', '# Samples', 'Result', 'Reference', 'Accuracy', 'Precision', 'Percentage', 'Tag', 'Clean'], tablefmt="fancy_grid"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Tool for benchmarking hyperparameters')
    parser.add_argument('-t', metavar='topologies', type=str, nargs='+', help='topologies to test', required=False)
    parser.add_argument('--from_history', action='store_true', help='Read the topology data from the history')
    parser.add_argument('-s', default='100000', type=int, help='number of samples')
    parser.add_argument('-c', default='4', help='number of cores')
    parser.add_argument('-v', default='0', type=int, help='Set verbosity: 0-10')
    parser.add_argument('-b', default='manual', help="benchmark name, in: 'manual' (select topologies by hand) or, "+
                                                      "'1loop', '2loop', '3loop', 'quick', 'long' which can be combined with a '+' sign.")
    parser.add_argument('--phase', default='both', choices=['real','imag','both'], help='the phase for the integration')
    parser.add_argument('--integrator', default='auto', choices=['auto','cuhre','vegas'], help='Specify the integrator to use') 
    parser.add_argument('--n_start', default='100000', help='n_start for vegas')
    parser.add_argument('--n_increase', default='100000', help='n_increase for vegas')
    parser.add_argument('--history_path', default='default', help='specify a JSON file path to store the result')
    args = parser.parse_args()

    samples = []

    _VERBOSITY = args.v    
    
    # A list of runs to be considered, each identified by a dictionary with the following entries:
    # {
    #   'topology'   : <topology_name>,
    #   'phase'      : <phase_to_consider>,
    #   'samples'    : <sample_size_to_consider>,
    #   'integrator' : <integrator_to_consider>,
    #   'n_start'    : <integrator_to_consider>, 
    #   'n_increase' : <integrator_to_consider>, 
    # }
    # This object is stored in an instance of a specialized list class "Benchmark" 
    defaults = {
        'phase'      : args.phase,
        'samples'    : args.s,
        'integrator' : args.integrator,
        'n_start'    : args.n_start,
        'n_increase' : args.n_increase
    }
    benchmark_runs = Benchmark(args.b, manual_specifications = [] if not args.t else [
        BenchmarkRun(t, **defaults) for t in args.t
    ])

    if args.from_history:
        historical_data = get_history(pjoin(file_path,"historical_benchmarks.json") if args.history_path=='default' else args.history_path)
        samples = [d for d in historical_data if d['topology'] in [r['topology'] for r in benchmark_runs]]
    else:
        pbar = tqdm([run['topology'] for run in benchmark_runs])
        print('')
        print("Now running a benchmark involving the following topologies:")
        print(">>> %s"%(', '.join(r['topology'] for r in benchmark_runs)))
        print('')
        for i_run, pbar_element in enumerate(pbar):
            pbar.set_description(pbar_element)
            result = benchmark_runs[i_run](n_cores=args.c)
            if _VERBOSITY>0: print("Result for topology '%s':"%benchmark_runs[i_run]['topology'])          
            if _VERBOSITY>0: render_data([result], benchmark_runs[i_run]['samples'])
            samples.append(result)

    if args.from_history or len(benchmark_runs) > 1:
        print("All results for: %s%s"%(
            '-t %s '%(' '.join(args.t)) if args.t else '',
            '-b %s '%args.b if args.b!='manual' else ''
        )) 
        render_data(samples, args.s, sort=True)

    # ask to save data
    if not args.from_history:
        if args.history_path.lower()=='default':
            if input("Do you want to save the new run? [y/N]: ") in ['y','Y']:
                save_to_history(samples, pjoin(file_path,"historical_benchmarks.json"))
        else:
            save_to_history(samples, pjoin(file_path,args.history_path))
