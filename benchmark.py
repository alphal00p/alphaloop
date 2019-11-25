#!/usr/bin/env python3
import subprocess
from tabulate import tabulate
from uncertainties import ufloat
import time
import yaml
import json
from tqdm import tqdm
import os
import shutil
import argparse
from math import sqrt
import sys
import glob 
import time
from pprint import pprint, pformat
from datetime import datetime

file_path = os.path.dirname(os.path.realpath( __file__ ))
pjoin = os.path.join

W  = '\033[0m'  # white (normal)
R  = '\033[31m' # red
G  = '\033[32m' # green
O  = '\033[33m' # orange
B  = '\033[34m' # blue
P  = '\033[35m' # purple
BOLD = '\033[1m'# bold

_VERBOSITY = 0
_TABLE_FORMAT = "fancy_grid"
_CONFIG_FILE_PATH = pjoin(file_path,"LTD", "hyperparameters.yaml")

# Specify a prefix so as to avoid collisions with other runs
_PREFIX = str(time.time()).replace('.','_')+'_'

_RUN_LOCALLY = True
_N_CORES = 4
_RUN_DIR = pjoin(file_path,'run_dir')
_WALL_TIME = 24
_ACCOUNT = 'eth5e'

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

    def __call__(self,n_cores=1,**opts):
        return self.get_rust_result(self['topology'], self['phase'], self['samples'], n_cores, self['integrator'],
            ltd_extra_args=[
                '--n_start', str(self['n_start']),
                '--n_increase', str(self['n_increase']),               
            ],
            **opts
        )

    @classmethod
    def get_rust_result(cls,topology, phase, num_samples, cores, integrator, ltd_extra_args=None, collect_only=False):

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
    
        no_analytical_found = False
        if analytical_result[0]==analytical_result[1]==0.:
            print("WARNING: Topology '%s' does not exist or does not specify an analytical result. The benchmark tool is meant to be used for topologies with a target result, so we will set the analytical result equal to the one obtained by RUST."%topology)
            no_analytical_found = True
        
        #cargo_options = ["cargo","run", "--release", "--bin", "ltd"]
        #if _VERBOSITY<=1:
        #    cargo_options.append("--quiet")
        #cargo_options.append("--")
        # Run the executable directly:
        ltd_exe = pjoin(file_path,"rust_backend","target","release","ltd")
        cargo_options = [ltd_exe,]
            
        prefix = pjoin(_RUN_DIR,'_'.join([_PREFIX,topology])+'_')
        ltd_options = [ "-t", topology, 
                       "--state_filename_prefix", prefix, 
                       "--log_file_prefix", prefix, 
                       "--res_file_prefix", prefix ]
        ltd_options.extend(["-s", str(num_samples)])
        ltd_options.extend(["-c", str(cores)])
        ltd_options.extend(["-f", pjoin(_RUN_DIR,'hyperparameters.yaml')])
        ltd_options.extend(["-p", pjoin(_RUN_DIR,'amplitudes.yaml')])        
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
            'cwd' : _RUN_DIR,
        }
        if not collect_only:
            # remove info from old runs
            try:
                os.remove(pjoin(_RUN_DIR,"%s%s_res.dat" % (prefix,topology) ))
            except:
                pass
            try:
                os.remove(pjoin(_RUN_DIR,"%s%s_state.dat" % (prefix,topology) ))
            except:
                pass
            if _RUN_LOCALLY:
                print("Now running topology %s"%topology)
                if _VERBOSITY>1: 
                    print("Running command: %s"%(' '.join(str(opt) for opt in ["cargo",]+cargo_options+["--",]+ltd_options)))
                    subprocess_options['stdout'] = subprocess.PIPE
       
                if _VERBOSITY>1:
                    rust_process = subprocess.Popen(
                        cargo_options+ltd_options, 
                        **subprocess_options
                    )
                    for line in iter(rust_process.stdout.readline, b''):
                        sys.stdout.write(line.decode())
                else:
                    output = subprocess.check_output(
                        cargo_options+ltd_options, 
                        **subprocess_options
                    ).strip().decode()
            else:
                print("Now launch job for topology %s"%topology)
                submission_script = open(pjoin(file_path,'submission_template.run'),'r').read()
                open(pjoin(_RUN_DIR,'submitter.run'),'w').write(submission_script%{
		    'job_name' : '%sjob'%os.path.basename(os.path.normpath(prefix)),
                    'n_hours' : _WALL_TIME,
                    'n_cpus_per_task' : int(cores),
                    'account': _ACCOUNT,
                    'output' : '%s/LTD_runs/logs/%sjob.out'%(os.environ['SCRATCH'], os.path.basename(os.path.normpath(prefix))),
                    'error' : '%s/LTD_runs/logs/%sjob.err'%(os.environ['SCRATCH'], os.path.basename(os.path.normpath(prefix))),
                    'executable_line' : ' '.join(cargo_options+ltd_options)
	        })
                subprocess.call(['sbatch','submitter.run'], cwd=_RUN_DIR)
                return None



        if not os.path.isfile(pjoin(_RUN_DIR,"%s%s_res.dat" % (prefix,topology) )):
            return None

        # read the output file
        with open(pjoin(_RUN_DIR,"%s%s_res.dat" % (prefix,topology) ), 'r') as f:
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

        if no_analytical_found:
            if phase == 'both':
                analytical_result = (rust_result['result'][0], rust_result['result'][1])
            elif phase == 'imag':
                analytical_result = (0.0, rust_result['result'][0])
            elif phase == 'real':
                analytical_result = (rust_result['result'][0], 0.0)

        result = {
            'revision': git_revision,
            'diff': git_diff == '',
            'num_samples': num_samples,
            'topology': topology,
            'result': tuple(integral_result),
            'error': tuple(error),
            'analytical_result': analytical_result,
        }
        return result

class BenchmarkRun1loop(BenchmarkRun):

    def __init__(self,*args,**opts):
        super(BenchmarkRun1loop, self).__init__(*args,**opts)
        self['phase'] = self.get('phase','both')        
        self['samples'] = self.get('samples',100*Units.M)
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

class BenchmarkRunHighloop(BenchmarkRun):

    def __init__(self,*args,**opts):
        super(BenchmarkRunHighloop, self).__init__(*args,**opts)
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

    def get_1loop(self, **opts):
        
        res = []
        for topology in self._ALL_1LOOP_TOPOLOGIES:
            res.append(BenchmarkRun1loop(topology, **opts))

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

    def get_short_Valentin(self):
        res = []

        # 1-loop topologies
        res.extend(self.get_1loop(samples=20*Units.M))

        # 2-loop topologies
        res.append(BenchmarkRun2loop("T2_6P_2L_Weinzierl_A", n_start=1000*Units.K, n_increase=10*Units.K, samples=100*Units.M))
        res.append(BenchmarkRun2loop("T2_6P_2L_Weinzierl_B", n_start=1000*Units.K, n_increase=10*Units.K, samples=100*Units.M))
        res.append(BenchmarkRun2loop("T2_6P_2L_Weinzierl_C", n_start=1000*Units.K, n_increase=10*Units.K, samples=100*Units.M))
        res.append(BenchmarkRun2loop("T2_6P_2L_Weinzierl_D", n_start=1000*Units.K, n_increase=10*Units.K, samples=100*Units.M))
        res.append(BenchmarkRun2loop("T2_6P_2L_Weinzierl_E", n_start=1000*Units.K, n_increase=10*Units.K, samples=100*Units.M))
        res.append(BenchmarkRun2loop("T2_6P_2L_Weinzierl_F", n_start=1000*Units.K, n_increase=10*Units.K, samples=100*Units.M))
        res.append(BenchmarkRun2loop("T3_DoubleBox_Weinzierl",n_start=1000*Units.K, n_increase=10*Units.K, samples=100*Units.M))

        # 3-loop topologies
        res.append(BenchmarkRun3loop("T4_TripleBox_Weinzierl",n_start=1000*Units.K, n_increase=10*Units.K, samples=100*Units.M))

        # 4-loop topologies
        res.append(BenchmarkRunHighloop("T4_Quadruple_Box_Weinzierl",n_start=1000*Units.K, n_increase=10*Units.K, samples=100*Units.M))

        return res

    def get_Valentin(self):
        res = []
        
        # 1-loop topologies
        res.extend(self.get_1loop())

        # 2-loop topologies
        res.append(BenchmarkRun2loop("T2_6P_2L_Weinzierl_A", n_start=1000*Units.K, n_increase=10*Units.K, samples=300*Units.M))
        res.append(BenchmarkRun2loop("T2_6P_2L_Weinzierl_B", n_start=1000*Units.K, n_increase=10*Units.K, samples=300*Units.M))
        res.append(BenchmarkRun2loop("T2_6P_2L_Weinzierl_C", n_start=1000*Units.K, n_increase=10*Units.K, samples=300*Units.M))
        res.append(BenchmarkRun2loop("T2_6P_2L_Weinzierl_D", n_start=1000*Units.K, n_increase=10*Units.K, samples=300*Units.M))        
        res.append(BenchmarkRun2loop("T2_6P_2L_Weinzierl_E", n_start=1000*Units.K, n_increase=10*Units.K, samples=300*Units.M))
        res.append(BenchmarkRun2loop("T2_6P_2L_Weinzierl_F", n_start=1000*Units.K, n_increase=10*Units.K, samples=300*Units.M))
        res.append(BenchmarkRun2loop("T3_DoubleBox_Weinzierl",n_start=1000*Units.K, n_increase=10*Units.K, samples=300*Units.M))
        
        # 3-loop topologies
        res.append(BenchmarkRun3loop("T4_TripleBox_Weinzierl",n_start=1000*Units.K, n_increase=10*Units.K, samples=300*Units.M))
 
        # 4-loop topologies
        res.append(BenchmarkRunHighloop("T4_Quadruple_Box_Weinzierl",n_start=1000*Units.K, n_increase=10*Units.K, samples=300*Units.M))

        return res

    def get_PS1PS2_1loop(self):
        
        res = []        
        PS1PS2_1loop = [

'1L_4P_PS1',
'1L_4P_PS1_massive',
'1L_4P_PS2',
'1L_4P_PS2_massive',

'1L_5P_PS1',
'1L_5P_PS1_massive',
'1L_5P_PS2',
'1L_5P_PS2_massive',

'1L_6P_PS1',
'1L_6P_PS1_massive',
'1L_6P_PS2',
'1L_6P_PS2_massive',

'1L_8P_PS1',
'1L_8P_PS1_massive',
'1L_8P_PS2',
'1L_8P_PS2_massive',

]
        for topo in PS1PS2_1loop:
            #res.append(BenchmarkRun1loop(topo, samples=50*Units.M))
            res.append(BenchmarkRun2loop(topo, n_start=1000*Units.K, n_increase=1000*Units.K, samples=300*Units.M))
        return res

    def get_PS1PS2_2loop(self):

        res = []        
        PS1PS2_2loop = [

'2L_4P_Ladder_PS1',
'2L_4P_Ladder_PS1_massive',
'2L_4P_Ladder_PS2',
'2L_4P_Ladder_PS2_massive',

'2L_5P_Planar_PS1',
'2L_5P_Planar_PS1_massive',
'2L_5P_Planar_PS2',
'2L_5P_Planar_PS2_massive',

'2L_6P_A_PS1',
'2L_6P_A_PS1_massive',
'2L_6P_A_PS2',
'2L_6P_A_PS2_massive',

'2L_6P_B_PS1',
'2L_6P_B_PS1_massive',
'2L_6P_B_PS2',
'2L_6P_B_PS2_massive',

'2L_6P_C_PS1',
'2L_6P_C_PS1_massive',
'2L_6P_C_PS2',
'2L_6P_C_PS2_massive',

'2L_6P_D_PS1',
'2L_6P_D_PS1_massive',
'2L_6P_D_PS2',
'2L_6P_D_PS2_massive',

'2L_6P_E_PS1',
'2L_6P_E_PS1_massive',
'2L_6P_E_PS2',
'2L_6P_E_PS2_massive',

'2L_6P_F_PS1',
'2L_6P_F_PS1_massive',
'2L_6P_F_PS2',
'2L_6P_F_PS2_massive',

#'2L_8P_PS1', # TOO HARD?
#'2L_8P_PS1_massive', # TOO HARD?
#'2L_8P_PS2', # TOO HARD?
#'2L_8P_PS2_massive', # TOO HARD?

]
        for topo in PS1PS2_2loop:
            res.append(BenchmarkRun2loop(topo, n_start=1000*Units.K, n_increase=1000*Units.K, samples=300*Units.M))
        return res

    def get_PS1PS2_3loop(self):

        res = []        
        PS1PS2_3loop = [

'3L_4P_Ladder_PS1',
'3L_4P_Ladder_PS1_massive',
'3L_4P_Ladder_PS2',
# '3L_4P_Ladder_PS2_massive', # FAILURE OF ECOX SOLVER

'3L_5P_Planar_PS1',
'3L_5P_Planar_PS1_massive',
'3L_5P_Planar_PS2',
# '3L_5P_Planar_PS2_massive', # FAILURE OF ECOX SOLVER


]
        for topo in PS1PS2_3loop:
            res.append(BenchmarkRun3loop(topo, n_start=1000*Units.K, n_increase=1000*Units.K, samples=300*Units.M))
        return res

    def get_PS1PS2_4loop(self):

        res = []        
        PS1PS2_4loop = [

'4L_4P_Ladder_PS1',
'4L_4P_Ladder_PS1_massive',
'4L_4P_Ladder_PS2',
'4L_4P_Ladder_PS2_massive',

# 'FISHNET_2x2_PS1', # TOO HARD? 
# 'FISHNET_2x2_PS1_massive', # TOO HARD?
# 'FISHNET_2x2_PS2', # TOO HARD?
# 'FISHNET_2x2_PS2_massive', # TOO HARD?

]
        for topo in PS1PS2_4loop:
            res.append(BenchmarkRunHighloop(topo, n_start=1000*Units.K, n_increase=1000*Units.K, samples=300*Units.M))
        return res


    def get_PS3_1loop(self):
        
        res = []        
        PS1PS2_1loop = [
'1L_4P_PS3',
'1L_4P_PS3_massive',
'1L_5P_PS3',
'1L_5P_PS3_massive',
'1L_6P_PS3',
'1L_6P_PS3_massive',
'1L_8P_PS3',
'1L_8P_PS3_massive',
]
        for topo in PS1PS2_1loop:
            #res.append(BenchmarkRun1loop(topo, samples=50*Units.M))
            res.append(BenchmarkRun2loop(topo, n_start=1000*Units.K, n_increase=1000*Units.K, samples=300*Units.M))
        return res

    def get_PS3_2loop(self):

        res = []        
        PS3_2loop = [ 
'TM1_bot',
'TM1_top',

'2L_4P_Ladder_PS3',
'2L_4P_Ladder_PS3_massive',
'2L_5P_Planar_PS3',
'2L_5P_Planar_PS3_massive',
'2L_6P_A_PS3',
'2L_6P_A_PS3_massive',
'2L_6P_B_PS3',
'2L_6P_B_PS3_massive',
'2L_6P_C_PS3',
'2L_6P_C_PS3_massive',
'2L_6P_D_PS3',
'2L_6P_D_PS3_massive',
'2L_6P_E_PS3',
'2L_6P_E_PS3_massive',
'2L_6P_F_PS3',
'2L_6P_F_PS3_massive',
'2L_8P_PS3',
'2L_8P_PS3_massive',
]
        for topo in PS3_2loop:
            res.append(BenchmarkRun2loop(topo, n_start=1000*Units.K, n_increase=1000*Units.K, samples=300*Units.M))
        return res

    def get_PS3_3loop(self):

        res = []        
        PS3_3loop = [

'3L_4P_Ladder_PS3',
'3L_4P_Ladder_PS3_massive',
'3L_5P_Planar_PS3',
'3L_5P_Planar_PS3_massive',


]
        for topo in PS3_3loop:
            res.append(BenchmarkRun3loop(topo, n_start=1000*Units.K, n_increase=1000*Units.K, samples=300*Units.M))
        return res

    def get_PS3_4loop(self):

        res = []        
        PS3_4loop = [

'4L_4P_Ladder_PS1',
'4L_4P_Ladder_PS1_massive',
'4L_4P_Ladder_PS2',
'4L_4P_Ladder_PS2_massive',

# 'FISHNET_2x2_PS1', # Too hard for cvxpy? 
# 'FISHNET_2x2_PS1_massive', # Too hard for cvxpy?
# 'FISHNET_2x2_PS2', # Too hard for cvxpy?
# 'FISHNET_2x2_PS2_massive', # Too hard for cvxpy?

]
        for topo in PS3_4loop:
            res.append(BenchmarkRunHighloop(topo, n_start=1000*Units.K, n_increase=1000*Units.K, samples=300*Units.M))
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
    t  = time.time()

    # get the hyperpamater file for this run
    hyperparam_resource = pjoin(_RUN_DIR,'hyperparameters.yaml')
    hyperparamaters = {}
    try:
        with open(hyperparam_resource, 'r') as f:
            hyperparamaters = yaml.safe_load(f)
    except:
        pass

    git_diff = subprocess.check_output(["git", "diff"]).strip().decode()

    historical_data = {t: {
            'hyperparameters': hyperparamaters,
            'diff': git_diff,
            'samples': samples,
        }
    }

    try:
        with open(output_path, 'r') as f:
            historical_data.update(json.load(f))
    except:
        pass

    with open(output_path, 'w') as f:
        json.dump(historical_data, f, indent=4)

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
                    "{:6}".format(precision) if precision is not None else 'N/A',
                    (R + '%.2g'%percentage + W if percentage > 1.0 else G + '%.2g'%percentage + W) if percentage is not None else 'N/A',
                    sample['revision'], sample['diff']]
            )
   
    print(tabulate(data, ['Topology', '# Samples', 'Result', 'Reference', 'Accuracy', 'Precision', 'Percentage', 'Tag', 'Clean'], tablefmt=_TABLE_FORMAT))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Tool for benchmarking hyperparameters')
    parser.add_argument('-t', metavar='topologies', type=str, nargs='+', help='topologies to test', required=False)
    parser.add_argument('--from_history', action='store_true', help='Read the topology data from the history')
    parser.add_argument('-s', default='100000', type=int, help='number of samples')
    parser.add_argument('-c', default='4', help='number of cores')
    parser.add_argument('--wall_time', default='24', type=int, help='Set wall time')
    parser.add_argument('-v', default='0', type=int, help='Set verbosity: 0-10')
    parser.add_argument('--table_format', default='fancy_grid', help="Chose the table render format, useful to disable some non-utf8 characters not supported by some terminals. Choose in: simple, plain, grid, fancy_grid, github, pipe, html, and more..." )
    parser.add_argument('-b', default='manual', help="benchmark name, in: 'manual' (select topologies by hand) or, "+
                                                      "'1loop', '2loop', '3loop', 'quick', 'long' which can be combined with a '+' sign.")
    parser.add_argument('--phase', default='both', choices=['real','imag','both'], help='the phase for the integration')
    parser.add_argument('--integrator', default='auto', choices=['auto','cuhre','vegas'], help='Specify the integrator to use') 
    parser.add_argument('--n_start', default='100000', help='n_start for vegas')
    parser.add_argument('--n_increase', default='100000', help='n_increase for vegas')
    parser.add_argument('--history_path', default='default', help='specify a JSON file path to store the result')
    parser.add_argument('--config_path', default=pjoin(file_path, 'LTD', 'hyperparameters.yaml'), help='specify a path to a hyperparameters.yaml file to consider.')
    parser.add_argument('--prefix', default='bm', help='Specify a prefix for the results.')
    parser.add_argument('--gather', action='store_true', help='Gather results.')
    parser.add_argument('--clean', action='store_true', help='Clean existing results.')
    parser.add_argument('--cluster', action='store_true', help='Launch jobs on cluster.')
    parser.add_argument('--save', action='store_true', help='Save results in json dump.')    
    parser.add_argument('--run_dir', default=pjoin(file_path, 'run_dir'), help='Specify the run directory.')    
    parser.add_argument('--account',default='eth5e', help='Cluster account to budget the job to.')
    parser.add_argument('--show_hyperparameters', default='0', choices=[0,1,2,3], type=int, help='level of hyperparameter printing in history mode')
    args = parser.parse_args()

    samples = []

    _VERBOSITY = args.v
    _WALL_TIME = args.wall_time
    _TABLE_FORMAT = args.table_format
    _CONFIG_FILE_PATH = args.config_path
    _RUN_DIR = pjoin(file_path,args.run_dir)
    _ACCOUNT = args.account

    if args.clean:
        print("Cleaning up directory %s."%_RUN_DIR)
        for f in ( glob.glob(pjoin(_RUN_DIR,'*_state.dat'))+
                   glob.glob(pjoin(_RUN_DIR,'*_res.dat'))+
                  glob.glob(pjoin(_RUN_DIR,'*.log')) ):
            os.remove(f)
        sys.exit(1)
        
    _PREFIX = args.prefix
    _RUN_LOCALLY = (not args.cluster)
    if _RUN_LOCALLY:
       _N_CORES = args.c
    else:
       _N_CORES = 36
    if not os.path.isdir(_RUN_DIR):
        print("Run directory %s not found."%_RUN_DIR)
        sys.exit(1)

    if _CONFIG_FILE_PATH != pjoin(_RUN_DIR,'hyperparameters.yaml'):
        shutil.copyfile(_CONFIG_FILE_PATH, pjoin(_RUN_DIR,'hyperparameters.yaml'))
    
    if not os.path.isfile(pjoin(_RUN_DIR,'amplitudes.yaml')):
        shutil.copyfile(pjoin(file_path,'LTD','amplitudes.yaml'), pjoin(_RUN_DIR,'amplitudes.yaml'))
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
        historical_data = get_history(
            pjoin(_RUN_DIR,"historical_benchmarks.json") if args.history_path=='default' else args.history_path)

        for t, v in historical_data.items():
            samples = [s for s in v['samples'] if s['topology'] in [r['topology'] for r in benchmark_runs]]
            if len(samples) > 0:
                print(BOLD + 'Run from {}'.format(datetime.fromtimestamp(float(t)).strftime("%a %d %H:%M:%S")) + W)
                
                if args.show_hyperparameters == 1:
                    # print some basic hyperparameters
                    dh = v['hyperparameters']['Deformation']
                    basic_params = [{
                        'fixed': {
                                'M_ij': dh['fixed']['M_ij'],
                                'include_normal_source': dh['fixed']['include_normal_source'],
                                'normalisation_of_subspace_components': dh['fixed']['normalisation_of_subspace_components'],
                                'normalisation_per_number_of_sources': dh['fixed']['normalisation_per_number_of_sources'],
                                'normalize_per_source': dh['fixed']['normalize_per_source'],
                                },
                        'normalize_on_E_surfaces_m': dh['normalize_on_E_surfaces_m'],
                        'overall_scaling_constant': dh['overall_scaling_constant'],
                        'scaling': {
                            'branch_cut_m':  dh['scaling']['branch_cut_m'],
                            'cut_propagator_check': dh['scaling']['cut_propagator_check'],
                            'expansion_check': dh['scaling']['expansion_check'],
                            'expansion_check_strategy': dh['scaling']['expansion_check_strategy'],
                            'expansion_threshold': dh['scaling']['expansion_threshold'],
                            'lambda': dh['scaling']['lambda'],
                            'non_cut_propagator_check':dh['scaling']['non_cut_propagator_check'],
                            'positive_cut_check': dh['scaling']['positive_cut_check'],
                            'source_branch_cut_m': dh['scaling']['source_branch_cut_m'],
                            'source_branch_cut_multiplier': dh['scaling']['source_branch_cut_multiplier'],
                            'source_branch_cut_threshold': dh['scaling']['source_branch_cut_threshold'],
                        }
                    }]
                    print(yaml.dump(basic_params, indent=4, sort_keys=True))
                elif args.show_hyperparameters == 2:
                    print(yaml.dump([v['hyperparameters']['Deformation']], indent=4, sort_keys=True))
                elif args.show_hyperparameters == 3:
                    print(yaml.dump([v['hyperparameters']], indent=4, sort_keys=True))

                render_data(samples, args.s, sort=True)
    else:
        print("Now running a benchmark involving the following topologies:")
        print(">>> %s"%(', '.join(r['topology'] for r in benchmark_runs)))
        print('')        
        if _RUN_LOCALLY and not args.gather: 
            pbar = tqdm([run['topology'] for run in benchmark_runs])
        else:
            pbar = [run['topology'] for run in benchmark_runs]
        for i_run, pbar_element in enumerate(pbar):
            if _RUN_LOCALLY and not args.gather: pbar.set_description(pbar_element)
            result = benchmark_runs[i_run](n_cores=_N_CORES,collect_only=args.gather)
            if result is None:
                continue
            if _VERBOSITY>0: print("Result for topology '%s':"%benchmark_runs[i_run]['topology'])          
            if _VERBOSITY>0: render_data([result], benchmark_runs[i_run]['samples'])
            samples.append(result)

        if (len(benchmark_runs) >= 1 and len(samples)>1) or (_VERBOSITY==0 and len(samples)==1):
            print("All results for: %s%s"%(
                '-t %s '%(' '.join(args.t)) if args.t else '',
                '-b %s '%args.b if args.b!='manual' else ''
            )) 
            render_data(samples, args.s, sort=True)

    # ask to save data
    if not args.from_history and args.save:
        if args.history_path.lower()=='default':
            #if input("Do you want to save the new run? [y/N]: ") in ['y','Y']:
            print("Saving result to file %s."%pjoin(_RUN_DIR,"historical_benchmarks.json"))                
            save_to_history(samples, pjoin(_RUN_DIR,"historical_benchmarks.json"))
        else:
            print("Saving result to file %s."%pjoin(file_path,args.history_path))            
            save_to_history(samples, pjoin(file_path,args.history_path))
