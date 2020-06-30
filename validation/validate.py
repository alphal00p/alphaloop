#!/usr/bin/python
"""
Run Validation for different processes
"""


import getopt
import glob
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml

pjoin = os.path.join

# Define the number of cores to use for each integration
CORES = 4

# This are personal path that have to be fixed
AL_PATH = "/home/andrea/BitBucket/alphaloop/"
MG_PATH = "/home/andrea/Programs/MG5_aMC_v2_7_2_py3"
VALIDATION_PATH = pjoin(AL_PATH, 'validation')

# Create Validation folder
Path(VALIDATION_PATH).mkdir(parents=True, exist_ok=True)

# Rust executable
rALPHA = pjoin(AL_PATH, 'rust_backend', 'target', 'debug', 'ltd')


def set_hyperparameters(process_name, sg_name, multi_channeling=False, workspace=None, multi_settings={}):
    '''Create the correct hyperparameter file'''
    path = pjoin(AL_PATH, 'LTD', 'hyperparameters.yaml')
    hyperparameters = yaml.load(open(path, 'r'), Loader=yaml.Loader)

    try:
        min_jpt = multi_settings['Selectors']['jet']['min_jpt']
    except KeyError:
        min_jpt = hyperparameters['Selectors']['jet']['min_jpt']

    # Set some default values for the hyperparameter
    # Set values for General
    hyperparameters['General']['multi_channeling'] = multi_channeling
    hyperparameters['General']['multi_channeling_including_massive_propagators'] = multi_channeling
    hyperparameters['General']['deformation_strategy'] = 'fixed'
    hyperparameters['General']['topology'] = sg_name
    hyperparameters['General']['res_file_prefix'] = "%s_" % process_name
    hyperparameters['General']['log_file_prefix'] =\
        "stats/%s_%s_%d_" % (process_name, sg_name, min_jpt)
    hyperparameters['General']['partial_fractioning_multiloop'] = True
    hyperparameters['General']['partial_fractioning_threshold'] = -1
    # Set values for Integrator
    hyperparameters['Integrator']['dashboard'] = False
    hyperparameters['Integrator']['integrator'] = 'vegas'
    # Set values for CrossSeciton
    print(hyperparameters['Selectors'].keys())
    hyperparameters['CrossSection']['numerator_source'] = 'FORM'
    hyperparameters['CrossSection']['picobarn'] = True

    # General settings can be parsed though the multi_settings valiable
    # It must have the same structure as hyperparameters
    # It's not necessary that it contains all the elements
    for name, settings in multi_settings.items():
        for setting, value in settings.items():
            if name not in hyperparameters or setting not in hyperparameters[name]:
                print("Unknown setting hyperparameters[{}][{}] does not exits!".format(
                    name, setting))
                sys.exit(1)
            if isinstance(hyperparameters[name][setting], dict):
                hyperparameters[name][setting].update(value)
            else:
                hyperparameters[name][setting] = value

    # Create the hyperpamater file
    output_path = pjoin(workspace, 'hyperparameters',
                        '%s_%s_%d.yaml' % (process_name, sg_name, min_jpt))
    with open(output_path, "w") as stream:
        yaml.dump(hyperparameters, stream)  # default_flow_style=None)
    return output_path


def run_super_graph(process_name, sg_name, aL_path_output, suffix='', multi_settings={}, workspace=None):
    ''' run supergraph for specific process '''
    # Create hyperparameter file
    hyperparameters_path = set_hyperparameters(
        process_name,
        sg_name,
        workspace=workspace,
        multi_settings=multi_settings)
    print(hyperparameters_path)
    # run integration
    rust_executable = [rALPHA]
    rust_executable += ['-f', '%s' % hyperparameters_path]
    rust_executable += ['--cross_section',
                        pjoin(aL_path_output, 'Rust_inputs', '%s.yaml' % sg_name)]
    rust_executable += ['-c', '%d' % CORES]

    # Store output
    #err_path = pjoin(WORKSPACE,'%s.err'%sg_name)
    #out_path = pjoin(WORKSPACE,'%s.out'%sg_name)

    r = subprocess.run(rust_executable,
                       cwd=workspace,
                       env=os.environ.update(
                           {'MG_NUMERATOR_PATH': '%s/' % aL_path_output}),
                       capture_output=True)
    if r.returncode != 0:
        print(r.stdout.decode('UTF-8'))
        print(r.stderr.decode('UTF-8'))


def get_result(filename):
    with open(filename, 'r') as stream:
        try:
            ss = ""
            for i in range(10):
                line = stream.readline()
                if '...' in line:
                    break
                ss += line
            # print(ss)
            result = yaml.safe_load(ss)
            return [result['neval'],
                    result['result'][0], result['error'][0],
                    result['result'][1], result['error'][1]]
        except yaml.YAMLError as exc:
            print(exc)


if __name__ == "__main__":
    # For the moment run for LO
    suffix = 'LO'

    # Extra settings
    hyper_settings = {'General': {},
                      'CrossSection': {},
                      'Deformation': {},
                      'Integrator': {},
                      'Observables': {},
                      'Selectors': {},
                      'Parameterization': {},
                      }

    # Integrator Settings
    hyper_settings['Integrator']['n_max'] = int(1e9)
    hyper_settings['Integrator']['n_new'] = int(1e5)
    hyper_settings['Integrator']['n_start'] = int(1e5)
    hyper_settings['Integrator']['n_increase'] = int(1e5)
    # Selector Settings
    hyper_settings['Selectors']['active_selectors'] = ['jet']
    hyper_settings['Selectors']['jet'] = {'min_jets': 1, 'min_jpt': 100}

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hp:', [
                                   'run', 'collect', 'process='])
    except getopt.GetoptError:
        print("\n>>> Try something like: %s --process=epem_a_ttx --run\n" %
              sys.argv[0])
        raise
    opts = {o: v for o, v in opts}
    # Parse arguments
    if opts.get('--process') is None:
        raise getopt.GetoptError(
            "\nMissing process definition: %s --process=epem_a_ttx run\n" % sys.argv[0])
    if opts.get('--run') is None and opts.get('--collect') is None:
        raise getopt.GetoptError("\nNeed to specify --run or --collect\n")
    process_name = opts['--process']

    # define alpha loop input folder
    aL_process_output = pjoin(MG_PATH, 'TEST_QGRAF_%s' % process_name)

    # load instructions
    instructions = yaml.load(
        open(pjoin(aL_process_output, 'Rust_inputs',
                   'all_QG_supergraphs.yaml'), 'r'),
        Loader=yaml.Loader)

    # Define workspace for the computation of all the diagrams
    WORKSPACE = pjoin(VALIDATION_PATH, '%s_%s' % (process_name, suffix))

    #############
    #   RUN
    #############
    if opts.get('--run') is not None:
        # Create workspace
        Path(WORKSPACE).mkdir(parents=True, exist_ok=True)
        Path(pjoin(WORKSPACE, 'stats')).mkdir(parents=True, exist_ok=True)
        Path(pjoin(WORKSPACE, 'hyperparameters')).mkdir(
            parents=True, exist_ok=True)
        sg_list = [topo['name'] for topo in instructions['topologies']]
        for sg_name in sg_list:
            run_super_graph(process_name,
                            sg_name,
                            aL_process_output, suffix='_%s' % suffix,
                            multi_settings=hyper_settings,
                            workspace=WORKSPACE)
    #############
    #  COLLECT
    #############
    if opts.get('--collect') is not None:
        sg_list = [topo['name'] for topo in instructions['topologies']]
        output_file = pjoin(WORKSPACE, "results.csv")
        results = []
        print(output_file)
        for sg in instructions['topologies']:
            filename = glob.glob("%s/*%s*.dat" % (WORKSPACE, sg['name']))[0]
            print("Extracting %s" % sg['name'], end="\n")

            data = [sg['name'], sg['multiplicity']]
            data.extend(get_result(filename))
            results += [data]
        print("")
        df = pd.DataFrame(results, columns=[
            "id", "multiplicity", "neval", "real", "real_err", "imag", "imag_err"])
        df.to_csv(output_file, sep=',', encoding='utf-8', index=False)
