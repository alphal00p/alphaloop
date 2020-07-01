#!/usr/bin/python3.8
"""
Run Validation for different processes
"""


import glob
import os
import subprocess
import sys
from argparse import ArgumentParser
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


def set_hyperparameters(process_name, sg_name, multi_channeling=False, workspace=None, no_jets=False, multi_settings={}):
    '''Create the correct hyperparameter file'''
    path = pjoin(AL_PATH, 'LTD', 'hyperparameters.yaml')
    hyperparameters = yaml.load(open(path, 'r'), Loader=yaml.Loader)

    try:
        min_jpt = multi_settings['Selectors']['jet']['min_jpt']
    except KeyError:
        min_jpt = hyperparameters['Selectors']['jet']['min_jpt']

    if no_jets:
        min_jpt = 'nocut'
    # Set some default values for the hyperparameter
    # Set values for General
    hyperparameters['General']['multi_channeling'] = multi_channeling
    hyperparameters['General']['multi_channeling_including_massive_propagators'] = multi_channeling
    hyperparameters['General']['deformation_strategy'] = 'fixed'
    hyperparameters['General']['topology'] = "%s_%s" % (sg_name, min_jpt)
    hyperparameters['General']['res_file_prefix'] = "%s_" % process_name
    hyperparameters['General']['log_file_prefix'] =\
        "stats/%s_%s_%s_" % (process_name, sg_name, min_jpt)
    hyperparameters['General']['partial_fractioning_multiloop'] = True
    hyperparameters['General']['partial_fractioning_threshold'] = -1
    # Set values for Integrator
    hyperparameters['Integrator']['dashboard'] = False
    hyperparameters['Integrator']['integrator'] = 'vegas'
    hyperparameters['Integrator']['internal_parallelization'] = True
    hyperparameters['Integrator']['n_vec'] = int(1e8)
    # Set values for CrossSection
    hyperparameters['CrossSection']['numerator_source'] = 'FORM'
    hyperparameters['CrossSection']['picobarns'] = True

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
                        '%s_%s_%s.yaml' % (process_name, sg_name, min_jpt))
    with open(output_path, "w") as stream:
        yaml.dump(hyperparameters, stream)  # default_flow_style=None)
    return output_path


def run_super_graph(process_name, sg_name, aL_path_output, suffix='', multi_settings={}, no_jets=False, workspace=None):
    ''' run supergraph for specific process '''
    # Create hyperparameter file
    hyperparameters_path = set_hyperparameters(
        process_name,
        sg_name,
        workspace=workspace,
        no_jets=no_jets,
        multi_settings=multi_settings)
    print(hyperparameters_path)
    # run integration
    rust_executable = [rALPHA]
    rust_executable += ['-f', '%s' % hyperparameters_path]
    rust_executable += ['--cross_section',
                        pjoin(aL_path_output, 'Rust_inputs', '%s.yaml' % sg_name)]
    rust_executable += ['-c', '%d' % CORES]

    # Store output
    err_path = pjoin(WORKSPACE, '%s.err' % sg_name)
    out_path = pjoin(WORKSPACE, '%s.out' % sg_name)

    r = subprocess.run(rust_executable,
                       cwd=workspace,
                       stdout=open(out_path, 'w'),
                       stderr=open(err_path, 'w'),
                       env=os.environ.update(
                           {'MG_NUMERATOR_PATH': '%s/' % aL_path_output}),
                       capture_output=False)
    if r.returncode != 0:
        print("\033[1;31mFAIL: see {} and {} for further details.\033[0m".format(
            err_path, out_path))
    #    print(r.stdout.decode('UTF-8'))
    #    print(r.stderr.decode('UTF-8'))


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

    # Extra settings
    hyper_settings = {'General': {},
                      'CrossSection': {},
                      'Deformation': {},
                      'Integrator': {},
                      'Observables': {},
                      'Selectors': {},
                      'Parameterization': {},
                      }

    # Parse arguments
    parser = ArgumentParser()
    parser.add_argument("--process", dest="process", help="process name")
    parser.add_argument("-g", "--generate", action="store_true", dest="generate", default=False,
                        help="Generate alphaLoop output")
    parser.add_argument("-r", "--run", action="store_true", dest="run", default=False,
                        help="run the integration for all SG")
    parser.add_argument("-c", "--collect", action="store_true", dest="collect", default=False,
                        help="Collect the reults from individual SQ in one CSV file")
    parser.add_argument("-v", "--validate", action="store_true", dest="validate", default=False,
                        help="Compare with bench results")
    parser.add_argument("--no_jets", action="store_true", dest="no_jets", default=False,
                        help="No jets in the final state")
    parser.add_argument("--min_jpt",  dest="min_jpt", default=100, type=int,
                        help="Jet cutoff")
    parser.add_argument("--n_max", dest="n_max", default=int(1e9), type=int,
                        help="max number of evaluation with Vegas")
    parser.add_argument("--cores", dest="cores", default=CORES, type=int,
                        help="number of cores used during integration")
    parser.add_argument("-@", dest="order", default='LO', type=str,
                        help="Perturabtive QCD order")
    args = parser.parse_args()
    if args.process is None:
        raise ValueError(
            "Missing process definition: %s --process=epem_a_ttx --run" % sys.argv[0])
    if not any(getattr(args, k) for k in ['run', 'collect', 'generate', 'validate']):
        raise ValueError(
            "\nNeed to specify at leas one of --run, --collect, --generate and --validate\n")

    process_name = args.process
    suffix = args.order
    CORES = args.cores
    # Integrator Settings
    hyper_settings['Integrator']['n_max'] = args.n_max
    hyper_settings['Integrator']['n_new'] = int(1e5)
    hyper_settings['Integrator']['n_start'] = int(1e5)
    hyper_settings['Integrator']['n_increase'] = int(1e5)
    # Selector Settings
    hyper_settings['Selectors']['active_selectors'] = [
    ] if args.no_jets else ['jet']
    hyper_settings['Selectors']['jet'] = {
        'min_jets': 1,
#        'dR': 0.4,
        'min_jpt': args.min_jpt}
    # define alpha loop input folder
    aL_process_output = pjoin(MG_PATH, 'TEST_QGRAF_%s_%s' %
                              (process_name, suffix))
    # Define workspace for the computation of all the diagrams
    WORKSPACE = pjoin(VALIDATION_PATH, '%s_%s' % (process_name, suffix))

    #############
    #  Generate
    #############
    if args.generate:
        MG5 = pjoin(MG_PATH, 'bin', 'mg5_aMC')
        CARD = pjoin(VALIDATION_PATH, 'cards', "%s_%s.aL" %
                     (process_name, suffix))

        r = subprocess.run([MG5, '--mode=alphaloop', CARD],
                           cwd=MG_PATH,
                           capture_output=False)

    # load instructions
    instructions = yaml.load(
        open(pjoin(aL_process_output, 'Rust_inputs',
                   'all_QG_supergraphs.yaml'), 'r'),
        Loader=yaml.Loader)
    #############
    #   RUN
    #############
    if args.run:
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
                            no_jets=args.no_jets,
                            workspace=WORKSPACE)
    #############
    #  COLLECT
    #############
    if args.collect:
        sg_list = [topo['name'] for topo in instructions['topologies']]
        output_file = pjoin(WORKSPACE, "%s_%s_results.csv" %
                            (process_name, args.min_jpt))
        results = []
        for sg in instructions['topologies']:
            filename = glob.glob("%s/*%s_%s*.dat" %
                                 (WORKSPACE, sg['name'], 'nocut' if args.no_jets else args.min_jpt))[0]
            print("Extracting %s: %s" % (sg['name'], filename), end="\n")

            data = [sg['name'], sg['multiplicity']]
            data.extend(get_result(filename))
            results += [data]
        print("")
        df = pd.DataFrame(results, columns=[
            "id", "multiplicity", "neval", "real", "real_err", "imag", "imag_err"])
        df.to_csv(output_file, sep=',', encoding='utf-8', index=False)

    #############
    #  VALIDATE
    #############
    if args.validate:
        # Bench result
        bench_results = pd.read_csv(
            pjoin(VALIDATION_PATH, 'bench', "%s_epem.csv" % suffix))
        # Filter by process
        bench_results = bench_results[bench_results['Process'] == process_name]
        mycut = 'nocut' if args.no_jets else str(args.min_jpt)
        # Filter by min_jpt
        bench_results = bench_results[ bench_results['min_ptj[GeV]'] == mycut]
        print("\033[1mBENCH:\033[0m\n",bench_results)

        # alphaLoop result
        aL_results = pd.read_csv(pjoin(WORKSPACE, "%s_%s_results.csv" %
                            (process_name, args.min_jpt)))
        print("\033[1maL SGs:\033[0m\n", aL_results)
        aL_total = aL_results[['real', 'real_err', 'imag','imag_err']].multiply(aL_results['multiplicity'],axis='index')
        print("\033[1maL Total:\033[0m\n", aL_total.sum())

