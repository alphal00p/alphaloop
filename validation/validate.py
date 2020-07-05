#!/usr/bin/python3.8
"""
Run Validation for different processes
"""


import datetime
import glob
import os
import re
import subprocess
import sys
from argparse import ArgumentParser
from pathlib import Path

import contextlib
import pandas as pd
import progressbar
import yaml

pjoin = os.path.join

class ValidateError(Exception):
    """ Error for the validation phase."""
    pass
# pandas float formatting
pd.options.display.float_format = '{:e}'.format

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


def set_hyperparameters(process_name, sg_name, workspace=None, min_jets=0, multi_settings={}):
    '''Create the correct hyperparameter file'''
    path = pjoin(AL_PATH, 'LTD', 'hyperparameters.yaml')
    hyperparameters = yaml.load(open(path, 'r'), Loader=yaml.Loader)

    try:
        min_jpt = multi_settings['Selectors']['jet']['min_jpt']
    except KeyError:
        min_jpt = hyperparameters['Selectors']['jet']['min_jpt']

    if min_jets == 0:
        min_jpt = 'nocut'
    # Set some default values for the hyperparameter
    # Set values for General
    hyperparameters['General']['deformation_strategy'] = 'none'
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
    hyperparameters['Integrator']['n_vec'] = int(1e4)
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


def run_super_graph(process_name, sg_name, aL_path_output, suffix='', multi_settings={}, min_jets=0, workspace=None, force=False):
    ''' run supergraph for specific process '''
    # Create hyperparameter file
    hyperparameters_path = set_hyperparameters(
        process_name,
        sg_name,
        workspace=workspace,
        min_jets=min_jets,
        multi_settings=multi_settings)
    # run integration
    rust_executable = [rALPHA]
    rust_executable += ['-f', '%s' % hyperparameters_path]
    rust_executable += ['--cross_section',
                        pjoin(aL_path_output, 'Rust_inputs', '%s.yaml' % sg_name)]
    rust_executable += ['-c', '%d' % CORES]

    # Store output
    err_path = pjoin(WORKSPACE, '%s' % os.path.basename(
        hyperparameters_path).replace('yaml', 'err'))
    out_path = pjoin(WORKSPACE, '%s' % os.path.basename(
        hyperparameters_path).replace('yaml', 'out'))
    log_path = pjoin(WORKSPACE, '%s' % os.path.basename(
        hyperparameters_path).replace('yaml', 'log'))

    log_info = {'EventInfo': '', 'IntegrandStatistics': '',
                'StartTime': None, 'LastUpdatedTime': None,
                'ElapsedTime': None, 'Cores': CORES}
    accepted = 0
    rejected = 0
    with progressbar.ProgressBar(prefix='%s | {variables.accepted}\u2713  {variables.rejected}\u2717, res: {variables.real_result} : ' % sg_name,
                                 max_value=multi_settings['Integrator']['n_max'] if force else progressbar.UnknownLength,
                                 variables={'total_samples': '0',
                                            'accepted': '0',
                                            'rejected': '0',
                                            'real_result': 'N/A'}
                                 ) as bar:
        r = subprocess.Popen(rust_executable,
                             cwd=workspace,
                             stdout=subprocess.PIPE,
                             stderr=open(err_path, 'w'),
                             env=os.environ.update(
                                 {'RUST_BACKTRACE': '1', 'MG_NUMERATOR_PATH': '%s/' % aL_path_output}),
                             bufsize=1,
                             universal_newlines=True)
        bar.update(0)
        log_info['StartTime'] = bar.start_time

        total_events = 0
        with open(out_path, 'w') as stdout:
            for line in r.stdout:
                if all(s in line for s in ['[1]', 'chisq', 'df']):
                    split = line.replace("(", "").split(" ")
                    value = float(split[split.index('+-')-1])
                    err = float(split[split.index('+-')+1])
                    df = float(split[split.index('df)\n')-1])
                    if df == 0:
                        chisq = 0
                    else:
                        chisq = float(split[split.index('\tchisq')+1])/df
                    bar.update(
                        real_result="{:e} +- {:.2e} (\u03C7\u00B2 {:.2f})".format(value, err, chisq))
                    stdout.write(line)
                elif 'IntegrandStatistics' in line:
                    log_info['IntegrandStatistics'] = parse_rust_dict(line)
                elif 'EventInfo' in line:
                    log_info['EventInfo'] = parse_rust_dict(line)
                    total_events = log_info['IntegrandStatistics']['total_samples']
                    accepted = log_info['EventInfo']['accepted_event_counter']
                    rejected = log_info['EventInfo']['rejected_event_counter']
                    bar.update(accepted='{:.2%}'.format(
                        accepted/(accepted+rejected)))
                    bar.update(rejected='{:.2%}'.format(
                        rejected/(accepted+rejected)))
                    # Sometimes Vegas exceed in the last iteration the n_max
                    # based on the definition of n_new, n_increase
                    if force and total_events > bar.max_value:
                        bar.max_value = total_events
                    bar.update(total_events)
                else:
                    stdout.write(line)
            log_info['LastUpdatedTime'] = bar.get_last_update_time()
            log_info['ElapsedTime'] = str(
                log_info['LastUpdatedTime']-log_info['StartTime'])
            r.wait()
            if r.returncode != 0:
                print("\033[1;31mFAIL: see {} and {} for further details.\033[0m".format(
                    err_path, out_path))

    with open(log_path, 'w' if force else 'a') as stdlog:
        yaml.dump([log_info], stdlog, default_flow_style=None)


def parse_rust_dict(rust_dict_string):
    """ Parse the printed output of rust into a python dict """
    rust_dict = re.search(r'\{(.*)\}', rust_dict_string).group(0)
    rust_dict = re.sub(
        r'Complex \{ re: (.*?), im: (.*?)\}', r'[\1, \2]', rust_dict)
    rust_dict = re.sub(
        r'Complex \{ re: (.*?), im: (.*?)\}', r'[\1, \2]', rust_dict)
    return eval(re.sub(r' ([^:^,]+):', r'"\1":', rust_dict))


def get_result(filename):
    """ Read the information in the .dat file """
    with open(filename, 'r') as stream:
        try:
            ss = ""
            line = stream.readline()
            while line and '...' not in line:
                ss += line
                line = stream.readline()
            result = yaml.safe_load(ss)
            # If only one output store it as real phase
            if len(result['result']) == 1:
                return [result['neval'],
                        result['result'][0], result['error'][0],
                        'N/A', 'N/A']
            else:
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
    parser.add_argument("--force", action="store_true", dest="force", default=False,
                        help="Allow to overwrite an existent alphaLoop output \
                            and regenerate the state files for the integration")
    parser.add_argument("-r", "--run", action="store_true", dest="run", default=False,
                        help="run the integration for all SG")
    parser.add_argument("-c", "--collect", action="store_true", dest="collect", default=False,
                        help="collect the reults from individual SQ in one CSV file")
    parser.add_argument("-v", "--validate", action="store_true", dest="validate", default=False,
                        help="compare with bench results")
    parser.add_argument("--refine", dest="refine", default=0, type=int,
                        help="run refinement REFINE number of times integrating the worst SG")
    parser.add_argument("--min_jets", dest="min_jets", default=0, type=int,
                        help="minimum number of jets in the dinal state")
    parser.add_argument("--min_jpt",  dest="min_jpt", default=100, type=int,
                        help="Jet cutoff")
    parser.add_argument("--n_max", dest="n_max", default=int(1e9), type=int,
                        help="max number of evaluation with Vegas")
    parser.add_argument("--multi_channeling", action="store_true", dest="multi_channeling", default=False,
                        help="enable multi-channeling during integration")
    parser.add_argument("--cores", dest="cores", default=CORES, type=int,
                        help="number of cores used during integration")
    parser.add_argument("-@", dest="order", default='LO', type=str,
                        help="perturabtive QCD order")
    parser.add_argument("--diag_name", dest="diag_name", default=None, type=str,
                        help="Integrate single SG: selecting by name")
    parser.add_argument("--diag_id", dest="diag_id", default=None, type=int,
                        help="Integrate single SG: selecting by position in the collected table")
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
    hyper_settings['Integrator']['integrated_phase'] = 'real'
    hyper_settings['Integrator']['state_filename_prefix'] = "%s_%s_" % (
        process_name, "nocut" if args.min_jets == 0 else args.min_jpt)
    hyper_settings['Integrator']['keep_state_file'] = True
    hyper_settings['Integrator']['load_from_state_file'] = True
    hyper_settings['Integrator']['n_max'] = args.n_max
    hyper_settings['Integrator']['n_new'] = int(1e4)
    hyper_settings['Integrator']['n_start'] = int(1e4)
    hyper_settings['Integrator']['n_increase'] = int(1e4)
    # Selector Settings
    hyper_settings['Selectors']['active_selectors'] = [
    ] if args.min_jets == 0 else ['jet']
    hyper_settings['Selectors']['jet'] = {
        'min_jets': args.min_jets,
        'dR': 0.4,
        'min_jpt': args.min_jpt}
    # Set multi-channeling strategy
    hyper_settings['General']['multi_channeling'] = args.multi_channeling
    hyper_settings['General']['multi_channeling_including_massive_propagators'] = args.multi_channeling

    # Define paths to executables, files and workspace
    aL_process_output = pjoin(MG_PATH, 'TEST_QGRAF_%s_%s' %
                              (process_name, suffix))
    CROSS_SECTION_SET = pjoin(aL_process_output, 'Rust_inputs',
                              'all_QG_supergraphs.yaml')
    WORKSPACE = pjoin(VALIDATION_PATH, '%s_%s' % (process_name, suffix))
    COLLECTION_PATH = pjoin(WORKSPACE, "%s_%s_results.csv" %
                            (process_name, 'nocut' if args.min_jets == 0 else args.min_jpt))
    MG5 = pjoin(MG_PATH, 'bin', 'mg5_aMC')
    CARD = pjoin(VALIDATION_PATH, 'cards', "%s_%s.aL" % (process_name, suffix))

    #############
    #  Generate
    #############
    if args.generate:
        if os.path.exists(aL_process_output) and not args.force:
            print(
                "\033[1;32;48mSkip generation use --force to regenerate!\033[0m")
        else:
            r = subprocess.run([MG5, '--mode=alphaloop', CARD],
                               cwd=MG_PATH,
                               capture_output=False)

    # load instructions
    try:
        instructions = yaml.load(
            open(CROSS_SECTION_SET, 'r'), Loader=yaml.Loader)
    except FileNotFoundError:
        print(
            "\033[1;31mFAIL: missing {}.\nHave you run with --generate?\033[0m".format(CROSS_SECTION_SET))
        raise

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
        for sg_id, sg_name in enumerate(sg_list):
            VEGAS_STATE_FILE = pjoin(WORKSPACE, "%s_%s_%s_state.dat" %
                                     (process_name, 'nocut' if args.min_jets == 0 else args.min_jpt, sg_name))
           
            if args.diag_id is not None:
                if sg_id != args.diag_id:
                    continue
            elif args.diag_name is not None:
                if sg_name != args.diag_name:
                    continue
            if args.force:
                with contextlib.suppress(FileNotFoundError):
                    os.remove(VEGAS_STATE_FILE)
            elif os.path.exists(VEGAS_STATE_FILE):
                print(
                    "\033[1;32;48mLoading state!\033[0m (Use --force to overwrite)")
            print(VEGAS_STATE_FILE, os.path.exists(VEGAS_STATE_FILE))
            run_super_graph(process_name,
                            sg_name,
                            aL_process_output, suffix='_%s' % suffix,
                            multi_settings=hyper_settings,
                            min_jets=args.min_jets,
                            workspace=WORKSPACE,
                            force=args.force or not os.path.exists(VEGAS_STATE_FILE))

    #############
    #  COLLECT
    #############
    if args.collect:
        sg_list = [topo['name'] for topo in instructions['topologies']]
        results = []
        for sg in instructions['topologies']:
            try:
                filename = glob.glob("%s/*%s_%s*.dat" %
                                     (WORKSPACE, sg['name'], 'nocut' if args.min_jets == 0 else args.min_jpt))[0]
            except IndexError:
                print(
                    "\033[1;31mFAIL: Cannot find .dat file for {0}.\nTry to generate it with --run --diag_name={0}\033[0m".format(sg['name']))
                raise

            #print("Extracting %s: %s" % (sg['name'], filename), end="\n")

            data = [sg['name'], sg['multiplicity']]
            data.extend(get_result(filename))
            eval_time = pd.Timedelta(0)
            for log in yaml.safe_load(open(filename.replace('_res.dat', '.log'), 'r')):
                eval_time += pd.to_timedelta(log['ElapsedTime'])
            data.extend([eval_time])
            results += [data]

        df = pd.DataFrame(results, columns=[
            "name", "multiplicity", "neval", "real", "real_err", "imag", "imag_err", 'eval_time'])
        df.to_csv(COLLECTION_PATH, sep=',', encoding='utf-8', index=False)

    #############
    #  Refine
    #############
    if args.refine > 0:
        CALL_BASE_ARGS = [arg % args.__dict__ for arg in ['python', os.path.basename(__file__),
                                                          "--process=%(process)s",
                                                          "-@%(order)s",
                                                          "--n_max=%(n_max)d",
                                                          "--cores=%(cores)d",
                                                          "--min_jets=%(min_jets)d",
                                                          "--min_jpt=%(min_jpt)d"]]
        if args.multi_channeling:
            CALL_BASE_ARGS += ["--multi_channeling"]
        
        # Collect
        r = subprocess.run(CALL_BASE_ARGS + ['-cv'])
        if r.returncode != 0:
            raise ValidateError()
        for ref_i in range(args.refine):
            # Extract worst diag
            worst_SG = pd.read_csv(COLLECTION_PATH).sort_values(by=['real_err'], ascending=False)['name'].values[0]
            print("\033[1;32;48mRefine (%d/%d): %s\033[0m"%(ref_i+1,args.refine,worst_SG))
            # Run with more points
            subprocess.run(CALL_BASE_ARGS + ['-rc', '--diag_name=%s'%worst_SG])
            if r.returncode != 0:
                raise ValidateError()
            ranked = pd.read_csv(COLLECTION_PATH).sort_values(by=['real_err'], ascending=False)
            print("\033[1;32;48mMoved down by %s positions!\033[0m"%(list(ranked['name']==worst_SG).index(True)))
            print("\033[1maL SGs ERROR SORT:\033[0m\n",ranked)
            print()
    
    #############
    #  VALIDATE
    #############
    import numpy as np
    if args.validate:
        print(
            "\033[1mProcess:\033[0m\n\tname=%(process)s, min_jets=%(min_jets)d, min_jpt=%(min_jpt)d" % args.__dict__)
        # Bench result
        bench_results = pd.read_csv(
            pjoin(VALIDATION_PATH, 'bench', "%s_epem.csv" % suffix))
        # Filter by process
        bench_results = bench_results[bench_results['Process'] == process_name]
        # Filter by min_jpt
        if args.min_jets == 0:
            bench_results = bench_results[bench_results['min_ptj[GeV]'].isnull(
            )]
        else:
            bench_results = bench_results[bench_results['min_ptj[GeV]']
                                          == args.min_jpt]
        print("\033[1mBENCH:\033[0m\n", bench_results)

        # alphaLoop result
        aL_results = pd.read_csv(COLLECTION_PATH)
        print("\033[1maL SGs:\033[0m\n", aL_results)
        aL_total_res = aL_results[['real', 'imag']]\
            .multiply(aL_results['multiplicity'], axis='index').sum()
        aL_total_err = np.sqrt((aL_results[['real_err', 'imag_err']]
                                .multiply(aL_results['multiplicity'], axis='index')**2).sum())
        print("\033[1maL SGs ERROR SORT:\033[0m\n",
              aL_results.sort_values(by=['real_err'], ascending=False))

        print("\033[1mCompare:\033[0m")
        diff = (aL_total_res['real'] - bench_results['cross-section[pb]'])
        rel = abs(diff/bench_results['cross-section[pb]'])
        mg5_prec = bench_results['MCerror[pb]'] / \
            abs(bench_results['cross-section[pb]'])
        aL_prec = [aL_total_err['real_err']/abs(aL_total_res['real'])]
        if not bench_results.empty:
            print("\t\033[1mMG res:\033[0m {:.5e} +/- {:.5e}".format(
                bench_results['cross-section[pb]'].values[0],
                bench_results['MCerror[pb]'].values[0]))
        print("\t\033[1maL res:\033[0m {:.5e} +/- {:.5e}".format(
            aL_total_res['real'], aL_total_err['real_err']))
        print()
        print("\t\033[1mMG precision:\033[0m", mg5_prec.values)
        print("\t\033[1maL precision:\033[0m", aL_prec)
        print("\t\033[1m|(MG-aL)/MG|:\033[0m", rel.values)
