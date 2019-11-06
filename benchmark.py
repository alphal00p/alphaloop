import subprocess
from tabulate import tabulate
from uncertainties import ufloat
import time
import yaml
import json
from tqdm import tqdm
import os
import argparse

W  = '\033[0m'  # white (normal)
R  = '\033[31m' # red
G  = '\033[32m' # green
O  = '\033[33m' # orange
B  = '\033[34m' # blue
P  = '\033[35m' # purple

def get_rust_result(topology, phase, num_samples, cores):
    # remove info from old runs
    try:
        os.remove("rust_backend/%s_res.dat" % topology)
    except:
        pass

    # TODO: check for errors
    output = subprocess.check_output(["cargo", "run", "--quiet", "--release", "--bin", "ltd", "--", "-t", 
        topology, "-s", str(num_samples), "-c", str(cores)], cwd="rust_backend").strip().decode()

    # read the output file
    with open("rust_backend/%s_res.dat" % topology, 'r') as f:
        result = yaml.safe_load(f)

    # get git revision
    git_revision = subprocess.check_output(["git", "describe", "--always"]).strip().decode()

    # get git diff
    git_diff = subprocess.check_output(["git", "diff"]).strip().decode()

    # get analytical result
    with open("LTD/topologies.yaml", 'r') as f:
        topologies = yaml.safe_load(f)
        for t in topologies:
            if t['name'] == topology:
                analytical_result = (t['analytical_result_real'], t['analytical_result_imag'])

    return {
        'revision': git_revision,
        'diff': git_diff,
        'num_samples': num_samples,
        'topology': topology,
        'result': result['result'],
        'error': result['error'],
        'analytical_result': analytical_result,
        'time': time.time()
    }

def get_history():
    historical_data = []

    try:
        with open("historical_benchmarks.json", 'r') as f:
            historical_data = json.load(f)
    except:
        pass

    return historical_data

def save_to_history(samples):
    historical_data = []
    try:
        with open("historical_benchmarks.json", 'r') as f:
            historical_data = json.load(f)
    except:
        pass

    historical_data.extend(samples)

    with open('historical_benchmarks.json', 'w') as f:
        json.dump(historical_data, f)

def get_score_for_sample(sample):
    scores = [None, None]

    for phase in [0, 1]:
        if sample['result'][phase] is not None:
            scores[phase] = abs(sample['analytical_result'][phase] - sample['result'][phase]) / sample['error'][phase]

    return scores

def render_data(samples):
    """ Render the data in a table. """
    data = []

    for sample in samples:
        score = get_score_for_sample(sample)

        for (phase, phase_name) in [(0, 'Real'), (1, 'Imag')]:
            if sample['result'][phase] is None:
                continue
            
            data.append(
                [sample['topology'] + ' ' + phase_name, "{:,}".format(int(sample['num_samples'])),
                    ufloat(sample['result'][phase], sample['error'][phase]), 
                    sample['analytical_result'][phase],
                    R + str(score[phase]) + W if score[phase] > 2.0 else G + str(score[phase]) + W,
                    sample['revision'], sample['diff'] == '']
            )

    print(tabulate(data, ['Topology', '# Samples', 'Result', 'Reference', 'Score', 'Tag', 'Clean'], tablefmt="fancy_grid"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Benchmark hyperparameters')
    parser.add_argument('topologies', metavar='topologies', type=str, nargs='+',
                        help='topologies to test')
    parser.add_argument('--from_history', action='store_true', help='Read the topology data from the history')
    parser.add_argument('-s', default='100000', help='number of samples')
    parser.add_argument('-c', default='4', help='number of cores')
    parser.add_argument('--phase', default='both', choices=['real','imag','both'], help='the phase for the integration')
    args = parser.parse_args()

    samples = []
    
    if args.from_history:
        historical_data = get_history()
        samples = [d for d in historical_data if d['topology'] in args.topologies]
    else:
        pbar = tqdm(args.topologies)
        for topology in pbar:
            pbar.set_description(topology)
            result = get_rust_result(topology, args.phase, args.s, args.c)
            samples.append(result)

    render_data(samples)

    # ask to save data
    if not args.from_history and input("Do you want to save the new run? [y/N]: ") in ['y','Y']:
       save_to_history(samples)
