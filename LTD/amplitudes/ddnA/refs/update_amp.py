#!/bin/python
import yaml

old_amplitude_file = "./amplitudes.yaml"
old_topolgy_file = "./topologies.yaml"


# UPDATE amplitude file
with open(old_amplitude_file, 'r') as f:
    amps = yaml.load(f, Loader=yaml.Loader)

for amp in amps:
    amp['amp_type'] = 'qqbarphotonsNLO'
    amp['n_loops'] = 1
    for diag in amp['diagrams']:
        diag['tensor_coefficients_split'] = []
        diag['denominators'] = [[0, p] for p in diag['denominators']]

with open("./amplitudes_v2.yaml", 'w') as stream:
    stream.write(yaml.dump(amps, Dumper=yaml.Dumper, default_flow_style=None))

# UPDATE topology file
with open(old_topolgy_file, 'r') as f:
    topos = yaml.load(f, Loader=yaml.Loader)

for top in topos:
    prop_n = 1
    top['maximum_ratio_expansion_threshold'] = -1.0
    for ll in top['loop_lines']:
        for prop in ll['propagators']:
            prop['name'] = "l%d"%prop_n
            prop_n += 1

with open("./topologies_v2.yaml", 'w') as stream:
    stream.write(yaml.dump(topos, Dumper=yaml.Dumper, default_flow_style=None))
