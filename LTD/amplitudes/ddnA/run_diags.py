#!/bin/python 
import os
import subprocess
import sys

import numpy as np
import yaml

pjoin = os.path.join

# Read complete runcard as dictionary
with open("./runcard_qqbarNa.yaml", 'r') as f:
    runcard = yaml.load(f, Loader=yaml.Loader)

# Read amplitudes
with open("./refs/amplitudes_v2.yaml", 'r') as f:
    amps = yaml.load(f, Loader=yaml.Loader)

# Get from runcard the location of the alphaloop folder
alphaloop_dir = runcard['directories']['alphaloop']
if True:
    sys.path.append(pjoin(alphaloop_dir, 'LTD/'))
    from amplitudes_rewrite.amp_exporter import AmpExporter


# run diagrams
for n_photon in [4]:
    # Get Polarizations
    with open("./qqbar%da_pols.yaml" % n_photon, 'w') as output:
        cmd = ['./get_polarizations.sh dd%dA' % n_photon]
        pols = subprocess.check_output(
            cmd, shell=True, cwd='./refs').decode()
        output.write(pols)
    with open("./qqbar%da_pols.yaml" % n_photon, 'r') as f:
        pols = yaml.load(f, Loader=yaml.Loader)

    # Update Runcard
    amp = next(filter(lambda x: x['name'] == "dd%dA" % n_photon, amps))
    runcard['amplitude']['name'] = 'qqbar%da' % n_photon
    runcard['amplitude']['diagram_list'] = './qqbar%da.json' % n_photon
    runcard['amplitude'].update(pols)
    runcard['amplitude']['external_data']['in_momenta'] = amp['ps'][:2]
    runcard['amplitude']['external_data']['out_momenta'] = [
        (-np.array(q)).tolist() for q in amp['ps'][2:-1]]

    runcard_path = './runcard_qqbar%da.yaml' % n_photon
    with open(runcard_path, 'w') as output:
        yaml.dump(runcard, output, default_flow_style=None)

    os.system('rm -rf qqbar%da' % n_photon)
    exporter = AmpExporter.from_runcard(runcard_path)
    exporter.export()
    os.system('rm %s' % runcard_path)
    os.system('rm qqbar%da_pols.yamnl' % n_photon)
