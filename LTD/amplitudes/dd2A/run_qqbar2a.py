import os
import sys

import yaml

pjoin = os.path.join


# Define path to runcard
runcard_path = "./runcard_qqbar2a.yaml"

# Get from runcard the location of the alphaloop folder
with open(runcard_path, "r") as f:
    runcard = yaml.load(f, Loader=yaml.Loader)
alphaloop_dir = runcard['directories']['alphaloop']

if True:
    sys.path.append(pjoin(alphaloop_dir, 'LTD/'))
    from amplitudes_rewrite.amp_exporter import AmpExporter


os.system('rm -rf qqbar2a ')
exporter = AmpExporter.from_runcard(runcard_path)
exporter.export()
