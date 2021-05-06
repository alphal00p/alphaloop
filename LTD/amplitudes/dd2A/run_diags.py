import os
import sys

import yaml

pjoin = os.path.join
# Import reference runcard
with open("./runcard_qqbar2a.yaml", 'r') as f:
    runcard = yaml.load(f, Loader=yaml.Loader)

# Get from runcard the location of the alphaloop folder
alphaloop_dir = runcard['directories']['alphaloop']
if True:
    sys.path.append(pjoin(alphaloop_dir, 'LTD/'))
    from amplitudes_rewrite.amp_exporter import AmpExporter

# run diagrams
#for diag_n in range(1, 11):
#for diag_n in [4,8,10]:
for diag_n in [2]:
    runcard['amplitude']['name'] = 'qqbar2a_%d' % diag_n
    runcard['amplitude']['diagram_list'] = './qqbar2a_%d.json' % diag_n

    runcard_path = './runcard_qqbar2a_%d.yaml' % diag_n
    with open(runcard_path, 'w') as output:
        yaml.dump(runcard, output)

    os.system('rm -rf qqbar2a_%d' % diag_n)
    exporter = AmpExporter.from_runcard(runcard_path)
    exporter.export()
    os.system('rm %s' % runcard_path)
