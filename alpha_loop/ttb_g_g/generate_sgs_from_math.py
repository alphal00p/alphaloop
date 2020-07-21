import shutil
import copy
import re
import math
import progressbar
import yaml
import sys
import os
import subprocess
from pathlib import Path


#########################################
# Define system
if True:
    pjoin = os.path.join
    root_path = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, os.path.normpath(pjoin(root_path, 'test_process')))
    sys.path.insert(0, os.path.normpath(pjoin(root_path, "../../",'alpha_loop')))
    sys.path.insert(0, os.path.normpath(pjoin(root_path, "../../")))
    import ttbaaLeftDiags as ttbaald
    import amplitudes as amp

TMP_OUTPUT = pjoin(root_path, 'TMPDIR')
FRM_OUTPUT = pjoin(root_path, 'FORM')
FRM_WORKSPACE = pjoin(FRM_OUTPUT,'workspace')

## Clear all folders
subprocess.run(['rm', '-r', TMP_OUTPUT, ], capture_output=True)
subprocess.run(['rm', '-r', FRM_OUTPUT, ], capture_output=True)

## Build output system
Path(FRM_OUTPUT).mkdir(parents=True, exist_ok=True)
Path(FRM_WORKSPACE).mkdir(parents=True, exist_ok=True)
Path(TMP_OUTPUT).mkdir(parents=True, exist_ok=True)

## BUILD SUPEGRPAPHS
left_diags = ttbaald.graphs
right_diags = []

for g in [left_diags[1]]:
    if 'analytic_num' not in g.keys():
        g['analytic_num'] = "1"
    right_diags+=[amp.to_right_diag(g,to_effective_vertex=True)]

count = 0
super_graphs =[]
for ld in left_diags:
    for rd in right_diags:        
        super_graphs+=[amp.sew_amp_diags(ld,rd)]

## PERFORM COLOR DECOMPOSITION
al_path=os.path.normpath(pjoin(root_path, "../../",'alpha_loop'))
all_c_decom_graphs = amp.perform_color_decomp(super_graphs,FRM_WORKSPACE,al_path)

for i,gs in enumerate(all_c_decom_graphs):
    amp.save_dict_to_file(gs,'ttbaa_colorStruc'+str(i+1))
