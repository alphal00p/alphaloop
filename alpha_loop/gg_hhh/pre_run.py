import yaml
import copy
import re
import math
import progressbar
import yaml
import sys
import os
import subprocess
from pathlib import Path
import numpy as np

RUST_PATH = "/home/armin/my_programs/pynloop/rust_backend"
CROSS_SEC_NAME ='ggHHHsgs'
do_prerun = True
pre_comp = 100000
num_rust_cores = 36
prec_goal = 1.e-4

if True:
    pjoin = os.path.join
    root_path = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, os.path.normpath(pjoin(root_path,'TMPDIR')))


cs_file=CROSS_SEC_NAME + '.yaml'
cs_file_PATH = os.path.normpath(pjoin(root_path,'TMPDIR'))

with open(cs_file_PATH+'/'+cs_file) as file:
    cs = yaml.load(file,Loader=yaml.FullLoader)


sys_cmd_exp = 'export MG_NUMERATOR_PATH='+root_path+'/'
sys_cmd_cd_rust = 'cd '+ RUST_PATH

if do_prerun:
    for tp in cs['topologies']:
        sys_cmd_run ='cargo run --bin ltd -- --cross_section '+ cs_file_PATH +'/'+ tp['name'] +'.yaml -c='+str(num_rust_cores) +' -s=' + str(pre_comp)
        os.system(sys_cmd_exp + '&&' + sys_cmd_cd_rust + '&&' + sys_cmd_run)

max_res_abs = 0.
for tp in cs['topologies']:
    with open(RUST_PATH+'/'+tp['name']+'_res'+'.dat') as file:
        res_abs = np.linalg.norm(np.array((yaml.load(file,Loader=yaml.FullLoader))['result']))
        tp['res_abs'] = copy.copy(res_abs)
        if max_res_abs < res_abs:
            max_res_abs = copy.copy(res_abs)

modified_cs = []
for tp in cs['topologies']:
    if np.absolute(tp['res_abs']/max_res_abs)>prec_goal:
        modified_cs+=[ {k:tp[k] for k in ('additional_LMBs','multiplicity','name') if k in tp}]
new_cs = {}
new_cs['name'] = CROSS_SEC_NAME + 'mod'
new_cs['topologies'] = modified_cs

with open(cs_file_PATH+'/modified_'+CROSS_SEC_NAME+'.yaml', 'w') as outfile:
    yaml.dump(new_cs, outfile)
