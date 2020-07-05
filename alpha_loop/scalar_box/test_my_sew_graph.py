import copy
import re
import math
import progressbar
import yaml
import sys
import os
import subprocess
from pathlib import Path

import logging
# Print logger coming from alphaLoop.FORM_processing
logger = logging.getLogger('alphaLoop.FORM_processing')
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('>> %(asctime)s - %(name)s: line %(lineno)d - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

pdg = {
    'd': 1, 'd~': -1,
    'u': 2, 'u~': -2,
    's': 3, 's~': -3,
    'c': 4, 'c~': -4,
    'b': 5, 'b~': -5,
    't': 6, 't~': -6,

    'e-': 11, 'e+': -11,

    'g': 21, 'a': 22,

    'h': 25
}

MG_PATH = "/home/armin/my_programs/MG5_aMC_v3_0_2_py3/"

if True:
    pjoin = os.path.join
    root_path = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, os.path.normpath(pjoin(root_path, "../../",'alpha_loop')))
    sys.path.insert(0, os.path.normpath(pjoin(root_path, "../../")))
    sys.path.insert(0, MG_PATH)
    import alpha_loop.interface as interface

    from LTD.ltd_utils import Colour
    from FORM_processing import FORM_processing_options, FORMSuperGraphList, FORMProcessor
    import models.model_reader as model_reader


def standalone_qgraph_files(qg_supergraphs_file, output_file):
    print('Create standalone QGRAF from: {}'.format(qg_supergraphs_file))

    with open(output_file, 'w') as f:
        f.write("graphs=[]\n")
        f.write("graph_names=[]\n")
        with open(qg_supergraphs_file, 'r') as stream:
            n_graph = 0
            for line in stream.read().splitlines():
                if 'graphs=[]' in line or 'graph_names=[]' in line:
                    continue
                if 'graphs.append' in line:
                    f.write('graph_names+=[\"SG_QG{}\"]\n'.format(n_graph))
                    n_graph += 1
                f.write(line+'\n')


# path to the python diagram output files
diagrams_python_source = None
# Path to UFO model to load.
model = 'aL_sm-no_widths'
# Process definition to consider.
process = 'g g > h h'
# perturbative order
perturbative_order = 'LO'
# Number of FORM cores
cores = 2
# Number of iterations for numerator optimization
optim_iter = 100
# Always generate the full energy polynomial
full_energy_poly = False
# Model restriction card to consider.
restrict_card = os.path.normpath(pjoin(root_path, "../../", 'models',
                      'aL_sm', 'restrict_no_widths.dat'))


####################################################
# Import model
####################################################
if cores > 1:
    FORM_processing_options['parallel'] = True
    FORM_processing_options['cores'] = cores
FORM_processing_options['extra-options'] = '-D OPTIMITERATIONS=' + \
    str(optim_iter)
if full_energy_poly:
    FORM_processing_options['extra-options'] += ' -D FULL_ENERGY_POLY=1'

cli = interface.alphaLoopInterface()

cli.do_import('model %s' % model)
computed_model = model_reader.ModelReader(cli._curr_model)
computed_model.set_parameters_and_couplings(restrict_card)
process_definition = cli.extract_process(process, proc_number=0)


# Generated the diagram file
#diag_path_qg = root_path+'/test_process_intefered/qqbaaNLO.py'
#diagrams_python_source_QG = 'all_QG_supergraphs.py'
#standalone_qgraph_files(diag_path_qg, diagrams_python_source_QG)
diagrams_python_source_QG = '/home/armin/my_programs/pynloop/alpha_loop/scalar_box/my_single_sg.py'

# Import
# Here I should include valentin PS-generator
external_momenta = {
    'q1':[5.0000000000000000e+02,  0.0000000000000000e+00,  0.0000000000000000e+00,  5.0000000000000000e+02],
    'q2':[5.0000000000000000e+02,  0.0000000000000000e+00,  0.0000000000000000e+00, -5.0000000000000000e+02],
    'q3':[5.0000000000000000e+02,  0.0000000000000000e+00,  0.0000000000000000e+00,  5.0000000000000000e+02],
    'q4':[5.0000000000000000e+02,  0.0000000000000000e+00,  0.0000000000000000e+00, -5.0000000000000000e+02]}
cut_momentum ={
    'c1':[4.9999999999999977e+02,  1.0740197658523471e+02,  4.3070555989194781e+02, -1.9321629357729731e+02],
    }

super_graphs_list = FORMSuperGraphList.from_dict(
    diagrams_python_source_QG, merge_isomorphic_graphs=False, verbose=True, model=computed_model,
    external_mom=external_momenta,fixed_cut_mom=cut_momentum,num_fixed_loops=1,is_scalar=True)

# Define number of jets and final state ids
n_jets = 0
final_state_particle_ids = []
for p in process.split('>')[-1].split(' '):
    if p == '':
        continue
    if abs(pdg[p]) in [1, 2, 3, 4, 21]:
        n_jets += 1
    else:
        final_state_particle_ids += [pdg[p]]
final_state_particle_ids = tuple(final_state_particle_ids)
print("\nfs_ids:{}, n_jets: {}\n".format(final_state_particle_ids, n_jets))

#########################################
# Define system
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

invalid_graphs = {'nocut': [], 'zero': []}


#########################################
# Run isomorphisms for QG
nocut_graphs_names = []
for graphs in super_graphs_list:
    for g in graphs:
        nocut_graphs_names += [g.name]

# Initialize FORMProcessor

form_processor = FORMProcessor(
    super_graphs_list, computed_model, process_definition)
form_processor.generate_squared_topology_files(
    TMP_OUTPUT, n_jets, final_state_particle_ids=final_state_particle_ids)


for graphs in form_processor.super_graphs_list:
    for g in graphs:
        nocut_graphs_names.remove(g.name)
for name in nocut_graphs_names:
    invalid_graphs['nocut'] += [name]

#########################################
# Build numerators
subprocess.run(['cp', os.path.normpath(pjoin(root_path, '../Templates/Source_make_opts')),
                os.path.normpath(pjoin(FRM_OUTPUT, 'make_opts'))])
subprocess.run(['cp', os.path.normpath(pjoin(
    root_path, '../Templates/FORM_output_makefile')), os.path.normpath(pjoin(FRM_OUTPUT, 'Makefile'))])


for graphs in form_processor.super_graphs_list:
        print("\t", [g.name for g in graphs])

#########################################
# Export mathematica drawings
diag_num = 0
for n, iso_graphs in enumerate(form_processor.super_graphs_list):
    for g in iso_graphs:
        with open(root_path+'/draws/Graph_{:06}.m'.format(diag_num), 'w') as f:
            pdf_file = 'Graph_{:06}.pdf'.format(diag_num)
            f.write(
                "GraphClass = If[$VersionNumber > 12, EdgeTaggedGraph, Graph];\n")
            f.write(
                "CreateEdge[u_,v_,t_]:=If[$VersionNumber > 12, DirectedEdge[u, v, t], DirectedEdge[u, v]];\n")
            f.write("aGraph=")
            f.write(g.get_mathematica_rendering_code(
                computed_model, FORM_id=n))
            f.write(
                ";\nExport[\""+pdf_file+"\", GraphicsGrid[{{aGraph}}], ImageSize -> {825.000000, 637.500000}];")
            diag_num += 1

#########################################
# Export numerator
form_processor.generate_numerator_functions(
    FRM_OUTPUT, output_format='c', workspace=FRM_WORKSPACE)

#########################################
# Check for zero numerators
for n, iso_graphs in enumerate(form_processor.super_graphs_list):
    if iso_graphs.is_zero:
        invalid_graphs['zero'] += [graph.name for graph in iso_graphs]




#########################################
# Print valid isomorphisms
print("\nInvalid Cuts: {}\n".format(invalid_graphs))