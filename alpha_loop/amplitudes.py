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
import importlib.util


###############################################################
#
# Define functions amplitude preparation (inteferences)
#
#############################################################


def append_loop_moms(graph):
    """gets the loop-momenta from a graph in qgraph format

    Args:
        graph ([dict]): [qgraf-format]

    Returns:
        [dict]: [qgraf-format]
    """
    if 'loop_momenta' in graph:
        return graph
    mom_str = ""
    for ee in graph['edges'].values():
        mom_str += ' ' + ee['momentum']
    p = re.compile(r'k\d*')
    loopmom = p.findall(mom_str)
    graph['loop_momenta'] = tuple(set(loopmom))
    return graph


def impose_mom_conservation(graph):
    ext_mom = graph['in_momenta']+graph['out_momenta']
    num_legs = 0
    for ee in graph['edges'].values():
        if (ee['type'] == 'out' or ee['type'] == 'in'):
            num_legs += 1
    if len(ext_mom) == num_legs:
        lhs1 = "-"+graph['out_momenta'][-1]  # -pLast
        lhs2 = graph['out_momenta'][-1]  # +pLast
        rhs1 = "-"+graph['in_momenta'][-1]
        rhs2 = graph['in_momenta'][-1]
        for m in graph['in_momenta'][:-1]:
            rhs1 += "-" + m
            rhs2 += "+" + m
        for m in graph['out_momenta'][:-1]:
            rhs1 += "+" + m
            rhs2 += "-" + m
        for ee in graph['edges'].values():
            ee['momentum'] = ee['momentum'].replace(lhs1, rhs1)
            ee['momentum'] = ee['momentum'].replace(lhs2, rhs2)
        for nn in graph['nodes'].values():
            # somehow momenta of in edges are not tuples but strings
            if isinstance(nn['momenta'], str):
                nn['momenta'] = nn['momenta'].replace(lhs1, rhs1)
                nn['momenta'] = nn['momenta'].replace(lhs2, rhs2)
            else:
                nn['momenta'] = tuple(x.replace(lhs1, rhs1)
                                      for x in nn['momenta'])
                nn['momenta'] = tuple(x.replace(lhs2, rhs2)
                                      for x in nn['momenta'])
        graph['out_momenta'] = tuple(x.replace(lhs2, rhs2)
                                     for x in graph['out_momenta'])
        graph['analytic_num'] = graph['analytic_num'].replace(lhs1, rhs1)
        graph['analytic_num'] = graph['analytic_num'].replace(lhs2, rhs2)
    return graph


def to_right_diag(left_diag, to_effective_vertex=False):
    """Translates left-diagram to right-diagram (complex-conjugation + flow reversal)

    Args:
        left_diag ([dict]): [left-diagram]
    """
    def rev_tup(tuples):
        new_tup = tuples[::-1]
        return new_tup

    # make edge labeling continous: since n1 > n2 scattering with n1!=n2 will introduce negative edge ids
    new_edge_ids = list(str(x) for x in sorted(
        list(range(1, len(left_diag['edges'])+1))))
    old_edge_ids = sorted(list((left_diag['edges']).keys()))
    
    for i, ne in enumerate(new_edge_ids):
        for nn in (left_diag['nodes']).values():
            nn['edge_ids'] = tuple(
                map(lambda x: ne if x == old_edge_ids[i] else x, nn['edge_ids']))
    for i, ee in enumerate(list(left_diag['edges'])):
        left_diag['edges'][new_edge_ids[i]
                           ] = left_diag['edges'].pop(old_edge_ids[i])
    for nn in (left_diag['nodes']).values():
        nn['edge_ids'] = tuple(int(x) for x in nn['edge_ids'])
    for ee in list((left_diag['edges']).keys()):
        left_diag['edges'][int(ee)] = left_diag['edges'].pop(ee)
    right_diag = copy.deepcopy(left_diag)

    # merge edges if we want effective vertex only
    if to_effective_vertex:
        ec = 0
        nc = 0
        # delete edges
        for ee in list((right_diag['edges']).keys()):
            if right_diag['edges'][ee]['type'] == 'virtual':
                del right_diag['edges'][ee]
            else:
                ec += 1
        # delete "internal" vertices
        for nn in list((right_diag['nodes']).keys()):
            if right_diag['nodes'][nn]['vertex_id'] == 0:
                del right_diag['nodes'][nn]
            else:
                nc += 1
        # relabel edge id's
        new_edge_ids = list(str(x) for x in sorted(
            list(range(1, len(right_diag['edges'])+1))))
        old_edge_ids = sorted(list((right_diag['edges']).keys()))

        for i, ne in enumerate(new_edge_ids):
            for nn in (right_diag['nodes']).values():
                nn['edge_ids'] = tuple(
                    map(lambda x: ne if x == old_edge_ids[i] else x, nn['edge_ids']))

        for i, ee in enumerate(list(right_diag['edges'])):
            right_diag['edges'][new_edge_ids[i]
                                ] = right_diag['edges'].pop(old_edge_ids[i])
        for nn in (right_diag['nodes']).values():
            nn['edge_ids'] = tuple(int(x) for x in nn['edge_ids'])
        for ee in list((right_diag['edges']).keys()):
            right_diag['edges'][int(ee)] = right_diag['edges'].pop(ee)

        # add internal vertex (fake vertex)
        right_diag['nodes'][nc+1] = {'edge_ids': tuple(x for x in list((right_diag['edges']).keys())),
            'vertex_id': 0,
            'momenta': tuple(map(lambda x: right_diag['edges'][x]['momentum'] if right_diag['edges'][x]['type'] == 'in' else '-'+str(right_diag['edges'][x]['momentum']), list((right_diag['edges']).keys()))),
            'indices': (), 'fake_vertex': True}

        for ee in (right_diag['edges']).values():
            ee['vertices'] = tuple(
                map(lambda x: nc+1 if x > nc else x, list((ee['vertices']))))

    # complex conjugation. will be performed in form
    right_diag['analytic_num'] = "hermconjugate(" + \
        right_diag['analytic_num']+")"
    if to_effective_vertex:
        right_diag['analytic_num'] = "hermconjugate(1)"
    node_offset = max(len(left_diag['nodes']), len(right_diag['nodes']))+1
    edge_offset = max(len(left_diag['edges']), len(right_diag['edges']))+1

    # reverse edge orientations but keep momentum flow
    for ee in list(right_diag['edges']):
        right_diag['edges'][ee]['vertices'] = tuple(
            x+node_offset for x in rev_tup(right_diag['edges'][ee]['vertices']))
        if right_diag['edges'][ee]['type'] == 'in':
            right_diag['edges'][ee]['type'] = 'out'
        elif right_diag['edges'][ee]['type'] == 'out':
            right_diag['edges'][ee]['type'] = 'in'

    for ee in list(right_diag['edges']):
        newEdgeKey = ee+edge_offset
        right_diag['edges'][ee]['name'] = right_diag['edges'][ee]['name']+'r'
        right_diag['edges'][newEdgeKey] = right_diag['edges'].pop(ee)

    for nn in list(right_diag['nodes']):
        right_diag['nodes'][nn]['edge_ids'] = tuple(
            x+edge_offset for x in right_diag['nodes'][nn]['edge_ids'])
        if right_diag['nodes'][nn]['vertex_id'] == -1:
            right_diag['nodes'][nn]['vertex_id'] = -2
        elif right_diag['nodes'][nn]['vertex_id'] == -2:
            right_diag['nodes'][nn]['vertex_id'] = -1
        else:
            momTup = ()
            for var in right_diag['nodes'][nn]['momenta']:
                mom = ""
                cc = 0
                for char in var:
                    if char == '+':
                        mom = mom + '-'
                    elif char == '-':
                        mom = mom + '+'
                    else:
                        if cc == 0:
                            mom = mom + "-" + char
                        else:
                            mom = mom + char
                    cc += 1
                momTup = momTup + (mom,)
            right_diag['nodes'][nn]['momenta'] = momTup
        right_diag['nodes'][nn+node_offset] = right_diag['nodes'].pop(nn)
    return right_diag


def sew_amp_diags(ld, rd):
    """sews left- and right-diagram to a cut forward scattering diagram

    Args:
        left_diag ([dict]): [dictionary containing the left-diagram]
        right_diag ([dict]): [dictionary containing the right diagram]
    """
    
    sewed_graph = {}
    left_diag = copy.deepcopy(ld)
    right_diag = copy.deepcopy(rd)
    # get loop-momenta
    left_diag = append_loop_moms(left_diag)
    right_diag = append_loop_moms(right_diag)

    # impose momentum conservation
    left_diag = impose_mom_conservation(left_diag)
    right_diag = impose_mom_conservation(right_diag)

    # check if loop-momenta need relabeling (they are not different: only relevant for squaring)
    offset = max(len(left_diag['loop_momenta']),
                 len(left_diag['loop_momenta']))
    if len(left_diag['loop_momenta'])+len(right_diag['loop_momenta']) != len(set(left_diag['loop_momenta']+right_diag['loop_momenta'])):

        for i in range(1, len(right_diag['loop_momenta'])+1):
            lhs = right_diag['loop_momenta'][i-1]
            rhs = 'k' + str(i+offset)
            for ee in right_diag['edges'].values():
                ee['momentum'] = ee['momentum'].replace(lhs, rhs)
            for nn in right_diag['nodes'].values():
                if isinstance(nn['momenta'], str):
                    nn['momenta'] = nn['momenta'].replace(lhs, rhs)
                else:
                    nn['momenta'] = tuple(x.replace(lhs, rhs)
                                          for x in nn['momenta'])
            right_diag['loop_momenta'] = tuple(
                x.replace(lhs, rhs) for x in right_diag['loop_momenta'])
            right_diag['analytic_num'] = right_diag['analytic_num'].replace(
                lhs, rhs)
            full_offset = offset+i
    else:
        full_offset = offset

    # fill external edges
    sewed_graph['edges'] = {}
    for ee in list(left_diag['edges']):
        if left_diag['edges'][ee]['type'] == 'in':
            sewed_graph['edges'][ee] = copy.deepcopy(left_diag['edges'][ee])
    for ee in list(right_diag['edges']):
        if right_diag['edges'][ee]['type'] == 'out':
            sewed_graph['edges'][ee] = copy.deepcopy(right_diag['edges'][ee])
    # fill vertual uncut edges
    for ee in list(left_diag['edges']):
        if left_diag['edges'][ee]['type'] == 'virtual':
            sewed_graph['edges'][ee] = copy.deepcopy(left_diag['edges'][ee])
    for ee in list(right_diag['edges']):
        if right_diag['edges'][ee]['type'] == 'virtual':
            sewed_graph['edges'][ee] = copy.deepcopy(right_diag['edges'][ee])
    # fill nodes
    sewed_graph['nodes'] = {}
    for nl in left_diag['nodes']:
        if left_diag['nodes'][nl]['vertex_id'] != -2:
            sewed_graph['nodes'][nl] = copy.deepcopy(left_diag['nodes'][nl])
    for nr in right_diag['nodes']:
        if right_diag['nodes'][nr]['vertex_id'] != -1:
            sewed_graph['nodes'][nr] = copy.deepcopy(right_diag['nodes'][nr])
    # fill cut-edges
    cc = 1
    for el in list(left_diag['edges']):
        if left_diag['edges'][el]['type'] == 'out':
            sewed = False

            for er in list(right_diag['edges']):

                if right_diag['edges'][er]['type'] == 'in' and right_diag['edges'][er]['momentum'] == left_diag['edges'][el]['momentum']:
                    key = copy.copy(el)
                    ec = {}
                    ec['type'] = "virtual"
                    ec['name'] = "c"+str(cc)
                    cc += 1
                    ec['momentum'] = left_diag['edges'][el]['momentum']
                    ec['vertices'] = (
                        left_diag['edges'][el]['vertices'][0], right_diag['edges'][er]['vertices'][1])
                    ec['PDG'] = copy.copy(left_diag['edges'][el]['PDG'])
                    ec['indices'] = copy.copy((left_diag['edges'][el]).get('indices', []))
                    ec['is_cut'] = True
                    full_edge=copy.deepcopy({key: copy.deepcopy(ec)})
                    sewed_graph['edges'].update(full_edge)
                    # relabel the edge_ids
                    for nn in sewed_graph['nodes'].values():
                        # test =False
                        # if er in nn['edge_ids']:
                        #     test =True
                        # if test:
                        #     print((nn['edge_ids'],er,el))
                        nn['edge_ids']=tuple(
                            map(lambda x: el if x == er else x, nn['edge_ids']))
                        # if test:
                        #     print((nn['edge_ids'],er,el))
                        #     print("************")
                    sewed=True
                    break
            if sewed == False:
                sys.exit("Could not sew graphs! Abort")

    sewed_graph['analytic_num']="(" + \
        left_diag['analytic_num']+")*("+right_diag['analytic_num']+")"
    sewed_graph['overall_factor']='1'
    # declare formally external momenta to now beeing loop-momenta
    for n, mom in enumerate(left_diag['out_momenta'][:-1]):
        rhs='k'+str(full_offset+n+1)
        for ee in sewed_graph['edges'].values():
            ee['momentum']=ee['momentum'].replace(mom, rhs)
        for nn in sewed_graph['nodes'].values():
            if isinstance(nn['momenta'], str):
                nn['momenta']=nn['momenta'].replace(mom, rhs)
            else:
                nn['momenta']=tuple(x.replace(mom, rhs)
                                      for x in nn['momenta'])
        sewed_graph['analytic_num']=sewed_graph['analytic_num'].replace(
            mom, rhs)

    
    # we need continous vertex labels for igraph
    old_vert=list(set(sewed_graph['nodes'].keys()))
    new_vert=["n"+str(i) for i in range(len(old_vert))]
    
    
    for i, vnew in enumerate(new_vert):
        for ee in sewed_graph['edges'].values():
            ee['vertices']=tuple(
                map(lambda x: vnew if x == old_vert[i] else x, ee['vertices']))
    for i, nn in enumerate(old_vert):
        sewed_graph['nodes'][new_vert[i]
                             ]=sewed_graph['nodes'].pop(nn)
    
    for i, vnew in enumerate(new_vert):
        for ee in sewed_graph['edges'].values():
            ee['vertices']=tuple(
                map(lambda x: i+1 if x == vnew else x, ee['vertices']))
    for i, vnew in enumerate(list(sewed_graph['nodes'])):
        sewed_graph['nodes'][i+1]=sewed_graph['nodes'].pop(vnew)


    # rename edges (ask ben why)
    # look at edge-map-lin
    # look at edge-name-map
    propCount=1
    in_count=0
    for ee in sewed_graph['edges'].values():
        if ee['type'] == 'virtual':
            ee['name']='p' + str(propCount)
            propCount += 1
        elif ee['type'] == 'in':
            ee['name']='q'+ee["momentum"][1:]
            in_count += 1
    for ee in sewed_graph['edges'].values():
        if ee['type'] == 'out':
            ee['name']='q'+str(int(ee["momentum"][1:])+in_count)
    

    #SET INDICES
    cc=1
    for eedict in sewed_graph['edges'].values():
        if eedict['type'] == 'virtual':
            eedict['indices'] = (cc,cc+1)
            cc+=2
        else:
            eedict['indices'] = (cc,)
            cc+=1

    # (nodes): assign indices and PDGs
    # delete indices in nodes
    for knn, nn in sewed_graph['nodes'].items():
        for key in list(nn):
            if key == 'indices' or key == 'PDGs':
                del nn[key]
            if key == 'fake_vertex':
                sewed_graph['effective_vertex_id']=knn
                del nn[key]
    # # assign new indices
    for nn , valdic in sewed_graph['nodes'].items():
        for eid in valdic['edge_ids']:
            
            vv = sewed_graph['edges'][eid]['vertices']
            idxs = sewed_graph['edges'][eid]['indices']
            
            if nn in vv:
                if vv[0] == nn:
                    valdic['indices'] = valdic.get('indices',()) + (idxs[0],)
                else:
                    valdic['indices'] = valdic.get('indices',())+ (idxs[-1],)
    sewed_graph["diag_set"] = left_diag.get("diag_set")
    sewed_graph["index_shift"] = left_diag.get("index_shift",0)
    return copy.deepcopy(sewed_graph)


def save_sg_dict_to_file(dictList, fileName='dict_dump'):
    pjoin=os.path.join
    root_path=os.path.dirname(os.path.realpath(__file__))

    f=open(pjoin(root_path, fileName+'.py'), 'w')
    f.write('graphs=[]\n')
    f.write('graph_names=[]\n')
    for cc, dic in enumerate(dictList):
        f.write('graph_names+=["SG_QG'+str(cc)+'"]\n')
        f.write('graphs.append(')
        f.write(str(dic))
        f.write(')\n')
    f.close()


class amplitude():
    """ import and manipulation of amplitudes
    """
    def __init__(self, runcard, output, src_dir):
        self.runcard=runcard
        self.output_dir=output
        self.dirs=['FORM', 'TEMPDIR', 'SGS_LIST', 'FORM/workspace']
        self.initialized=False
        self.root_dir=''
        self.super_graphs=intefered_sgs_list()
        self.form_files=["ChisholmIdentities.prc",
            "coltracebased.h", "amp_numerator.frm", "color_basis.frm"]
        self.src_dir=src_dir
        self.name=""
        self.from_math=True
        self.interfence_type=None
        self.math_sg="./"
        self.perform_color=True
        self.color_per_graph=False
        self.external_data={}

    def create_dir(self, path):
        try:
                os.mkdir(path)
        except OSError:
                print("Creation of the directory %s failed" % path)
        else:
                print("Successfully created the directory %s " % path)


    def set_up_working_dir(self, abs_path_do_working_dir, clean_up=True):
        self.root_dir=abs_path_do_working_dir
        Path(self.root_dir).mkdir(parents=True, exist_ok=True)

        Path(os.path.join(self.root_dir, 'Cards')).mkdir(
            parents=True, exist_ok=True)
        # TODO add param_card in Source

        if not os.path.exists(self.root_dir):
            print("create root directory and add templates")
            self.create_dir(self.root_dir)
        # create directory tree
        for ddir in self.dirs:
            if clean_up is True and os.path.exists(os.path.join(self.root_dir, ddir)):
                shutil.rmtree(os.path.join(self.root_dir, ddir))
            Path(os.path.join(self.root_dir, ddir)).mkdir(
                parents=True, exist_ok=True)
        # copy form files
        self.form_workspace=os.path.join(self.root_dir, 'FORM', 'workspace')
        for ff in self.form_files:
            if not os.path.exists(os.path.join(self.form_workspace, ff)):
                shutil.copy(os.path.join(self.src_dir, ff),
                            os.path.join(self.form_workspace, ff))








    def create_intefered_amp(self):
        # here we need to include the python qgraf-output
        if self.from_math:
            mod_name='math_supergraph'

            spec=importlib.util.spec_from_file_location(mod_name, self.math_sg)
            math_supergraph=importlib.util.module_from_spec(spec)
            spec.loader.exec_module(math_supergraph)
            # BUILD SUPEGRPAPHS
            left_diags=math_supergraph.graphs
            right_diags=[]
            for i,g in enumerate(left_diags):
                if 'analytic_num' not in g.keys():
                    g['analytic_num']="1"
                if 'diag_set' not in g.keys():
                    g['diag_set'] = "diag_set_xxx" + str(i)
                if self.interfence_type == "effective_vertex":
                    right_diags += [to_right_diag(g, to_effective_vertex=True)]
                else:
                    right_diags += [to_right_diag(g,
                                                  to_effective_vertex=False)]


            super_graphs=[]
            if self.interfence_type == "effective_vertex":
                for cc, ld in enumerate(left_diags):
                    super_graphs += [sew_amp_diags(ld, right_diags[cc])]
            else:
                for ld in left_diags:
                    for rd in right_diags:
                        super_graphs += [sew_amp_diags(ld, rd)]


            self.math_sg=super_graphs

    def def_process(self, proc_dict):
            if isinstance(proc_dict, dict):

                self.name=(proc_dict['process_specification']).pop(
                    'name', 'new_amplitude')
                self.from_math=(proc_dict['process_specification']).pop(
                    'generate_from_math', True)
                self.interfence_type=(proc_dict['process_specification']).pop(
                    'interference_type', "effective_vertex")
                self.math_sg=(proc_dict['process_specification']).pop(
                    'math_sg', './')
                self.perform_color=(proc_dict['process_specification']).pop(
                    'perform_color', True)
                self.color_per_graph=(proc_dict['process_specification']).pop(
                    'color_per_graph', False)
                self.external_data=(proc_dict['process_specification']).get(
                    'external_data', {})

            else:
                sys.exit("type has to be a dictionary")
            if self.interfence_type != "effective_vertex":
                sys.exit(
                    "So far only the interference type `effective_vertex` is supported within the LTD-code")





    def initialize_process(self, clean_up=True):
        path_to_runcard=self.runcard
        abs_path_do_working_dir=os.path.abspath(self.output_dir)
        if not self.initialized:
            if abs_path_do_working_dir == '' and self.root_dir == '':

                dir_yaml=os.path.dirname(path_to_runcard)
                self.set_up_working_dir(os.path.join(
                    dir_yaml, "process_dir"), clean_up=clean_up)
                with open(path_to_runcard) as f:
                    self.def_process(yaml.load(f, Loader=yaml.FullLoader))
            elif abs_path_do_working_dir != '' and self.root_dir == '':
                self.set_up_working_dir(os.path.join(
                    abs_path_do_working_dir, "process_dir"), clean_up=clean_up)
                with open(path_to_runcard) as f:
                    self.def_process(yaml.load(f, Loader=yaml.FullLoader))

            self.create_intefered_amp()
            self.initialized=True
            save_sg_dict_to_file(self.math_sg, os.path.join(
                self.root_dir, 'TEMPDIR', self.name+'_sgs_intefered'))
            out_dicts=[]

            _MANDATORY_EXTERNAL_DATA=['in_momenta', 'out_momenta', 'spinor_v',
                'spinor_vbar', 'spinor_u', 'spinor_ubar', 'pol',  'cpol',  'n_in', 'n_out']
            for entry in _MANDATORY_EXTERNAL_DATA:
                if entry not in self.external_data.keys():
                    print("Missing entry in external_data:", entry)
                    sys.exit("Extend the run-card.yaml")
            for entry in _MANDATORY_EXTERNAL_DATA[2:-2]:
                if len(self.external_data.get(entry)) > 0:
                    for elem in self.external_data.get(entry):
                        for el in elem:
                            if len(el) < 2:
                                print(
                                    entry, "is supposed to be complex. Its components are of the form [Re, Im] ")
                                sys.exit("Extend the run-card.yaml")



            # perform color decomposition
            if self.perform_color:
                self.super_graphs=intefered_sgs_list(self.math_sg)
                self.super_graphs.perform_color_decomp(
                    self.form_workspace, CPERGRAPH=self.color_per_graph)
                # dump to dict

                for i, gs in enumerate(self.super_graphs.sgs_list):
                    save_sg_dict_to_file(gs, os.path.join(
                        self.root_dir, 'SGS_LIST', self.name+'_color_struc'+str(i+1)))
                    out_dicts += [os.path.join(self.root_dir,
                                               'SGS_LIST', self.name+'_color_struc'+str(i+1))]
        return out_dicts
############################################################################################################

class intefered_sgs_list():
    """ Inteferred supergraphs
    """
    def __init__(self, sgs_list=[{}]):
        self.sgs_list=sgs_list
        self.color_decomp=False

    def generate_numerator_functions(self, workspace, graph_dict, SIGID, colpergraph):
        """ Use form to perform the color separations"""

        FORM_vars={}
        pjoin=os.path.join
        n_incoming=sum(
            [1 for edge in (graph_dict['edges']).values() if edge['type'] == 'in'])
        FORM_vars['NINITIALMOMENTA']=n_incoming
        n_loops=len((graph_dict['edges']).keys()) - \
            len((graph_dict['nodes']).keys()) + 1
        FORM_vars['NFINALMOMENTA']=n_loops
        FORM_vars['SGID']=SIGID
        FORM_vars['CPERGRAPH']=colpergraph
        FORM_vars['INDSHIFT'] = graph_dict.get("index_shift",0)
        i_graph=int(FORM_vars['SGID'])

        form_input=graph_dict['analytic_num']
        # verify input
        if form_input.find(".")!=-1:
            sys.exit("ERROR IN NUMERATOR: The numerator contains . , which is not allowed. Write scalar-products as sp(mom1+mom2+...,momN+...)")

        
        for cc in re.findall("x+[A-Za-z]",form_input):
            if not cc[1].isupper():
                err = "ERROR IN NUMERATOR: "+cc+" has to have an capital letter as the second character"
                sys.exit(err)

        
        selected_workspace=workspace
        FORM_source=pjoin(selected_workspace, 'color_basis.frm')

        with open(pjoin(selected_workspace, 'input_%d.h' % i_graph), 'w') as f:
            f.write('L F = {};'.format(form_input))
        print(' '.join([
            'form',
        ] +
            ['-D %s=%s' % (k, v) for k, v in FORM_vars.items()] +
            [FORM_source, ]
        ))

        r=subprocess.run(' '.join([
            'form',
        ] +
            ['-D %s=%s' % (k, v) for k, v in FORM_vars.items()] +
            [FORM_source, ]
        ),
            shell=True,
            cwd=selected_workspace,
            capture_output=True)
        if r.returncode != 0:
            raise "FORM processing failed with error:\n%s" % (
                r.stdout.decode('UTF-8'))


    def perform_color_decomp(self, FORM_WS, CPERGRAPH=0):
        from pathlib import Path
        pjoin=os.path.join
        all_color=[]

        for i, g in enumerate(self.sgs_list):
            if CPERGRAPH:
                self.generate_numerator_functions(FORM_WS, g, i, colpergraph=1)
            else:
                self.generate_numerator_functions(FORM_WS, g, i, colpergraph=0)
            f=open(pjoin(FORM_WS, 'SG_'+str(i)+'_color_decomp.txt'), 'r')
            color=copy.deepcopy(eval((f.read()).replace('\n', '')))
            all_color += list(color.keys())
            f.close()
            g['color_structures']=color
        all_color=list(dict.fromkeys(all_color))

        graphs_color_sep=[]
        for c in all_color:
            struc=[]
            for g in self.sgs_list:
                if c in (g['color_structures']).keys():
                    tmp_g=copy.deepcopy(g)
                    tmp_g['analytic_num']=copy.copy(g['color_structures'][c])
                    del tmp_g['color_structures']
                    tmp_g['color_struc']=c
                    struc += [tmp_g]
            graphs_color_sep += [struc]
        self.sgs_list=graphs_color_sep
        self.color_decomp=True
