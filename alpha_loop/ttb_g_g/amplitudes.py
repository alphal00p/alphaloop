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
        for ee in list((right_diag['edges']).keys()):
            if right_diag['edges'][ee]['type'] == 'virtual':
                del right_diag['edges'][ee]
            else:
                ec += 1
        for nn in list((right_diag['nodes']).keys()):
            if right_diag['nodes'][nn]['vertex_id'] == 0:
                del right_diag['nodes'][nn]
            else:
                nc += 1
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
        # add internal vertex
        right_diag['nodes'][nc+1] = {'edge_ids': tuple(x for x in list((right_diag['edges']).keys())), 'vertex_id': 0,
                                     'momenta': tuple(map(lambda x: right_diag['edges'][x]['momentum'] if right_diag['edges'][x]['type'] == 'in' else '-'+str(right_diag['edges'][x]['momentum']), list((right_diag['edges']).keys())))}

        for ee in (right_diag['edges']).values():
            ee['vertices'] = tuple(
                map(lambda x: nc+1 if x > nc else x, list((ee['vertices']))))

    # complex conjugation. will be perforemed in form
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

    # check if loop-momenta need relabeling (they are not different)
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
                    ec['indices'] = copy.copy(
                        left_diag['edges'][el]['indices'])
                    full_edge = copy.deepcopy({key: copy.deepcopy(ec)})
                    sewed_graph['edges'].update(full_edge)
                    # relabel the edge_ids
                    for nn in sewed_graph['nodes'].values():
                        # test =False
                        # if er in nn['edge_ids']:
                        #     test =True
                        # if test:
                        #     print((nn['edge_ids'],er,el))
                        nn['edge_ids'] = tuple(
                            map(lambda x: el if x == er else x, nn['edge_ids']))
                        # if test:
                        #     print((nn['edge_ids'],er,el))
                        #     print("************")
                    sewed = True
                    break
            if sewed == False:
                sys.exit("Could not sew graphs! Abort")

    sewed_graph['analytic_num'] = "(" + \
        left_diag['analytic_num']+")*("+right_diag['analytic_num']+")"
    sewed_graph['overall_factor'] = '1'
    # declare formally external momenta to now beeing loop-momenta
    for n, mom in enumerate(left_diag['out_momenta'][:-1]):
        rhs = 'k'+str(full_offset+n+1)
        for ee in sewed_graph['edges'].values():
            ee['momentum'] = ee['momentum'].replace(mom, rhs)
        for nn in sewed_graph['nodes'].values():
            if isinstance(nn['momenta'], str):
                nn['momenta'] = nn['momenta'].replace(mom, rhs)
            else:
                nn['momenta'] = tuple(x.replace(mom, rhs)
                                      for x in nn['momenta'])
        sewed_graph['analytic_num'] = sewed_graph['analytic_num'].replace(
            mom, rhs)

    # we need continous vertex labels for igraph
    new_vert = sorted(list(range(1, len(sewed_graph['nodes'])+1)))
    # old keys will always be larger than new keys
    old_vert = sorted(list(sewed_graph['nodes'].keys()))
    for i, vnew in enumerate(new_vert):
        for ee in sewed_graph['edges'].values():
            ee['vertices'] = tuple(
                map(lambda x: vnew if x == old_vert[i] else x, ee['vertices']))
    for i, nn in enumerate(list(sewed_graph['nodes'])):
        sewed_graph['nodes'][new_vert[i]
                             ] = sewed_graph['nodes'].pop(old_vert[i])

    # rename edges (ask ben why)
    # look at edge-map-lin
    # look at edge-name-map
    propCount = 1
    extCount = 1
    for ee in sewed_graph['edges'].values():
        if ee['type'] == 'virtual':
            ee['name'] = 'q' + str(propCount)
            propCount += 1
        elif ee['type'] == 'in':
            ee['name'] = 'p'+str(extCount)
            extCount += 1
    for ee in sewed_graph['edges'].values():
        if ee['type'] == 'out':
            ee['name'] = 'p'+str(extCount)
            extCount += 1
    return copy.deepcopy(sewed_graph)


def save_dict_to_file(dictList, fileName='intefered_diags'):
    pjoin = os.path.join
    root_path = os.path.dirname(os.path.realpath(__file__))

    f = open(pjoin(root_path, fileName+'.py'), 'w')
    f.write('graphs=[]\n')
    f.write('graph_names=[]\n')
    for cc, dic in enumerate(dictList):
        f.write('graph_names+=["SG_QG'+str(cc)+'"]\n')
        f.write('graphs.append(')
        f.write(str(dic))
        f.write(')\n')
    f.close()


def generate_numerator_functions(workspace, plugin_path, graph_dict, SIGID):
    """ Use form to perform the color separations"""

    FORM_vars = {}
    pjoin = os.path.join
    n_incoming = sum(
        [1 for edge in (graph_dict['edges']).values() if edge['type'] == 'in'])
    FORM_vars['NINITIALMOMENTA'] = n_incoming
    n_loops = len((graph_dict['edges']).keys()) - \
        len((graph_dict['nodes']).keys()) + 1
    FORM_vars['NFINALMOMENTA'] = n_loops
    FORM_vars['SGID'] = SIGID
    i_graph = int(FORM_vars['SGID'])

    form_input = graph_dict['analytic_num']

    selected_workspace = workspace
    shutil.copy(pjoin(plugin_path, "numerator.frm"),
                pjoin(selected_workspace, 'numerator.frm'))
    shutil.copy(pjoin(plugin_path, "diacolor.h"),
                pjoin(selected_workspace, 'diacolor.h'))
    shutil.copy(pjoin(plugin_path, "color_basis.frm"),
                pjoin(selected_workspace, 'color_basis.frm'))
    shutil.copy(pjoin(plugin_path, "amp_numerator.frm"),
                pjoin(selected_workspace, 'amp_numerator.frm'))
    shutil.copy(pjoin(plugin_path, "ChisholmIdentities.prc"),
                pjoin(selected_workspace, "ChisholmIdentities.prc"))

    FORM_source = pjoin(selected_workspace, 'color_basis.frm')

    with open(pjoin(selected_workspace, 'input_%d.h' % i_graph), 'w') as f:
        f.write('L F = {};'.format(form_input))

    r = subprocess.run(' '.join([
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


def perform_color_decomp(graphs, FORM_WS, plugin_path):
    from pathlib import Path
    pjoin = os.path.join
    all_color =[]
    for i, g in enumerate(graphs):
        generate_numerator_functions(FORM_WS, plugin_path, g, i)
        f = open(pjoin(FORM_WS, 'SG_'+str(i)+'_color_decomp.txt'), 'r')
        color = copy.deepcopy(eval((f.read()).replace('\n', '')))
        all_color += list(color.keys())
        f.close()
        g['color_structures'] = color
    all_color = list(dict.fromkeys(all_color))
    graphs_color_sep = []
    for c in all_color:
        struc = []
        for g in graphs:
            if c in (g['color_structures']).keys():
                tmp_g = copy.deepcopy(g)
                tmp_g['analytic_num'] = copy.copy(g['color_structures'][c])
                del tmp_g['color_structures']
                tmp_g['color_struc'] = c
                struc += [tmp_g]
        graphs_color_sep += [struc]
    return graphs_color_sep
