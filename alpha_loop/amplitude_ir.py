#!/usr/bin/env python3
import sys
import copy

def to_right_diag(left_diag):
    """Translates left-diagram to right-diagram (complex-conjugation + flow reversal)

    Args:
        left_diag ([dict]): [left-diagram]
    """
    def rev_tup(tuples): 
        new_tup = tuples[::-1] 
        return new_tup

    right_diag=left_diag
    # complex conjugation. will be perforemed in form
    right_diag['analytic_num']="conjugate("+right_diag['analytic_num']+")"

    
    node_offset = len(left_diag['nodes'])

    # reverse edge orientations but keep momentum flow
    for ee in list(right_diag['edges']):
        right_diag['edges'][ee]['vertices']=tuple(x+node_offset for x in rev_tup(right_diag['edges'][ee]['vertices']))
        if right_diag['edges'][ee]['type']=='in':
            right_diag['edges'][ee]['type']='out'
        elif right_diag['edges'][ee]['type']=='out':
            right_diag['edges'][ee]['type']='in'
        newEdgeKey= tuple(x+node_offset for x in rev_tup(ee[:2]))+(-ee[2],)
        right_diag['edges'][ee]['name']=right_diag['edges'][ee]['name']+'r'
        right_diag['edges'][newEdgeKey]=right_diag['edges'].pop(ee)

    for nn in list(right_diag['nodes']):
        if right_diag['nodes'][nn]['vertex_id']==-1:
            right_diag['nodes'][nn]['vertex_id']=-2
        elif right_diag['nodes'][nn]['vertex_id']==-2:
            right_diag['nodes'][nn]['vertex_id']=-1
        else:
            momTup =()
            for var in right_diag['nodes'][nn]['momenta']:
                mom=""
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
                    cc+=1
                momTup = momTup + (mom,)  
            right_diag['nodes'][nn]['momenta']=momTup
        right_diag['nodes'][nn+node_offset]=right_diag['nodes'].pop(nn)       
    return right_diag

def sew_amp_diags(left_diag,right_diag):
    """sews left- and right-diagram to a cut forward scattering diagram

    Args:
        left_diag ([dict]): [dictionary containing the left-diagram]
        right_diag ([dict]): [dictionary containing the right diagram]
    """
    sewed_graph ={}
    sewed_graph["left_diagram"]=left_diag
    sewed_graph["right_diagram"]=right_diag
    sewed_graph["squared_diagram"]={}
    
    # fill external edges
    sewed_graph['squared_diagram']['edges']={}
    for ee in list(left_diag['edges']):
        if left_diag['edges'][ee]['type'] == 'in':
            sewed_graph['squared_diagram']['edges'][ee]=left_diag['edges'][ee]
    for ee in list(right_diag['edges']):
        if right_diag['edges'][ee]['type'] == 'out':
            sewed_graph['squared_diagram']['edges'][ee]=right_diag['edges'][ee]
    # fill vertual uncut edges
    for ee in list(left_diag['edges']):
        if left_diag['edges'][ee]['type'] == 'virtual':
            sewed_graph['squared_diagram']['edges'][ee]=left_diag['edges'][ee]
    for ee in list(right_diag['edges']):
        if right_diag['edges'][ee]['type'] == 'virtual':
            sewed_graph['squared_diagram']['edges'][ee]=right_diag['edges'][ee]
    # fill cut-edges
    cc=1
    for el in list(left_diag['edges']):
        if left_diag['edges'][el]['type'] == 'out':
            sewed =False
            for er in list(right_diag['edges']):
                if right_diag['edges'][er]['type'] == 'in' and right_diag['edges'][er]['momentum']==left_diag['edges'][el]['momentum']:
                    ec=copy.copy((el[0],er[1],el[2]))
                    ec = {}
                    ec['type']="cut"
                    ec['name']="c"+str(cc)
                    cc+=1
                    ec['momentum']=left_diag['edges'][el]['momentum']
                    ec['vertices']=(el[0],er[1])
                    sewed_graph['squared_diagram']['edges'].update(ec)
                    sewed=True
                    break
            if sewed==False:
                  sys.exit("Could not sew graphs! Abort")
    sewed_graph["squared_diagram"]['analytic_num']="("+left_diag['analytic_num']+")*("+right_diag['analytic_num']+")"
    return sewed_graph

import pandas as pd
sys.path.append('/home/armin/my_programs/pynloop/alpha_loop/test_process')
import qqbaaNLO
mygraphs=qqbaaNLO.graphs
for g in mygraphs:
    if 'analytic_num' not in g.keys():
        g['analytic_num'] = "1"
    # indices and PDGs on edges and vertices have no meaning for the scalar integrals we consider
    # lets drop them, before they are misintepreted
    for ee in list(g['edges']):
        del g['edges'][ee]['indices']
        del g['edges'][ee]['PDG']
    for nn in list(g['nodes']):
        del g['nodes'][nn]['indices']
        del g['nodes'][nn]['PDGs']


ld=mygraphs[0]
rd=copy.deepcopy(ld)
rd = to_right_diag(rd)

print(ld['edges'])
print("\r\n")
print(rd['edges'])
print("\r\n")

gg=sew_amp_diags(ld,rd)
print(gg['squared_diagram']['edges'])
