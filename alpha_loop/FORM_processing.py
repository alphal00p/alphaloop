#!/usr/bin/env python3

import copy
import logging
import os
from pathlib import Path
from pprint import pprint, pformat
import math
import igraph
import time
import numpy as np
from itertools import combinations_with_replacement
from collections import OrderedDict

import progressbar
from itertools import chain, product
import sys
import subprocess
import argparse
import shutil
pjoin = os.path.join

if __name__ == "__main__":
    root_path = os.path.dirname(os.path.realpath( __file__ ))
    sys.path.insert(0, pjoin(root_path,os.path.pardir))
    sys.path.insert(0, pjoin(root_path,os.path.pardir,os.path.pardir,os.path.pardir))
    if 'MG5DIR' in os.environ:
        sys.path.insert(0, os.environ['MG5DIR'])
    else:
        print("\033[91mYou are using ./FORM_processing.py in standalone, it is recommended then "+
              "that you define the environment variable 'MG5DIR' pointing to the root directory of MG5aMC in your system.\033[0m")
        sys.exit(1)

import alpha_loop.utils as utils
import re

import madgraph.core.base_objects as base_objects
import madgraph.various.misc as misc
import madgraph.iolibs.file_writers as writers
from madgraph import MadGraph5Error, InvalidCmd, MG5DIR
import models.model_reader as model_reader
import multiprocessing

import LTD.squared_topologies
import LTD.ltd_utils
import LTD.partial_fractioning

logger = logging.getLogger('alphaLoop.FORM_processing')

if __name__ == "__main__":
    logging.basicConfig()
    logger.setLevel(logging.INFO)

plugin_path = os.path.dirname(os.path.realpath( __file__ ))


FORM_processing_options = {
    'FORM_path': 'form', 
    'tFORM_path': 'tform', 
    # Define the extra aguments for the compilation
    'compilation-options': [],
    'cores': multiprocessing.cpu_count(), 
    'extra-options': {'OPTIMITERATIONS': 1000},
    # If None, only consider the LMB originally chosen.
    # If positive and equal to N, consider the first N LMB from the list of LMB automatically generated
    # If negative consider all possible LMBs.
    'number_of_lmbs' : None,
    # If None, the reference LMB will be the one originally chosen.
    # If positive and equal to N, the Nth LMB will be used for the reference implementation of the supergraph.
    'reference_lmb' : None,
    'FORM_call_sig_id_offset_for_additional_lmb' : 1000000,
    'generate_integrated_UV_CTs' : True,
    'generate_renormalisation_graphs' : True,
    'UV_min_dod_to_subtract' : 0,
    'selected_epsilon_UV_order' : 0,
    # Select how to include the finite part of the renormalisation. Possible values are:
    # a) 'together' : all contributions are ketp.
    # b) 'only' : only the contribution from the finite part of the renormalisation is kept.
    # c) 'removed' : the finitie part of the renormalisation is removed. 
    'renormalisation_finite_terms' : 'together'
}

# Can switch to tmpdir() if necessary at some point
FORM_workspace = pjoin(plugin_path,'FORM_workspace')
Path(FORM_workspace).mkdir(parents=True, exist_ok=True)
resources_to_link = ['diacolor.h']
for resource in resources_to_link:
    if not os.path.exists(pjoin(FORM_workspace,resource)):
        utils.ln(pjoin(plugin_path,resource),starting_dir=FORM_workspace)

class FormProcessingError(MadGraph5Error):
    """ Error for the FORM processing phase."""
    pass

class FORMSuperGraph(object):
    """ Simplified SuperGraph object with only the necessary information to have it processed by FORM."""

    _FORM_Feynman_rules_conventions = {
        # SSS
        ( 1, 1, 1 ): (0, 1, 2),
        # SSSS
        ( 1, 1, 1 ): (0, 1, 2),
        # GGG
        ( 3, 3, 3 ): (0, 1, 2),
        # GGGG
        ( 3, 3, 3, 3 ): (0, 1, 2, 3),
        # FxFV
        ( -2, 2, 3 ): (0, 2, 1),
        # GHxGHV
        ( -1, 1, 3 ): (0, 2, 1),
        # FxFS
        ( -2, 1, 2 ): (0, 1, 2),
        # SSS
        ( 1, 1, 1 ): (0, 1, 2),
        # SSSS
        ( 1, 1, 1, 1 ): (0, 1, 2),
        # VxVV (e.g. W+ W- a )
        ( -3, 3, 3 ): (1, 0, 2),
        # Fx F
        (-2, 2): (0, 1),
        # V V
        (-3,3): (0,1),
        (3,3): (0,1),
        # S S
        (1,1): (0,1),
        # S V V
        (1,3,3): (0,1,2),
        # S V V V
        (1,3,3,3): (0,1,2,3),
        # S V V V V
        (1,3,3,3,3): (0,1,2,3,4),
        # S S V V
        (1,1,3,3): (0,1,2,3),
        # S S V V V
        (1,1,3,3,3): (0,1,2,3,4),
        # S S S V V
        (1,1,1,3,3): (0,1,2,3,4),
        # S S S V V V
        (1,1,1,3,3,3): (0,1,2,3,4,5)
    }

    _include_momentum_routing_in_rendering=False
    _include_edge_name_in_rendering=True
    _rendering_size = (1.0*(11.0*60),1.0*(8.5*60)) # 1.0 prefactor should be about 1 landscape A4 format per graph
    # Choose graph layout strategy. Interesting options are in comment.
    _graph_layout_strategy = '{"PackingLayout"->"ClosestPacking"}' 
    #_graph_layout_strategy = 'GraphLayout -> "SpringElectricalEmbedding"' 
    #_graph_layout_strategy = 'GraphLayout -> "SpringEmbedding"' 
    #_graph_layout_strategy = 'GraphLayout -> {"LayeredEmbedding", "Orientation" -> Left, "RootVertex" -> "I1"}' 
    #_graph_layout_strategy = 'GraphLayout -> {"LayeredDigraphEmbedding", "Orientation" -> Left}'
    
    # '{"SpringEmbedding"}' gives interesting results too.

    def __init__(self, *args,
        call_identifier=None,
        name=None,
        edges=None,
        nodes=None,
        overall_factor="1",
        multiplicity=1,
        benchmark_result=0.0,
    ):
        """ initialize a FORM SuperGraph from several options."""

        self.is_zero = False
        self.edges = edges
        self.nodes = nodes
        self.overall_factor = overall_factor
        self.multiplicity = multiplicity
        # Give the possibility of specifying a benchmark result to the output yaml file for future ref.
        self.benchmark_result = benchmark_result
        # A hashable call signature
        self.call_identifier = call_identifier
        if name is None:
            self.name = str(self.call_identifier)
        else:
            self.name = name

        self.squared_topology = None
        self.replacement_rules = None
        self.integrand_info = {}

        # Store copies of self for different choices of LMB to be used for cross-check.
        self.additional_lmbs = []

    def filter_valid_cuts(self, cuts):
        """ 
        filter graph base on a list of allowed cuts whose entires are defined as:
            ([allowd pdg], n) : to ensure to cut "n" times particles contained 
                                in the pdg list
            ('any', n)        : cut "n" edges that could be anything 
                                (usefull for extra real radiations)
        example: e+e- > ggh @NLO => [([21],2), ([25],1), ('any',1)]
        """
        g = igraph.Graph()
        g.add_vertices(len(self.nodes))
        undirected_edges = set()

        for e in self.edges.values():
            undirected_edges.add(tuple(sorted(e['vertices'])))

        cut_edges = [[] for _ in range(len(cuts))]
        incoming_vertices = []
        outgoing_vertices = []
        edge_colors = {tuple(sorted([x - 1 for x in ue])) : [] for ue in undirected_edges}
        for ue in undirected_edges:
            multiple = 0
            for e in self.edges.values():
                if tuple(sorted(e['vertices'])) == ue:
                    multiple += 1
                    if e['type'] == 'virtual':
                        for ci, c in enumerate(cuts):
                            if c[0] == 'any':
                                cut_edges[ci] += [tuple(sorted([x - 1 for x in ue]))]
                            else:
                                if abs(e['PDG']) in c[0]:
                                    cut_edges[ci] += [tuple(sorted([x - 1 for x in ue]))]
                        edge_colors[tuple(sorted([x - 1 for x in ue]))] += [abs(e['PDG'])]
                    elif e['type'] == 'in':
                        incoming_vertices += [min(e['vertices'])-1]
                    elif e['type'] == 'out':
                        outgoing_vertices += [min(e['vertices'])-1]
            g.add_edges([tuple(sorted([x - 1 for x in ue]))]*multiple)

        take_cuts = []
        cut_colors = []
        n_optional = 0
        valid_cut = False
        for ci, cut in enumerate(cuts):
            cut_colors += [cut[0]] * cut[1]
            if cut[0] == 'any':
                n_optional += cut[1]
                take_cuts += [cut_edges[ci]+[None]] * cut[1]
            else:
                take_cuts += [cut_edges[ci]] * cut[1]
        #print(take_cuts)
        #print(cut_colors)
        #print(incoming_vertices)
        #print(outgoing_vertices)
        invalid_cuts = []
        count_checks = 0
        #print(cuts)
        for cut_edges in product(*take_cuts):
            # When a valid is cut we know we need to keep this graph
            if valid_cut:
                break
            
            # Check if cut has to be dropped based on previous failed attempts
            if any(all(cut_edges.count(c) >= veto_c.count(c) for c in set(veto_c)) for veto_c in invalid_cuts):
                continue
            count_checks += 1
            
            # Check valid color cuts
            ec = copy.deepcopy(edge_colors)
            for ce, cc in zip(cut_edges,cut_colors):
                if cc == 'any':
                    continue
                for color in cc:
                    if color in ec[ce]:
                        ec[ce].remove(color)
                        break
                else:
                    break
            else:
                # Allow for Pure Virtual corrections
                virtual_loops = cut_edges.count(None)
                for _ in range(virtual_loops):
                    cut_edges = list(cut_edges)
                    cut_edges.remove(None)
                    cut_edges = tuple(cut_edges)

                # Apply set of cuts
                gtmp = g.copy()
                for ci in range(len(cut_edges)):
                    if not gtmp.are_connected(*cut_edges[ci]):
                        # update invalid cuts
                        new_veto = cut_edges[:ci+1]
                        invalid_cuts = list(filter(lambda veto_c: not all(veto_c.count(c) >= new_veto.count(c) for c in set(new_veto)),invalid_cuts))
                        invalid_cuts += [new_veto]
                        break
                    gtmp.delete_edges(gtmp.get_eid(*cut_edges[ci]))
                    if not gtmp.is_connected():
                        if ci+1 < len(cut_edges):
                            if ci+1 < len(cut_edges) - n_optional + virtual_loops:
                                # update invalid cuts
                                new_veto = cut_edges[:ci+1]
                                invalid_cuts = list(filter(lambda veto_c: not all(veto_c.count(c) >= new_veto.count(c) for c in set(new_veto)),invalid_cuts))
                                invalid_cuts += [new_veto]
                            break
                        else:
                            # check that the vertices are correctly connected 
                            # to the left and right paths
                            # Meaning two vertices involved in a cut should belong do
                            # the opposite disconnected graphs
                            if any(len(gtmp.get_shortest_paths(c[0],to=c[1])[0])>0 for c in cut_edges):
                                break
                            # check that incoming and outgoing are on separate disconnected graphs
                            if any(len(gtmp.get_shortest_paths(incoming_vertices[0],to=c)[0])==0 for c in incoming_vertices[1:]) or\
                               any(len(gtmp.get_shortest_paths(incoming_vertices[0],to=c)[0])>0 for c in outgoing_vertices):
                                break
                            valid_cut = True
        return valid_cut


    def get_mathematica_rendering_code(self, model, FORM_id=None, lmb_id=None):
        """ Generate mathematica expression for drawing this graph."""

        repl_dict = {
            'edge_font_size'     : 10,
            'vetex_font_size'    : 10,
            'width'              : self._rendering_size[0],
            'height'             : self._rendering_size[1],
            'graph_layout_strategy' : self._graph_layout_strategy
        }
        if self.call_identifier and all(k in self.call_identifier for k in ['proc_id','left_diagram_id','right_diagram_id']):
            graph_name='MG: %s'%('P%(proc_id)dL%(left_diagram_id)dR%(right_diagram_id)d'%self.call_identifier)
        else:
            graph_name='MG: %s'%self.name
        if FORM_id is not None:
            graph_name +=' | FORM: #%d'%FORM_id
        if lmb_id is not None:
            graph_name +=' | LMB: #%d'%lmb_id
        repl_dict['graph_name'] = graph_name

        # Special name rendering rules
        def get_part_name(pdg):
            quark_names = {1:"d",2:"u",3:"s",4:"c",5:"b",6:"t"}
            if pdg==11:
                return r"\!\(\*SuperscriptBox[\(e\), \(+\)]\)"
            elif pdg==-11:
                return r"\!\(\*SuperscriptBox[\(e\), \(-\)]\)"
            elif pdg==22:
                return r"\[Gamma]"
            elif pdg in [-1,-2,-3,-4,-5,-6]:
                return r"\!\(\*OverscriptBox[\(%s\), \(_\)]\)"%quark_names[abs(pdg)]
            else:
                return model.get_particle(pdg).get_name()

        node_key_to_node_name = {}
        def get_node_label(node_key):
            if node_key in node_key_to_node_name:
                return node_key_to_node_name[node_key]
            node = self.nodes[node_key]
            if 'renormalisation_vertex_n_loops' in node:
                label = 'UV%dL'%node['renormalisation_vertex_n_loops']
                label += '%dP'%node['renormalisation_vertex_n_shrunk_edges']
                # If that label is not unique, then prefix it the actual key.
                if label in node_key_to_node_name.values():
                    label = '%s_%s'%(str(node_key),label)
            else:
                label = str(node_key)

            node_key_to_node_name[node_key] = label
            return label

        # Generate edge list
        edge_template = """Labeled[Style[CreateEdge["%(in_node)s","%(out_node)s",%(edge_key)d]%(edge_style)s,%(edge_color)s,Thickness[%(thickness)f]],"%(edge_label)s"]"""
        all_edge_definitions = []
        all_edge_shape_functions = []
        for edge_key, edge_data in self.edges.items():
            if not isinstance(edge_key, tuple):
                edge_key = (*edge_data['vertices'], edge_key)
            edge_repl_dict = {}
            # is_LMB = ('name' in edge_data and 'LMB' in str(edge_data['name']).upper())
            abs_sig = ( [abs(s) for s in edge_data['signature'][0]], [abs(s) for s in edge_data['signature'][1]])
            is_LMB = (sum(abs_sig[0]) == 1 and sum(abs_sig[1]) == 0)
            if is_LMB:
                edge_repl_dict['thickness'] = 0.005
            else:
                edge_repl_dict['thickness'] = 0.002
            edge_repl_dict['in_node'] = get_node_label(edge_key[0])
            edge_repl_dict['out_node'] = get_node_label(edge_key[1])
            edge_repl_dict['edge_key'] = edge_key[2]
            edge_repl_dict['arrow_style'] = 'Arrow' if not is_LMB else 'HalfFilledDoubleArrow'
            edge_repl_dict['arrow_size'] = 0.015 if not is_LMB else 0.025
            all_edge_shape_functions.append(
                'CreateEdge["%(in_node)s","%(out_node)s",%(edge_key)d]->GraphElementData["%(arrow_style)s", "ArrowSize" -> %(arrow_size)f]'%edge_repl_dict
            )
            if 'name' in edge_data and 'CUT' in str(edge_data['name']).upper():
                edge_repl_dict['edge_style'] = ",Dashed"
            else:
                edge_repl_dict['edge_style'] = ""
            if edge_data['PDG'] in [-1,-2,-3,-4,-5,1,2,3,4,5]:
                color = "Cyan"
            elif edge_data['PDG'] in [-6,6]:
                color = "Blue"
            elif edge_data['PDG'] in [21,]:
                color = "Red"
            elif edge_data['PDG'] in [82,-82]:
                color = "Pink"
            elif edge_data['PDG'] in [25,]:
                color = "Green"
            else:
                color = "Gray"
            edge_repl_dict['edge_color'] = color
            edge_label_pieces = ['psi' if edge_data['PDG'] == 1337 else get_part_name(edge_data['PDG']),]
            if 'name' in edge_data and self._include_edge_name_in_rendering:
                edge_label_pieces.append(edge_data['name'])
            if self._include_momentum_routing_in_rendering:
                edge_label_pieces.append(edge_data['momentum'])
            edge_label = "|".join(edge_label_pieces)
            edge_repl_dict['edge_label'] = edge_label
            all_edge_definitions.append(edge_template%edge_repl_dict)

        repl_dict['edge_lists'] = ',\n'.join(all_edge_definitions)
        repl_dict['edge_shape_definitions'] = ',\n'.join(all_edge_shape_functions)
        return \
"""Labeled[GraphClass[{
%(edge_lists)s
},
EdgeShapeFunction -> {
%(edge_shape_definitions)s
},
EdgeLabelStyle -> Directive[FontFamily -> "CMU Typewriter Text", FontSize -> %(edge_font_size)d, Bold],
VertexLabelStyle -> Directive[FontFamily -> "CMU Typewriter Text", FontSize -> %(vetex_font_size)d, Bold],
VertexSize -> Large,
VertexLabels -> Placed[Automatic,Center],
GraphLayout -> %(graph_layout_strategy)s,
ImageSize -> {%(width)f, %(height)f}
],"%(graph_name)s"]"""%repl_dict
    
    def draw(self, model, output_dir,FORM_id=None, lmb_id=None):
        """ Outputs the mathematica code for rendering this FORMSuperGraph."""
        
        if FORM_id is not None:
            file_name = 'Graph_%04d'%FORM_id
        else:
            file_name = 'Graph_%s'%self.name
        
        if lmb_id is not None:
            file_name += '_LMB_%04d'%lmb_id

        MM_code = \
"""GraphClass = If[$VersionNumber > 12, EdgeTaggedGraph, Graph];
CreateEdge[u_,v_,t_]:=If[$VersionNumber > 12, DirectedEdge[u, v, t], DirectedEdge[u, v]];
aGraph=%s;
"""%self.get_mathematica_rendering_code(model,FORM_id=FORM_id, lmb_id=lmb_id)
        # Export to PDF in landscape format. One graph per page for now.
        # The 1.2 multiplier accounts for margins
        MM_code += 'Export["%s.pdf", GraphicsGrid[{{aGraph}}], ImageSize -> {%f, %f}];'%(
                    file_name,1.25*self._rendering_size[0],1.25*self._rendering_size[1])
        open(pjoin(output_dir,'%s.m'%file_name),'w').write(MM_code)

    def generate_numerator_form_input(self, additional_overall_factor='', only_algebra=False):
        # create the input file for FORM
        form_diag = self.overall_factor+additional_overall_factor
        for node in self.nodes.values():
            if node['vertex_id'] < 0:
                continue
        
            form_diag += '*\nvx({},{},{})'.format(
                ','.join(str(p) for p in node['PDGs']),
                ','.join(node['momenta']),
                ','.join(str(i) for i in node['indices']),
            )

        for edge in self.edges.values():
            form_diag += '*\nprop({},{},{},{})'.format(
                edge['PDG'],
                edge['type'],
                edge['momentum'],
                ','.join(str(i) for i in edge['indices']),
            )

        if only_algebra:
            return form_diag

        if self.replacement_rules is None:
            raise AssertionError("No energy configurations specified for numerator: run the denominator generation first")

        # now add all the replacement rules
        form_diag += self.replacement_rules

        return form_diag

    def get_topo_generator(self, specified_LMB=None):
        """ Returns a topology generator for that FORMSuperGraph."""

        topo_edges = copy.deepcopy(self.edges)

        original_LMB = {}
        external_edges = []
        other_edges = []
        edge_name_to_key = {}
        # Relabel edges according to alphaLoop conventions:
        for edge_key, edge_data in topo_edges.items():
            # Fix for QGRAF pipeline
            if not isinstance(edge_key, tuple):
                edge_key = (*edge_data['vertices'], edge_key) 
          
            if edge_data['type'] == 'virtual':
                if not edge_data['name'].startswith('p'):
                    edge_data['name'] = 'p%s'%edge_data['name']
                other_edges.append((edge_data['name'],edge_key[0],edge_key[1]))
            else:
                if not edge_data['name'].startswith('q'):
                    edge_data['name'] = 'q%s'%edge_data['name'][1:]
                external_edges.append((edge_data['name'],edge_data['vertices'][0],edge_data['vertices'][1]))
            edge_name_to_key[edge_data['name']]=edge_key

            # Test if it is a defining edge of the lmb
            abs_sig = ( [abs(s) for s in edge_data['signature'][0]], [abs(s) for s in edge_data['signature'][1]])
            if sum(abs_sig[0]) == 1 and sum(abs_sig[1]) == 0:
                original_LMB[abs_sig[0].index(1)]=edge_data['name']

        topo_edges = external_edges+other_edges

        # Set the LMB to a sorted one
        original_LMB = sorted(list(original_LMB.items()),key=lambda e: e[0])
        assert(all(oLMBe[0]==i for i,oLMBe in enumerate(original_LMB)))
        original_LMB = [oLMBe[1] for oLMBe in original_LMB]

        topo_generator = LTD.ltd_utils.TopologyGenerator(topo_edges)
        topo_generator.generate_momentum_flow( loop_momenta = (original_LMB if specified_LMB is None else specified_LMB) )
        original_LMB = [edge_name_to_key[oLMBe] for oLMBe in original_LMB]

        return topo_generator, edge_name_to_key, original_LMB

    def generate_additional_LMBs(self):
        """ Depending on the FORM options 'number_of_lmbs' and 'reference_lmb', this function fills in the attribute 
        additional_lmbs of this class."""

        if FORM_processing_options['number_of_lmbs'] is None:
            return
    
        topo_generator, edge_name_to_key, original_LMB = self.get_topo_generator()
        edge_key_to_name = {v:k for k,v in edge_name_to_key.items()}

        all_lmbs = topo_generator.loop_momentum_bases()
        all_lmbs= [ tuple([edge_name_to_key[topo_generator.edge_map_lin[e][0]] for e in lmb]) for lmb in all_lmbs]
        
        # Then overwrite the reference LMB if the user requested it
        if FORM_processing_options['reference_lmb'] is not None:
            original_LMB = all_lmbs[(FORM_processing_options['reference_lmb']-1)%len(all_lmbs)]
            # Regenerate the topology with this new overwritten LMB
            topo_generator, _, _ = self.get_topo_generator(specified_LMB=[ edge_key_to_name[e_key] for e_key in original_LMB ] )
            # And adjust all signatures (incl. the string momenta assignment accordingly)
            signatures = topo_generator.get_signature_map()
            signatures = { edge_name_to_key[edge_name]: [sig[0],
                    [ i+o for i,o in zip(sig[1][:len(sig[1])//2],sig[1][len(sig[1])//2:]) ]
                ] for edge_name, sig in signatures.items() }
            for edge_key, edge_data in self.edges.items():
                # Fix for QGRAF pipeline
                if not isinstance(edge_key, tuple):
                    edge_key = (*edge_data['vertices'], edge_key) 
          
                edge_data['signature'] = signatures[edge_key]
                edge_data['momentum'] = FORMSuperGraph.momenta_decomposition_to_string(edge_data['signature'], set_outgoing_equal_to_incoming=False)
            for node_key, node_data in self.nodes.items():
                node_data['momenta'] = tuple([                    
                        FORMSuperGraph.momenta_decomposition_to_string(
                            [ 
                                [ s*(1 if self.edges[e_key]['vertices'][1]==node_key else -1) for s in self.edges[e_key]['signature'][0] ],
                                [ s*(1 if self.edges[e_key]['vertices'][1]==node_key else -1) for s in self.edges[e_key]['signature'][1] ],
                            ],
                            set_outgoing_equal_to_incoming=False)
                    for e_key in node_data['edge_ids']])

        original_lmb_signatures = topo_generator.get_signature_map()
        original_lmb_signatures = { edge_name_to_key[edge_name]: [sig[0],
                    [ i+o for i,o in zip(sig[1][:len(sig[1])//2],sig[1][len(sig[1])//2:]) ]
                ] for edge_name, sig in original_lmb_signatures.items() }

        # Now generate copies of this supergraph with different LMBs
        additional_lmbs_SGs = []
        for i_lmb, lmb in enumerate(all_lmbs):
            if lmb==original_LMB:
                continue
            if FORM_processing_options['number_of_lmbs']>=0 and len(additional_lmbs_SGs)==FORM_processing_options['number_of_lmbs']:
                break

            # Generate the topology with this additional LMB
            other_lmb_topo_generator, _, _ = self.get_topo_generator(specified_LMB=[ edge_key_to_name[e_key] for e_key in lmb ])
            # And adjust all signatures (incl. the string momenta assignment accordingly)
            other_lmb_signatures = other_lmb_topo_generator.get_signature_map()
            other_lmb_signatures = { edge_name_to_key[edge_name]: [sig[0],
                    [ i+o for i,o in zip(sig[1][:len(sig[1])//2],sig[1][len(sig[1])//2:]) ]
                ] for edge_name, sig in other_lmb_signatures.items() }

            # Compute the affine transformation to go from the original LMB to this additional LMB
            affine_transfo_other_lmb_to_original = []
            for defining_edge_key in lmb:
                other_sig = original_lmb_signatures[defining_edge_key]
                affine_transfo_other_lmb_to_original.append( [ list(other_sig[0]), list(other_sig[1]) ]  )
            affine_transfo_original_lmb_to_other_lmb = []
            for defining_edge_key in original_LMB:
                other_sig = other_lmb_signatures[defining_edge_key]
                affine_transfo_original_lmb_to_other_lmb.append( [ list(other_sig[0]), list(other_sig[1]) ]  )

            other_LMB_super_graph = copy.deepcopy(self)
            # Flag these additional supergraphs as *NOT* original representative by setting their 'additional_lmbs' to an integer instead
            # which identifies which additional LMB it corresponds to.
            other_LMB_super_graph.additional_lmbs = i_lmb+1
            for edge_key, edge_data in other_LMB_super_graph.edges.items():
                # Fix for QGRAF pipeline
                if not isinstance(edge_key, tuple):
                    edge_key = (*edge_data['vertices'], edge_key) 
          
                edge_data['signature'] = other_lmb_signatures[edge_key]
                edge_data['momentum'] = FORMSuperGraph.momenta_decomposition_to_string(edge_data['signature'], set_outgoing_equal_to_incoming=False)
            for node_key, node_data in other_LMB_super_graph.nodes.items():
                node_data['momenta'] = tuple([
                    FORMSuperGraph.momenta_decomposition_to_string(
                        [ 
                            [ s*(1 if other_LMB_super_graph.edges[e_key]['vertices'][1]==node_key else -1) for s in other_LMB_super_graph.edges[e_key]['signature'][0] ],
                            [ s*(1 if other_LMB_super_graph.edges[e_key]['vertices'][1]==node_key else -1) for s in other_LMB_super_graph.edges[e_key]['signature'][1] ],
                        ],
                        set_outgoing_equal_to_incoming=False)
                     for e_key in node_data['edge_ids']])

            additional_lmbs_SGs.append( (i_lmb+1, affine_transfo_other_lmb_to_original, affine_transfo_original_lmb_to_other_lmb, other_LMB_super_graph ) )

        self.additional_lmbs = additional_lmbs_SGs

    def derive_signatures(self):
        n_incoming = sum([1 for edge in self.edges.values() if edge['type'] == 'in'])
        n_loops = len(self.edges) - len(self.nodes) + 1

        # parse momentum
        p = re.compile(r'(^|\+|-)(k|p)(\d*)')
        for edge in self.edges.values():
            parsed = [e.groups() for e in p.finditer(edge['momentum'])]
            signature = ([0 for _ in range(n_loops)], [0 for _ in range(n_incoming)])
            for pa in parsed:
                if pa[1] == 'p':
                    signature[1][int(pa[2]) - 1] = 1 if pa[0] == '' or pa[0] == '+' else -1
                else:
                    signature[0][int(pa[2]) - 1] = 1 if pa[0] == '' or pa[0] == '+' else -1

            edge['signature'] = signature

    def impose_signatures(self):

        for eid, e in self.edges.items():
            e['momentum'] = FORMSuperGraph.momenta_decomposition_to_string(e['signature'], set_outgoing_equal_to_incoming=False)
            neg_mom =  FORMSuperGraph.momenta_decomposition_to_string([[-s for s in sp] for sp in e['signature']], set_outgoing_equal_to_incoming=False)
            for i, vi in enumerate(e['vertices']):
                e_index = self.nodes[vi]['edge_ids'].index(eid)
                mom = list(self.nodes[vi]['momenta'])
                mom[e_index] = neg_mom if i == 0 else e['momentum']
                self.nodes[vi]['momenta'] = mom

    @classmethod
    def sort_edges(cls, model, edges_to_sort):
        """ Sort the edges adjacent to a node so as to map the conventions
        from the Feynmna rules encoded in Form."""

        # The only relevant properties are the spin and wether an edge
        # is a particle or antiparticle.
        edge_particles = [model.get_particle(e['PDG']) for e in edges_to_sort]
        identities = [ ( p.get('spin') * (1 if p.get('is_part') else -1) , position ) for position, p in enumerate(edge_particles) ]
        identities.sort(key=lambda el: el[0])
        canonical_identifier = tuple([identity for identity, position in identities])

        if canonical_identifier not in cls._FORM_Feynman_rules_conventions:
            raise FormProcessingError("Conventions for FORM Feynman rules of signature {} not specifed.".format(canonical_identifier))
        
        new_position = cls._FORM_Feynman_rules_conventions[canonical_identifier]

        return [ edges_to_sort[identities[position][1]] for position in new_position ]

    @classmethod
    def momenta_decomposition_to_string(cls, momenta_decomposition, set_outgoing_equal_to_incoming=True):
        """ Turns ((1,0,0,-1),(1,1)) into 'k1-k4+p1+p2'"""

        res = ""
        first=True
        # The outgoing momenta are set element-wise equal to the incoming ones.
        if set_outgoing_equal_to_incoming:
            momenta_decomposition = [
                momenta_decomposition[0],
                [inp+outp for inp, outp in zip(
                    momenta_decomposition[1][:len(momenta_decomposition[1])//2],
                    momenta_decomposition[1][len(momenta_decomposition[1])//2:]
                )]
            ]
        # Fuse outgoing and incoming
        for symbol, mom_decomposition in zip(('k','p'),momenta_decomposition):
            for i_k, wgt in enumerate(mom_decomposition):
                if wgt!=0:
                    if first:
                        if wgt<0:
                            res+="-"
                        first=False
                    else:
                        if wgt<0:
                            res+="-"
                        else:
                            res+="+"
                    if abs(wgt)!=1:
                        res+="%d*"%abs(wgt)
                    res+="%s%d"%(symbol,(i_k+1))

        return res

    @classmethod
    def from_LTD2SuperGraph(cls, LTD2_super_graph):
        """ Creates a FORMSuperGraph from an LTD2_super_graph."""

        # First make sure the momentum routing of the LTD2_super_graph
        # is set:
        LTD2_super_graph.set_momentum_routing()

        overall_factor = "1"

        if LTD2_super_graph.symmetry_factor is not None:
            multiplicity = LTD2_super_graph.symmetry_factor
        else:
            multiplicity = 1

        # Let us just match the overall phase picked for the MG num:
        overall_phase = complex(-1.0,0.0)**len(LTD2_super_graph.cuts)
        if overall_phase.imag != 0:
            raise FormProcessingError("No support for overall complex phase yet (Ben: how do we put a complex number in FORM? ^^)")
        else:
            overall_factor += '*%d'%int(overall_phase.real)

        fermion_factor = LTD2_super_graph.diag_left_of_cut.fermion_factor*LTD2_super_graph.diag_right_of_cut.fermion_factor
        overall_factor += '*%d'%fermion_factor

        model = LTD2_super_graph.model

        # Let us also include a factor -1 for each closed ghost loop.
        # This is automatically done in MadGraph already by considering the "wavefunction" of external scalar ghosts to be sqrt(i).
        # Said differently, we want a factor -1 for each pair of two ghosts in the final state
        n_ghosts = len([ 1 for c in LTD2_super_graph.cuts if model.get_particle(LTD2_super_graph.graph.edges[c[1]]['pdg']).get('ghost') ])
        assert(n_ghosts%2==0)
        if n_ghosts > 0:
            overall_factor += '*%d'%(-1**(n_ghosts//2))

        local_graph = copy.deepcopy(LTD2_super_graph.graph)
        
        # Collect the outer-most nodes
        outer_in_nodes = []
        outer_out_nodes = []
        for node_key, node_data in local_graph.nodes.items():
            # Assign the corresponding edge
            if node_key.startswith('I') or node_key.startswith('O'):
                for edge_key, edge_data in local_graph.edges.items():
                    if edge_key[0]==node_key or edge_key[1]==node_key:
                        node_data['edge_ids'] = tuple([edge_key,])
                        break
                if 'edge_ids' not in node_data:
                    node_data['edge_ids'] = tuple([])
            if node_key.startswith('I'):
                outer_in_nodes.append(node_key)
                node_data['momenta'] = ('p%d'%int(node_key[1:]),)
                node_data['vertex_id'] = -1 # Incoming node at assigned vertex ID = -1
            if node_key.startswith('O'):
                outer_out_nodes.append(node_key)
                node_data['momenta'] = ('p%d'%int(node_key[1:]),)
                node_data['vertex_id'] = -2 # Incoming node at assigned vertex ID = -2

        curr_index = 0
        # Assign indices to the outermost nodes
        for node_key in outer_in_nodes:
            this_index = int(node_key[1:])
            curr_index = max(curr_index,this_index)
            local_graph.nodes[node_key]['indices'] = (this_index,)

        n_incoming=int(curr_index)
        for node_key in outer_out_nodes:
            this_index = n_incoming+int(node_key[1:])
            curr_index = max(curr_index,this_index)
            local_graph.nodes[node_key]['indices'] = (this_index,)

        # assign indices to all edges now
        for edge_key, edge_data in local_graph.edges.items():
            edge_indices = []
            found_external_node = False
            for node_key in [edge_key[0],edge_key[1]]:
                if node_key in outer_in_nodes+outer_out_nodes:
                    edge_indices = [ local_graph.nodes[node_key]['indices'][0], ]
                    found_external_node = True
                    local_graph.nodes[node_key]['PDGs'] = (edge_data['pdg'],)
                else:
                    curr_index += 1
                    if not found_external_node:
                        edge_indices.append(curr_index)
            edge_data['indices'] = tuple(edge_indices)
            if edge_key[0] in outer_in_nodes or edge_key[1] in outer_in_nodes:
                edge_data['type'] = 'in'
#                loop_momenta = tuple([0 for _ in range(LTD2_super_graph.n_loops)])
#                external_momenta = [0 for _ in range(len(outer_in_nodes))]
#                if edge_key[0] in outer_in_nodes:
#                    external_momenta[int(edge_key[0][1:])-1] = 1
#                elif edge_key[1] in outer_in_nodes:
#                    external_momenta[int(edge_key[1][1:])-1] = -1
#                edge_data['momentum'] = (loop_momenta, tuple(external_momenta))
            elif edge_key[0] in outer_out_nodes or edge_key[1] in outer_out_nodes:
                edge_data['type'] = 'out'
#                loop_momenta = tuple([0 for _ in range(LTD2_super_graph.n_loops)])
#                external_momenta = [0 for _ in range(len(outer_in_nodes))]
#                if edge_key[1] in outer_out_nodes:
#                    external_momenta[int(edge_key[1][1:])-1] = 1
#                elif edge_key[0] in outer_out_nodes:
#                    external_momenta[int(edge_key[0][1:])-1] = -1
#                edge_data['momentum'] = (loop_momenta, tuple(external_momenta))
            else:
                edge_data['type'] = 'virtual'
            edge_data['PDG'] = edge_data['pdg']
            del edge_data['pdg']
            # Keep the name (used for rendering only)
#            del edge_data['name']

        # assign node extra information
        for node_key, node_data in local_graph.nodes.items():
            # Skip external nodes that have already been treated
            if node_key in outer_in_nodes+outer_out_nodes:
                continue
            # First collect all adjacent edges:
            adjacent_in_edges = [
                dict(list(in_edge_data.items())+[('key',(u,v,c))]) for u,v,c,in_edge_data in 
                local_graph.in_edges(node_key,data=True,keys=True)
            ]
            adjacent_out_edges = [
                dict(list(out_edge_data.items())+[('key',(u,v,c))]) for u,v,c,out_edge_data in 
                local_graph.out_edges(node_key,data=True,keys=True)
            ]
            # The direction only matters in so far as we must flip the
            # momentum carried by edges that are outgoing as well as
            # change their PDG from particle to antiparticle
            for data in adjacent_in_edges:
                if len(data['indices'])==2:
                    data['index_for_this_node'] =  data['indices'][1]
                else:
                    # This is an external edge with a single index.
                    data['index_for_this_node'] =  data['indices'][0]
            for data in adjacent_out_edges:
                data['index_for_this_node'] =  data['indices'][0]
                data['momentum'] = (
                    tuple([-p for p in data['momentum'][0]]),
                    tuple([-p for p in data['momentum'][1]])
                )
                data['PDG'] = model.get_particle(data['PDG']).get_anti_pdg_code()

            # We must now sort the edges according to the special rules for the vertex Feynman rules
            all_adjacent_edges = cls.sort_edges(model,adjacent_in_edges+adjacent_out_edges)
            node_data['PDGs'] = tuple([e['PDG'] for e in all_adjacent_edges])
            node_data['indices'] = tuple([e['index_for_this_node'] for e in all_adjacent_edges])
            node_data['momenta'] = tuple([cls.momenta_decomposition_to_string(e['momentum']) for e in all_adjacent_edges])
            node_data['edge_ids'] = tuple([e['key'] for e in all_adjacent_edges])

            # Example printout information about the vertex
            #misc.sprint("Vertex of node %s:\n%s"%(
            #        str(node_key),pformat(model.get_interaction(node_data['vertex_id']))))

        # Finally overwrite the edge momentum so as to be a string
        for edge_key, edge_data in local_graph.edges.items():
            edge_data['momentum'] =cls.momenta_decomposition_to_string(edge_data['momentum'])
            edge_data['vertices'] =(edge_key[0],edge_key[1])
            # In FORM conventions the fermion are going together with their flow, so we need
            # to flip the order of their fundamental/antifundamental indices so that FORM
            # builds the correct propagator. 
            # NO LONGER NEEDED: This now directly done right in FORM.
#            particle = model.get_particle(edge_data['PDG'])
#            if len(edge_data['indices'])>1 and (particle.get('spin')%2==0 or particle.get('ghost')):
#                if not particle.get('is_part'):
#                    edge_data['indices'] = tuple([edge_data['indices'][1],edge_data['indices'][0]])

        graph_name = 'P%(proc_id)dL%(left_diagram_id)dR%(right_diagram_id)d'%LTD2_super_graph.call_signature
        form_super_graph =  cls(
            name = graph_name,
            call_identifier=LTD2_super_graph.call_signature,
            edges = dict(local_graph.edges),
            nodes = dict(local_graph.nodes),
            overall_factor = overall_factor,
            multiplicity = multiplicity
        )
        #misc.sprint(graph_name)
        #misc.sprint(pformat(dict(LTD2_super_graph.graph.edges.items())))
        #misc.sprint(pformat(dict(local_graph.edges.items())))

        return form_super_graph

    @classmethod
    def from_dict(cls, file_path):
        """ Creates a FORMSuperGraph from a Python dict file path."""

        # TODO: Creates an instance from a Python dict dump.
        pass

    def to_dict(self, file_path=None):
        """ Outputs the FORMSuperGraph self to a Python dict file path."""

        dict_to_dump = {
            'edges' : self.edges,
            'nodes' : self.nodes,
            'overall_factor' : self.overall_factor,
            'multiplicity' : self.multiplicity
        }
        if file_path:
            open(file_path,'w').write(pformat(dict_to_dump))
        else:
            return dict_to_dump

    def generate_ltd_integrand(self, topo, workspace, numerator_call):
        """Construct a table of integrand descriptors"""
        integrand_body = ''
        max_diag_id = 0

        for cut, cut_loop_topos in zip(topo.cuts, topo.cut_diagrams):
            for diag_set, loop_diag_set in zip(cut['diagram_sets'], cut_loop_topos):
                energy_map, energies, constants, shift_map = [], [], [], []
                prop_mom_in_lmb = []
                signature_offset = 0
                loops = 0 # LTD loop count

                prop_id = {}
                counter = 0
                for di, diag_info in enumerate(loop_diag_set['diagram_info']):
                    for li, l in enumerate(diag_info['graph'].loop_lines):
                        if all(s == 0 for s in l.signature):
                            continue
                        for pi, p in enumerate(l.propagators):
                            prop_id[(di, li, pi)] = counter
                            counter += 1

                res = []
                propagators = []
                #ltd_signatures = {}
                propcount = 1
                for di, diag_info in enumerate(loop_diag_set['diagram_info']):
                    for li, l in enumerate(diag_info['graph'].loop_lines):
                        is_constant = all(s == 0 for s in l.signature)

                        #ltd_signatures[(di, li)] = [0]* l.signature

                        for p in l.propagators:
                            # contruct the momentum in the LMB, using that LTD
                            lmp = np.array([0]*topo.topo.n_loops)
                            for s, v in zip(l.signature, diag_info['loop_momentum_map']):
                                lmp += s * np.array(v[0])

                            # transport shift to LMB
                            shift = np.array([0]*topo.topo.n_loops)
                            extshift = np.array(p.parametric_shift[1])
                            for s, c in zip(p.parametric_shift[0], cut['cuts']):
                                shift += s * np.array(c['signature'][0])
                                extshift += s * np.array(c['signature'][1])
                            extshift = np.array(list(extshift[:len(extshift)//2]) + [0]*(len(extshift)//2)) +\
                                np.array(list(extshift[len(extshift)//2:]) + [0]*(len(extshift)//2))

                            totalmom = self.momenta_decomposition_to_string((lmp + shift, extshift), True)

                            # TODO: recycle energy computations since Es will appear more than once
                            if not is_constant:
                                prop_mom_in_lmb.append((lmp, shift, extshift))
                                energy_map.append(p.m_squared if not p.uv else 'mUV*mUV')
                                energies.append(totalmom)
                                shift_map.append(list(shift) + list(extshift))

                            for _ in range(p.power):
                                if is_constant:
                                    constants.append((totalmom, p.m_squared if not p.uv else 'mUV*mUV'))

                    diagres = []
                    # enumerate all cut options
                    for co in [x for css in diag_info['graph'].ltd_cut_structure for x in 
                            product(*[[(cs, (li, i)) for i in range(len(l.propagators))] for li, (cs, l) in enumerate(zip(css, diag_info['graph'].loop_lines)) if cs != 0])]:
                        ltd_closure_factor = int(np.prod([s for (s, _) in co]))
                        
                        # construct the cut basis to LTD loop momentum basis mapping, used to substitute the numerator
                        mat = [diag_info['graph'].loop_lines[li].signature for (_, (li, _)) in co]
                        if mat == []:
                            nmi = []
                        else:
                            nmi = np.linalg.inv(np.array(mat).transpose())

                        m = []
                        ltdenergy = ['ltd{0},{1}E{0}'.format(prop_id[(di, li, pi)], '+' if cut_sign == 1 else '-') for (cut_sign, (li, pi)) in co]
                        for i, r in enumerate(nmi):
                            mm = []
                            for (c, (_, (li, pi))) in zip(r, co):
                                if c != 0:
                                    ext = self.momenta_decomposition_to_string((prop_mom_in_lmb[prop_id[(di,li, pi)]][1], prop_mom_in_lmb[prop_id[(di, li, pi)]][2]), True)
                                    if ext == '':
                                        ext = '0'
                                    
                                    mm.append('{}ltd{}{}energies({})'.format('+' if c == 1 else '-', prop_id[(di, li, pi)], '-' if c == 1 else '+', ext))
                            m += ['c{}'.format(i + len(cut['cuts']) + signature_offset), '+'.join(mm)]

                        r = []
                        der = []
                        for li, l in enumerate(diag_info['graph'].loop_lines):
                            if all(s == 0 for s in l.signature):
                                continue

                            energy = []
                            energy_full = ''
                            sig_map = nmi.dot(l.signature)
                            for (sig_sign, (cut_sign, (lci, pci))) in zip(sig_map, co):
                                if sig_sign != 0:
                                    momp = self.momenta_decomposition_to_string((prop_mom_in_lmb[prop_id[(di, lci, pci)]][1], prop_mom_in_lmb[prop_id[(di, lci, pci)]][2]), True)
                                    energy.append('{},ltd{},{}energies({})'.format(int(sig_sign), prop_id[(di, lci, pci)], 
                                        '+' if -sig_sign == 1 else '-', 
                                        '0' if momp == '' else momp))
                                    # the full energy including the cut sign
                                    energy_full += '{}E{}{}energies({})'.format('+' if int(cut_sign * sig_sign) == 1 else '-', prop_id[(di, lci, pci)], 
                                        '+' if -int(sig_sign) == 1 else '-', 
                                        '0' if momp == '' else momp)

                            for pi, p in enumerate(l.propagators):
                                # TODO: recycle propagator ids if the functional form is the same!
                                powmod = '' if p.power == 1 else '^' +  str(p.power)
                                if (1, (li, pi)) in co:
                                    propagators.append('2*E{0}'.format(prop_id[(di, li, pi)]))
                                    r.append('prop({0},1,ltd{1},0,E{1}){2}'.format(propcount, prop_id[(di, li, pi)], powmod))
                                    if p.power > 1:
                                        # add derivative prescription
                                        der.append('ltd{},{}'.format(prop_id[(di,li, pi)], p.power - 1))
                                elif (-1, (li, pi)) in co:
                                    propagators.append('-2*E{0}'.format(prop_id[(di, li, pi)]))
                                    r.append('prop({0},-1,ltd{1},0,-E{1}){2}'.format(propcount, prop_id[(di, li, pi)], powmod))
                                    if p.power > 1:
                                        der.append('ltd{},{}'.format(prop_id[(di, li, pi)], p.power - 1))
                                else:
                                    momp = self.momenta_decomposition_to_string((prop_mom_in_lmb[prop_id[(di, li, pi)]][1], prop_mom_in_lmb[prop_id[(di, li, pi)]][2]), True)
                                    r.append('prop({0},{1},E{2}+energies({3})){4}*prop({5},{1},-E{2}+energies({3})){4}'.format(propcount, ','.join(energy), prop_id[(di, li, pi)],
                                        '0' if momp == '' else momp, powmod, propcount + 1))

                                    propagators.append('{}+E{}+energies({})'.format(energy_full, prop_id[(di, li, pi)], '0' if momp == '' else momp))
                                    propagators.append('{}-E{}+energies({})'.format(energy_full, prop_id[(di, li, pi)], '0' if momp == '' else momp))
                                    propcount += 1
                                propcount += 1
                        diagres.append('{}*{}{}{}({})'.format(
                            ltd_closure_factor,
                            'ltdcbtolmb({})*'.format(','.join(m)) if len(m) > 0 else '',
                            'ltdenergy({})*'.format(','.join(ltdenergy)) if len(ltdenergy) > 0 else '',
                            'der({})*'.format(','.join(der))  if len(der) > 0 else '',
                            '*'.join(r) if len(r) > 0 else '1'))

                    res.append('\n\t\t\t+'.join(diagres))
                    loops += diag_info['graph'].n_loops
                    signature_offset += diag_info['graph'].n_loops

                res = '\n\t\t*'.join(['({})'.format(l) for l in res])

                self.integrand_info[diag_set['id']] = (energy_map, constants, loops)
                max_diag_id = max(max_diag_id, diag_set['id'])
                integrand_body += 'Fill ltdtopo({}) = (-1)^{}*constants({})*\n\tellipsoids({})*\n\tallenergies({})*(\n\t\t{}\n);\n'.format(diag_set['id'],
                    loops, ','.join(c[0] for c in constants), ','.join(propagators), ','.join(energies), res)

        with open(pjoin(workspace, 'ltdtable_{}.h'.format(numerator_call)), 'w') as f:
            f.write("""
Auto S invd, E, shift, ltd;
S r, s;
CF a, num, ncmd, conf1, replace, energies, ellipsoids, ltdcbtolmb, ltdenergy, constants;
NF allenergies;
Set invdset: invd0,...,invd400;
CTable ltdtopo(0:{});

{}
""".format(max_diag_id, integrand_body))

    def generate_integrand(self, topo, workspace, numerator_call, progress_bar=None):
        """Construct a table of integrand descriptors"""
        integrand_body = ''
        max_diag_id = 0

        for cut, cut_loop_topos in zip(topo.cuts, topo.cut_diagrams):
            for diag_set, loop_diag_set in zip(cut['diagram_sets'], cut_loop_topos):
                signatures, n_props, energy_map, energies, constants, shift_map = [], [], [], [], [], []
                signature_offset = 0
                total_ltd_loops = topo.topo.n_loops - len(cut['cuts']) + 1
                for diag_info in loop_diag_set['diagram_info']:
                    for l in diag_info['graph'].loop_lines:
                        is_constant = all(s == 0 for s in l.signature)
                        if not is_constant:
                            signatures.append([0]*signature_offset + list(l.signature) + [0]*(total_ltd_loops
                                    - signature_offset - diag_info['graph'].n_loops)) # padded LTD signatures
                            n_props.append(sum(p.power for p in (l.propagators)))
                        for p in l.propagators:
                            # contruct the momentum in the LMB, using that LTD
                            lmp = np.array([0]*topo.topo.n_loops)
                            for s, v in zip(l.signature, diag_info['loop_momentum_map']):
                                lmp += s * np.array(v[0])

                            # transport shift to LMB
                            shift = np.array([0]*topo.topo.n_loops)
                            extshift = np.array(p.parametric_shift[1])
                            for s, c in zip(p.parametric_shift[0], cut['cuts']):
                                shift += s * np.array(c['signature'][0])
                                extshift += s * np.array(c['signature'][1])
                            extshift = np.array(list(extshift[:len(extshift)//2]) + [0]*(len(extshift)//2)) +\
                                np.array(list(extshift[len(extshift)//2:]) + [0]*(len(extshift)//2))

                            totalmom = self.momenta_decomposition_to_string((lmp + shift, extshift), True)

                            for _ in range(p.power):
                                # TODO: recycle energy computations when there are duplicate edges
                                if not is_constant:
                                    energy_map.append((p.m_squared if p.m_squared != 0 else 'small_mass_sq') if not p.uv else 'mUV*mUV')
                                    energies.append(totalmom)
                                    shift_map.append(list(shift) + list(extshift))
                                else:
                                    constants.append((totalmom, (p.m_squared if p.m_squared != 0 else 'small_mass_sq') if not p.uv else 'mUV*mUV'))
                    signature_offset += diag_info['graph'].n_loops

                if len(signatures) == 0:
                    # no loop dependence for this cut
                    res = '\t1\n'
                    resden = ''
                else:
                    pf = LTD.partial_fractioning.PartialFractioning(n_props, signatures,
                                            name=str(diag_set['id']), shift_map=np.array(shift_map).T,
                                            n_sg_loops=topo.topo.n_loops, ltd_index=len(cut['cuts']) - 1,
                                            progress_bar = progress_bar)
                    pf.shifts_to_externals()
                    res = pf.to_FORM()
                    res = '\n'.join(['\t' + l for l in res.split('\n')])
                    resden = ','.join(pf.den_library)

                self.integrand_info[diag_set['id']] = (energy_map, constants, total_ltd_loops)
                max_diag_id = max(max_diag_id, diag_set['id'])
                integrand_body += 'Fill pftopo({}) = constants({})*\nallenergies({})*\nellipsoids({})*(\n{});\n'.format(diag_set['id'],
                    ','.join(c[0] for c in constants), ','.join(energies), resden, res)

        with open(pjoin(workspace, 'pftable_{}.h'.format(numerator_call)), 'w') as f:
            f.write("""
Auto S invd, E, shift;
S r, s;
CF a, num, ncmd, conf1, ellipsoids, allenergies, replace, constants;
NF energies;
CTable pftopo(0:{});

{}
""".format(max_diag_id, integrand_body))


    def get_edge_scaling(self, pdg):
        # all scalings that deviate from -2
        scalings = {1: -1, 2: -1, 3: -1, 4: -1, 5: -1, 6: -1, 11: -1, 12: -1, 13: -1}
        return scalings[abs(pdg)] if abs(pdg) in scalings else -2

    def get_node_scaling(self, pdgs):
        # only the triple gluon vertex and the ghost gluon vertex have a non-zero scaling
        if pdgs == (21, 21, 21) or pdgs == (-82, 21, 82):
            return 1
        else:
            return 0

    def generate_squared_topology_files(self, root_output_path, model, process_definition, n_jets, numerator_call, final_state_particle_ids=(),jet_ids=None, write_yaml=True, bar=None,
        integrand_type=None, workspace=None):
        if workspace is None:
            workspace = pjoin(root_output_path, os.pardir, 'workspace')

        if bar:
            bar.update(i_lmb='1')
            max_lmb = 1
            if isinstance(self.additional_lmbs, list):
                max_lmb += len(self.additional_lmbs)
            bar.update(max_lmb='%d'%max_lmb)

        if self.is_zero:
            return False

        # Relabel edges according to alphaLoop conventions:
        for edge_key, edge_data in self.edges.items():
            edge_data['name'] = 'p' + edge_data['name'] if edge_data['type'] == 'virtual' else 'q' + edge_data['name'][1:]

        # TODO: sort such that the first 4 entries are external (it seems to happen by chance now every time)
        edge_map_lin = [(e['name'], e['vertices'][0], e['vertices'][1]) for e in self.edges.values()]
        assert(e[0] != 'q' or int(e[1:]) < 5 for e in edge_map_lin)

        particle_ids = { e['name']: e['PDG'] for e in self.edges.values() }
        particle_masses = {e['name']: 0. if e['PDG'] == 1337 else model['parameter_dict'][model.get_particle(e['PDG']).get('mass')].real for e in self.edges.values()}

        num_incoming = sum(1 for e in edge_map_lin if e[0][0] == 'q') // 2

        if num_incoming == 1:
            external_momenta = {'q1': [500., 0., 0., 0.], 'q2': [500., 0., 0., 0.]}
            #external_momenta = {'q1': [1., 0., 0., 0.], 'q2': [1., 0., 0., 0.]}
            p = np.array(external_momenta['q1'])
        else:
            external_momenta = {'q1': [500., 0., 0., 500.], 'q2': [500., 0., 0., -500.], 'q3': [500., 0., 0., 500.], 'q4': [500., 0., 0., -500.]}
            #external_momenta = {'q1': [1., 0., 0., 1.], 'q2': [1., 0., 0., -1.], 'q3': [1., 0., 0., 1.], 'q4': [1., 0., 0., -1.]}
            p = np.array(external_momenta['q1']) + np.array(external_momenta['q2'])

        loop_momenta = []
        n_loops = len(self.edges) - len(self.nodes) + 1
        for loop_var in range(n_loops):
            lm = next((ee['name'], ee['signature'][0][loop_var]) for ee in self.edges.values() if all(s == 0 for s in ee['signature'][1]) and \
                sum(abs(s) for s in ee['signature'][0]) == 1 and ee['signature'][0][loop_var] == 1)
            loop_momenta.append(lm)

        if isinstance(self.additional_lmbs, list):
            call_signature_ID = numerator_call
        else:
            call_signature_ID = self.additional_lmbs*FORM_processing_options['FORM_call_sig_id_offset_for_additional_lmb']+numerator_call

        topo = LTD.squared_topologies.SquaredTopologyGenerator(edge_map_lin,
            self.name, ['q1', 'q2'][:num_incoming], n_jets, external_momenta,
            loop_momenta_names=tuple([l for l,s in loop_momenta]),
            loop_momenta_signs=tuple([s for l,s in loop_momenta]),
            particle_ids=particle_ids,
            masses=particle_masses,
            final_state_particle_ids=final_state_particle_ids,
            jet_ids=jet_ids,
            overall_numerator=1.0,
            numerator_structure={},
            FORM_numerator={'call_signature': {'id': call_signature_ID}},
            FORM_integrand={'call_signature': {'id': call_signature_ID}},
            edge_weights={e['name']: self.get_edge_scaling(e['PDG']) for e in self.edges.values()},
            vertex_weights={nv: self.get_node_scaling(n['PDGs']) for nv, n in self.nodes.items()},
            generation_options=FORM_processing_options,
            analytic_result=(self.benchmark_result if hasattr(self,"benchmark_result") else None)
        )
        # check if cut is possible
        if len(topo.cuts) == 0:
            logger.info("No cuts for graph {}".format(self.name))
            return False

        self.generate_replacement_rules(topo)
        self.squared_topology = topo

        # Also generate the squared topology yaml files of all additional LMB topologies used for cross-check
        if isinstance(self.additional_lmbs, list):
            for i_lmb, (_,_,_,other_lmb_supergraph) in enumerate(self.additional_lmbs):
                if bar:
                    bar.update(i_lmb='%d'%(i_lmb+2))
                other_lmb_supergraph.generate_squared_topology_files(root_output_path, model, process_definition, n_jets, numerator_call, 
                        final_state_particle_ids=final_state_particle_ids,jet_ids=jet_ids, write_yaml=write_yaml,workspace=workspace,
                        bar=bar, integrand_type=integrand_type)

        if integrand_type is not None:
            if integrand_type == "LTD":
                self.generate_ltd_integrand(topo, workspace, call_signature_ID)
            if integrand_type == "PF":
                self.generate_integrand(topo, workspace, call_signature_ID, progress_bar = bar)

        if write_yaml:
            if isinstance(self.additional_lmbs, int):
                topo.export(pjoin(root_output_path, "%s_LMB%d.yaml"%(self.name,self.additional_lmbs)))
            else:
                topo.export(pjoin(root_output_path, "%s.yaml"%self.name))

        return True


    def generate_replacement_rules(self, topo):
        # collect the transformations of the bubble
        configurations = []
        bubble_to_cut = OrderedDict()
        for cut, cut_loop_topos in zip(topo.cuts, topo.cut_diagrams):
            for diag_set, loop_diag_set in zip(cut['diagram_sets'], cut_loop_topos):
                trans = ['1']
                diag_set_uv_conf = []

                for diag_info, loop_diag_info in zip(diag_set['diagram_info'], loop_diag_set['diagram_info']):
                    der_edge = None
                    if diag_info['derivative'] is not None:
                        trans.append('1/2') # add a factor 1/2 since the bubble will appear in two cuts

                        # if the cutkosky cut has a negative sign, we derive in -p^0 instead of p^0.
                        # here we compensate for this sign
                        trans.append(str(next(c for c in cut['cuts'] if c['edge'] == diag_info['derivative'][0])['sign']))

                        der_edge = diag_info['derivative'][1]
                        ext_mom = next(ee for ee in self.edges.values() if ee['name'] == diag_info['derivative'][0])['momentum']
                        der_mom = next(ee for ee in self.edges.values() if ee['name'] == diag_info['derivative'][1])['momentum']
                        ext_sig = next(ee for ee in self.edges.values() if ee['name'] == diag_info['derivative'][0])['signature']
                        der_sig = next(ee for ee in self.edges.values() if ee['name'] == diag_info['derivative'][1])['signature']

                        if diag_info['bubble_momenta'] not in bubble_to_cut:
                            bubble_to_cut[diag_info['bubble_momenta']] = (diag_info['derivative'][0], set())

                        if diag_info['derivative'][0] == diag_info['derivative'][1]:
                            # numerator derivative
                            index = next(i for i, bub in enumerate(bubble_to_cut.keys()) if bub == diag_info['bubble_momenta'])
                            trans.append('der(pbubble' + str(index) + ')')
                        else:
                            # check if we pick up a sign change due to the external momentum flowing in the opposite direction
                            signs = [se * sc for se, sc in zip(der_sig[0] + der_sig[1], ext_sig[0] + ext_sig[1]) if se * sc != 0]
                            assert(len(set(signs)) == 1)
                            trans.append('-2*{}({})'.format('1*' if signs[0] == 1 else '-1*', der_mom))
                            bubble_info = bubble_to_cut[diag_info['bubble_momenta']]
                            bubble_to_cut[diag_info['bubble_momenta']] = (bubble_info[0], bubble_info[1] | {der_edge})

                    # translate the UV forest
                    uv_subgraphs = []

                    if diag_info['uv_info'] is None:
                        continue

                    uv_info = diag_info['uv_info']

                    uv_sig = '*'.join(['t{}^{}'.format(x, uv_info['derived_loop_lines'].count(x)) for x in range(len(uv_info['loop_lines']))])

                    # get the external momenta without sign
                    # note: the signatures are unaffected by the bubble treatment from before
                    external_momenta = set(m for e in uv_info['external_edges'] 
                        for m in self.momenta_decomposition_to_string(next(ee for ee in self.edges.values() if ee['name'] == e)['signature'], False)
                                .replace('-', '+').split('+') if m != '')
                    
                    uv_props = []
                    for i, uv_ll_props in enumerate(uv_info['loop_lines']):
                        # loop through all parametric shifts of all the propagators in this loop line
                        # find the associated LTD loop line
                        ll_sig, propagators = next(ll for ll in loop_diag_info['uv_loop_lines'] if set(p[0] for p in ll[1]).issuperset(set(uv_ll_props)))

                        # the parametric shift is given in terms of external momenta of the subgraph
                        # translate the signature and param_shift to momenta of the supergraph
                        loop_mom_sig = ''
                        for s, lmm in zip(ll_sig, loop_diag_info['loop_momentum_map']):
                            if s != 0:
                                loop_mom_sig += '{}({})'.format('+' if s == 1 else '-', self.momenta_decomposition_to_string(lmm, False))

                        for (edge_name, param_shift) in propagators:
                            ext_mom_sig = ''

                            if all(s == 0 for s in param_shift[1]):
                                continue

                            for (ext_index, s) in enumerate(param_shift[1]):
                                if s != 0:
                                    ext_mom = diag_info['graph'].edge_map_lin[diag_info['graph'].ext[ext_index]][0]
                                    ext_edge = next(ee for ee in self.edges.values() if ee['name'] == ext_mom)
                                    ext_mom_sig += '{}({})'.format('+' if s == 1 else '-',
                                        self.momenta_decomposition_to_string(ext_edge['signature'], False))

                            # the edge may have a raised power due to the bubble derivative
                            power = 2 if edge_name == diag_info['derivative_edge'] else 1
                            for _ in range(power):
                                uv_props.append('uvprop({},t{},{})'.format(loop_mom_sig, i, ext_mom_sig))
                    # it could be that there are no propagators with external momentum dependence when pinching duplicate edges
                    if uv_props == []:
                        uv_sig = 't0^0'
                        uv_props = ['1']

                    if diag_set['integrated_ct']:
                        uv_props.append('integratedctflag')

                    uv_conf = 'uvconf({},{},{},{})'.format(uv_sig, uv_info['taylor_order'], ','.join(external_momenta), '*'.join(uv_props))
                    uv_subgraphs.append('subgraph({}{},{})'.format(uv_info['graph_index'], 
                    (',' if len(uv_info['subgraph_indices']) > 0 else '') + ','.join(str(si) for si in uv_info['subgraph_indices']),
                    uv_conf))

                    if uv_subgraphs != []:
                        diag_set_uv_conf.append('*'.join(uv_subgraphs))

                # construct the map from the lmb to the cmb
                cmb_map = []
                n_loops = len(self.edges) - len(self.nodes) + 1
                for i in range(n_loops):
                    s = ''
                    for cc, cs in enumerate(loop_diag_set['cb_to_lmb'][i * n_loops:i * n_loops+n_loops]):
                        if cs != 0:
                            s += '{}c{}'.format('+' if cs == 1 else '-', cc + 1)
                            if cc < len(cut['cuts'][:-1]):
                                d = self.momenta_decomposition_to_string(([0] * n_loops, cut['cuts'][cc]['signature'][1]), False)
                                if d != '':
                                    # note the sign inversion
                                    s += '{}({})'.format('-' if cs == 1 else '+', d)
                    cmb_map.append(('k' + str(i + 1), s))

                cmb_map = 'cmb({})'.format(','.join(d for c in cmb_map for d in c))
                # store which momenta are LTD momenta
                conf = 'c0' if n_loops == len(cut['cuts'][:-1]) else ','.join('c' + str(i + 1) for i in range(len(cut['cuts'][:-1]), n_loops) )

                conf = 'conf({},{},{},{})'.format(diag_set['id'], cmb_map, conf, '*'.join(trans))
                if diag_set_uv_conf != []:
                    conf += '*\n    uv({})'.format('*'.join(diag_set_uv_conf))

                configurations.append(conf)
        
        # replace external momentum in all bubble vertices and edges with a dummy momentum
        # TODO: there is an awful lot of string manipulation here...
        mommap = []
        for i, ((bubble_edges, bubble_ext_edges), (cut_edge, bubble_derivative_edges)) in enumerate(bubble_to_cut.items()):
            ce = next(ee for ee in self.edges.values() if ee['name'] == cut_edge)

            if len(bubble_derivative_edges) == 0:
                # the bubble is a single vertex: update it
                bubble_v = next(v for v in self.nodes.values() if set(self.edges[eid]['name'] for eid in v['edge_ids']) == set(bubble_ext_edges))
                bubble_v['momenta'] = ('pbubble' + str(i) if bubble_v['momenta'][0] == ce['momentum'] else '-pbubble' + str(i),
                    'pbubble' + str(i) if bubble_v['momenta'][1] == ce['momentum'] else '-pbubble' + str(i))

            for edge in bubble_derivative_edges:
                e = next(ee for ee in self.edges.values() if ee['name'] == edge)

                signs = [se * sc for se, sc in zip(e['signature'][0] + e['signature'][1], ce['signature'][0] + ce['signature'][1]) if se * sc != 0]
                assert(len(set(signs)) == 1)
                e['momentum'] += '{}(pbubble{}-({}))'.format('+' if signs[0] == 1 else '-', i, ce['momentum'])

                # update the vertices attached to this edge
                ext_mom_vars = [emp for emp in ce['momentum'].replace('-', '+').split('+') if emp != '']
                ext_mom_signs = ['-' if ('-' + emp) in ce['momentum'] else '+' for emp in ext_mom_vars]
                for vi in e['vertices']:
                    v = self.nodes[vi]
                    momenta = list(v['momenta'])
                    for mi, (m, eid) in enumerate(zip(momenta, v['edge_ids'])):
                        # only update the momentum of the external (cut) edge or the current edge
                        if self.edges[eid]['name'] == edge or self.edges[eid]['name'] not in bubble_edges:
                            # find a momentum that is shared between the external momentum and the edge to determine the sign
                            (ext_var_in_mom, ext_sign_in_mom) = next((emp, ems) for emp, ems in zip(ext_mom_vars, ext_mom_signs) if emp in m)
                            sign = ('-' if ext_sign_in_mom == '+' else '+') if ('-' + ext_var_in_mom) in m else ('+' if ext_sign_in_mom == '+' else '-')
                            momenta[mi] += '{}(pbubble{}-({}))'.format(sign, i, ce['momentum'])
                    v['momenta'] = tuple(momenta)

            mommap.append('subs(pbubble{},{})'.format(i, ce['momentum']))

        self.replacement_rules = ''
        if len(mommap) > 0:
            self.replacement_rules = '*\n' + '*\n'.join(mommap)

        self.replacement_rules += '*\nconfigurations(\n  +{})'.format('\n  +'.join(configurations))


class FORMSuperGraphIsomorphicList(list):
    """ Container class for a list of FROMSuperGraph with the same denominator structure"""
    
    def __init__(self, graph_list):
        """ Instantiates a list of FORMSuperGraphs from a list of either
        FORMSuperGraph instances or LTD2SuperGraph instances."""

        for g in graph_list:
            if isinstance(g,FORMSuperGraph):
                self.append(g)
            else:
                self.append(FORMSuperGraph.from_LTD2SuperGraph(g))

    def generate_numerator_form_input(self, additional_overall_factor='',):
        return '+'.join(g.generate_numerator_form_input(additional_overall_factor) for g in self)

    def generate_numerator_functions(self, additional_overall_factor='', output_format='c', workspace=None, FORM_vars=None, active_graph=None,process_definition=None):
        """ Use form to plugin Feynman Rules and process the numerator algebra so as
        to generate a low-level routine in file_path that encodes the numerator of this supergraph."""

        _MANDATORY_FORM_VARIABLES = ['SGID','NINITIALMOMENTA','NFINALMOMENTA','SELECTEDEPSILONORDER','UVRENORMFINITEPOWERTODISCARD']

        if FORM_vars is None:
            raise FormProcessingError("FORM_vars must be supplied when calling generate_numerator_functions.")
        FORM_vars = dict(FORM_vars)

        if active_graph is None:
            characteristic_super_graph = self[0]
        else:
            characteristic_super_graph = active_graph

        if 'NINITIALMOMENTA' not in FORM_vars:
            n_incoming = sum([1 for edge in characteristic_super_graph.edges.values() if edge['type'] == 'in'])
            FORM_vars['NINITIALMOMENTA'] = n_incoming
        if 'NFINALMOMENTA' not in FORM_vars:
            n_loops = len(characteristic_super_graph.edges) - len(characteristic_super_graph.nodes) + 1
            FORM_vars['NFINALMOMENTA'] = n_loops

        if FORM_vars is None or not all(opt in FORM_vars for opt in _MANDATORY_FORM_VARIABLES):
            raise FormProcessingError("The following variables must be supplied to FORM: %s"%str(_MANDATORY_FORM_VARIABLES))

        FORM_vars.update(FORM_processing_options['extra-options'])

        i_graph = int(FORM_vars['SGID'])

        # write the form input to a file
        if active_graph is None:
            form_input = self.generate_numerator_form_input(additional_overall_factor)
        else:
            form_input = characteristic_super_graph.generate_numerator_form_input(additional_overall_factor)

        if workspace is None:
            selected_workspace = FORM_workspace
            FORM_source = pjoin(plugin_path,"numerator.frm")
        else:
            selected_workspace = workspace
            shutil.copy(pjoin(plugin_path,"numerator.frm"),pjoin(selected_workspace,'numerator.frm'))
            shutil.copy(pjoin(plugin_path,"diacolor.h"),pjoin(selected_workspace,'diacolor.h'))
            FORM_source = pjoin(selected_workspace,'numerator.frm')

        with open(pjoin(selected_workspace,'input_%d.h'%i_graph), 'w') as f:
            f.write('L F = {};'.format(form_input))

        FORM_cmd = ' '.join([
                FORM_processing_options["FORM_path"],
                ]+
                [ '-D %s=%s'%(k,v) for k,v in FORM_vars.items() ] +
                [ '-M', '-l'] +
                [ FORM_source, ]
        )
        r = subprocess.run(FORM_cmd,
            shell=True,
            cwd=selected_workspace,
            capture_output=True)
        # TODO understand why FORM sometimes returns return_code 1 event though it apparently ran through fine
        if r.returncode != 0 and not os.path.isfile(pjoin(selected_workspace,'out_%d.proto_c'%i_graph)):
            raise FormProcessingError("FORM processing failed with error:\n%s\nFORM command to reproduce:\ncd %s; %s"%(
                    r.stdout.decode('UTF-8'),
                    selected_workspace, FORM_cmd
            ))

        # return the code for the numerators
        if not os.path.isfile(pjoin(selected_workspace,'out_%d.proto_c'%i_graph)):
            raise FormProcessingError(
                    ( "FORM failed to produce an output for super graph ID=%d. Output file not found at '%s'."%
                                                                    (i_graph,pjoin(selected_workspace,'out_%d.proto_c'%i_graph)))+
                "\nFORM command to reproduce:\ncd %s; %s"%(selected_workspace,FORM_cmd)
            )

        with open(pjoin(selected_workspace,'out_%d.proto_c'%i_graph), 'r') as f:
            num_code = f.read()

        return num_code

    def to_dict(self, file_path=None):
        """ Store that into a dict."""
        to_dump = [g.to_dict() for g in self]
        if file_path:
            open(file_path,'w').write(pformat(to_dump))
        else:
            return to_dump
    
    def multiplicity_factor(self,iso_id, workspace, form_source):
        output_match = re.compile(r'isoF=(.*?);')
        multiplicity = 0
        reference = self[0].generate_numerator_form_input('', only_algebra=True)
        FORM_vars = {}
        FORM_vars['SGID'] = iso_id
        FORM_vars['ID0'] = 0
        for i_graph, g in enumerate(self):
            mapped = g.generate_numerator_form_input('', only_algebra=True)
            with open(pjoin(workspace,'iso_check_{}_{}_{}.frm'.format(iso_id, 0, i_graph+1)), 'w') as f:
                FORM_vars['IDn'] = i_graph+1
                f.write("L F1 = %s;\n"%reference)
                f.write("L F2 = %s;\n"%mapped)
            
            r = subprocess.run(' '.join([
                FORM_processing_options["tFORM_path"],
                ]+
                ['-w%d' % FORM_processing_options['cores']]+
                [ '-D %s=%s'%(k,v) for k,v in FORM_vars.items() ]+
                [ form_source, ]
            ),
            shell=True,
            cwd=workspace,
            capture_output=True)
            if r.returncode != 0:
                raise FormProcessingError("FORM processing failed with error:\n%s"%(r.stdout.decode('UTF-8')))

            output = r.stdout.decode('UTF-8').replace(' ','').replace('\n','')
            factor = re.sub(r'rat\(([-0-9]+),([-0-9]+)\)', r'(\1)/(\2)', output_match.findall(output)[0])

            if "rat" in factor:
                raise FormProcessingError("Multiplicity not found: {} / {} is not rational (iso_check_%(SGID)d_%(ID0)d_%(IDn)d)".format(self[0].name,g.name )%FORM_vars)
            else:
                multiplicity += eval(factor)

            #logger.info("{} = {} * {}".format(self[0].name, factor, g.name ))
        return multiplicity    

    def generate_squared_topology_files(self, root_output_path, model, process_definition, n_jets, numerator_call, final_state_particle_ids=(), jet_ids=None, workspace=None, bar=None,
            integrand_type=None):
        for i, g in enumerate(self):
            # Now we generate the squared topology only for the first isomorphic graph
            # to obtain the replacement rules for the bubble.
            # The other elements of the isomorphic set are going to contribute only at the 
            # numerator 
            if i==0:
                r = g.generate_squared_topology_files(root_output_path, model, process_definition, n_jets, numerator_call, 
                                final_state_particle_ids, jet_ids=jet_ids, write_yaml=i==0, workspace=workspace, bar=bar, integrand_type=integrand_type)
            else:
                g.replacement_rules = self[0].replacement_rules
        #print(r)
        return r

class FORMSuperGraphList(list):
    """ Container class for a list of FORMSuperGraphIsomorphicList."""

    extension_names = {'c': 'c'}

    def __init__(self, graph_list, name='SGL'):
        """ Instantiates a list of FORMSuperGraphs from a list of either
        FORMSuperGraphIsomorphicList instances, FORMSuperGraph, or LTD2SuperGraph instances."""

        self.name = name

        for g in graph_list:
            if isinstance(g, FORMSuperGraph):
                self.append(FORMSuperGraphIsomorphicList([g]))
            elif isinstance(g, FORMSuperGraphIsomorphicList):
                self.append(g)
            else:
                self.append(FORMSuperGraphIsomorphicList([FORMSuperGraph.from_LTD2SuperGraph(g)]))

    @classmethod
    def from_squared_topology(cls, edge_map_lin, name, incoming_momentum_names, model, loop_momenta_names=None, particle_ids={},benchmark_result=0.0):
        vertices = [v for e in edge_map_lin for v in e[1:]]

        topo_generator = LTD.ltd_utils.TopologyGenerator(edge_map_lin, {})
        topo_generator.generate_momentum_flow(loop_momenta_names)
        sig = topo_generator.get_signature_map()

        edges = {i: {
            'PDG': 1337 if e[0] not in particle_ids else particle_ids[e[0]],
            'indices': (1 + i * 2,) if vertices.count(e[1]) == 1 or vertices.count(e[2]) == 1 else (1 + i * 2, 1 + i * 2 + 1),
            'signature': [sig[e[0]][0],
                    [ i+o for i,o in zip(sig[e[0]] [1][:len(sig[e[0]] [1])//2], sig[e[0]] [1][len(sig[e[0]] [1])//2:])]],
            'name': e[0],
            'type': 'in' if e[0] in incoming_momentum_names else ('out' if vertices.count(e[1]) == 1 or vertices.count(e[2]) == 1 else 'virtual'),
            'vertices': tuple(e[1:]),
            } for i, e in enumerate(edge_map_lin)
        }

        nodes = {v: {
            'PDGs': tuple(e['PDG'] if v == e['vertices'][1] or e['PDG'] == 1337 else model.get_particle(e['PDG']).get_anti_pdg_code() for e in edges.values() if v in e['vertices']),
            'edge_ids': tuple(ei for ei, e in edges.items() if v in e['vertices']),
            'indices': tuple(e['indices'][0] if v == e['indices'][0] or len(e['indices']) == 1 else e['indices'][1] for e in edges.values() if v in e['vertices']),
            'momenta': tuple('DUMMY' for e in edges.values() if v in e['vertices']),
            'vertex_id': -1 if vertices.count(v) == 1 else 100,
            } for v in set(vertices)
        }

        for eid, e in edges.items():
            e['momentum'] = FORMSuperGraph.momenta_decomposition_to_string(e['signature'], set_outgoing_equal_to_incoming=False)
            neg_mom =  FORMSuperGraph.momenta_decomposition_to_string([[-s for s in sp] for sp in e['signature']], set_outgoing_equal_to_incoming=False)

            # TODO: why was this here?
            #if len(e['vertices']) == 1:
            #    continue

            for i, vi in enumerate(e['vertices']):
                e_index = nodes[vi]['edge_ids'].index(eid)
                mom = list(nodes[vi]['momenta'])
                mom[e_index] = neg_mom if i == 0 else e['momentum']
                nodes[vi]['momenta'] = mom

        # set the correct particle ordering for all the edges
        for n in nodes.values():
            if len(n['PDGs']) > 1 and all(nn != 1337 for nn in n['PDGs']):
                edge_order = FORMSuperGraph.sort_edges(model, [{'PDG': pdg, 'index': i} for i, pdg in enumerate(n['PDGs'])])
                for g in ('PDGs', 'indices', 'momenta', 'edge_ids'):
                    n[g] = tuple(n[g][eo['index']] for eo in edge_order)


        form_graph = FORMSuperGraph(name=name, edges = edges, nodes=nodes,
            overall_factor='-1', # there is an overall factor of -1 wrt Forcer, presumable due to the Wick rotation?
            multiplicity = 1,
            benchmark_result=benchmark_result
        )
        return FORMSuperGraphList([FORMSuperGraphIsomorphicList([form_graph])], name=name + '_set')


    @classmethod
    def from_dict(cls, dict_file_path, first=None, merge_isomorphic_graphs=False, verbose=False, model = None, workspace=None, cuts=None):
        """ Creates a FORMSuperGraph list from a dict file path."""
        from pathlib import Path
        p = Path(dict_file_path)
        sys.path.insert(0, str(p.parent))
        m = __import__(p.stem)

        logger.info("Imported {} supergraphs.".format(len(m.graphs)))

        # Filter specific graphs by name 
        #filter_graphs = ['SG_QG3','SG_QG4']
        #m.graphs = [ g for (g,name) in zip(m.graphs, m.graph_names) if name in filter_graphs]
        #m.graph_names = ['SG_MG3','SG_QG4']
        #m.graph_names = [name for name in m.graph_names if name in filter_graphs ]

        # Now convert the vertex names to be integers according to QGRAF format:
        for i, g in enumerate(m.graphs):
            if isinstance(list(g['nodes'].keys())[0], str):
                new_nodes ={}
                node_names={}
                for i_node, (n_key, n) in enumerate(g['nodes'].items()):
                    node_names[n_key]=i_node+1
                    new_nodes[i_node+1] = n
                for n in g['nodes'].values():
                    n['edge_ids'] = tuple([ (node_names[edge_key[0]],node_names[edge_key[1]],edge_key[2]) for edge_key in n['edge_ids'] ])
                new_edges = {}
                for edge_key, edge in g['edges'].items():
                    edge['vertices']=(node_names[edge_key[0]],node_names[edge_key[1]])
                    new_edges[(node_names[edge_key[0]],node_names[edge_key[1]],edge_key[2])]=edge
                g['nodes'] = new_nodes
                g['edges'] = new_edges

        # Check for the rooting of the loop momenta from the QGRAF OUTPUT
        # If necessary flip the flow direction
        for name, g in zip(m.graph_names, m.graphs):
            n_loops = len(g['edges']) - len(g['nodes']) + 1
            flows = []
            for lms in {e['momentum'] for e in g['edges'].values() if 'p' not in e['momentum'] and e['momentum'].count('k') == 1}:
                lm = lms.replace('-', '')
                assert(re.match('^k[0-9]+$', lm))
                if lm in flows or '-'+lm in flows:
                    continue
                flows += [lms]
            assert(len(flows) == n_loops)
            
            for lms in flows:
                if '-' not in lms:
                    continue
                lm = lms.replace('-', '')
                if verbose:
                    logger.info("QGraf remap loop momentum {} for graph {}: {} -> {}".format(lm, name, lm, lms))
                subs = [('^-'+lm, 'TMP'),('\+'+lm, '-TMP'), ('-'+lm, '+TMP'), ('^'+lm, '-TMP'), ('TMP', lm)]
                for e in g['edges'].values():
                    for sub in subs:
                        e['momentum'] = re.sub(*sub, e['momentum'])
                for n in g['nodes'].values():
                    # skip if external
                    if n['vertex_id'] < 0:
                        continue
                    moms = list(n['momenta'])
                    for sub in subs:
                        for n_mom in range(len(moms)):
                            moms[n_mom] = re.sub(*sub, moms[n_mom])
                    n['momenta'] = tuple(moms)
        
        full_graph_list = []
        for i, g in enumerate(m.graphs):
            if hasattr(m,'graph_names'):
                graph_name=m.graph_names[i]
            else:
                graph_name=p.stem + '_' + str(i)
            # convert to FORM supergraph
            form_graph = FORMSuperGraph(name=graph_name, edges = g['edges'], nodes=g['nodes'], 
                        overall_factor=g['overall_factor'], multiplicity = g.get('multiplicity',1) )
            form_graph.derive_signatures()
            full_graph_list.append(form_graph)

        # Adjust edge orientation of all repeated edges so that they have identical signatures
        for g in full_graph_list:

            # First detect edges that need flipping
            all_signatures = {}
            for e_key, e in g.edges.items():
                pos_signature = tuple(tuple(s for s in sp) for sp in e['signature'])
                neg_signature = tuple(tuple(-s for s in sp) for sp in e['signature'])
                if pos_signature in all_signatures:
                    all_signatures[pos_signature].append((e_key,+1))
                elif neg_signature in all_signatures:
                    all_signatures[neg_signature].append((e_key,-1))
                else:
                    all_signatures[pos_signature] = [(e_key,+1)]

            for sig, edges in all_signatures.items():
                for edge_key, sign in edges:
                    if sign==+1:
                        continue
                    # For now we only allow to flip particles that are not fermionic otherwise 
                    # special care must be taken when flipping an edge and normally fermionic lines
                    # should already have been aligned at this stage.
                    flipped_part = model.get_particle(g.edges[edge_key]['PDG'])
                    if flipped_part.is_fermion():
                        raise FormProcessingError("A *fermionic* repeated edge with oppoosite signatures was found. This should happen for bosons only.")

                    g.edges[edge_key]['PDG'] = flipped_part.get_anti_pdg_code()
                    g.edges[edge_key]['signature'] = [[-s for s in sp] for sp in g.edges[edge_key]['signature']]
                    g.edges[edge_key]['vertices'] = [
                        g.edges[edge_key]['vertices'][1],
                        g.edges[edge_key]['vertices'][0],
                    ]

            # Now adust the string momenta of edges and nodes accordingly.            
            g.impose_signatures()


        if first is not None:
            logger.info("Taking first {} supergraphs.".format(first))
            full_graph_list = full_graph_list[:first]

        if not cuts is None:
            graph_filtered = {'DUMP':[], 'KEEP':[]}
            with progressbar.ProgressBar(prefix='Filter SG with valid cuts: {variables.keep}\u2713  {variables.drop}\u2717 : ', max_value=len(full_graph_list),variables={'keep': '0', 'drop':'0'}) as bar:
                for sgid, graph in enumerate(full_graph_list):
                    if graph.filter_valid_cuts(cuts):
                        graph_filtered["KEEP"] += [graph]
                    else: 
                        graph_filtered["DUMP"] += [graph]
                    bar.update(sgid+1)
                    bar.update(keep=len(graph_filtered['KEEP']))
                    bar.update(drop=len(graph_filtered['DUMP']))
            #for k,v in graph_filtered.items():
            #    print(k,":")
            #    print("\t",np.array([g.name for g in v]))
            full_graph_list[:] = graph_filtered['KEEP']
            logger.info("\033[1mRemoved {} graphs with no valid cuts\033[0m".format(len(graph_filtered['DUMP'])))

        iso_groups = []

        import time
        import sympy as sp
        # group all isomorphic graphs
        n_externals = max(len([1 for e in graph.edges.values() if e['type']=='in' or e['type']=='out']) for graph in full_graph_list)
        if model is None:
            pdg_primes = {pdg : sp.prime(i + n_externals + 1) for i, pdg in enumerate([1,2,3,4,5,6,11,12,13,21,22,25,82])}
        else:
            pdg_primes = {pdg : sp.prime(i + n_externals + 1) for i, pdg in enumerate([p['pdg_code'] for p in model['particles']])}
        
        for graph in full_graph_list:
            g = igraph.Graph()
            g.add_vertices(len(graph.nodes))
            undirected_edges = set()

            for e in graph.edges.values():
                undirected_edges.add(tuple(sorted(e['vertices'])))

            edge_colors = []
            ext_id = 0
            for ue in undirected_edges:
                e_color = 1
                for e in graph.edges.values():
                    if tuple(sorted(e['vertices'])) == ue:
                        #TODO: Reserve the first #externals primes for external edges
                        # with the current implementation it still allows swap 
                        # of final/initial state edges
                        if e['type'] == 'in':
                            ext_id += 1
                            e_color *= sp.prime(ext_id)
                        else:
                            e_color *= pdg_primes[abs(e['PDG'])]
                edge_colors.append(e_color)
                g.add_edges([tuple(sorted([x - 1 for x in ue]))])
            
            for (ref_graph, ref_colors), graph_list in iso_groups:
                v_maps = ref_graph.get_isomorphisms_vf2(g, edge_color1=ref_colors, edge_color2=edge_colors)
                if v_maps != [] and merge_isomorphic_graphs:
                    v_map = v_maps[0]
                    if verbose: 
                        logger.info("\033[1m{} <- {}\033[0m".format(graph_list[0].name,graph.name))
                    
                    # map the signature from the new graph to the reference graph using the vertex map
                    edge_map = [next(re.index for re in ref_graph.es if (re.source, re.target) == (v_map[e.source], v_map[e.target]) or \
                        (re.source, re.target) == (v_map[e.target], v_map[e.source])) for e in g.es]

                    # go through all loop momenta and make the matrix
                    # and the shifts
                    mat = []
                    shifts = []
                    selected_edges = []
                    n_loops = len(graph_list[0].edges) - len(graph_list[0].nodes) + 1
                    for loop_var in range(n_loops):
                        # find the edge that has (+-)k`loop_var`
                        # TODO: Fix for higher degeneracy
                        orig_edge = next(ee for ee in graph_list[0].edges.values() if all(s == 0 for s in ee['signature'][1]) and \
                            sum(abs(s) for s in ee['signature'][0]) == 1 and ee['signature'][0][loop_var] != 0)
                        orig_edge_in_ref = next(ee for ee in ref_graph.es if tuple(sorted(orig_edge['vertices'])) == (ee.source + 1, ee.target + 1))

                        # now get the corresponding edge(s) in the new graph
                        en = g.es[edge_map.index(orig_edge_in_ref.index)]
                        new_edges = [ee for ee in graph.edges.values() \
                            if (ee['vertices'] == (en.source + 1, en.target + 1) or ee['vertices'] == (en.target + 1, en.source + 1)) and abs(ee['PDG'])== abs(orig_edge['PDG']) and ee not in selected_edges]
                       
                        # in the case of a multi-edge, we simply take the one with the shortest momentum string
                        new_edge = sorted(new_edges, key=lambda x: len(x['momentum']))[0]
                        selected_edges.append(new_edge)
                        new_edge['name'] += 'LMB'
                        
                        orientation_factor = 1
                        if model.get_particle(abs(new_edge['PDG']))['self_antipart']:
                            if v_map[new_edge['vertices'][0]-1]+1 != orig_edge['vertices'][0]:
                                orientation_factor = -1
                        else:
                            match_direction = v_map[new_edge['vertices'][0]-1]+1 == orig_edge['vertices'][0]
                            match_particle = new_edge['PDG'] == orig_edge['PDG']
                            if not match_direction:# != match_particle:
                                orientation_factor = -1

                        shifts.append([-orientation_factor * s * orig_edge['signature'][0][loop_var] for s in new_edge['signature'][1]])
                        mat.append([orientation_factor * s * orig_edge['signature'][0][loop_var] for s in new_edge['signature'][0]])

                    # Map all signatures and momentum labels of the new graph
                    shifts = np.array(shifts).T
                    mat = np.array(np.matrix(mat).T.I)
                    for e in graph.edges.values():
                        # skip external edges
                        if e["type"] != "virtual":
                            continue

                        new_sig = mat.dot(np.array(e['signature'][0]))
                        new_shift = e['signature'][1]+shifts.dot(new_sig)
                        
                        e['signature'] = (list(new_sig), list(new_shift))
                        e['momentum'] = graph.momenta_decomposition_to_string(e['signature'], set_outgoing_equal_to_incoming=False)
                    
                    if verbose:
                        parser = re.compile(r'(^|\+|-)(k|p)(\d*)')
                        n_incoming = sum([1 for edge in graph.edges.values() if edge['type'] == 'in'])
                        n_loops = len(graph.edges) - len(graph.nodes) + 1
                        old_momenta =['k%d'%(i+1) for i in range(n_loops)] 
                        new_momenta =[] 
                        for momentum in old_momenta:
                            parsed = [e.groups() for e in parser.finditer(momentum)]
                            signature = ([0 for _ in range(n_loops)], [0 for _ in range(n_incoming)])
                            for pa in parsed:
                                if pa[1] == 'p':
                                    signature[1][int(pa[2]) - 1] = 1 if pa[0] == '' or pa[0] == '+' else -1
                                else:
                                    signature[0][int(pa[2]) - 1] = 1 if pa[0] == '' or pa[0] == '+' else -1
                            new_sig = mat.dot(np.array(signature[0]))
                            new_shift = signature[1]+shifts.dot(new_sig)
                            new_momenta.append(graph.momenta_decomposition_to_string((list(new_sig), list(new_shift)), set_outgoing_equal_to_incoming=False))
                        logger.info("Map: {} -> {}".format(tuple(old_momenta), tuple(new_momenta)))
                    
                    
                    # Map vertices
                    for v in graph.nodes.values():
                        if v['vertex_id']<0:
                            continue
                        # skip external edge        
                        n_incoming = sum([1 for edge in graph.edges.values() if edge['type'] == 'in'])
                        n_loops = len(graph.edges) - len(graph.nodes) + 1

                        # parse momentum
                        parser = re.compile(r'(^|\+|-)(k|p)(\d*)')
                        new_momenta = []
                        for momentum in v['momenta']:
                            parsed = [e.groups() for e in parser.finditer(momentum)]
                            signature = ([0 for _ in range(n_loops)], [0 for _ in range(n_incoming)])
                            for pa in parsed:
                                if pa[1] == 'p':
                                    signature[1][int(pa[2]) - 1] = 1 if pa[0] == '' or pa[0] == '+' else -1
                                else:
                                    signature[0][int(pa[2]) - 1] = 1 if pa[0] == '' or pa[0] == '+' else -1
                            new_sig = mat.dot(np.array(signature[0]))
                            new_shift = signature[1]+shifts.dot(new_sig)
                            new_momenta.append(graph.momenta_decomposition_to_string((list(new_sig), list(new_shift)), set_outgoing_equal_to_incoming=False))
                        v['momenta'] = tuple(new_momenta)
                    graph_list.append(graph)
                    break
            else:
                iso_groups.append(((g, edge_colors), [graph]))

        logger.info("\033[1m{} unique supergraphs\033[m".format(len(iso_groups)))

        FORM_iso_sg_list = FORMSuperGraphList([FORMSuperGraphIsomorphicList(iso_group) for _, iso_group in iso_groups], name=p.stem)

        if workspace is not None:
            selected_workspace = workspace
            shutil.copy(pjoin(plugin_path,"multiplicity.frm"),pjoin(selected_workspace,'multiplicity.frm'))
            shutil.copy(pjoin(plugin_path,"numerator.frm"),pjoin(selected_workspace,'numerator.frm'))
            shutil.copy(pjoin(plugin_path,"diacolor.h"),pjoin(selected_workspace,'diacolor.h'))
            FORM_source = pjoin(selected_workspace,'multiplicity.frm')

            with progressbar.ProgressBar(
                prefix = 'Processing Isomorphic sets (ISO #{variables.iso_id}/%d, Zeros: {variables.zero_multiplicity_count})'%len(iso_groups),
                max_value=len(iso_groups),
                variables = {'iso_id' : '0', 'zero_multiplicity_count' : '0'}
                ) as bar:
            
                zero_multiplicity_count = 0
                for iso_id, iso_graphs in enumerate(FORM_iso_sg_list):
                    bar.update(iso_id='%d'%iso_id, zero_multiplicity_count='%d'%zero_multiplicity_count)
                    multiplicity = iso_graphs.multiplicity_factor(iso_id, selected_workspace, FORM_source)
                    if multiplicity == 0:
                        zero_multiplicity_count += 1
                        iso_graphs[:] = []
                    else:
                        iso_graphs[:] = iso_graphs[:1]
                        iso_graphs[0].multiplicity = multiplicity
                    bar.update(iso_id + 1)
                bar.update(zero_multiplicity_count='%d'%zero_multiplicity_count)

                for iso_id in reversed(range(len(FORM_iso_sg_list))):
                    if FORM_iso_sg_list[iso_id] == []:
                        del FORM_iso_sg_list[iso_id]
                
            if zero_multiplicity_count > 0 :
                logger.info("%d isomorphic sets have 0 multiplicity -> they are dropped!"%zero_multiplicity_count)
            
            #for iso_id, iso_graphs in enumerate(FORM_iso_sg_list):
            #    g = iso_graphs[0]
            #    print("#{}\t{:10}: {}".format(iso_id, g.name, g.multiplicity))


        # Generate additional copies for different LMB of the representative graph of each isomorphic set.
        # Note that depending on the option FORM_processing['representative_lmb'], the function below can also modify the LMB of the representative sg.
        for FORM_iso_sg in FORM_iso_sg_list:
            FORM_iso_sg[0].generate_additional_LMBs()

        return FORM_iso_sg_list

    def to_dict(self, file_path):
        """ Outputs the FORMSuperGraph list to a Python dict file path."""

        # TODO: Dump to Python dict
        pass

    def generate_integrand_functions(self, root_output_path, additional_overall_factor='',
                                    params={}, output_format='c', workspace=None, header="", integrand_type=None,process_definition=None):
        header_map = {'header': header}
        """ Generates optimised source code for the graph numerator in several
        files rooted in the specified root_output_path."""

        if len(self)==0:
            raise FormProcessingError("Cannot generat numerators for an empty list of supergraphs.")

        if output_format not in self.extension_names:
            raise FormProcessingError("This FORMSuperGraphList instance requires at least one entry for generating numerators.")

        # add all numerators in one file and write the headers
        numerator_header = """#include <tgmath.h>
#include <quadmath.h>
#include <signal.h>
#include "%(header)snumerator.h"
"""

        var_pattern = re.compile(r'Z\d*_')
        input_pattern = re.compile(r'lm(\d*)')
        conf_exp = re.compile(r'conf\(([^)]*)\)\n')
        return_exp = re.compile(r'return ([^;]*);\n')
        float_pattern = re.compile(r'((\d+\.\d*)|(\.\d+))')

        # TODO: multiprocess this loop
        all_numerator_ids = []
        with progressbar.ProgressBar(
            prefix = 'Processing numerators (graph #{variables.i_graph}/%d, LMB #{variables.i_lmb}/{variables.max_lmb}) with FORM ({variables.timing} ms / supergraph) : '%len(self),
            max_value=len(self),
            variables = {'timing' : '0', 'i_graph' : '0', 'i_lmb': '0', 'max_lmb': '0'}
            ) as bar:
            total_time = 0.

            for i_graph, graph in enumerate(self):
                graph.is_zero = True

                graphs_to_process = []
                if isinstance(graph[0].additional_lmbs,list) and graph[0].additional_lmbs != []:
                    graphs_to_process.append( (0,i_graph,graph[0]) )
                    graphs_to_process.extend([(g.additional_lmbs,i_graph,g) for _,_,_,g in graph[0].additional_lmbs])
                else:
                    # By setting the active graph to None we will then sum overall members of the iso set.
                    graphs_to_process.append( (0,i_graph, None) )
                bar.update(max_lmb='%d'%len(graphs_to_process))
                for i_lmb, i_g, active_graph in graphs_to_process:
                    bar.update(i_graph='%d'%(i_graph+1), i_lmb='%d'%(i_lmb+1))
                    i = i_lmb*FORM_processing_options['FORM_call_sig_id_offset_for_additional_lmb']+i_g
                    all_numerator_ids.append(i)
                    time_before = time.time()

                    with open(pjoin(root_output_path, 'workspace', 'out_integrand_{}.proto_c'.format(i))) as f:
                        num = f.read()

                    total_time += time.time()-time_before
                    num = num.replace('i_', 'I')
                    num = input_pattern.sub(r'lm[\1]', num)
                    num = num.replace('\nZ', '\n\tZ') # nicer indentation

                    confs = []
                    integrand_main_code = ''
                    conf_secs = num.split('#CONF')
                    for conf_sec in conf_secs[1:]:
                        conf_sec = conf_sec.replace("#CONF\n", '')

                        # parse the configuration id
                        conf = list(conf_exp.finditer(conf_sec))[0].groups()[0].split(',')
                        conf_id = int(conf[0])
                        confs.append(conf_id)

                        conf_sec = conf_exp.sub('', conf_sec)

                        # parse the constants
                        mom_map, constants, loops = graph[0].integrand_info[conf_id]
                        const_secs = conf_sec.split('#CONSTANTS')[1]
                        const_secs = const_secs.replace('\n', '')
                        const_code = '*'.join('\t({}-{})'.format(e, mass) for e, (_, mass) in zip(const_secs.split(','), constants) if e != '')

                        # parse the energies
                        energy_secs = conf_sec.split('#ENERGIES')[1]
                        energy_secs = energy_secs.replace('\n', '')
                        energy_code = '\n'.join('\tdouble complex E{} = sqrt({}+{});'.format(i, e, mass) for i, (e, mass) in enumerate(zip(energy_secs.split(','), mom_map)) if e != '')
                        if integrand_type == "PF":
                            const_code = const_code + ('*' if len(const_code) > 0 and len(mom_map) > 0 else '') + '*'.join('2.*E{}'.format(i) for i in range(len(mom_map)))
                        if const_code == '':
                            const_code = '1'

                        # parse the denominators
                        denom_secs = conf_sec.split('#ELLIPSOIDS')[1]
                        denom_secs = denom_secs.replace('\n', '')
                        denom_code = '\n'.join('\tdouble complex invd{} = 1./({});'.format(i, d) for i, d in enumerate(denom_secs.split(',')) if d != '')

                        conf_sec = conf_sec.split('#ELLIPSOIDS')[-1]
                        returnval = list(return_exp.finditer(conf_sec))[0].groups()[0]
                        conf_sec = return_exp.sub('return pow(2.*pi*I,{})/({})*({});\n'.format(loops, const_code, returnval), conf_sec)

                        # collect all temporary variables
                        temp_vars = list(sorted(set(var_pattern.findall(conf_sec))))

                        if len(temp_vars) > 0:
                            graph.is_zero = False

                        main_code = '{}\n{}\n{}'.format(energy_code, denom_code, conf_sec)
                        main_code = main_code.replace('logmUV', 'log(mUV*mUV)').replace('logmu' , 'log(mu*mu)').replace('logmt' , 'log(mass_t*mass_t)')
                        integrand_main_code += '\nstatic inline double complex %(header)sevaluate_{}_{}(double complex lm[], double complex params[]) {{\n\t{}\n{}}}'.format(i, conf_id,
                            'double complex {};'.format(','.join(temp_vars)) if len(temp_vars) > 0 else '', main_code
                        )

                        main_code_f128 = main_code.replace('pow', 'cpowq').replace('sqrt', 'csqrtq').replace('log', 'clogq').replace('pi', 'M_PIq').replace('double complex', '__complex128')
                        main_code_f128 = float_pattern.sub(r'\1q', main_code_f128)
                        integrand_main_code += '\n' + '\nstatic inline __complex128 %(header)sevaluate_{}_{}_f128(__complex128 lm[], __complex128 params[]) {{\n\t{}\n{}}}'.format(i, conf_id,
                            '__complex128 {};'.format(','.join(temp_vars)) if len(temp_vars) > 0 else '', main_code_f128
                        )

                    integrand_main_code += \
"""
double complex %(header)sevaluate_{}(double complex lm[], double complex params[], int conf) {{
   switch(conf) {{
{}
    }}
}}
""".format(i,
        '\n'.join(
        ['\t\tcase {}: return %(header)sevaluate_{}_{}(lm, params);'.format(conf, i, conf) for conf in sorted(confs)] +
        (['\t\tdefault: return 0.;']) # ['\t\tdefault: raise(SIGABRT);'] if not graph.is_zero else 
        ))


                    integrand_main_code += \
"""
__complex128 %(header)sevaluate_{}_f128(__complex128 lm[], __complex128 params[], int conf) {{
   switch(conf) {{
{}
    }}
}}
""".format(i,
        '\n'.join(
        ['\t\tcase {}: return %(header)sevaluate_{}_{}_f128(lm, params);'.format(conf, i, conf) for conf in sorted(confs)] +
        (['\t\tdefault: return 0.;']) # ['\t\tdefault: raise(SIGABRT);'] if not graph.is_zero else 
        ))

                    bar.update(timing='%d'%int((total_time/float(i_graph+1))*1000.0))
                    bar.update(i_graph+1)

                    writers.CPPWriter(pjoin(root_output_path, '%(header)sintegrand{}.c'%header_map).format(i)).write((numerator_header + integrand_main_code)%header_map)

        numerator_code = \
"""
#include <tgmath.h>
#include <quadmath.h>
#include <signal.h>

{}
{}

double complex %(header)sevaluate(double complex lm[], double complex params[], int diag, int conf) {{
    switch(diag) {{
{}
    }}
}}

__complex128 %(header)sevaluate_f128(__complex128 lm[], __complex128 params[], int diag, int conf) {{
    switch(diag) {{
{}
    }}
}}
""".format(
    '\n'.join('double complex %(header)sevaluate_{}(double complex[], double complex[], int conf);'.format(i) for i in all_numerator_ids),
    '\n'.join('__complex128 %(header)sevaluate_{}_f128(__complex128[], __complex128[], int conf);'.format(i) for i in all_numerator_ids),
    '\n'.join(
    ['\t\tcase {}: return %(header)sevaluate_{}(lm, params, conf);'.format(i, i) for i in all_numerator_ids]+
    ['\t\tdefault: raise(SIGABRT);']
    ),
    '\n'.join(
    ['\t\tcase {}: return %(header)sevaluate_{}_f128(lm, params, conf);'.format(i, i) for i in all_numerator_ids]+
    ['\t\tdefault: raise(SIGABRT);']
    )
    )
        writers.CPPWriter(pjoin(root_output_path, '%(header)sintegrand.c'%header_map)).write(numerator_code%header_map)


    def generate_numerator_functions(self, root_output_path, additional_overall_factor='', params={}, 
                                    output_format='c', workspace=None, header="", integrand_type=None,process_definition=None):
        header_map = {'header': header}
        """ Generates optimised source code for the graph numerator in several
        files rooted in the specified root_output_path."""

        if len(self)==0:
            raise FormProcessingError("Cannot generat numerators for an empty list of supergraphs.")

        if output_format not in self.extension_names:
            raise FormProcessingError("This FORMSuperGraphList instance requires at least one entry for generating numerators.")

        # add all numerators in one file and write the headers
        numerator_header = """#include <tgmath.h>
#include <quadmath.h>
#include <signal.h>
#include "%(header)snumerator.h"
"""

        var_pattern = re.compile(r'Z\d*_')
        input_pattern = re.compile(r'lm(\d*)')
        energy_exp = re.compile(r'f\(([^)]*)\)\n')
        split_number = re.compile(r'\\\n\s*')
        return_exp = re.compile(r'return ([^;]*);\n')
        float_pattern = re.compile(r'((\d+\.\d*)|(\.\d+))')

        FORM_vars={
            'SELECTEDEPSILONORDER':'%d'%FORM_processing_options['selected_epsilon_UV_order']
        }

        if FORM_processing_options['renormalisation_finite_terms']=='together':
            # Keep all terms, so set the discared power to 0
            FORM_vars['UVRENORMFINITEPOWERTODISCARD'] = 0
        elif FORM_processing_options['renormalisation_finite_terms']=='only':
            # Discard all terms not prefixed with UVRenormFinite, so set discarded power to -1
            FORM_vars['UVRENORMFINITEPOWERTODISCARD'] = -1
        elif FORM_processing_options['renormalisation_finite_terms']=='removed':
            # Discard all terms prefixed with UVRenormFinite, so set discarded power to 1
            FORM_vars['UVRENORMFINITEPOWERTODISCARD'] = 1
        else:
            raise FormProcessingError("The FORM processing option 'renormalisation_finite_terms' "+
                                      "can only take the following value: 'together', 'only' or 'removed', but not '%s'."%FORM_processing_options['renormalisation_finite_terms'])

        _MANDATORY_FORM_VARIABLES = ['SGID','NINITIALMOMENTA','NFINALMOMENTA','SELECTEDEPSILONORDER','UVRENORMFINITEPOWERTODISCARD']

        if integrand_type is not None:
            FORM_vars['INTEGRAND'] = integrand_type

        # TODO: multiprocess this loop
        max_buffer_size = 0
        all_numerator_ids = []
        with progressbar.ProgressBar(
            prefix = 'Processing numerators (graph #{variables.i_graph}/%d, LMB #{variables.i_lmb}/{variables.max_lmb}) with FORM ({variables.timing} ms / supergraph) : '%len(self),
            max_value=len(self),
            variables = {'timing' : '0', 'i_graph' : '0', 'i_lmb': '0', 'max_lmb': '0'}
            ) as bar:
            total_time = 0.

            for i_graph, graph in enumerate(self):
                graph.is_zero = True

                graphs_to_process = []
                if isinstance(graph[0].additional_lmbs,list) and graph[0].additional_lmbs != []:
                    graphs_to_process.append( (0,i_graph,graph[0]) )
                    graphs_to_process.extend([(g.additional_lmbs,i_graph,g) for _,_,_,g in graph[0].additional_lmbs])
                else:
                    # By setting the active graph to None we will then sum overall members of the iso set.
                    graphs_to_process.append( (0,i_graph, None) )
                bar.update(max_lmb='%d'%len(graphs_to_process))
                for i_lmb, i_g, active_graph in graphs_to_process:
                    bar.update(i_graph='%d'%(i_graph+1), i_lmb='%d'%(i_lmb+1))
                    i = i_lmb*FORM_processing_options['FORM_call_sig_id_offset_for_additional_lmb']+i_g
                    all_numerator_ids.append(i)
                    FORM_vars['SGID']='%d'%i
                    time_before = time.time()
                    num = graph.generate_numerator_functions(additional_overall_factor,
                                    workspace=workspace, FORM_vars=FORM_vars, active_graph=active_graph,process_definition=process_definition)
                    total_time += time.time()-time_before
                    num = num.replace('i_', 'I')
                    num = input_pattern.sub(r'lm[\1]', num)
                    num = num.replace('\nZ', '\n\tZ') # nicer indentation

                    confs = []
                    numerator_main_code = ''
                    conf_secs = num.split('#CONF')
                    for conf_sec in conf_secs[1:]:
                        conf_sec = conf_sec.replace("#CONF\n", '')

                        # collect all temporary variables
                        temp_vars = list(sorted(set(var_pattern.findall(conf_sec))))

                        # parse the configuration id
                        conf = list(energy_exp.finditer(conf_sec))[0].groups()[0].split(',')
                        conf_id = conf[0]
                        
                        ltd_vars = [v for v in conf[1:] if v != 'c0']
                        max_rank = 8 # FIXME: determine somehow
                        # create the dense polynomial in the LTD energies
                        numerator_pows = [j for i in range(max_rank + 1) for j in combinations_with_replacement(range(len(ltd_vars)), i)]

                        mono_secs = conf_sec.split('#NEWMONOMIAL')

                        num_body = energy_exp.sub('', '\n'.join(mono_secs[0].split('\n')[2:]))
                        
                        rank = 0
                        max_index = 0
                        return_statements = []
                        for mono_sec in mono_secs[1:]:
                            mono_sec = mono_sec.replace('#NEWMONOMIAL\n', '')
                            mono_sec = split_number.sub('', mono_sec)
                            # parse monomial powers and get the index in the polynomial in the LTD basis
                            pows = list(energy_exp.finditer(mono_sec))[0].groups()[0].split(',')
                            pows = tuple(sorted(ltd_vars.index(r) for r in pows if r != 'c0'))
                            rank = max(rank, len(pows))
                            index = numerator_pows.index(pows)
                            max_index = max(index, max_index)
                            max_buffer_size = max(index, max_buffer_size)
                            mono_sec = energy_exp.sub('', mono_sec)
                            returnval = list(return_exp.finditer(mono_sec))[0].groups()[0]
                            return_statements.append((index, return_exp.sub('\tout[{}] = {};'.format(index, returnval), mono_sec)))

                        for r in sorted(return_statements, key=lambda x: x[0]):
                            num_body += r[1]

                        num_body = num_body.replace('logmUV', 'log(mUV*mUV)').replace('logmu' , 'log(mu*mu)').replace('logmt' , 'log(mass_t*mass_t)')

                        confs.append((conf_id, rank))

                        if len(temp_vars) > 0:
                            graph.is_zero = False

                        numerator_main_code += '\n// polynomial in {}'.format(','.join(conf[1:]))
                        numerator_main_code += '\nstatic inline int %(header)sevaluate_{}_{}(double complex lm[], double complex params[], double complex* out) {{\n\t{}'.format(i, conf_id,
                            'double complex {};'.format(','.join(temp_vars)) if len(temp_vars) > 0 else ''
                        ) + num_body + '\n\treturn {}; \n}}\n'.format(max_index + 1)

                        numerator_main_code_f128 = num_body.replace('pow', 'cpowq').replace('sqrt', 'csqrtq').replace('log', 'clogq').replace('pi', 'M_PIq').replace('double complex', '__complex128')
                        numerator_main_code_f128 = float_pattern.sub(r'\1q', numerator_main_code_f128)
                        numerator_main_code += '\n' + '\nstatic inline int %(header)sevaluate_{}_{}_f128(__complex128 lm[], __complex128 params[], __complex128* out) {{\n\t{}\n{}\n\treturn {};\n}}\n'.format(i, conf_id,
                            '__complex128 {};'.format(','.join(temp_vars)) if len(temp_vars) > 0 else '', numerator_main_code_f128, max_index + 1
                        )

                    numerator_main_code += \
"""
int %(header)sevaluate_{}(double complex lm[], double complex params[], int conf, double complex* out) {{
   switch(conf) {{
{}
    }}
}}

int %(header)sevaluate_{}_f128(__complex128 lm[], __complex128 params[], int conf, __complex128* out) {{
    switch(conf) {{
{}
    }}
}}

int %(header)sget_rank_{}(int conf) {{
   switch(conf) {{
{}
    }}
}}
""".format(i,
        '\n'.join(
        ['\t\tcase {}: return %(header)sevaluate_{}_{}(lm, params, out);'.format(conf, i, conf) for conf, _ in sorted(confs)] +
        (['\t\tdefault: out[0] = 0.; return 1;']) # ['\t\tdefault: raise(SIGABRT);'] if not graph.is_zero else 
        ),
        i,
        '\n'.join(
        ['\t\tcase {}: return %(header)sevaluate_{}_{}_f128(lm, params, out);'.format(conf, i, conf) for conf, _ in sorted(confs)] +
        (['\t\tdefault: out[0] = 0.; return 1;']) # ['\t\tdefault: raise(SIGABRT);'] if not graph.is_zero else 
        ),
        i,
        '\n'.join(
        ['\t\tcase {}: return {};'.format(conf, rank) for conf, rank in sorted(confs)] +
        (['\t\tdefault: return 0;']) # ['\t\tdefault: raise(SIGABRT);'] if not graph.is_zero else 
        ))

                    bar.update(timing='%d'%int((total_time/float(i_graph+1))*1000.0))
                    bar.update(i_graph+1)

                    writers.CPPWriter(pjoin(root_output_path, '%(header)snumerator{}.c'%header_map).format(i)).write((numerator_header + numerator_main_code)%header_map)

        numerator_code = \
"""
#include <tgmath.h>
#include <quadmath.h>
#include <signal.h>

{}
{}
{}

int %(header)sevaluate(double complex lm[], double complex params[], int diag, int conf, double complex* out) {{
    switch(diag) {{
{}
    }}
}}

int %(header)sevaluate_f128(__complex128 lm[], __complex128 params[], int diag, int conf, __complex128* out) {{
    switch(diag) {{
{}
    }}
}}

int %(header)sget_buffer_size() {{
    return {};
}}

int %(header)sget_rank(int diag, int conf) {{
    switch(diag) {{
{}
    }}
}}
""".format(
    '\n'.join('int %(header)sevaluate_{}(double complex[], double complex[], int conf, double complex*);'.format(i) for i in all_numerator_ids),
    '\n'.join('int %(header)sevaluate_{}_f128(__complex128[], __complex128[], int conf, __complex128*);'.format(i) for i in all_numerator_ids),
    '\n'.join('int %(header)sget_rank_{}(int conf);'.format(i) for i in all_numerator_ids),
    '\n'.join(
    ['\t\tcase {}: return %(header)sevaluate_{}(lm, params, conf, out);'.format(i, i) for i in all_numerator_ids]+
    ['\t\tdefault: raise(SIGABRT);']
    ),
    '\n'.join(
    ['\t\tcase {}: return %(header)sevaluate_{}_f128(lm, params, conf, out);'.format(i, i) for i in all_numerator_ids]+
    ['\t\tdefault: raise(SIGABRT);']
    ), 
    (max_buffer_size + 1) * 2,
    '\n'.join(
    ['\t\tcase {}: return %(header)sget_rank_{}(conf);'.format(i, i) for i in all_numerator_ids]+
    ['\t\tdefault: raise(SIGABRT);']
    )
    )
        writers.CPPWriter(pjoin(root_output_path, '%(header)snumerator.c'%header_map)).write(numerator_code%header_map)


        params = copy.deepcopy(params)
        params['mUV'] = 'params[0]'
        params['mu'] = 'params[1]'
        params['gs'] = 'params[2]'
        params['small_mass_sq'] = 'params[3]'

        header_code = \
"""
#ifndef NUM_H
#define NUM_H

{}

#endif
""".format('\n'.join('#define  {} {}'.format(k, v) for k, v in params.items()))

        writers.CPPWriter(pjoin(root_output_path, '%(header)snumerator.h'%header_map)).write(header_code%header_map)

        if integrand_type is not None:
            self.generate_integrand_functions(root_output_path, additional_overall_factor='', 
                            params={}, output_format='c', workspace=None, integrand_type=integrand_type,process_definition=process_definition)

    def generate_squared_topology_files(self, root_output_path, model, process_definition, n_jets, final_state_particle_ids=(), jet_ids=None, filter_non_contributing_graphs=True, workspace=None,
        integrand_type=None):
        if workspace is None:
            workspace = pjoin(root_output_path, os.pardir, 'workspace')
        topo_collection = {
            'name': self.name,
            'topologies': []
        }

        # Set sensible jet_ids if none
        if jet_ids is None:
            jet_ids=tuple(list(range(1,6))+list(range(-1,-6,-1))+[21,82,-82,])

        contributing_supergraphs = []

        non_zero_graph = 0
        with progressbar.ProgressBar(prefix='Generating squared topology files (graph #{variables.i_graph}/%d, LMB #{variables.i_lmb}/{variables.max_lmb}, PF #{variables.PF_config}, {variables.timing} ms / supergraph) : '%len(self), 
                max_value=len(self), variables = {'timing' : '0', 'i_graph' : '0', 'i_lmb': '0', 'max_lmb' : '0', 'PF_config': 'N/A'} ) as bar:
            total_time=0.0
            for i, g in enumerate(self):
                time_before = time.time()
                bar.update(i_graph='%d'%(i+1))
                if g.generate_squared_topology_files(root_output_path, model, process_definition, n_jets, numerator_call=non_zero_graph, 
                                                            final_state_particle_ids=final_state_particle_ids,jet_ids=jet_ids, 
                                                            workspace=workspace, bar=bar, integrand_type=integrand_type):
                    topo_collection['topologies'].append({
                        'name': g[0].name,
                        # Let us not put it there but in the topology itself
                        #'benchmark_result': g[0].benchmark_result,
                        'multiplicity': g[0].multiplicity
                        ,'additional_LMBs': [
                            {
                                'name' : '%s_LMB%d'%(other_supergraph.name,other_supergraph.additional_lmbs),
                                'defining_lmb_to_this_lmb' : original_lmb_to_other_lmb_affine_transfo,
                                'this_lmb_to_defining_lmb' : other_lmb_to_original_lmb_affine_transfo
                            }
                            for i_lmb, other_lmb_to_original_lmb_affine_transfo, 
                                original_lmb_to_other_lmb_affine_transfo, other_supergraph in g[0].additional_lmbs
                        ]

                    })
                    non_zero_graph += 1
                    contributing_supergraphs.append(g)

                total_time += time.time()-time_before
                bar.update(timing='%d'%int((total_time/float(i+1))*1000.0))
                bar.update(i+1)
        
        if filter_non_contributing_graphs:
            self[:] = contributing_supergraphs

        if FORM_processing_options['generate_renormalisation_graphs']:
            renormalization_graphs = self.generate_renormalization_graphs(model,process_definition)
        else:
            renormalization_graphs = []
        self.extend([FORMSuperGraphIsomorphicList([g]) for g in renormalization_graphs])

        with progressbar.ProgressBar(prefix='Generating renormalization squared topology files (graph #{variables.i_graph}/%d, LMB #{variables.i_lmb}/{variables.max_lmb}, PF #{variables.PF_config}, {variables.timing} ms / supergraph) : '%len(renormalization_graphs), 
                max_value=len(renormalization_graphs), variables = {'timing' : '0', 'i_graph' : '0', 'i_lmb': '0', 'max_lmb' : '0', 'PF_config': 'N/A'} ) as bar:
            total_time=0.0
            for i, g in enumerate(renormalization_graphs):
                time_before = time.time()
                bar.update(i_graph='%d'%(i+1))
                if g.generate_squared_topology_files(root_output_path, model, process_definition, n_jets, numerator_call=non_zero_graph, 
                                                            final_state_particle_ids=final_state_particle_ids,jet_ids=jet_ids,
                                                            workspace=workspace, bar=bar, integrand_type=integrand_type):
                    topo_collection['topologies'].append({
                        'name': g.name,
                        'multiplicity': g.multiplicity,
                        'additional_LMBs': []
                    })
                    non_zero_graph += 1
                    contributing_supergraphs.append(g)

                total_time += time.time()-time_before
                bar.update(timing='%d'%int((total_time/float(i+1))*1000.0))
                bar.update(i+1)

        try:
            import yaml
            from yaml import Loader, Dumper
        except ImportError:
            raise BaseException("Install yaml python module in order to import topologies from yaml.")

        open(pjoin(root_output_path, self.name + '.yaml'), 'w').write(yaml.dump(topo_collection, Dumper=Dumper))

    def get_renormalization_vertex(self, in_pdgs, loop_count, model, process_definition, is_external_bubble=False):
        
        # This overall factor is an overall fudge which will need to be understood together with how exactly we must
        # inherit the overall phase of the reference loop SG that this vertex is supposed to renormalised. 
        overall_factor = '(-1)'

        if loop_count != 1:
            raise FormProcessingError("Renormalisation is currently implemented only at one-loop.")
            
        if any( model.get_particle(pdg).get('ghost') for pdg in in_pdgs):
            # For now we absorbed all renormalisation constants into the gluon contribution so that we set the ghost 
            # renormalisation conditions to zero. This implies that the local UV pole cancellation check must be done
            # by combining the ghost and gluon supergraph with appropriate routing.
            return '0'

        # WARNING: All the things below are a SHAME and must be generalised as soon as we're done with the publication....

        # All renormalisation 2-point vertices are for now coded up directly in numerator.frm
        if len(in_pdgs)==2: #and not is_external_bubble:
            return overall_factor

        # No renormalisation for shrunk loop-induced vertices involving only gluons and the Higgs
        if all(pdg in [21,25] for pdg in in_pdgs):
            return '0'

        hardcoded_mass_parameters = {
            6   : 'mass_t',
            -6  : 'mass_t',
        }
        
        hardcoded_log_quark_mass = {
            6   : 'logmt',
            -6  : 'logmt'
        }

        def get_particle_mass(pdg):
            if pdg in hardcoded_mass_parameters:
                return hardcoded_mass_parameters[pdg]
            return model.get_particle(pdg).get('mass')

        pdgs = sorted([abs(pdg) for pdg in in_pdgs],reverse=True)
        quark_pdgs = range(1,7)

        active_quarks = [q_pdg for q_pdg in quark_pdgs if q_pdg not in process_definition.get('forbidden_particles')]
        n_massless_quarks = len([1 for q_pdg in active_quarks if get_particle_mass(q_pdg).upper()=='ZERO' ])
        massive_quark_pdgs = [ q_pdg for q_pdg in active_quarks if get_particle_mass(q_pdg).upper()!='ZERO' ]

        symbol_replacement_dict = {
            'C_F' : '(4/3)',
            'T_F' : '(1/2)',
            'C_A' : '3',
            'pi'  : 'pi',
            'gs'  : 'gs',
            'ep'  : 'ep',
            'mu_r': 'mu_r',
            'n_f' : '%d'%n_massless_quarks,
            'UVRenormFINITE_marker': 'UVRenormFINITE'
        }

        delta_Z_massless_quark = '%(C_F)s*%(gs)s^2/16/%(pi)s^2*(-1/%(ep)s)'
        delta_Z_massive_quark = '%(C_F)s*%(gs)s^2/16/%(pi)s^2*((-1/%(ep)s) + %(UVRenormFINITE_marker)s*(-4-3*(logmu - %(log_quark_mass)s)) )'

        # Note the many overall factors (1/2) come from the fact that the expressions are derived from the alpha couplings (i.e. g^2/4\pi)
        # and not for the coupling g directly, so that a factor 1/2 is warranted.

        delta_Z_gluon = []
        # Add the gluon and ghost contribution to the gluon wavefunction renormalisation
        delta_Z_gluon.append('( (1/2)*%(C_A)s*%(gs)s^2/48/%(pi)s^2 * ( 5/%(ep)s ) )'%symbol_replacement_dict)
        # Add the massless quark contributions to the gluon wavefunction renormalisation
        delta_Z_gluon.append('( (1/2)*%(n_f)s*%(T_F)s*%(gs)s^2/48/%(pi)s^2 * ( -4/%(ep)s ) )'%symbol_replacement_dict)
        # Add the massive quark contribution to the gluon wavefunction renormalisation
        for q_pdg in massive_quark_pdgs:

            if q_pdg not in hardcoded_log_quark_mass:
                raise FormProcessingError( 
                    ("There is not hard-coded symbol for the logarithm of the mass of particle with PDG=%d.\n"%q_pdg)+
                    "Are you sure you did not to specify forbidden particles in your process definition, using the / x y z syntax?")

            delta_Z_gluon.append('( (1/2)*%(T_F)s*%(gs)s^2/48/%(pi)s^2 * ( -4/%(ep)s + %(UVRenormFINITE_marker)s*(-4*(logmu - %(log_quark_mass)s)) ) )'%(
                dict(symbol_replacement_dict,**{
                    'log_quark_mass':hardcoded_log_quark_mass[q_pdg],
                })
            ))
        # Combine all terms for the gluon wavefunction renormalisation
        delta_Z_gluon = '(%s)'%('+'.join(delta_Z_gluon))

        delta_g_s = []
        # Add the gluon and ghost contribution to the strong coupling renormalisation
        delta_g_s.append('( (1/2)*%(C_A)s*%(gs)s^2/48/%(pi)s^2 * ( -11/%(ep)s ) )'%symbol_replacement_dict)
        # Add the massless quark contributions to the strong coupling renormalisation
        delta_g_s.append('( (1/2)*%(n_f)s*%(T_F)s*%(gs)s^2/48/%(pi)s^2 * ( 4/%(ep)s ) )'%symbol_replacement_dict)
        # Add the massive quark contribution to the strong coupling renormalisation
        for q_pdg in massive_quark_pdgs:
            delta_g_s.append('( (1/2)*%(T_F)s*%(gs)s^2/48/%(pi)s^2 * ( 4/%(ep)s + %(UVRenormFINITE_marker)s*(4*(logmu - %(log_quark_mass)s)) ) )'%(
                dict(symbol_replacement_dict,**{
                    'log_quark_mass':hardcoded_log_quark_mass[q_pdg],
                })
            ))
        # Combine all terms for the strong coupling renormalisation
        delta_g_s = '(%s)'%('+'.join(delta_g_s))

        # Combine all terms for the yukawa renormalisation
        # WARNING: He consider here only contribution from the top and bottom as these are the only fermions we consider as
        # having a yukawa interaction in the model, but we do allow for them to be massless and adjust the renormalisation accordingly.
        potential_yukawa_quark_pdgs = [5,6]
        active_yukawa_quarks = [q_pdg for q_pdg in potential_yukawa_quark_pdgs if q_pdg not in process_definition.get('forbidden_particles')]
        n_massless_yukawa_quarks = len([1 for q_pdg in active_yukawa_quarks if get_particle_mass(q_pdg).upper()=='ZERO' ])
        massive_yukawa_quark_pdgs = [ q_pdg for q_pdg in active_yukawa_quarks if get_particle_mass(q_pdg).upper()!='ZERO' ]

        delta_yukawa = []
        # Add the massless quark contributions to the yukawa renormalisation
        delta_yukawa.append('( -(1/2)*%(n_f_yukawa)s*2*%(C_F)s*%(gs)s^2/16/%(pi)s^2 * ( 3/%(ep)s ) )'%(
            dict(symbol_replacement_dict,**{
                'n_f_yukawa':'%d'%n_massless_yukawa_quarks,
            })
        ))
        # Add the massive quark contribution to the yukawa renormalisation
        for q_pdg in massive_yukawa_quark_pdgs:
            delta_yukawa.append('( - (1/2)*2*%(C_F)s*%(gs)s^2/16/%(pi)s^2 * ( 3/%(ep)s + %(UVRenormFINITE_marker)s*(4 + 3*(logmu - %(log_quark_mass)s)) ) )'%(
                dict(symbol_replacement_dict,**{
                    'log_quark_mass':hardcoded_log_quark_mass[q_pdg],
                })
            ))
        # Combine all terms for the yukawa renormalisation
        delta_yukawa = '(%s)'%('+'.join(delta_yukawa))

        ########################################################
        # Now build the specified UV renormalisation vertices
        ########################################################

#       Self-energies are for now coded up direclty in numerator.frm

        # quark self-energy
#        if len(pdgs) == 2 and pdgs[0] in quark_pdgs and pdgs[1]==pdgs[0]:
#            quark_mass = get_particle_mass(pdgs[0])
#            if quark_mass=='ZERO':
#                res = ('(-(%s))'%delta_Z_massless_quark)%symbol_replacement_dict
#                res = res
#                return res
#            else:
#                res = ('(-(%s))'%delta_Z_massive_quark)%(
#                    dict(symbol_replacement_dict,**{
#                        'quark_mass':quark_mass,
#                        'log_quark_mass':hardcoded_log_quark_mass[pdgs[0]],                            
#                    }))
#                return res

        # gluon self-energy
#        if len(pdgs) == 2 and pdgs[0]==21 and pdgs[1]==pdgs[0]:
#            return '(-(%s))'%delta_Z_gluon

        # a/z qq vertex
        if len(pdgs) == 3 and (pdgs[0] in [22,23]) and (pdgs[1] in quark_pdgs) and (pdgs[2]==pdgs[1]):
            loop_multiplicity_factor = 1
            quark_mass = get_particle_mass(pdgs[1])
            res = '%s*'%overall_factor
            if quark_mass=='ZERO':
                res += ('(+(2)*(1/2)*(%s))'%delta_Z_massless_quark)%symbol_replacement_dict
            else:
                res += ('(+(2)*(1/2)*(%s))'%delta_Z_massive_quark)%(
                    dict(symbol_replacement_dict,**{
                        'quark_mass':quark_mass,
                        'log_quark_mass':hardcoded_log_quark_mass[pdgs[1]],
                    }))
            return '((1/%d)*%s)'%(loop_multiplicity_factor,res)

        # hqq vertex
        if len(pdgs) == 3 and (pdgs[0] in [25,]) and (pdgs[1] in quark_pdgs) and (pdgs[2]==pdgs[1]):
            loop_multiplicity_factor = 1
            res = [ '(+%s)'%delta_yukawa, ]
            quark_mass = get_particle_mass(pdgs[1])
            if quark_mass=='ZERO':
                res.append( ('(+(2)*(1/2)*(%s))'%delta_Z_massless_quark)%symbol_replacement_dict )
            else:
                res.append( ('(+(2)*(1/2)*(%s))'%delta_Z_massive_quark)%(
                    dict(symbol_replacement_dict,**{
                        'quark_mass':quark_mass,
                        'log_quark_mass':hardcoded_log_quark_mass[pdgs[1]],
                    })) )
            return '((1/%d)*%s*(%s))'%(loop_multiplicity_factor,overall_factor,('+'.join(res)))

        # gqq vertex
        if len(pdgs) == 3 and (pdgs[0] in [21,]) and (pdgs[1] in quark_pdgs) and (pdgs[2]==pdgs[1]):
            loop_multiplicity_factor = 2
            res = [ 
                '(+%s)'%delta_g_s, 
                '(+%s)'%delta_Z_gluon
            ]
            quark_mass = get_particle_mass(pdgs[1])
            if quark_mass=='ZERO':
                res.append( ('(+(2)*(1/2)*(%s))'%delta_Z_massless_quark)%symbol_replacement_dict )
            else:
                res.append( ('(+(2)*(1/2)*(%s))'%delta_Z_massive_quark)%(
                    dict(symbol_replacement_dict,**{
                        'quark_mass':quark_mass,
                        'log_quark_mass':hardcoded_log_quark_mass[pdgs[1]],
                    })) )
            return '((1/%d)*%s*(%s))'%(loop_multiplicity_factor, overall_factor,('+'.join(res)))

        # ggg vertex
        if len(pdgs) == 3 and all(pdg==21 for pdg in pdgs):
            # 2 per n_f, 2 per ghost and 1 for the gluon
            loop_multiplicity_factor = 1 + 2 + 2 * n_massless_quarks
            res = [ 
                '(+%s)'%delta_g_s,
                '(+3*(%s))'%delta_Z_gluon,
            ]
            return '((1/%d)*%s*(%s))'%(loop_multiplicity_factor, overall_factor,('+'.join(res)))

        # gggg vertex
        if len(pdgs) == 4 and all(pdg==21 for pdg in pdgs):
            # 6 per n_f, 6 per ghost and 3 for the gluon
            loop_multiplicity_factor = 3 + 6 + 6 * n_massless_quarks
            res = [ 
                '(+2*(%s))'%delta_g_s,
                '(+4*(%s))'%delta_Z_gluon,
            ]
            return '((1/%d)*%s*(%s))'%(loop_multiplicity_factor,overall_factor,('+'.join(res)))

        return None

    @classmethod
    def shrink_edges(cls, edges, nodes, edges_to_shrink, bubble_external_shrinking_step=False):
        subgraph_pdgs = set()
        for ek in list(edges.keys()):
            ee = edges[ek]
            if ee['name'] not in edges_to_shrink:
                continue
            subgraph_pdgs.add(ee['PDG'])

            # remove the edge from the vertices
            for na in ee['vertices']:
                node = nodes[na]
                edge_index = node['edge_ids'].index(ek)
                for g in ('PDGs', 'indices', 'momenta', 'edge_ids'):
                    node[g] = tuple(ind for i, ind in enumerate(node[g]) if i != edge_index)

            # fuse vertices
            if ee['vertices'][0] != ee['vertices'][1]:
                node = nodes[ee['vertices'][0]]
                # START TEMPORARY HACK (Just for the node rendering purposes)
                if not bubble_external_shrinking_step:
                    # Mark this node as being shrunk. Eventually there will be a more general way of handling
                    # multiple effective shrunk renormalisation nodes. For now, I just flag it here as one-loop
                    node['renormalisation_vertex_n_loops'] = 1
                    node['renormalisation_vertex_n_shrunk_edges'] = len(edges_to_shrink)
                # END TEMPORARY HACK
                for na in ee['vertices'][1:]:
                    n = nodes[na]
                    for g in ('PDGs', 'indices', 'momenta', 'edge_ids'):
                        node[g] = tuple(list(node[g]) + list(n[g]))
                    # START TEMPORARY HACK (Just for the node rendering purposes)
                    if bubble_external_shrinking_step:
                        for g in ['renormalisation_vertex_n_loops', 'renormalisation_vertex_n_shrunk_edges']:
                            if g not in node and g in n:
                                node[g] = n[g]
                    # END TEMPORARY HACK
                    del nodes[na]

                vert_to_replace = ee['vertices'][1]
                vert_to_replace_with = ee['vertices'][0]
                for e in edges.values():
                    e['vertices'] = tuple(vert_to_replace_with if ind == vert_to_replace else ind for ind in e['vertices'])

            del edges[ek]

        return subgraph_pdgs

    def generate_renormalization_graphs(self, model, process_definition):
        """Generate all required renormalization graphs for the graphs in the supergraph list"""

        import sympy as sp
        n_externals = max(len([1 for e in graph[0].edges.values() if e['type']=='in' or e['type']=='out']) for graph in self)
        if model is None:
            pdg_primes = {pdg : sp.prime(i + n_externals + 1) for i, pdg in enumerate([1,2,3,4,5,6,11,12,13,21,22,25,82])}
        else:
            pdg_primes = {pdg : sp.prime(i + n_externals + 1) for i, pdg in enumerate([p['pdg_code'] for p in model['particles']])}

        renormalization_graphs = []
        for i_sg, gs in enumerate(self):
            squared_topology = gs[0].squared_topology
            for i_cut, cut_info in enumerate(squared_topology.cuts):
                for i_diag_set, diag_set in enumerate(cut_info['diagram_sets']):
                    if len(diag_set['uv_propagators']) == 0:
                        continue

                    # create the renormalization graph per diagram set
                    nodes = copy.deepcopy(gs[0].nodes)
                    edges = copy.deepcopy(gs[0].edges)
                    vertex_factors = []
                    uv_effective_vertex_edges = []

                    # BEGIN PIECE OF MULTIPLICITY HACK
                    has_pure_external_gluon_effective_vertex = False
                    # END PIECE OF MULTIPLICITY HACK

                    # Compute the bubble edges so as to be able to determine if the bubble to computed the renormalisation vertex for is "internal" or "external".
                    bubble_edges = set()
                    for c in cut_info['cuts']:
                        sig = next(ee['signature'] for ee in edges.values() if ee['name'] == c['edge'])
                        inv_sig = tuple([-s for s in x] for x in sig)
                        bubble_edges |= set(ee['name'] for ee in edges.values() if (ee['signature'] == sig or ee['signature'] == inv_sig) and ee['name'] != c['edge'])

                    for diag_info in diag_set['diagram_info']:
                        if diag_info['uv_vertices'] is not None and len(diag_info['uv_vertices']) > 0:
                            graph = diag_info['graph']
                            # look up every uv vertex in the table
                            for v, l in diag_info['uv_vertices']:
                                edges_on_vertex = [e[0] for e in graph.edge_map_lin if v in e[1:]]
                                uv_effective_vertex_edges.append((edges_on_vertex, l))
                                edge_pdgs = tuple(next(ee for ee in edges.values() if ee['name'] == e)['PDG'] for e in edges_on_vertex)
                                # BEGIN PIECE OF MULTIPLICITY HACK
                                if not has_pure_external_gluon_effective_vertex:
                                    has_pure_external_gluon_effective_vertex = set(edge_pdgs)==set([21,])
                                # END PIECE OF MULTIPLICITY HACK
                                is_external_bubble = ( (len(edges_on_vertex)==2) and any( (e in bubble_edges) for e in edges_on_vertex) )
                                vertex_contrib = self.get_renormalization_vertex(edge_pdgs, l, model, process_definition, is_external_bubble=is_external_bubble)
                                if vertex_contrib is None:
                                    logger.warning("WARNING: unknown renormalization vertex {} at {} loops".format(edge_pdgs, l))
                                    #raise AssertionError("Unknown renormalization vertex {} at {} loops".format(edge_pdgs, l))
                                    # It can often happen, for example for (but not only) Loop-Induced processes, that there must be
                                    # not renormalisation vertex associated to an integrated UV CT. For instance g g H H.
                                    # Ideally one should code those explicitly in the function `get_renormalization_vertex`, but for 
                                    # now it is also OK to simply let it return None with the above warning and use `0` as the numerator
                                    # for this renormalisation supergraph, so that it hopefully gets removed automatically at the time 
                                    # of explicitly building its numerator.
                                    vertex_factors.append('0')
                                else:
                                    vertex_factors.append('(' + vertex_contrib + ')')

                    # shrink the UV subgraph in the edge list
                    subgraph_pdgs = self.shrink_edges(edges, nodes, diag_set['uv_propagators'],bubble_external_shrinking_step=False)

                    # BEGIN PIECE OF MULTIPLICITY HACK
                    is_made_of_gluons_only = set(subgraph_pdgs)==set([21,])
                    # END PIECE OF MULTIPLICITY HACK

                    # Recomputed the bubble edges now the that UV subgraph has been shrunk, they will be removed at the end
                    bubble_edges = set()
                    for c in cut_info['cuts']:
                        sig = next(ee['signature'] for ee in edges.values() if ee['name'] == c['edge'])
                        inv_sig = tuple([-s for s in x] for x in sig)
                        bubble_edges |= set(ee['name'] for ee in edges.values() if (ee['signature'] == sig or ee['signature'] == inv_sig) and ee['name'] != c['edge'])

                    # compute the proper multiplicity by dividing out the symmetry factor of the UV divergent subgraphs
                    # TODO: make multi-loop and multi-uv subgraph compatible
                    # In the new paradigm this is done automatically as it is divided out in the corresponding renormalisation counterterms.
                    multiplicity = gs[0].multiplicity
#                    if len(set([21, 22]) & subgraph_pdgs) == 0:
#                        multiplicity /= 2

                    # set the correct particle ordering for all the edges
                    for n in nodes.values():
                        if len(n['PDGs']) > 1:
                            edge_order = FORMSuperGraph.sort_edges(model, [{'PDG': pdg, 'index': i} for i, pdg in enumerate(n['PDGs'])])
                            for g in ('PDGs', 'indices', 'momenta', 'edge_ids'):
                                n[g] = tuple(n[g][eo['index']] for eo in edge_order)

                    # check if the remaining graph is seen before
                    # NOTE: this must be exactly the same procedure as the graph isomorphic filtering,
                    # to obtain the correct multiplicity
                    # TODO: The multiplicity factor must be computed accordingly to what was already done
                    # at LO for the original loop SG, that is including the color factors as well.
                    g = igraph.Graph()
                    vertex_map = list(nodes.keys())
                    g.add_vertices(len(vertex_map))
                    undirected_edges = set()

                    for e in edges.values():
                        undirected_edges.add(tuple(sorted(vertex_map.index(v) for v in e['vertices'])))

                    edge_colors = []
                    ext_id = 0
                    for ue in undirected_edges:
                        e_color = 1
                        for e in edges.values():
                            if tuple(sorted(vertex_map.index(v) for v in e['vertices'])) == ue:
                                # find the associated edge to get the PDG
                                edge = next(ee for ee in edges.values() if ee['name'] == e['name'])
                                if edge['type'] == 'in':
                                    ext_id += 1
                                    e_color *= sp.prime(ext_id)
                                else:
                                    e_color *= pdg_primes[abs(edge['PDG'])]
                        edge_colors.append(e_color)
                        g.add_edges([ue])

                    v_colors = [0]*len(vertex_map)
                    for v, l in uv_effective_vertex_edges:
                        # find the vertex that is shrunk and give it a color
                        n = next(ni for ni, n in nodes.items() if set(v) == set(edges[e]['name'] for e in n['edge_ids']))
                        v_colors[vertex_map.index(n)] = l

                    for (ref_graph, ref_v_colors, ref_e_colors), ref_ren_graph in renormalization_graphs:

                        # BEGIN PIECE OF MULTIPLICITY HACK
                        # If this loop super graph is a multi-gluon one but does not contain only gluons in its loop content
                        # then we skip it here immediately as we want to instead use the loop graphs with gluon only as a 
                        # representative so that we can inherit the right symmetry factor
#                        if has_pure_external_gluon_effective_vertex and not is_made_of_gluons_only:
#                            break
                        # END PIECE OF MULTIPLICITY HACK

                        if ref_graph.get_isomorphisms_vf2(g,
                            color1=ref_v_colors,
                            color2=v_colors,
                            edge_color1=ref_e_colors,
                            edge_color2=edge_colors) != []:

                            # Accumulate the multiplicity in the renormalisation SG since the right quantity
                            # will be divided out in the definition of the renormalisation vertex
                            # However we must not aggregate diagram sets from the same supergraph and same uv propgators
                            if not any(
                                (   ( self[o_sg][0].name == gs[0].name ) and 
                                        ( set(self[o_sg][0].squared_topology.cuts[o_cut]['diagram_sets'][o_diag_set]['uv_propagators']) ==
                                          set(diag_set['uv_propagators']) )
                                )  for o_sg,o_cut,o_diag_set in ref_ren_graph.matching_loop_subdiags):
                                    ref_ren_graph.multiplicity += multiplicity
                                    ref_ren_graph.matching_loop_subdiags.append((i_sg,i_cut,i_diag_set))
                                    ref_ren_graph.name += '_%sDS%d'%(gs[0].name,diag_set['id'])

                            # BEGIN PIECE OF MULTIPLICITY HACK
#                            if not has_pure_external_gluon_effective_vertex:
#                                if multiplicity != ref_ren_graph.multiplicity:
#                                    raise FormProcessingError(
#                                        "Two loop SG mapped to the same renormalisation graph have "+
#                                        "different multiplicities even though it is expected that they should not.")
                            # END PIECE OF MULTIPLICITY HACK
                            
                            break
                    else:
                        # remove all bubble edges, since they will cancel with the effective vertex
                        #subgraph_pdgs = self.shrink_edges(edges, nodes, bubble_edges, bubble_external_shrinking_step=True)
                        
                        # set the correct particle ordering for all the edges
                        for n in nodes.values():
                            if len(n['PDGs']) > 1:
                                edge_order = FORMSuperGraph.sort_edges(model, [{'PDG': pdg, 'index': i} for i, pdg in enumerate(n['PDGs'])])
                                for n_keys in ('PDGs', 'indices', 'momenta', 'edge_ids'):
                                    n[n_keys] = tuple(n[n_keys][eo['index']] for eo in edge_order)

                        # WARNING: Double-check that it is indeed always the right thing to do to inherit the overall factor
                        # in the renormalisation graph from the reference one.

                        # add a new supergraph for this renormalization component
                        form_graph = FORMSuperGraph(name='renorm_{}DS{}'.format(gs[0].name,diag_set['id']),
                                edges=edges, nodes=nodes,
                                overall_factor='(%s)*(%s)'%(
                                    gs[0].overall_factor, ('*'.join(vertex_factors))
                                ), multiplicity=multiplicity)

                        # Add a container with all SG names mapped to this renormalisation graph
                        form_graph.matching_loop_subdiags = [(i_sg,i_cut,i_diag_set)]

                        # set a basis
                        topo_generator = LTD.ltd_utils.TopologyGenerator([(e['name'], e['vertices'][0], e['vertices'][1]) for e in edges.values()])
                        topo_generator.generate_momentum_flow()
                        sig = topo_generator.get_signature_map()

                        for eid, e in edges.items():
                            e['signature'] = [sig[e['name']][0],
                                [ i+o for i,o in zip(sig[e['name']] [1][:len(sig[e['name']] [1])//2],sig[e['name']] [1][len(sig[e['name']] [1])//2:]) ]]
                            e['momentum'] = FORMSuperGraph.momenta_decomposition_to_string(e['signature'], set_outgoing_equal_to_incoming=False)

                            if len(e['vertices']) == 1:
                                continue

                            for i, vi in enumerate(e['vertices']):
                                e_index = nodes[vi]['edge_ids'].index(eid)
                                mom = list(nodes[vi]['momenta'])
                                mom[e_index] = '-({})'.format(e['momentum']) if i == 0 else e['momentum']
                                nodes[vi]['momenta'] = mom

                        renormalization_graphs.append(((g, v_colors, edge_colors), form_graph))

        return [g for _, g in renormalization_graphs]

    def produce_output(self):
        form_processor = FORMProcessor(self, computed_model, process_definition)
        TMP_OUTPUT = pjoin(root_path, 'TEST_' + self.name.split('_')[0])
        Path(TMP_OUTPUT).mkdir(parents=True, exist_ok=True)
        TMP_workspace = pjoin(TMP_OUTPUT, 'FORM', 'workspace')
        Path(TMP_workspace).mkdir(parents=True, exist_ok=True)
        TMP_FORM = pjoin(TMP_OUTPUT, 'FORM')
        #FORM_processing_options['compilation-options'] += ['-e', "OPTIMIZATION_LVL=2"]
        form_processor.generate_squared_topology_files(TMP_OUTPUT, 0,
                    workspace=TMP_workspace,
                    integrand_type='PF')
        form_processor.generate_numerator_functions(TMP_FORM, output_format='c',
                    workspace=TMP_workspace,
                    integrand_type='PF')
        # copy the makefile, stripping the first line
        source_file = open('Templates/FORM_output_makefile', 'r')
        source_file.readline()
        with open(pjoin(TMP_FORM, 'Makefile'), 'w') as target_file:
            shutil.copyfileobj(source_file, target_file)
        Path(pjoin(TMP_OUTPUT, 'lib')).mkdir(parents=True, exist_ok=True)
        FORMProcessor.compile(TMP_FORM)

        drawings_output_path = pjoin(TMP_OUTPUT, 'Drawings')
        Path(drawings_output_path).mkdir(parents=True, exist_ok=True)
        shutil.copy(pjoin(plugin_path, 'Templates','Drawings_makefile'),
                    pjoin(drawings_output_path,'Makefile'))
        form_processor.draw(drawings_output_path)


class FORMProcessor(object):
    """ A class for taking care of the processing of a list of FORMSuperGraphList.
    Useful because many aspects common to all supergraphs and function do not belong to FORMSuperGraphList.
    """

    def __init__(self, super_graphs_list, model, process_definition):
        """ Specify aditional information such as the model that is useful for FORM processing."""
        self.super_graphs_list = super_graphs_list
        self.model = model
        self.process_definition = process_definition
        if isinstance(self.process_definition, base_objects.ProcessDefinition):
            all_processes = list(proc for proc in self.process_definition)
            self.repr_process = all_processes[0]
        else:
            self.repr_process = self.process_definition

        if not FORM_processing_options['generate_integrated_UV_CTs']:
            logger.warning('%s\n\nGeneration of integrated UV CTs is disabled per user request. Physical results will be incorrect.\n\n%s'%(utils.bcolors.RED,utils.bcolors.ENDC))
        if not FORM_processing_options['generate_renormalisation_graphs']:
            logger.warning('%s\n\nGeneration of renormalisation contributions is disabled per user request. Physical results will be incorrect.\n\n%s'%(utils.bcolors.RED,utils.bcolors.ENDC))
        if FORM_processing_options['UV_min_dod_to_subtract']>0:
            logger.warning('%s\n\nAs per user request, not all UV divegences will be locally subtracted. Numerical integration may be divergent in the UV.\n\n%s'%(utils.bcolors.RED,utils.bcolors.ENDC))
        if FORM_processing_options['UV_min_dod_to_subtract']<0:
            logger.warning('%s\n\nAs per user request, subleading UV divegences will also be locally subtracted. Very complicated integrands may result from this.\n\n%s'%(utils.bcolors.RED,utils.bcolors.ENDC))
        if FORM_processing_options['selected_epsilon_UV_order']!=0:
            logger.warning('%s\n\nAs per user request, the selected epsilon order to be exported will be %d. This must be for pole cancellation check only. \n\n%s'%(
                                                        utils.bcolors.RED, FORM_processing_options['selected_epsilon_UV_order'], utils.bcolors.ENDC))
        if FORM_processing_options['renormalisation_finite_terms']!='together':
            logger.warning(("%s\n\nAs per user request, the finite part of the renormalisation counterterms have been elected to be included "+
                            "as '%s' (and not the default 'together'). Results will be incorrect if not post-processed.\n\n%s")%(
                utils.bcolors.RED, FORM_processing_options['renormalisation_finite_terms'], utils.bcolors.ENDC))

    def draw(self, output_dir):
        """ For now simply one Mathematica script per supergraph."""

        for i_graph, super_graphs in enumerate(self.super_graphs_list):
            super_graphs[0].draw(self.model, output_dir, FORM_id=i_graph)
            # Draw supergraphs for additional LMBs
            if isinstance(super_graphs[0].additional_lmbs,list):
                for i_lmb,_,_,sg in super_graphs[0].additional_lmbs:
                    sg.draw(self.model, output_dir, FORM_id=i_graph, lmb_id=i_lmb)

    def generate_numerator_functions(self, root_output_path, output_format='c',workspace=None, header="", integrand_type=None, include_hel_avg_factor=None):
        assert(header in ['MG', 'QG', ''])

        params = {
            'mass_t': self.model['parameter_dict'][self.model.get_particle(6).get('mass')].real,
            'gs': self.model['parameter_dict']['G'].real,
            'ge': math.sqrt(4. * math.pi / self.model['parameter_dict']['aEWM1'].real),
            'gy': self.model['parameter_dict']['mdl_yt'].real / math.sqrt(2.),
            'ghhh': 6. * self.model['parameter_dict']['mdl_lam'].real,
            'pi': 'M_PI',
        }

        if include_hel_avg_factor is None:
            helicity_averaging_factor = 1
            for leg in self.repr_process.get('legs'):
                # Skip final states
                if leg.get('state') is True:
                    continue

                helicity_averaging_factor *= len(self.model.get_particle(leg.get('id')).get_helicity_states())
            helicity_averaging_factor = "/" + str(helicity_averaging_factor)
            additional_overall_factor = helicity_averaging_factor
        else:
            additional_overall_factor = '*(%d)'%(int(include_hel_avg_factor))

        return self.super_graphs_list.generate_numerator_functions(
            root_output_path,
            output_format=output_format,
            additional_overall_factor=additional_overall_factor,
            params=params,workspace=workspace, header=header,
            integrand_type=integrand_type,
            process_definition=self.process_definition
        )

    @classmethod
    def compile(cls, root_output_path, arg=[]):

        if os.path.isfile(pjoin(root_output_path,'Makefile')):
            try:
                logger.info("Now compiling FORM-generated numerators with options: %s ..."%(' '.join(FORM_processing_options['compilation-options'])))
                misc.compile(arg=FORM_processing_options['compilation-options'] ,cwd=root_output_path,mode='cpp', nb_core=FORM_processing_options["cores"])
            except MadGraph5Error as e:
                logger.info("%sCompilation of FORM-generated numerator failed:\n%s%s"%(
                    utils.bcolors.RED,str(e),utils.bcolors.ENDC))
        else:
            logger.warning(("\n%sYou are running FORM_processing directly from the __main__ of FORM_processing.py.\n"+
                           "You will thus need to compile numerators.c manually.%s")%(utils.bcolors.GREEN, utils.bcolors.ENDC))

    def generate_squared_topology_files(self, root_output_path, n_jets, final_state_particle_ids=(), jet_ids=None, filter_non_contributing_graphs=True, workspace=None,
        integrand_type=None):
        self.super_graphs_list.generate_squared_topology_files(
            root_output_path, self.model, self.process_definition, n_jets, final_state_particle_ids, jet_ids=jet_ids, filter_non_contributing_graphs=filter_non_contributing_graphs, workspace=workspace,
            integrand_type=integrand_type
        )


if __name__ == "__main__":
   
    parser = argparse.ArgumentParser(description='Generate numerators with FORM and yaml for Rust.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--diagrams_python_source', default='None', type=str,
                        help='path to the python diagram output files.')
    parser.add_argument('--model', default='sm-no_widths', type=str,
                        help='Path to UFO model to load.')
    parser.add_argument('--process', default='e+ e- > a > d d~', type=str,
                        help='Process definition to consider.')
    parser.add_argument('--cores', default=1, type=int,
                        help='Number of FORM cores')
    parser.add_argument('--optim_iter', default=100, type=int,
                        help='Number of iterations for numerator optimization')
    parser.add_argument('--restrict_card',
        default=pjoin(os.environ['MG5DIR'],'models','sm','restrict_no_widths.dat'), 
                        type=str, help='Model restriction card to consider.')
    args = parser.parse_args()

    if args.cores > 1:
        FORM_processing_options['cores'] = args.cores
    FORM_processing_options['extra-options'] = {'OPTIMITERATIONS': args.optim_iter}

    import alpha_loop.interface as interface
    cli = interface.alphaLoopInterface()

    cli.do_import('model %s'%args.model)
    computed_model = model_reader.ModelReader(cli._curr_model)
    computed_model.set_parameters_and_couplings(args.restrict_card)        
    process_definition=cli.extract_process(args.process, proc_number=0)

    # result is -2 Zeta[3] 3 Pi/(16 Pi^2)^3 = -5.75396*10^-6
    mercedes = FORMSuperGraphList.from_squared_topology([('q1', 0, 1), ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 6),
            ('p4', 6, 5), ('p5', 5, 1), ('p6', 2, 4), ('p7', 3, 4), ('p8', 4, 5), ('q2', 6, 7)], "Mercedes", ['q1'], computed_model,
            loop_momenta_names=('p1', 'p2', 'p3')).produce_output()

    # result is -5 Zeta[5] 4 Pi/(16 Pi^2)^4 = -1.04773*10^-7
    FORMSuperGraphList.from_squared_topology([('q1', 0, 1), ('p1', 1, 2), ('p2', 2, 7), ('p3', 7, 3), ('p4', 3, 6),
        ('p5', 6, 5), ('p6', 5, 1), ('p7', 2, 4), ('p8', 3, 4), ('p9', 4, 5), ('p10', 7, 4), ('q2', 6, 8)],
        "DoubleMercedes", ['q1'], computed_model, loop_momenta_names=('p1', 'p2', 'p3', 'p4')).produce_output()

    # result is -5 /2 Zeta[5] 4 Pi/(16 Pi^2)^4 = -5.23865e-08
    FORMSuperGraphList.from_squared_topology([('q1', 0, 1), ('p1', 1, 3), ('p2', 1, 4), ('p3',2, 3), ('p4', 2, 5),
        ('p5', 3, 6), ('p6',4, 7), ('p7', 4, 8), ('p8', 5, 7), ('p9', 5, 8), ('p10', 6, 7), ('p11',6, 8), ('q2', 2, 9)],
        "STF4L17", ['q1'], computed_model).produce_output()
