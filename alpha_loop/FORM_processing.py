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
from itertools import chain
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

logger = logging.getLogger('alphaLoop.FORM_processing')

if __name__ == "__main__":
    logging.basicConfig()
    logger.setLevel(logging.INFO)

plugin_path = os.path.dirname(os.path.realpath( __file__ ))


FORM_processing_options = {
    'FORM_path': 'form', 
    'cores': multiprocessing.cpu_count(), 
    'extra-options': '-D OPTIMITERATIONS=1000',
    # If None, only consider the LMB originally chosen.
    # If positive and equal to N, consider the first N LMB from the list of LMB automatically generated
    # If negative consider all possible LMBs.
    'number_of_lmbs' : None,
    # If None, the reference LMB will be the one originally chosen.
    # If positive and equal to N, the Nth LMB will be used for the reference implementation of the supergraph.
    'reference_lmb' : None,
    'FORM_call_sig_id_offset_for_additional_lmb' : 1000000
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
    }

    _include_momentum_routing_in_rendering=False
    _include_edge_name_in_rendering=False
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
    ):
        """ initialize a FORM SuperGraph from several options."""

        self.is_zero = False
        self.edges = edges
        self.nodes = nodes
        self.overall_factor = overall_factor
        self.multiplicity = multiplicity
        # A hashable call signature
        self.call_identifier = call_identifier
        if name is None:
            self.name = str(self.call_identifier)
        else:
            self.name = name

        self.replacement_rules = None

        # Store copies of self for different choices of LMB to be used for cross-check.
        self.additional_lmbs = []

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
            edge_repl_dict['in_node'] = str(edge_key[0])
            edge_repl_dict['out_node'] = str(edge_key[1])
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
            edge_label_pieces = [get_part_name(edge_data['PDG']),]
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

    def generate_squared_topology_files(self, root_output_path, model, n_jets, numerator_call, final_state_particle_ids=(),jet_ids=None, write_yaml=True, bar=None):

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
        particle_masses = {e['name']: model['parameter_dict'][model.get_particle(e['PDG']).get('mass')].real for e in self.edges.values()}

        num_incoming = sum(1 for e in edge_map_lin if e[0][0] == 'q') // 2

        if num_incoming == 1:
            external_momenta = {'q1': [500., 0., 0., 0.], 'q2': [500., 0., 0., 0.]}
            #external_momenta = {'q1': [1., 0., 0., 0.], 'q2': [1., 0., 0., 0.]}
            p = np.array(external_momenta['q1'])
        else:
            external_momenta = {'q1': [500., 0., 0., 500.], 'q2': [500., 0., 0., -500.], 'q3': [500., 0., 0., 500.], 'q4': [500., 0., 0., -500.]}
            #external_momenta = {'q1': [1., 0., 0., 1.], 'q2': [1., 0., 0., -1.], 'q3': [1., 0., 0., 1.], 'q4': [1., 0., 0., -1.]}
            p = np.array(external_momenta['q1']) + np.array(external_momenta['q2'])

        # compute mUV
        FORM_processing_options['mUV'] = 2 * math.sqrt(p[0]**2 - p[1]**2 - p[2]**2 - p[3]**2)
        FORM_processing_options['logmUV'] = math.log(FORM_processing_options['mUV']**2)

        loop_momenta = []
        n_loops = len(self.edges) - len(self.nodes) + 1
        for loop_var in range(n_loops):
            lm = next((ee['name'], ee['signature'][0][loop_var]) for ee in self.edges.values() if all(s == 0 for s in ee['signature'][1]) and \
                sum(abs(s) for s in ee['signature'][0]) == 1 and ee['signature'][0][loop_var] != 0)
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
            edge_weights={e['name']: self.get_edge_scaling(e['PDG']) for e in self.edges.values()},
            vertex_weights={nv: self.get_node_scaling(n['PDGs']) for nv, n in self.nodes.items()},
        )
        # check if cut is possible
        if len(topo.cuts) == 0:
            logger.info("No cuts for graph {}".format(self.name))
            return False

        self.generate_replacement_rules(topo)

        # Also generate the squared topology yaml files of all additional LMB topologies used for cross-check
        if isinstance(self.additional_lmbs, list):
            for i_lmb, (_,_,_,other_lmb_supergraph) in enumerate(self.additional_lmbs):
                if bar:
                    bar.update(i_lmb='%d'%(i_lmb+2))
                other_lmb_supergraph.generate_squared_topology_files(root_output_path, model, n_jets, numerator_call, 
                                    final_state_particle_ids=final_state_particle_ids,jet_ids=jet_ids, write_yaml=write_yaml)

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

                        # numerator derivative
                        if ext_mom == der_mom:
                            index = next(i for i, bub in enumerate(bubble_to_cut.keys()) if bub == diag_info['bubble_momenta'])
                            trans.append('der(pbubble' + str(index) + ')')
                        else:
                            # check if we pick up a sign change due to the external momentum flowing in the opposite direction
                            signs = [se * sc for se, sc in zip(der_sig[0] + der_sig[1], ext_sig[0] + ext_sig[1]) if se * sc != 0]
                            assert(len(set(signs)) == 1)
                            trans.append('-2*({}{})'.format('+' if signs[0] == 1 else '-', der_mom))
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
                            power = 2 if edge_name == der_edge else 1
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
        for i, (bubble_edges, (cut_edge, bubble_derivative_edges)) in enumerate(bubble_to_cut.items()):
            ce = next(ee for ee in self.edges.values() if ee['name'] == cut_edge)
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

    def generate_numerator_functions(self, additional_overall_factor='', output_format='c', workspace=None, FORM_vars=None, active_graph=None):
        """ Use form to plugin Feynman Rules and process the numerator algebra so as
        to generate a low-level routine in file_path that encodes the numerator of this supergraph."""

        _MANDATORY_FORM_VARIABLES = ['SGID','NINITIALMOMENTA','NFINALMOMENTA']

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

        r = subprocess.run(' '.join([
                FORM_processing_options["FORM_path"],
                ]+
                [ '-D %s=%s'%(k,v) for k,v in FORM_vars.items() ]+
                [ FORM_source, ]
            ),
            shell=True,
            cwd=selected_workspace,
            capture_output=True)
        if r.returncode != 0:
            raise FormProcessingError("FORM processing failed with error:\n%s"%(r.stdout.decode('UTF-8')))

        # return the code for the numerators
        if not os.path.isfile(pjoin(selected_workspace,'out_%d.proto_c'%i_graph)):
            raise FormProcessingError("FORM failed to produce an output for super graph ID=%d. Output file not found at '%s'."%
                                                                    (i_graph,pjoin(selected_workspace,'out_%d.proto_c'%i_graph)))

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
        output_match = re.compile("isoF=[-10]+;")
        factor_match = re.compile("[-10]+")
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
                FORM_processing_options["FORM_path"],
                ]+
                [ '-D %s=%s'%(k,v) for k,v in FORM_vars.items() ]+
                [ form_source, ]
            ),
            shell=True,
            cwd=workspace,
            capture_output=True)
            if r.returncode != 0:
                raise FormProcessingError("FORM processing failed with error:\n%s"%(r.stdout.decode('UTF-8')))

            output = r.stdout.decode('UTF-8').replace(' ','').replace('\n','')
            factor = int(factor_match.findall(output_match.findall(output)[0])[0])

            if factor == 0:
                raise FormProcessingError("Multiplicity not found: {} =/= (+/-) * {}. (iso_check_%(SGID)d_%(ID0)d_%(IDn)d)".format(self[0].name,g.name )%FORM_vars)
            elif factor == 10:
                multiplicity = 0
                continue
            elif factor == -1 or factor == 1:
                multiplicity += factor
            else:
                raise FormProcessingError("Unknown isoF for multiplicity factor : usiF={}".format(factor))

            #logger.info("{} = ({:+d}) * {}".format(self[0].name, factor, g.name ))
        return multiplicity    

    def generate_squared_topology_files(self, root_output_path, model, n_jets, numerator_call, final_state_particle_ids=(), jet_ids=None, bar=None ):
        for i, g in enumerate(self):
            # Now we generate the squared topology only for the first isomorphic graph
            # to obtain the replacement rules for the bubble.
            # The other elements of the isomorphic set are going to contribute only at the 
            # numerator 
            if i==0:
                r = g.generate_squared_topology_files(root_output_path, model, n_jets, numerator_call, 
                                final_state_particle_ids, jet_ids=jet_ids, write_yaml=i==0, bar=bar)
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
    def from_dict(cls, dict_file_path, first=None, merge_isomorphic_graphs=False, verbose=False, model = None, workspace=None):
        """ Creates a FORMSuperGraph list from a dict file path."""
        from pathlib import Path
        p = Path(dict_file_path)
        sys.path.insert(0, str(p.parent))
        m = __import__(p.stem)

        logger.info("Imported {} supergraphs.".format(len(m.graphs)))

        # Filter specific graphs by name 
        #filter_graphs = ['SG_QG8','SG_QG9']
        #m.graphs = [ g for (g,name) in zip(m.graphs, m.graph_names) if name in filter_graphs]
        #m.graph_names = ['SG_MG8','SG_QG9']
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

        # Check for the rooting of the loop momenta form the QGRAF OUTPUT
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

        if first is not None:
            logger.info("Taking first {} supergraphs.".format(first))
            full_graph_list = full_graph_list[:first]

        
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

        logger.info("{} unique supergraphs".format(len(iso_groups)))

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

    def generate_numerator_functions(self, root_output_path, additional_overall_factor='', params={}, output_format='c', workspace=None, header=""):
        header_map = {'header': header}
        """ Generates optimised source code for the graph numerator in several
        files rooted in the specified root_output_path."""

        if len(self)==0:
            raise FormProcessingError("Cannot generat numerators for an empty list of supergraphs.")

        if output_format not in self.extension_names:
            raise FormProcessingError("This FORMSuperGraphList instance requires at least one entry for generating numerators.")

        # add all numerators in one file and write the headers
        numerator_header = """#include <math.h>
#include <complex.h>
#include <signal.h>
#include "%(header)snumerator.h"
"""

        var_pattern = re.compile(r'Z\d*_')
        input_pattern = re.compile(r'lm(\d*)')
        energy_exp = re.compile(r'f\(([^)]*)\)\n')
        split_number = re.compile(r'\\\n\s*')
        return_exp = re.compile(r'return ([^;]*);\n')

        FORM_vars={}

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
                    num = graph.generate_numerator_functions(additional_overall_factor,workspace=workspace, FORM_vars=FORM_vars, active_graph=active_graph)
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

                        confs.append((conf_id, rank))

                        if len(temp_vars) > 0:
                            graph.is_zero = False

                        # TODO: Check that static inline vs inline induces no regression!
                        numerator_main_code += '\n// polynomial in {}'.format(','.join(conf[1:]))
                        numerator_main_code += '\nstatic inline int %(header)sevaluate_{}_{}(double complex lm[], double complex* out) {{\n\t{}'.format(i, conf_id,
                            'double complex {};'.format(','.join(temp_vars)) if len(temp_vars) > 0 else ''
                        ) + num_body + '\n\treturn {}; \n}}\n'.format(max_index + 1)


                    numerator_main_code += \
"""
int %(header)sevaluate_{}(double complex lm[], int conf, double complex* out) {{
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
        ['\t\tcase {}: return %(header)sevaluate_{}_{}(lm, out);'.format(conf, i, conf) for conf, _ in sorted(confs)] +
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
#include <complex.h>
#include <signal.h>

{}
{}

int %(header)sevaluate(double complex lm[], int diag, int conf, double complex* out) {{
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
    '\n'.join('int %(header)sevaluate_{}(double complex[], int conf, double complex*);'.format(i) for i in all_numerator_ids),
    '\n'.join('int %(header)sget_rank_{}(int conf);'.format(i) for i in all_numerator_ids),
    '\n'.join(
    ['\t\tcase {}: return %(header)sevaluate_{}(lm, conf, out);'.format(i, i) for i in all_numerator_ids]+
    ['\t\tdefault: raise(SIGABRT);']
    ), 
    (max_buffer_size + 1) * 2,
    '\n'.join(
    ['\t\tcase {}: return %(header)sget_rank_{}(conf);'.format(i, i) for i in all_numerator_ids]+
    ['\t\tdefault: raise(SIGABRT);']
    )
    )
        writers.CPPWriter(pjoin(root_output_path, '%(header)snumerator.c'%header_map)).write(numerator_code%header_map)


        params['mUV'] = FORM_processing_options['mUV']
        params['logmUV'] = FORM_processing_options['logmUV']
        header_code = \
"""
#ifndef NUM_H
#define NUM_H

{}

#endif
""".format('\n'.join('#define  {} {}'.format(k, v) for k, v in params.items()))

        writers.CPPWriter(pjoin(root_output_path, '%(header)snumerator.h'%header_map)).write(header_code%header_map)

    def generate_squared_topology_files(self, root_output_path, model, n_jets, final_state_particle_ids=(), jet_ids=None, filter_non_contributing_graphs=True):
        topo_collection = {
            'name': self.name,
            'topologies': []
        }

        # Set sensible jet_ids if none
        if jet_ids is None:
            jet_ids=tuple(list(range(1,6))+list(range(-1,-6,-1))+[21,82,-82,])

        contributing_supergraphs = []

        with progressbar.ProgressBar(prefix='Generating squared topology files (graph #{variables.i_graph}/%d, LMB #{variables.i_lmb}/{variables.max_lmb}, {variables.timing} ms / supergraph) : '%len(self), 
                max_value=len(self), variables = {'timing' : '0', 'i_graph' : '0', 'i_lmb': '0', 'max_lmb' : '0'} ) as bar:
            non_zero_graph = 0
            total_time=0.0
            for i, g in enumerate(self):
                time_before = time.time()
                bar.update(i_graph='%d'%(i+1))
                if g.generate_squared_topology_files(root_output_path, model, n_jets, numerator_call=non_zero_graph, 
                                                            final_state_particle_ids=final_state_particle_ids,jet_ids=jet_ids, bar=bar):
                    topo_collection['topologies'].append({
                        'name': g[0].name,
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


        try:
            import yaml
            from yaml import Loader, Dumper
        except ImportError:
            raise BaseException("Install yaml python module in order to import topologies from yaml.")

        open(pjoin(root_output_path, self.name + '.yaml'), 'w').write(yaml.dump(topo_collection, Dumper=Dumper))

        if filter_non_contributing_graphs:
            self[:] = contributing_supergraphs

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

    def draw(self, output_dir):
        """ For now simply one Mathematica script per supergraph."""

        for i_graph, super_graphs in enumerate(self.super_graphs_list):
            super_graphs[0].draw(self.model, output_dir, FORM_id=i_graph)
            # Draw supergraphs for additional LMBs
            if isinstance(super_graphs[0].additional_lmbs,list):
                for i_lmb,_,_,sg in super_graphs[0].additional_lmbs:
                    sg.draw(self.model, output_dir, FORM_id=i_graph, lmb_id=i_lmb)

    def generate_numerator_functions(self, root_output_path, output_format='c',workspace=None, header=""):
        assert(header in ['MG', 'QG', ''])

        params = {
            'mass_t': self.model['parameter_dict'][self.model.get_particle(6).get('mass')].real,
            'gs': self.model['parameter_dict']['G'].real,
            'ge': math.sqrt(4. * math.pi / self.model['parameter_dict']['aEWM1'].real),
            'gy': self.model['parameter_dict']['mdl_yt'].real / math.sqrt(2.),
            'ghhh': 6. * self.model['parameter_dict']['mdl_lam'].real,
        }

        helicity_averaging_factor = 1
        for leg in self.repr_process.get('legs'):
            # Skip final states
            if leg.get('state') is True:
                continue

            helicity_averaging_factor *= len(self.model.get_particle(leg.get('id')).get_helicity_states())
        helicity_averaging_factor = "/" + str(helicity_averaging_factor)

        additional_overall_factor = helicity_averaging_factor
        return self.super_graphs_list.generate_numerator_functions(
            root_output_path, output_format=output_format,
            additional_overall_factor=additional_overall_factor,
            params=params,workspace=workspace, header=header)

    @classmethod
    def compile(cls, root_output_path):

        if os.path.isfile(pjoin(root_output_path,'Makefile')):
            try:
                logger.info("Now compiling FORM-generated numerators...")
                misc.compile(cwd=root_output_path,mode='cpp', nb_core=FORM_processing_options["cores"])
            except MadGraph5Error as e:
                logger.info("%sCompilation of FORM-generated numerator failed:\n%s%s"%(
                    utils.bcolors.RED,str(e),utils.bcolors.ENDC))
        else:
            logger.warning(("\n%sYou are running FORM_processing directly from the __main__ of FORM_processing.py.\n"+
                           "You will thus need to compile numerators.c manually.%s")%(utils.bcolors.GREEN, utils.bcolors.ENDC))

    def generate_squared_topology_files(self, root_output_path, n_jets, final_state_particle_ids=(), jet_ids=None, filter_non_contributing_graphs=True):
        self.super_graphs_list.generate_squared_topology_files(
            root_output_path, self.model, n_jets, final_state_particle_ids, jet_ids=jet_ids, filter_non_contributing_graphs=filter_non_contributing_graphs
        )


if __name__ == "__main__":
   
    parser = argparse.ArgumentParser(description='Generate numerators with FORM and yaml for Rust.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('diagrams_python_source', default=None, type=str,
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
    FORM_processing_options['extra-options'] = '-D OPTIMITERATIONS=' + str(args.optim_iter)

    import alpha_loop.interface as interface
    cli = interface.alphaLoopInterface()

    cli.do_import('model %s'%args.model)
    computed_model = model_reader.ModelReader(cli._curr_model)
    computed_model.set_parameters_and_couplings(args.restrict_card)        
    process_definition=cli.extract_process(args.process, proc_number=0)

    super_graph_list = FORMSuperGraphList.from_dict(args.diagrams_python_source)
    form_processor = FORMProcessor(super_graph_list, computed_model, process_definition)
    TMP_OUTPUT = pjoin(root_path,'TMPDIR')
    Path(TMP_OUTPUT).mkdir(parents=True, exist_ok=True)
    form_processor.generate_squared_topology_files(TMP_OUTPUT, 0, final_state_particle_ids=(6, 6, 25, 25))
    form_processor.generate_numerator_functions(TMP_OUTPUT, output_format='c')
