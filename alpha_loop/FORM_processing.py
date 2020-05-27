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

import progressbar
from itertools import chain
import sys
import subprocess
import argparse
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

import madgraph.various.misc as misc
import madgraph.iolibs.file_writers as writers
from madgraph import MadGraph5Error, InvalidCmd, MG5DIR
import models.model_reader as model_reader

import LTD.squared_topologies
import LTD.ltd_utils

logger = logging.getLogger('alphaLoop.FORM_processing')

if __name__ == "__main__":
    logging.basicConfig()
    logger.setLevel(logging.INFO)

plugin_path = os.path.dirname(os.path.realpath( __file__ ))


FORM_processing_options = {'FORM_path': 'form', 'TFORM_path': 'tform', 'parallel': False, 'cores': '-w1', 'extra-options': '-D OPTIMITERATIONS=100'}

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
        # FxFV
        ( -2, 1, 2 ): (0, 1, 2),
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
    ):
        """ initialize a FORM SuperGraph from several options."""

        self.is_zero = False
        self.edges = edges
        self.nodes = nodes
        self.overall_factor = overall_factor
        # A hashable call signature
        self.call_identifier = call_identifier
        if name is None:
            self.name = str(self.call_identifier)
        else:
            self.name = name

    def get_mathematica_rendering_code(self, model, FORM_id=None):
        """ Generate mathematica expression for drawing this graph."""

        repl_dict = {
            'edge_font_size'     : 10,
            'vetex_font_size'    : 10,
            'width'              : self._rendering_size[0],
            'height'             : self._rendering_size[1],
            'graph_layout_strategy' : self._graph_layout_strategy
        }
        graph_name='MG: %s'%self.name
        if FORM_id is not None:
            graph_name +=' | FORM: #%d'%FORM_id
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
            edge_repl_dict = {}
            is_LMB = ('name' in edge_data and 'LMB' in str(edge_data['name']).upper())
            if is_LMB:
                edge_repl_dict['thickness'] = 0.004
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
            elif edge_data['PDG'] in [82,]:
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
    
    def draw(self, model, output_dir,FORM_id=None):
        """ Outputs the mathematica code for rendering this FORMSuperGraph."""
        
        if FORM_id is not None:
            file_name = 'Graph_%04d'%FORM_id
        else:
            file_name = 'Graph_%s'%self.name

        MM_code = \
"""GraphClass = If[$VersionNumber > 12, EdgeTaggedGraph, Graph];
CreateEdge[u_,v_,t_]:=If[$VersionNumber > 12, DirectedEdge[u, v, t], DirectedEdge[u, v]];
aGraph=%s;
"""%self.get_mathematica_rendering_code(model,FORM_id=FORM_id)
        # Export to PDF in landscape format. One graph per page for now.
        # The 1.2 multiplier accounts for margins
        MM_code += 'Export["%s.pdf", GraphicsGrid[{{aGraph}}], ImageSize -> {%f, %f}];'%(
                    file_name,1.25*self._rendering_size[0],1.25*self._rendering_size[1])
        open(pjoin(output_dir,'%s.m'%file_name),'w').write(MM_code)

    def generate_numerator_form_input(self, additional_overall_factor=''):
        # create the input file for FORM
        form_diag = self.overall_factor + additional_overall_factor
        for node in self.nodes.values():
            if node['vertex_id'] < 0:
                continue
        
            form_diag += '*vx({},{},{})'.format(
                ','.join(str(p) for p in node['PDGs']),
                ','.join(node['momenta']),
                ','.join(str(i) for i in node['indices']),
            )

        for edge in self.edges.values():
            form_diag += '*prop({},{},{},{})'.format(
                edge['PDG'],
                edge['type'],
                edge['momentum'],
                ','.join(str(i) for i in edge['indices']),
            )

        return form_diag


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

        # TODO: correctly assign the overall factor.
        overall_factor = "1"

        model = LTD2_super_graph.model
        local_graph = copy.deepcopy(LTD2_super_graph.graph)

        # Collect the outer-most nodes
        outer_in_nodes = []
        outer_out_nodes = []
        for node_key, node_data in local_graph.nodes.items():
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
                dict(in_edge_data) for u,v,in_edge_data in 
                local_graph.in_edges(node_key,data=True)
            ]
            adjacent_out_edges = [
                dict(out_edge_data) for u,v,out_edge_data in 
                local_graph.out_edges(node_key,data=True)
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

            # Example printout information about the vertex
            #misc.sprint("Vertex of node %s:\n%s"%(
            #        str(node_key),pformat(model.get_interaction(node_data['vertex_id']))))

        # Finally overwrite the edge momentum so as to be a string
        for edge_key, edge_data in local_graph.edges.items():
            edge_data['momentum'] =cls.momenta_decomposition_to_string(edge_data['momentum'])
            particle = model.get_particle(edge_data['PDG'])
            # In the above conventions the fermion are going against their flow, so we need
            # to flip the order of their fundamental/antifundamental indices so that FORM
            # builds the correct propagator. 
            if len(edge_data['indices'])>1 and particle.get('spin')==2:
                if particle.get('is_part'):
                    edge_data['indices'] = tuple([edge_data['indices'][1],edge_data['indices'][0]])

        graph_name = 'P%(proc_id)dL%(left_diagram_id)dR%(right_diagram_id)d'%LTD2_super_graph.call_signature
        form_super_graph =  cls(
            name = graph_name,
            call_identifier=LTD2_super_graph.call_signature,
            edges = dict(local_graph.edges),
            nodes = dict(local_graph.nodes),
            overall_factor = overall_factor,
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
            'overall_factor' : self.overall_factor
        }
        if file_path:
            open(file_path,'w').write(pformat(dict_to_dump))
        else:
            return dict_to_dump

    def generate_squared_topology_files(self, root_output_path, n_jets, numerator_call, final_state_particle_ids=()):
        if self.is_zero:
            return False

        # the first 4 entries are the external momenta
        edge_map_lin = [(e['name'] if e['type'] == 'virtual' else 'q' + e['name'][1:], e['vertices'][0], e['vertices'][1]) for e in self.edges.values()]
        assert(e[0] != 'q' or int(e[1:]) < 5 for e in edge_map_lin)

        particle_ids = {e['name'] if e['type'] == 'virtual' else 'q' + e['name'][1:]: abs(e['PDG']) for e in self.edges.values()}

        external_momenta = {'q1': [1., 0., 0., 1.], 'q2': [1., 0., 0., -1.], 'q3': [1., 0., 0., 1.], 'q4': [1., 0., 0., -1.]}
        num_incoming = sum(1 for e in edge_map_lin if e[0][0] == 'q') // 2

        loop_momenta = []
        n_loops = len(self.edges) - len(self.nodes) + 1
        for loop_var in range(n_loops):
            # FIXME: what if the edge is -k?
            lm = next(ee['name'] for ee in self.edges.values() if all(s == 0 for s in ee['signature'][1]) and \
                sum(abs(s) for s in ee['signature'][0]) == 1 and ee['signature'][0][loop_var] != 0)
            loop_momenta.append(lm)

        topo = LTD.squared_topologies.SquaredTopologyGenerator(edge_map_lin,
            self.name, ['q1', 'q2'][:num_incoming], n_jets, external_momenta,
            loop_momenta_names=tuple(loop_momenta),
            particle_ids=particle_ids,
            final_state_particle_ids=final_state_particle_ids,
            overall_numerator=1.0,
            numerator_structure={},
            FORM_numerator={'call_signature': {'id': numerator_call}})

        # check if cut is possible
        if len(topo.cuts) == 0:
            logger.info("No cuts for graph {}".format(self.name))
            return False

        topo.export(pjoin(root_output_path, self.name + ".yaml"))

        return True


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

    def generate_numerator_functions(self, additional_overall_factor='', output_format='c'):
        """ Use form to plugin Feynman Rules and process the numerator algebra so as
        to generate a low-level routine in file_path that encodes the numerator of this supergraph."""

        # write the form input to a file
        form_input = self.generate_numerator_form_input(additional_overall_factor)

        with open(pjoin(FORM_workspace,'input.h'), 'w') as f:
            f.write('L F = {};'.format(form_input))

        r = subprocess.run([
                FORM_processing_options["TFORM_path"] if FORM_processing_options["parallel"] else FORM_processing_options["FORM_path"],
                FORM_processing_options["cores"],
                pjoin(plugin_path,"numerator.frm")
            ],
            cwd=FORM_workspace,
            capture_output=True)
        if r.returncode != 0:
            raise FormProcessingError("FORM processing failed with error:\n%s"%(r.stdout.decode('UTF-8')))

        # return the code for the numerators
        with open(pjoin(FORM_workspace,'out.proto_c'), 'r') as f:
            num_code = f.read()

        return num_code

    def to_dict(self, file_path=None):
        """ Store that into a dict."""
        to_dump = [g.to_dict() for g in self]
        if file_path:
            open(file_path,'w').write(pformat(to_dump))
        else:
            return to_dump

    def generate_squared_topology_files(self, root_output_path, n_jets, numerator_call, final_state_particle_ids=() ):
        # only generate a yaml file for the reference graph of the isomorphic set
        return self[0].generate_squared_topology_files(root_output_path, n_jets, numerator_call, final_state_particle_ids)

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
    def from_dict(cls, dict_file_path, first=None, merge_isomorphic_graphs=True):
        """ Creates a FORMSuperGraph list from a dict file path."""
        from pathlib import Path
        p = Path(dict_file_path)
        sys.path.insert(0, str(p.parent))
        m = __import__(p.stem)

        logger.info("Imported {} supergraphs.".format(len(m.graphs)))

        graph_list = []
        for i, g in enumerate(m.graphs):
            # convert to FORM supergraph
            form_graph = FORMSuperGraph(name=p.stem + '_' + str(i), edges = g['edges'], nodes=g['nodes'], overall_factor=g['overall_factor'])
            form_graph.derive_signatures()
            graph_list.append(form_graph)

        if first is not None:
            logger.info("Taking first {} supergraphs.".format(first))
            graph_list = graph_list[:first]

        iso_groups = []

        import time
        import sympy as sp

        # group all isomorphic graphs
        pdg_primes = {pdg : sp.prime(i + 1) for i, pdg in enumerate([1,2,3,4,5,6,11,12,13,21,22,25,82])}
        for graph in graph_list:
            g = igraph.Graph()
            g.add_vertices(len(graph.nodes))
            undirected_edges = set()

            for e in graph.edges.values():
                undirected_edges.add(tuple(sorted(e['vertices'])))

            edge_colors = []
            for ue in undirected_edges:
                e_color = 1
                for e in graph.edges.values():
                    if tuple(sorted(e['vertices'])) == ue:
                        e_color *= pdg_primes[abs(e['PDG'])]
                edge_colors.append(e_color)
                g.add_edges([tuple(sorted([x - 1 for x in ue]))])

            for (ref_graph, ref_colors), graph_list in iso_groups:
                v_maps = ref_graph.get_isomorphisms_vf2(g, edge_color1=ref_colors, edge_color2=edge_colors)
                if v_maps != [] and merge_isomorphic_graphs:
                    # map the signature from the new graph to the reference graph using the vertex map
                    edge_map = [next(re.index for re in ref_graph.es if (re.source, re.target) == (v_maps[0][e.source], v_maps[0][e.target]) or \
                        (re.source, re.target) == (v_maps[0][e.target], v_maps[0][e.source])) for e in g.es]

                    # go through all loop momenta and make the matrix
                    mat = []
                    n_loops = len(graph_list[0].edges) - len(graph_list[0].nodes) + 1
                    for loop_var in range(n_loops):
                        # find the edge that has (+-)k`loop_var`
                        orig_edge = next(ee for ee in graph_list[0].edges.values() if all(s == 0 for s in ee['signature'][1]) and \
                            sum(abs(s) for s in ee['signature'][0]) == 1 and ee['signature'][0][loop_var] != 0)
                        orig_edge_in_ref = next(ee for ee in ref_graph.es if tuple(sorted(orig_edge['vertices'])) == (ee.source + 1, ee.target + 1))

                        # now get the corresponding edge(s) in the new graph
                        en = g.es[edge_map.index(orig_edge_in_ref.index)]
                        new_edges = [ee for ee in graph.edges.values() if ee['vertices'] == (en.source + 1, en.target + 1) or ee['vertices'] == (en.target + 1, en.source + 1)]
                        # in the case of a multi-edge, we simply take the one with the shortest momentum string
                        new_edge = sorted(new_edges, key=lambda x: len(x['momentum']))[0]
                        mat.append([s * orig_edge['signature'][0][loop_var] for s in new_edge['signature'][0]])

                    # now map all signatures and momentum labels of the new graph
                    mat = np.array(mat).T
                    for e in graph.edges.values():
                        # skip external edges
                        if e["type"] != "virtual":
                            continue

                        new_sig = mat.dot(np.array(e['signature'][0]))
                        e['signature'] = (list(new_sig), e['signature'][1])
                        e['momentum'] = graph.momenta_decomposition_to_string(e['signature'], set_outgoing_equal_to_incoming=False)

                        # TODO: replace this by a loop over nodes once the nodes store the edge information
                        v1 = graph.nodes[e['vertices'][0]]
                        edge_index_in_v = v1['indices'].index(e['indices'][0]) if e['indices'][0] in v1['indices'] else v1['indices'].index(e['indices'][1])
                        moms = list(v1['momenta'])
                        if e['indices'][0] in v1['indices']:
                            moms[edge_index_in_v] = e['momentum']
                        else:
                            sig = ([-s for s in e['signature'][0]], [-s for s in e['signature'][1]])
                            new_momentum = graph.momenta_decomposition_to_string(sig, set_outgoing_equal_to_incoming=False)
                            moms[edge_index_in_v] = new_momentum
                        v1['momenta'] = tuple(moms)

                    graph_list.append(graph)
                    break
            else:
                iso_groups.append(((g, edge_colors), [graph]))

        logger.info("{} unique supergraphs".format(len(iso_groups)))

        return FORMSuperGraphList([FORMSuperGraphIsomorphicList(iso_group) for _, iso_group in iso_groups], name=p.stem)

    def to_dict(self, file_path):
        """ Outputs the FORMSuperGraph list to a Python dict file path."""

        # TODO: Dump to Python dict
        pass

    def generate_numerator_functions(self, root_output_path, additional_overall_factor='', params={}, output_format='c'):
        """ Generates optimised source code for the graph numerator in several
        files rooted in the specified root_output_path."""

        if len(self)==0:
            raise FormProcessingError("Unsupported output format for numerator functions: '%s'"%output_format)

        if output_format not in self.extension_names:
            raise FormProcessingError("This FORMSuperGraphList instance requires at least one entry for generating numerators.")

        # add all numerators in one file and write the headers
        numerator_code = """#include <math.h>
#include <complex.h>
#include <signal.h>

{}
""".format('\n'.join('const double {} = {};'.format(k, v) for k, v in params.items()))

        pattern = re.compile(r'Z(()\d*)_')
        input_pattern = re.compile(r'lm(\d*)')

        graphs_zero = True
        with progressbar.ProgressBar(
            prefix = 'Processing numerators with FORM ({variables.timing} ms / supergraph) : ',
            max_value=len(self),
            variables = {'timing' : '0'}
            ) as bar:
            total_time = 0.
            for i, graph in enumerate(self):
                time_before = time.time()
                num = graph.generate_numerator_functions(additional_overall_factor)
                total_time += time.time()-time_before
                num = num.replace('i_', 'I')
                num = input_pattern.sub(r'lm[\1]', num)
                num = num.replace('\nZ', '\n\tZ') # nicer indentation
                max_intermediate_variable = max((int(index.groups()[0]) for index in pattern.finditer(num)), default=0)

                if max_intermediate_variable > 0:
                    graph.is_zero = False

                numerator_code += '\ndouble complex evaluate_{}(double complex lm[]) {{\n\t{}\n'.format(i,
                    'double complex {};'.format(','.join('Z' + str(i) + '_' for i in range(1, max_intermediate_variable + 1))) if max_intermediate_variable > 0 else ''
                ) + num + '}\n'

                bar.update(timing='%d'%int((total_time/float(i+1))*1000.0))
                bar.update(i+1)

        if graphs_zero:
            self[0].is_zero = True

        numerator_code += \
"""
double complex evaluate(double complex lm[], int i) {{
    switch(i) {{
{}
    }}
}}
""".format('\n'.join(
    ['\t\tcase {}: return evaluate_{}(lm);'.format(i, i) for i in range(len(self))]+
    ['\t\tdefault: raise(SIGABRT);']
    ))
        writers.CPPWriter(pjoin(root_output_path, 'numerator.c')).write(numerator_code)

    def generate_squared_topology_files(self, root_output_path, n_jets, final_state_particle_ids=()):
        topo_collection = {
            'name': self.name,
            'topologies': []
        }

        with progressbar.ProgressBar(prefix='Generating squared topology files : ', max_value=len(self)) as bar:
            for i, g in enumerate(self):
                if g.generate_squared_topology_files(root_output_path, n_jets, numerator_call=i, final_state_particle_ids=final_state_particle_ids):
                    topo_collection['topologies'].append({
                        'name': g[0].name,
                        'multiplicity': 1
                    })

                bar.update(i+1)

        try:
            import yaml
            from yaml import Loader, Dumper
        except ImportError:
            raise BaseException("Install yaml python module in order to import topologies from yaml.")

        open(pjoin(root_output_path, self.name + '.yaml'), 'w').write(yaml.dump(topo_collection, Dumper=Dumper))

class FORMProcessor(object):
    """ A class for taking care of the processing of a list of FORMSuperGraphList.
    Useful because many aspects common to all supergraphs and function do not belong to FORMSuperGraphList.
    """

    def __init__(self, super_graphs_list, model, process_definition):
        """ Specify aditional information such as the model that is useful for FORM processing."""
        self.super_graphs_list = super_graphs_list
        self.model = model
        self.process_definition = process_definition

    def draw(self, output_dir):
        """ For now simply one Mathematica script per supergraph."""

        for i_graph, super_graphs in enumerate(self.super_graphs_list):
            super_graphs[0].draw(self.model, output_dir, FORM_id=i_graph)

    def generate_numerator_functions(self, root_output_path, output_format='c'):
        params = {
            'mass_t': self.model['parameter_dict'][self.model.get_particle(6).get('mass')].real,
            'gs': self.model['parameter_dict']['G'].real,
            'ge': math.sqrt(4. * math.pi / self.model['parameter_dict']['aEWM1'].real),
            'gy': self.model['parameter_dict']['mdl_yt'].real / math.sqrt(2.),
        }

        helicity_averaging_factor = 1
        for leg in self.process_definition.get('legs'):
            # Skip final states
            if leg.get('state') is True:
                continue

            helicity_averaging_factor *= len(self.model.get_particle(leg.get('id')).get_helicity_states())
        helicity_averaging_factor = "/" + str(helicity_averaging_factor)

        # there is another 1/4 difference between FORM and MG that is unexplained
        additional_overall_factor = helicity_averaging_factor + '/4'

        return self.super_graphs_list.generate_numerator_functions(
            root_output_path, output_format=output_format,
            additional_overall_factor=additional_overall_factor,
            params=params)

    @classmethod
    def compile(cls, root_output_path):

        if os.path.isfile(pjoin(root_output_path,'Makefile')):
            try:
                logger.info("Now compiling FORM-generated numerators...")
                misc.compile(cwd=root_output_path,mode='cpp')
            except MadGraph5Error as e:
                logger.info("%sCompilation of FORM-generated numerator failed:\n%s%s"%(
                    utils.bcolors.RED,str(e),utils.bcolors.ENDC))
        else:
            logger.warning(("\n%sYou are running FORM_processing directly from the __main__ of FORM_processing.py.\n"+
                           "You will thus need to compile numerators.c manually.%s")%(utils.bcolors.GREEN, utils.bcolors.ENDC))

    def generate_squared_topology_files(self, root_output_path, n_jets, final_state_particle_ids=()):
        self.super_graphs_list.generate_squared_topology_files(root_output_path, n_jets, final_state_particle_ids)


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
        FORM_processing_options['parallel'] = True
        FORM_processing_options['cores'] = '-w' + str(args.cores)
    FORM_processing_options['extra-options'] = '-D OPTIMITERATIONS=' + str(args.optim_iter)
     
    import alpha_loop.interface as interface
    cli = interface.alphaLoopInterface()

    cli.do_import('model %s'%args.model)
    computed_model = model_reader.ModelReader(cli._curr_model)
    computed_model.set_parameters_and_couplings(args.restrict_card)        
    process_definition=cli.extract_process(args.process, proc_number=0)

    super_graph_list = FORMSuperGraphList.from_dict(args.diagrams_python_source)
    form_processor = FORMProcessor(super_graph_list, computed_model, process_definition)
    form_processor.generate_numerator_functions('./lib', output_format='c')
    form_processor.generate_squared_topology_files('./RUST_input', 2, final_state_particle_ids=(6, 6, 25))
