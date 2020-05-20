#!/usr/bin/env python3

import copy
import logging
import os
from pathlib import Path
from pprint import pprint, pformat

from itertools import chain
import sys
import subprocess
pjoin = os.path.join

if __name__ == "__main__":
    root_path = os.path.dirname(os.path.realpath( __file__ ))
    sys.path.insert(0, pjoin(root_path,os.path.pardir))
    sys.path.insert(0, pjoin(root_path,os.path.pardir,os.path.pardir,os.path.pardir))
import alpha_loop.utils as utils
import re

import madgraph.various.misc as misc
import madgraph.iolibs.file_writers as writers
from madgraph import MadGraph5Error, InvalidCmd, MG5DIR
import models.model_reader as model_reader

logger = logging.getLogger('alphaLoop.FORM_processing')

plugin_path = os.path.dirname(os.path.realpath( __file__ ))


FORM_processing_options = {'FORM_path': 'form'}

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
        # VxVV (e.g. W+ W- a )
        ( -3, 3, 3 ): (1, 0, 2),
    }

    def __init__(self, *args,
        call_identifier=None,
        name=None,
        edges=None,
        nodes=None,
        overall_factor="1",
    ):
        """ initialize a FORM SuperGraph from several options."""

        self.edges = edges
        self.nodes = nodes
        self.overall_factor = overall_factor
        # A hashable call signature
        self.call_identifier = call_identifier
        if name is None:
            self.name = str(self.call_identifier)
        else:
            self.name = name

    def generate_numerator_form_input(self):
        # make sure all indices are positive by adding an offset
        index_offset = 0
        #for node in self.nodes.values():
        #    index_offset = min(index_offset, min(node['indices']))
        #index_offset = -index_offset + 1

        # create the input file for FORM
        form_diag = self.overall_factor
        for node in self.nodes.values():
            if node['vertex_id'] < 0:
                continue
        
            form_diag += '*vx({},{},{})'.format(
                ','.join(str(p) for p in node['PDGs']),
                ','.join(node['momenta']),
                ','.join(str(i + index_offset) for i in node['indices']),
            )

        for edge in self.edges.values():
            form_diag += '*prop({},{},{},{})'.format(
                edge['PDG'],
                edge['type'],
                edge['momentum'],
                ','.join(str(i + index_offset) for i in edge['indices']),
            )

        return form_diag

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
            raise FormProcessingError("Conventions for FORM Feynman rules of signature '%s' not specifed."%canonical_identifier)
        
        new_position = cls._FORM_Feynman_rules_conventions[canonical_identifier]

        return [ edges_to_sort[identities[position][1]] for position in new_position ]

    @classmethod
    def momenta_decomposition_to_string(cls, momenta_decomposition):
        """ Turns ((1,0,0,-1),(1,1)) into 'k1-k4+p1+p2'"""

        res = ""
        first=True
        # The outgoing momenta are set element-wise equal to the incoming ones.
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
            del edge_data['name']

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

        form_super_graph =  cls(
            name = LTD2_super_graph.name,
            call_identifier=LTD2_super_graph.call_signature,
            edges = dict(local_graph.edges),
            nodes = dict(local_graph.nodes),
            overall_factor = overall_factor,
        )

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

    def generate_numerator_form_input(self):
        return '+'.join(g.generate_numerator_form_input() for g in self)

    def generate_numerator_functions(self, output_format='c'):
        """ Use form to plugin Feynman Rules and process the numerator algebra so as
        to generate a low-level routine in file_path that encodes the numerator of this supergraph."""

        # write the form input to a file
        form_input = self.generate_numerator_form_input()

        with open(pjoin(FORM_workspace,'input.h'), 'w') as f:
            f.write('L F = {};'.format(form_input))

        # TODO specify cores in the input
        r = subprocess.run([
                FORM_processing_options["FORM_path"], 
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

class FORMSuperGraphList(list):
    """ Container class for a list of FORMSuperGraphIsomorphicList."""

    extension_names = {'c': 'c'}

    def __init__(self, graph_list):
        """ Instantiates a list of FORMSuperGraphs from a list of either
        FORMSuperGraphIsomorphicList instances, FORMSuperGraph, or LTD2SuperGraph instances."""

        for g in graph_list:
            if isinstance(g, FORMSuperGraph):
                self.append(FORMSuperGraphIsomorphicList([g]))
            elif isinstance(g, FORMSuperGraphIsomorphicList):
                self.append(g)
            else:
                self.append(FORMSuperGraphIsomorphicList([FORMSuperGraph.from_LTD2SuperGraph(g)]))

    @classmethod
    def from_dict(cls, dict_file_path):
        """ Creates a FORMSuperGraph list from a dict file path."""
        from pathlib import Path
        p = Path(dict_file_path)
        sys.path.insert(0, str(p.parent))
        m = __import__(p.stem)

        logger.info("Imported {} supergraphs.".format(len(m.graphs)))

        graph_list = []
        for i, g in enumerate(m.graphs):
            # convert to FORM supergraph
            form_graph = FORMSuperGraph(name='Graph' + str(i), edges = g['edges'], nodes=g['nodes'], overall_factor=g['overall_factor'])
            graph_list.append(form_graph)

        # TODO: put isomorphic graphs in the same group
        iso_groups = [[g] for g in graph_list]

        return FORMSuperGraphList([FORMSuperGraphIsomorphicList(iso_group) for iso_group in iso_groups])

    def to_dict(self, file_path):
        """ Outputs the FORMSuperGraph list to a Python dict file path."""

        # TODO: Dump to Python dict
        pass

    def generate_numerator_functions(self, root_output_path, output_format='c'):
        """ Generates optimised source code for the graph numerator in several
        files rooted in the specified root_output_path."""

        if output_format not in self.extension_names:
            raise FormProcessingError("Unsupported output format for numerator functions: '%s'"%output_format)

        # add all numerators in one file and write the headers
        numerator_code = """#include <math.h>
            // TODO: define all masses!

            double evaluate(double k[], int i){
                switch(i) {
        """ + '\n'.join('\t\tcase {}: return evaluate_{}(k);'.format(i, i) for i in range(len(self))) + \
        """            
                }
            }
        """

        pattern = re.compile(r'Z(()\d*)_')

        for i, graph in enumerate(self):
            num = graph.generate_numerator_functions()

            max_intermediate_variable = max(int(index.groups()[0]) for index in pattern.finditer(num))

            numerator_code += 'double evaluate_{}(double k[]){{\n\tdouble {};\n'.format(i,
                ','.join('Z' + str(i) for i in range(1,max_intermediate_variable + 1))
            ) + num + '\n}'

        writers.CPPWriter(pjoin(root_output_path, 'numerator.c')).write(numerator_code)
        if os.path.isfile(pjoin(root_output_path,'Makefile')):
            try:
                misc.compile(cwd=root_output_path,mode='cpp')
            except MadGraph5Error as e:
                logger.info("%sCompilation of FORM-generated numerator failed:\n%s%s"%(
                    utils.bcolors.RED,str(e),utils.bcolors.ENDC))
        else:
            logger.warning(("\n%sYou are running FORM_processing directly from the __main__ of FORM_processing.py.\n"+
                           "You will thus need to compile numerators.c manually.%s")%(utils.bcolors.GREEN, utils.bcolors.ENDC))

class FORMProcessor(object):
    """ A class for taking care of the processing of a list of FORMSuperGraphList.
    Useful because many aspects common to all supergraphs and function do not belong to FORMSuperGraphList.
    """

    def __init__(self, super_graphs_list, model, process_definition):
        """ Specify aditional information such as the model that is useful for FORM processing."""
        self.super_graphs_list = super_graphs_list
        self.model = model
        self.process_definition = process_definition

        ############
        # FOR BEN: 
        # ##########

        # example of how to access various parameters
        PDG_code = 6
        misc.sprint('m_top_quark=',self.model['parameter_dict'][self.model.get_particle(PDG_code).get('mass')].real)
        PDG_code = 1
        misc.sprint('m_down_quark=',self.model['parameter_dict'][self.model.get_particle(PDG_code).get('mass')].real)
        PDG_code = 6
        misc.sprint('width_top_quark=',self.model['parameter_dict'][self.model.get_particle(PDG_code).get('width')].real)
        # And some various other useful coupling parameters:
        misc.sprint('gs=',self.model['parameter_dict']['G'].real)
        misc.sprint('alpha_EW=',1./self.model['parameter_dict']['aEWM1'].real)
        misc.sprint('y_t=',self.model['parameter_dict']['mdl_yt'].real)
        # Display all parameters as follows:
        misc.sprint(self.model['parameter_dict'].keys())

    def generate_numerator_functions(self, root_output_path, output_format='c'):

        return self.super_graphs_list.generate_numerator_functions(
            root_output_path, output_format=output_format)

if __name__ == "__main__":
    logging.info("TODO: process arguments and load an externally provided list of yaml dump files.")

    if len(sys.argv)==1 or not os.path.isfile(sys.argv[1]):
        raise FormProcessingError("Incorrect arguments: %s"%str(sys.argv))
    super_graph_path = sys.argv[1]

    import alpha_loop.interface as interface
    cli = interface.alphaLoopInterface()

    cli.do_import('model PLUGIN/alphaloop/models/aL_sm-no_widths')
    computed_model = model_reader.ModelReader(cli._curr_model)
    computed_model.set_parameters_and_couplings(pjoin(
        plugin_path,os.path.pardir,'models','aL_sm','restrict_no_widths.dat'))        
    process_definition=cli.extract_process('e+ e- > a > d d~', proc_number=0)

    super_graph_list = FORMSuperGraphList.from_dict(super_graph_path)
    form_processor = FORMProcessor(super_graph_list, computed_model, process_definition)
    form_processor.generate_numerator_functions('.', output_format='c')
