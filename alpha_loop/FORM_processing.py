import copy
import logging
import os
from pprint import pprint, pformat
try:
    import madgraph.various.misc as misc
except:
    pass
from itertools import chain
import sys

logger = logging.getLogger('alphaLoop.FORM_processing')

pjoin = os.path.join

class FormProcessingError(Exception):
    """ Error for the FORM processing phase."""
    pass

class FORMSuperGraph(object):
    """ Simplified SuperGraph object with only the necessary information to have it processed by FORM."""

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
        for node in self.nodes.values():
            index_offset = min(index_offset, min(node['indices']))
        index_offset = -index_offset + 1

        # create the input file for FORM
        form_diag = self.overall_factor
        for node_id, node in self.nodes.items():
            if node_id < 0:
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
    
    def generate_numerator_functions(self, file_path, output_format='c'):
        """ Use form to plugin Feynman Rules and process the numerator algebra so as
        to generate a low-level routine in file_path that encodes the numerator of this supergraph."""

        # TODO Ben
        #open(file_path,'w').write("TODO")

    @classmethod
    def from_LTD2SuperGraph(cls, LTD2_super_graph):
        """ Creates a FORMSuperGraph from an LTD2_super_graph."""

        # TODO: correctly assign the overall factor.
        overall_factor = "1"

        model = LTD2_super_graph.model
        all_edges=dict(copy.deepcopy(LTD2_super_graph.graph.edges))
        misc.sprint(all_edges)
        all_nodes=dict(copy.deepcopy(LTD2_super_graph.graph.nodes))
        misc.sprint(all_nodes)
        # Example printout information about the vertex
        for node_key, node_info in all_nodes.items():
            if node_info['vertex_id']>0:
                misc.sprint("Vertex of node %s:\n%s"%(
                    str(node_key),pformat(model.get_interaction(node_info['vertex_id']))))
        raise FormProcessingError("The rest of the conversion of an LTD2SuperGraph to "+
            "a FORMSuperGraph still remains to be done.")

        return cls(
            name = LTD2_super_graph.name,
            call_identifier=LTD2_super_graph.call_signature,
            edges = all_edges,
            nodes = all_nodes,
            overall_factor = overall_factor,
        )

    @classmethod
    def to_yaml(self, yaml_file_path):
        """ Outputs the FORMSuperGraph self to a yaml file path."""

        # TODO: Dump to yaml
        pass


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
        form_input = ''
        for g in self:
            form_input += g.generate_numerator_form_input()

        return form_input

    def generate_numerator_functions(self, file_path, output_format='c'):
        """ Use form to plugin Feynman Rules and process the numerator algebra so as
        to generate a low-level routine in file_path that encodes the numerator of this supergraph."""


        form_input = self.generate_numerator_form_input()

        print(form_input)

        # TODO Ben
        #open(file_path,'w').write("TODO")

        
class FORMSuperGraphList(list):
    """ Container class for a list of FORMSuperGraphIsomorphicList."""

    extension_names = {'rust': 'rs', 'c': 'c'}

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

        print("Imported {} graphs".format(len(m.graphs)))

        graph_list = []
        for i, g in enumerate(m.graphs):
            # convert to FORM supergraph
            form_graph = FORMSuperGraph(name='Graph' + str(i), edges = g['edges'], nodes=g['nodes'], overall_factor=g['overall_factor'])
            graph_list.append(form_graph)

        # TODO: put isomorphic graphs in the same group
        iso_groups = [[g] for g in graph_list]

        return FORMSuperGraphList([FORMSuperGraphIsomorphicList(iso_group) for iso_group in iso_groups])

    @classmethod
    def to_yaml(self, yaml_file_path):
        """ Outputs the FORMSuperGraph list to a yaml file path."""

        # TODO: Dump to yaml
        pass

    def generate_numerator_functions(self, root_output_path, output_format='rust'):
        """ Generates optimised source code for the graph numerator in several
        files rooted in the specified root_output_path."""

        if output_format not in self.extension_names:
            raise FormProcessingError("Unsupported output format for numerator functions: '%s'"%output_format)

        for i, graph in enumerate(self):
            graph.generate_numerator_functions(
                pjoin(root_output_path,'%s.%s'%('G' + str(i), self.extension_names[output_format]))
            )

if __name__ == "__main__":
    logging.info("TODO: process arguments and load an externally provided list of yaml dump files.")

    super_graph_list = FORMSuperGraphList.from_dict(sys.argv[1])
    super_graph_list.generate_numerator_functions('.', output_format='c')
