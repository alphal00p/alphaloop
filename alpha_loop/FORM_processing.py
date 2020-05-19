import copy
import logging
import os
from pprint import pprint, pformat
try:
    import madgraph.various.misc as misc
except:
    pass
from itertools import chain

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
    
    def generate_numerator_functions(self, file_path, output_format='rust'):
        """ Use form to plugin Feynman Rules and process the numerator algebra so as
        to generate a low-level routine in file_path that encodes the numerator of this supergraph."""

        # TODO Ben
        open(file_path,'w').write("TODO")

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
    def from_yaml(cls, yaml_file_path):
        """ Creates a FORMSuperGraph from a yaml file path."""

        # TODO: Creates an instance from a yaml dump.
        pass

    @classmethod
    def to_yaml(self, yaml_file_path):
        """ Outputs the FORMSuperGraph self to a yaml file path."""

        # TODO: Dump to yaml
        pass
        
class FORMSuperGraphList(list):
    """ Container class for a list of FROMSuperGraph."""

    extension_names = {'rust','rs'}

    def __init__(self, graph_list):
        """ Instantiates a list of FORMSuperGraphs from a list of either
        FORMSuperGraph instances or LTD2SuperGraph instances."""

        for g in graph_list:
            if isinstance(g,FORMSuperGraph):
                self.append(g)
            else:
                self.append(FORMSuperGraph.from_LTD2SuperGraph(g))

    @classmethod
    def from_yaml(cls, yaml_file_path):
        """ Creates a FORMSuperGraph list from a yaml file path."""

        # TODO: Creates an instance from a yaml dump.
        pass

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

        for graph in self:
            graph.generate_numerator_functions(
                pjoin(root_output_path,'%s.%s'%(graph.name,self.extension_names[output_format]))
            )

if __name__ == "__main__":
    logging.info("TODO: process arguments and load an externally provided list of yaml dump files.")