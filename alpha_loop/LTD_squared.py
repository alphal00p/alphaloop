#####################################################
#                                                   #
#  Source file of the alphaLoop MG5aMC plugin.      #
#                                                   #
#####################################################

import networkx as nx
import copy

class LTD2DiagramList(list):
    """ Class for storing a list of LTD2Diagrams."""

    def __init__(self, helas_matrix_element):
        """ Initialise this list of LTD2Diagrams from an instance of HelasMatrixElement."""

        self.helas_matrix_element = helas_matrix_element
        self.optimization = 1
        if len([wf for wf in self.helas_matrix_element.get_all_wavefunctions() if wf.get('number') == 1]) > 1:
            self.optimization = 0
        self.process = helas_matrix_element.get('processes')[0]
        self.model = self.process.get('model')

        self.wf_dict = {}
        self.vx_list = []

        for diag in self.helas_matrix_element.get('diagrams'):
            self.append(LTD2Diagram(diag, self.wf_dict, self.vx_list, self.process, self.model, self.optimization))

class LTD2Diagram(object):
    """ Class for storing a diagram, analoguous to base_amplitude.Diagram but with more topological information
    necessary for treating the corresponding contribution in the context of an LTD^2 computation."""

    def __init__(self, helas_diagram, wf_dict, vx_list, process, model, optimization):
        """ Initialise this LTD2Diagram from an instance of HelasDiagrams."""

        self.helas_diagram = helas_diagram
        self.process = process
        self.model = model
        self.optimization = optimization
        self.wf_dict = wf_dict
        self.vx_list = vx_list

        # The amplidue.Diagram representation of this helas_diagram
        self.base_diagram = self.helas_diagram.get('amplitudes')[0].get_base_diagram(self.wf_dict, self.vx_list, self.optimization)
        self.base_diagram.calculate_orders(self.model)

        self.helas_amplitude_numbers = tuple([amp.get('number') for amp in helas_diagram.get('amplitudes')])

        # Now build the networkX directed graph corresponding to this diagram
        # We also store the (u,v) pair of edges corresponding to the initial and final states.
        self.initial_state_edges = {}
        self.final_state_edges = {}
        self.graph, self.n_internal_nodes = self.build_graph()
        #self.draw()

    def draw(self):
        """ For debugging: use matplotlib to draw this diagram. """
        
        nx.draw_networkx(self.graph)
        import matplotlib.pyplot as plt
        plt.show()

    def build_graph(self):
        """ Build the networkx graph representation of this diagram."""

        graph = nx.DiGraph()

        initial_legs = tuple([leg for leg in self.process.get('legs') if leg.get('state') == False])
        final_legs = tuple([leg for leg in self.process.get('legs') if leg.get('state') == True])
        final_leg_numbers = set([l.get('number') for l in final_legs])

        # Add a node for all external legs
        leg_number_to_node_number = {}
        for leg in initial_legs:
            leg_number_to_node_number[leg.get('number')] = 'I%d'%leg.get('number')
            graph.add_node(leg_number_to_node_number[leg.get('number')], **{
                'vertex_id' : -1 # stands for incoming leg
                }
            )
        for leg in final_legs:
            leg_number_to_node_number[leg.get('number')] = 'O%d'%leg.get('number')
            graph.add_node(leg_number_to_node_number[leg.get('number')], **{
                'vertex_id' : -2 # stands for outgoing leg
                }
            )

        next_node_number = 0
        for vertex in self.base_diagram.get('vertices'):
            next_node_number += 1
            graph.add_node('L%d'%next_node_number, **{
                'vertex_id' : vertex.get('id'),
                }
            )
            leg_numbers_processed=[]
            for leg in vertex.get('legs'):
                if leg.get('number') in leg_numbers_processed:
                    continue
                leg_numbers_processed.append(leg.get('number'))
                if leg.get('number') in final_leg_numbers:
                    u_edge, v_edge = 'L%d'%next_node_number, leg_number_to_node_number[leg.get('number')]
                else:
                    u_edge, v_edge = leg_number_to_node_number[leg.get('number')], 'L%d'%next_node_number
                if leg.get('number') not in self.initial_state_edges and leg.get('number') not in self.final_state_edges:
                    if leg.get('state') == False:
                        self.initial_state_edges[leg.get('number')] = (u_edge, v_edge)
                    else:
                        self.final_state_edges[leg.get('number')] = (u_edge, v_edge)
                graph.add_edge( u_edge, v_edge , **{
                    'pdg' : leg.get('id'),
                    }
                )
                leg_number_to_node_number[leg.get('number')] = 'L%d'%next_node_number

        #print(self.base_diagram)
        #print(self.helas_amplitude_numbers)
        return graph, next_node_number

    def get_complex_conjugate(self):
        """ Returns a deep copy of self corresponding to flipping the orientation of all edges, and
        potentially eventually computing the additional factors of `i`."""

        complex_conjugated_LTD2_diagram = copy.deepcopy(self)
        complex_conjugated_LTD2_diagram.graph = complex_conjugated_LTD2_diagram.graph.reverse(copy=False)

        # Also flip the orientation of the external edges.
        for leg_number in complex_conjugated_LTD2_diagram.initial_state_edges:
            complex_conjugated_LTD2_diagram.initial_state_edges[leg_number] = (
                complex_conjugated_LTD2_diagram.initial_state_edges[leg_number][1].replace('L','R'),
                complex_conjugated_LTD2_diagram.initial_state_edges[leg_number][0].replace('I','O'),
            )
        for leg_number in complex_conjugated_LTD2_diagram.final_state_edges:
            complex_conjugated_LTD2_diagram.final_state_edges[leg_number] = (
                complex_conjugated_LTD2_diagram.final_state_edges[leg_number][1].replace('O','I'),
                complex_conjugated_LTD2_diagram.final_state_edges[leg_number][0].replace('L','R'),
            )
        node_relabeling = {}
        for number in self.initial_state_edges:
            node_relabeling['I%d'%number] = 'O%d'%number
            # Set the vertex ID flag to be an outgoing leg
            complex_conjugated_LTD2_diagram.graph.nodes['I%d'%number]['vertex_id']=-2
        for number in self.final_state_edges:
            node_relabeling['O%d'%number] = 'I%d'%number
            # Set the vertex ID flag to be an incoming leg
            complex_conjugated_LTD2_diagram.graph.nodes['O%d'%number]['vertex_id']=-1
        for node_number in range(1,self.n_internal_nodes+1):
            node_relabeling['L%d'%node_number] = 'R%d'%node_number

        #complex_conjugated_LTD2_diagram.draw()
        complex_conjugated_LTD2_diagram.graph = nx.relabel_nodes(complex_conjugated_LTD2_diagram.graph, node_relabeling, copy=True)

        return complex_conjugated_LTD2_diagram

class SuperGraphList(list):
    """ Class for storing a list of SuperGraph instances."""

    def __init__(self, LTD2_diagram_list):
        """ Instantiate a list of all super graphs from a list of LTD2 diagrams."""

        # First build all possible super graphs 
        # TODO there are certainly optimisations to consider here. This step will be
        # very slow for large number of diagrams (100+)
        all_super_graphs = []
        for i_diag_left, LTD2_diag_left in enumerate(LTD2_diagram_list):
            for LTD2_diag_right in LTD2_diagram_list[i_diag_left:]:
                all_super_graphs.append(SuperGraph(LTD2_diag_left,LTD2_diag_right))

        # Then filter isomorphic ones
        for super_graph in all_super_graphs:
            if not any(super_graph.is_isomorphic_to(graph_in_basis) for graph_in_basis in self):
                self.append(super_graph)
        
        # print("This process contains %d individual different super-graphs."%len(self))

class SuperGraph(object):
    """ Class representing a super graph in the LTD^2 formalism"""

    def __init__(self, LTD2_diagram_left, LTD2_diagram_right):
        """ Instantiate a super graph from two LTD2Diagram instances sitting respectively
        to the left and right of the Cutkosky cut."""

        self.diag_left_of_cut = copy.deepcopy(LTD2_diagram_left)
        self.diag_right_of_cut = LTD2_diagram_right.get_complex_conjugate()

        self.graph, self.cuts = self.sew_graphs(self.diag_left_of_cut, self.diag_right_of_cut)

    def draw(self):
        """ For debugging: use matplotlib to draw this super graph. """
        
        nx.draw_networkx(self.graph)
        import matplotlib.pyplot as plt
        plt.show()

    @classmethod
    def sew_graphs(cls, diag_left_of_cut, diag_right_of_cut):
        """ Combine two diagrams into one graph that corresponds to the super graph."""

        # for each final state leg, remove the corresponding edges and substitute them by one
        # edge directly connecting the two corresponding nodes.
        # Eventually, for ISR we may consider doing the same thing with the initial states.
        assert(diag_left_of_cut.final_state_edges.keys()==diag_right_of_cut.final_state_edges.keys())
        assert(diag_left_of_cut.initial_state_edges.keys()==diag_right_of_cut.initial_state_edges.keys())

        left_graph = diag_left_of_cut.graph
        right_graph = diag_right_of_cut.graph

        sewing_edges = {}
        for final_state_leg_number in diag_left_of_cut.final_state_edges:
            left_graph_edge = diag_left_of_cut.final_state_edges[final_state_leg_number]
            left_node = left_graph_edge[0]
            edge_attributes = dict(left_graph.edges[left_graph_edge])
            left_graph.remove_edge(left_graph_edge[0],left_graph_edge[1])
            right_graph_edge = diag_right_of_cut.final_state_edges[final_state_leg_number]
            right_node = right_graph_edge[1]
            right_graph.remove_edge(right_graph_edge[0],right_graph_edge[1])
            # Now add back to the whole graph, which we take to be the left graph, the edge
            # directly connecting the two nodes that use to link to the external leg_number final_state_leg_number
            sewing_edges[final_state_leg_number]= ((left_node, right_node),edge_attributes)

        # Now name edges
        initial_state_identifier=0
        for initial_state_leg_number in sorted(diag_left_of_cut.initial_state_edges.keys()):
            initial_state_identifier += 1
            left_graph.edges[diag_left_of_cut.initial_state_edges[initial_state_leg_number]]['name'] = 'q%d'%initial_state_identifier
        for initial_state_leg_number in sorted(diag_right_of_cut.initial_state_edges.keys()):
            initial_state_identifier += 1
            right_graph.edges[diag_right_of_cut.initial_state_edges[initial_state_leg_number]]['name'] = 'q%d'%initial_state_identifier

        internal_momenta_id= 0
        for edge in left_graph.edges:
            if 'name' not in left_graph.edges[edge]:
                internal_momenta_id += 1
                left_graph.edges[edge]['name'] = 'Lp%d'%internal_momenta_id
        internal_momenta_id= 0
        for edge in right_graph.edges:
            if 'name' not in right_graph.edges[edge]:
                internal_momenta_id += 1
                right_graph.edges[edge]['name'] = 'Rp%d'%internal_momenta_id

        # Now absorb all remaining edges of the right graph into the left one (the edges relabeling has already been 
        # applied so this is safe)
        sewed_graph = nx.compose(left_graph,right_graph)

        # And finally add back the sewing edges
        for leg_number in sorted(sewing_edges.keys()):
            (sewing_edge_u, sewing_edge_v), edge_attributes = sewing_edges[leg_number]
            # name the cut edge as `cut<i>` where <i> is the leg number being cut.
            edge_attributes['name']='cut%d'%leg_number
            sewed_graph.add_edge(sewing_edge_u, sewing_edge_v, **edge_attributes)
        
        return sewed_graph, tuple([
            (leg_number, sewing_edges[leg_number][0]) 
            for leg_number in sorted(sewing_edges.keys())])
    
    def is_isomorphic_to(self, other_super_graph):
        """ Uses networkx to decide if the two graphs are isomorphic."""
         
        return nx.is_isomorphic(self.graph, other_super_graph.graph,
            edge_match=lambda e1,e2: e1['pdg']==e2['pdg'],
            node_match=lambda n1,n2: n1['vertex_id']==n2['vertex_id']
        )
