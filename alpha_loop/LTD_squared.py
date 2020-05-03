#####################################################
#                                                   #
#  Source file of the alphaLoop MG5aMC plugin.      #
#                                                   #
#####################################################

import networkx as nx

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

class LTD2Diagram(list):
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
        self.graph = self.build_graph()
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
        next_node_number = 0
        for leg in initial_legs:
            leg_number_to_node_number[leg.get('number')] = leg.get('number')
            next_node_number=max(next_node_number,leg.get('number'))
            graph.add_node(leg.get('number'), attr_dict={
                'vertex_id' : -10 # stands for incoming leg
                }
            )
        for leg in final_legs:
            leg_number_to_node_number[leg.get('number')] = leg.get('number')
            next_node_number=max(next_node_number,leg.get('number'))
            graph.add_node(leg.get('number'), attr_dict={
                'vertex_id' : -20 # stands for outgoing leg
                }
            )

        for vertex in self.base_diagram.get('vertices'):
            next_node_number += 1
            graph.add_node(next_node_number, attr_dict={
                'vertex_id' : vertex.get('id'),
                }
            )
            leg_numbers_processed=[]
            for leg in vertex.get('legs'):
                if leg.get('number') in leg_numbers_processed:
                    continue
                leg_numbers_processed.append(leg.get('number'))
                if leg.get('number') in final_leg_numbers:
                    u_edge, v_edge = next_node_number, leg_number_to_node_number[leg.get('number')]
                else:
                    u_edge, v_edge = leg_number_to_node_number[leg.get('number')], next_node_number
                graph.add_edge( u_edge, v_edge , attr_dict={
                    'pdg' : leg.get('id'),
                    }
                )
                leg_number_to_node_number[leg.get('number')] = next_node_number

        #print(self.base_diagram)
        #print(self.helas_amplitude_numbers)
        return graph