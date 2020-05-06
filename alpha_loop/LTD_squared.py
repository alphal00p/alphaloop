#####################################################
#                                                   #
#  Source file of the alphaLoop MG5aMC plugin.      #
#                                                   #
#####################################################

import networkx as nx
import progressbar
import copy
import logging
import madgraph.various.misc as misc

import LTD.squared_topologies as squared_topology_processor
import alpha_loop.utils as utils

logger = logging.getLogger('alphaLoop.LTD2Processing')

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

        graph = nx.MultiDiGraph()

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
                        self.initial_state_edges[leg.get('number')] = (u_edge, v_edge, 0)
                    else:
                        self.final_state_edges[leg.get('number')] = (u_edge, v_edge, 0)
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
                0
            )
        for leg_number in complex_conjugated_LTD2_diagram.final_state_edges:
            complex_conjugated_LTD2_diagram.final_state_edges[leg_number] = (
                complex_conjugated_LTD2_diagram.final_state_edges[leg_number][1].replace('O','I'),
                complex_conjugated_LTD2_diagram.final_state_edges[leg_number][0].replace('L','R'),
                0
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

    def __init__(self, LTD2_diagram_list, proc_number):
        """ Instantiate a list of all super graphs from a list of LTD2 diagrams."""

        self.proc_number = proc_number

        # First build all possible super graphs 
        # TODO there are certainly optimisations to consider here. This step will be
        # very slow for large number of diagrams (100+)
        all_super_graphs = []
        logger.info("Building supergraphs from %d x %d diagramatic combinations..."%(
                                                len(LTD2_diagram_list),len(LTD2_diagram_list)))
        n_super_graphs_ignored = 0
        with progressbar.ProgressBar(
            prefix = 'Generating non-unique supergraphs : ',
            max_value=int((len(LTD2_diagram_list)*(len(LTD2_diagram_list)+1))/2)) as bar:
            i_bar = 0
            for i_diag_left, LTD2_diag_left in enumerate(LTD2_diagram_list):
                for i_diag_right, LTD2_diag_right in enumerate(LTD2_diagram_list[i_diag_left:]):
                    new_super_graph = SuperGraph(LTD2_diag_left,LTD2_diag_right, 
                        {   'proc_id' : proc_number, 
                            'left_diagram_id' : i_diag_left+1, 
                            'right_diagram_id' : i_diag_left+i_diag_right+1,
                        }
                    )
                    if not new_super_graph.should_be_considered():
                        n_super_graphs_ignored += 1
                        continue
                    all_super_graphs.append(new_super_graph)
                    i_bar += 1
                    bar.update(i_bar)

        if n_super_graphs_ignored>0:
            logger.critical("%salphaLoop removed a total of %d supergraphs because their kind is not supported yet. Results will be wrong!%s"%
                                (utils.bcolors.RED,n_super_graphs_ignored,utils.bcolors.ENDC))

        # Then filter isomorphic ones
        logger.info("Detecting isomorphisms between %d non-unique super-graphs..."%len(all_super_graphs))
        with progressbar.ProgressBar(
            prefix = 'Filtering supergraphs ({variables.n_unique_super_graphs_sofar} unique found so far) : ',
            max_value=len(all_super_graphs),
            variables = {'n_unique_super_graphs_sofar' : '0'}
            ) as bar:
            for i_graph, super_graph in enumerate(all_super_graphs):
                if not any(super_graph.is_isomorphic_to(graph_in_basis) for graph_in_basis in self):
                    self.append(super_graph)
                bar.update(n_unique_super_graphs_sofar='%d'%len(self))
                bar.update(i_graph+1)
        
        # print("This process contains %d individual different super-graphs."%len(self))

class SuperGraph(object):
    """ Class representing a super graph in the LTD^2 formalism"""

    def __init__(self, LTD2_diagram_left, LTD2_diagram_right, call_signature, name=None):
        """ Instantiate a super graph from two LTD2Diagram instances sitting respectively
        to the left and right of the Cutkosky cut."""

        # Typically a dictionary witht he following entries:
        # {'proc_id' : <i>, 'left_diagram_id' : <i>, 'right_diagram_id' : <i> }
        self.call_signature = call_signature

        self.diag_left_of_cut = copy.deepcopy(LTD2_diagram_left)
        self.diag_right_of_cut = LTD2_diagram_right.get_complex_conjugate()

        # The attributes below will be set by the function sew_graphs.
        self.graph = None
        self.cuts = None
        self.external_incoming_momenta = None
        self.sew_graphs(self.diag_left_of_cut, self.diag_right_of_cut)

        self.name = name

    def set_name(self, name):
        """ Define an accessor for the name attribute to make it clear that it will typically be set after 
        the generation of the supergraph."""
        self.name = name

    def should_be_considered(self):
        """ Place here all rules regarding whether this supergraph should be 
        considered without our squared LTD framework."""

        return True
        # TODO remove self-energies as we will be implementing them from
        # a separate procedure involving an indpendent process generation
        # featuring two-point vertices.

    def generate_yaml_input_file(self, file_path, model, alphaLoop_options):
        """ Generate the yaml input file for the rust_backend, fully specifying this squared topology."""

        local_DEBUG = False
        if local_DEBUG: misc.sprint("Now processing topology '%s'."%self.name)
        # Nodes are labelled like I<i>, O<i>, L<i> or R<i>.
        # We will map those integers by offsets of the letter
        letter_offsets = {'I':1000,'O':2000,'L':8000, 'R':9000}
        def cast_node_to_int(node_name):
            return letter_offsets[node_name[0]]+int(node_name[1:])

        def get_edges_list():
            edges_list = []
            for edge, edge_info in self.graph.edges.items():
                edges_list.append(
                    ( edge_info['name'], cast_node_to_int(edge[0]),cast_node_to_int(edge[1]) )
                )
            if local_DEBUG: misc.sprint(edges_list)
            return edges_list

        def get_incoming_momenta_names():
            incoming_momenta_names = [ self.graph.edges[edge[1]]['name'] for edge in self.external_incoming_momenta ]
            if local_DEBUG: misc.sprint(incoming_momenta_names)
            return incoming_momenta_names

        def get_external_momenta_assignment():
            # TODO This must go, it does not make much sense for it to be assigned in the context of using MG numerators.
            return {'q1': [1000., 0., 0., 1000.], 'q2': [1000., 0., 0., -1000.], 'q3': [1000., 0., 0., 1000.], 'q4': [1000., 0., 0., -1000.]}

        def get_loop_momenta_names():
            loop_momenta_names = [ self.graph.edges[edge[1]]['name'] for edge in self.cuts[:-1] ]
            if local_DEBUG: misc.sprint(loop_momenta_names)
            return loop_momenta_names

        def get_final_state_particle_ids():
            final_state_particle_ids = []
            all_external_pdgs = [ self.graph.edges[edge[1]]['pdg'] for edge in self.cuts ]
            for pdg in all_external_pdgs:
                if pdg not in alphaLoop_options['_jet_PDGs']:
                    final_state_particle_ids.append(pdg)
            if local_DEBUG: misc.sprint(final_state_particle_ids)
            return tuple(final_state_particle_ids)

        def get_njets_in_observable_process():
            njets_in_observable_process = max( (len(self.cuts)-len(get_final_state_particle_ids()))-alphaLoop_options['perturbative_order'].count('N'), 0)
            if local_DEBUG: misc.sprint(njets_in_observable_process)
            return njets_in_observable_process

        def get_particle_ids():
            particle_ids = { edge_info['name'] : edge_info['pdg'] for 
                edge, edge_info in self.graph.edges.items()
            }
            if local_DEBUG: misc.sprint(particle_ids)
            return particle_ids

        squared_topology = squared_topology_processor.SquaredTopologyGenerator(
            get_edges_list(),
            self.name, 
            get_incoming_momenta_names(), 
            get_njets_in_observable_process(),
            # Below needs to be removed when interfacing to MG5aMC numerators
            get_external_momenta_assignment(),
            loop_momenta_names=get_loop_momenta_names(),
            final_state_particle_ids=get_final_state_particle_ids(),
            particle_ids=get_particle_ids(),
            MG_numerator = {
                # The call signature is typically a dictionary with the following format
                # {'proc_id' : <i>, 'left_diagram_id' : <i>, 'right_diagram_id' : <i> }
                'call_signature' : self.call_signature
            },
            # The numerator specifications below are of no use in the context of using MG numerators.
            overall_numerator=1.0,
            numerator_structure={}
        )
        squared_topology.export(file_path)

    def draw(self):
        """ For debugging: use matplotlib to draw this super graph. """
        
        nx.draw_networkx(self.graph)
        import matplotlib.pyplot as plt
        plt.show()

    def sew_graphs(self, diag_left_of_cut, diag_right_of_cut):
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
            left_graph.remove_edge(*left_graph_edge)
            right_graph_edge = diag_right_of_cut.final_state_edges[final_state_leg_number]
            right_node = right_graph_edge[1]
            right_graph.remove_edge(*right_graph_edge)
            # Now add back to the whole graph, which we take to be the left graph, the edge
            # directly connecting the two nodes that use to link to the external leg_number final_state_leg_number
            # Use the final_state_number as the key for the multiDiGraph edge so as to support multi-edges.
            # We are guaranteed that the final state leg number can serve as a unique key.
            sewing_edges[final_state_leg_number]= ((left_node, right_node, final_state_leg_number),edge_attributes)

        # Now name edges
        initial_state_identifier=0
        external_incoming_momenta = []
        for initial_state_leg_number in sorted(diag_left_of_cut.initial_state_edges.keys()):
            initial_state_identifier += 1
            left_graph.edges[diag_left_of_cut.initial_state_edges[initial_state_leg_number]]['name'] = 'q%d'%initial_state_identifier
            external_incoming_momenta.append((initial_state_leg_number,diag_left_of_cut.initial_state_edges[initial_state_leg_number]))
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
        cuts = []
        for leg_number in sorted(sewing_edges.keys()):
            (sewing_edge_u, sewing_edge_v, edge_key), edge_attributes = sewing_edges[leg_number]
            # name the cut edge as `cut<i>` where <i> is the leg number being cut.
            edge_attributes['name']='cut%d'%leg_number
            cuts.append((leg_number, (sewing_edge_u, sewing_edge_v, edge_key)))
            sewed_graph.add_edge(sewing_edge_u, sewing_edge_v, key=edge_key, **edge_attributes)
        
        self.graph = sewed_graph
        self.cuts = tuple(cuts)
        self.external_incoming_momenta = tuple(external_incoming_momenta)

        return sewed_graph, 
    
    def is_isomorphic_to(self, other_super_graph):
        """ Uses networkx to decide if the two graphs are isomorphic."""

        def edge_match_function(e1,e2):
            # This function needs tu support multi-edges
            # For now all we do is making sure the set of PDGs of 
            # all edges connecting the two nodes match.
            return set(e['pdg'] for e in e1.values()) == \
                   set(e['pdg'] for e in e2.values())
        return nx.is_isomorphic(self.graph, other_super_graph.graph,
            edge_match=edge_match_function,
            node_match=lambda n1,n2: n1['vertex_id']==n2['vertex_id']
        )