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
from itertools import chain
from madgraph import MadGraph5Error, InvalidCmd, MG5DIR

import LTD.squared_topologies as squared_topology_processor
import alpha_loop.utils as utils

# For the self-energy treatment we have to generate additional species of
# particles and interactions which we will differentiate using this global offset.
self_energy_global_id_offset = 100000

logger = logging.getLogger('alphaLoop.LTD2Processing')

class LTD2Error(MadGraph5Error):
    """ Error for the alphaLoop LTD2 treatment."""
    pass

class LTD2Diagram(object):
    """ Class for storing a diagram, analoguous to base_amplitude.Diagram but with more topological information
    necessary for treating the corresponding contribution in the context of an LTD^2 computation."""

    def __init__(self, helas_diagram, wf_dict, vx_list, process, model, optimization, alphaLoop_options):
        """ Initialise this LTD2Diagram from an instance of HelasDiagrams."""

        self.helas_diagram = helas_diagram
        self.alphaLoop_options = alphaLoop_options
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

    def get_descendants(self, seed_node, veto_edges=None, veto_edge_func=None):
        """ Get all nodes reachable from seed_node, excluding edges passing function veto_edge_func."""

        if veto_edge_func is None:
            veto_edge_func = lambda e: False
        if veto_edges is None:
            veto_edges = []
        
        descendants=set([seed_node,])
        for (u,v,c,edge_info) in chain(
                self.graph.in_edges(seed_node,data=True,keys=True),
                self.graph.out_edges(seed_node,data=True,keys=True)
            ):
            edge_key=(u,v,c)
            if edge_key in veto_edges or veto_edge_func(edge_info):
                continue
            veto_edges.append(edge_key)
            other_node = edge_key[1] if edge_key[0]==seed_node else edge_key[0]
            descendants.add(other_node)
            if not (other_node.startswith('I') or other_node.startswith('O')):
                descendants |= self.get_descendants(other_node,veto_edges,veto_edge_func)
        
        return descendants

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

    def get_copy(self):
        """ Returns a copy of self with only relevant attributed deep-copied."""
        self_copy = copy.copy(self)
        self_copy.initial_state_edges = copy.deepcopy(self.initial_state_edges)
        self_copy.final_state_edges = copy.deepcopy(self.final_state_edges)
        self_copy.graph = copy.deepcopy(self.graph)
        self_copy.n_internal_nodes = copy.deepcopy(self.n_internal_nodes)

        return self_copy

    def get_complex_conjugate(self):
        """ Returns a deep copy of self corresponding to flipping the orientation of all edges, and
        potentially eventually computing the additional factors of `i`."""

        complex_conjugated_LTD2_diagram = self.get_copy()
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

class SelfEnergyLTD2Diagram(LTD2Diagram):
    """ Class for storing a diagram, analoguous to base_amplitude.Diagram but with more topological information
    necessary for treating the corresponding self-energy contribution in the context of an LTD^2 computation."""

    def __init__(self, *args, **opts):
        """ Initialise this self-energy LTD2Diagram from an instance of HelasDiagrams."""

        super(SelfEnergyLTD2Diagram, self).__init__(*args, **opts)

class LTD2DiagramList(list):
    """ Class for storing a list of LTD2Diagrams."""

    LTD2_diagram_class = LTD2Diagram

    def __init__(self, helas_matrix_element, alphaLoop_options):
        """ Initialise this list of LTD2Diagrams from an instance of HelasMatrixElement."""

        self.helas_matrix_element = helas_matrix_element
        self.alphaLoop_options = alphaLoop_options
        self.optimization = 1
        if len([wf for wf in self.helas_matrix_element.get_all_wavefunctions() if wf.get('number') == 1]) > 1:
            self.optimization = 0
        self.process = helas_matrix_element.get('processes')[0]
        self.model = self.process.get('model')

        self.wf_dict = {}
        self.vx_list = []

        for diag in self.helas_matrix_element.get('diagrams'):
            self.append(self.LTD2_diagram_class(diag, self.wf_dict, self.vx_list, 
                self.process, self.model, self.optimization, self.alphaLoop_options))

class SelfEnergyLTD2DiagramList(LTD2DiagramList):
    """ Class for handling/generating a list of SelfEnergyLTD2Diagrams."""

    LTD2_diagram_class = SelfEnergyLTD2Diagram

    def __init__(self, *args, **opts):
        super(SelfEnergyLTD2DiagramList,self).__init__(*args, **opts)

class SuperGraph(object):
    """ Class representing a super graph in the LTD^2 formalism"""

    def __init__(self, LTD2_diagram_left, LTD2_diagram_right, call_signature, name=None):
        """ Instantiate a super graph from two LTD2Diagram instances sitting respectively
        to the left and right of the Cutkosky cut."""

        # Typically a dictionary witht he following entries:
        # {'proc_id' : <i>, 'left_diagram_id' : <i>, 'right_diagram_id' : <i> }
        self.call_signature = call_signature

        self.model = LTD2_diagram_left.model
        self.process = LTD2_diagram_left.process
        self.alphaLoop_options = LTD2_diagram_left.alphaLoop_options
        self.diag_left_of_cut = LTD2_diagram_left.get_copy()
        self.diag_right_of_cut = LTD2_diagram_right.get_copy()

        # The attributes below will be set by the function sew_graphs.
        self.graph = None
        self.cuts = None
        self.external_incoming_momenta = None
        self.sew_graphs(
            # Use copies to avoid border effects from edges removal
            self.diag_left_of_cut.get_copy(), 
            self.diag_right_of_cut.get_complex_conjugate()
        )

        self.self_energy_structure = None

        self.name = name

    def get_self_energy_structures(self):
        """ Computes the self-energy structure.
        An empty dictionary denotes the lack of any self-energy
        """

        if self.self_energy_structure is not None:
            return self.self_energy_structure
        
        all_external_pdgs = {l.get('id') for l in self.process.get('legs')}

        descendants = {}
        for graph in [self.diag_left_of_cut.graph,self.diag_right_of_cut.graph]:
            for edge, edge_info in graph.edges.items():
                # Skip external edges
                if (edge[0][:1] in ['I','O']) or (edge[1][:1] in ['I','O']):
                    continue
                # Skip any edge corresponding to a particle that is not a jet
                # and not in the external states of the process
                if not (edge_info['pdg'] in self.alphaLoop_options['_jet_PDGs'] or 
                    edge_info['pdg'] in all_external_pdgs):
                    continue
                this_edge_descendants = frozenset(d for d in nx.descendants(graph,edge[1]) if d[:1] in ['I','O'])
                if this_edge_descendants in descendants:
                    descendants[this_edge_descendants].append(edge)
                else:
                    descendants[this_edge_descendants] = [ edge, ]

        # Now keep all group descendants which are associated to more than one edge.
        # These correspond to self-energy raised propagators.
        descendants = {k:v for k,v in descendants.items() if len(v)>1}

        self.self_energy_structure = descendants
        if any(len(v)>2 for v in self.self_energy_structure.values()):
            raise LTD2Error("The self-energy detection logic must be wrong.")

        return self.self_energy_structure

    def set_name(self, name):
        """ Define an accessor for the name attribute to make it clear that it will typically be set after 
        the generation of the supergraph."""
        self.name = name

    def should_be_considered_before_isomorphism(self):
        """ Place here all rules regarding whether this supergraph should be 
        considered without our squared LTD framework, but only if this filtering must take place
        before checking for isomorphisms between supergraphs."""

        # By default, all selection will be performed *after* checking for isomorphisms.
        return True

    def should_be_considered(self):
        """ Place here all rules regarding whether this supergraph should be 
        considered without our squared LTD framework."""

        if not self.alphaLoop_options['include_self_energies_from_squared_amplitudes']:
            if len(self.get_self_energy_structures()) > 0:
                return False

        return True

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
        
        # Remove all nodes associated to no edges (resulting from the sewing operation)
        sewed_graph.remove_nodes_from(list(nx.isolates(sewed_graph)))

        self.graph = sewed_graph
        self.cuts = tuple(cuts)
        self.external_incoming_momenta = tuple(external_incoming_momenta)
    
    def is_isomorphic_to(self, other_super_graph):
        """ Uses networkx to decide if the two graphs are isomorphic."""

        def edge_match_function(e1,e2):
            # This function needs tu support multi-edges
            # For now all we do is making sure the set of PDGs of 
            # all edges connecting the two nodes match.
            # TODO: Fix ambiguity with fermion flow and part / antipart
            # For now consider two edges equal whenever the *abs* of PDGs matches.
            return set(abs(e['pdg']) for e in e1.values()) == \
                   set(abs(e['pdg']) for e in e2.values())
        return nx.is_isomorphic(self.graph, other_super_graph.graph,
            edge_match=edge_match_function,
            node_match=lambda n1,n2: n1['vertex_id']==n2['vertex_id']
        )

    def get_subgraphs_info(self):
        """ Return information about subgraphs. For non self-energy supergraphs, this is actually trivial."""
        
        return [{
            'id' : 0,
            'left_edges': [ e_info['name'] for e_key, e_info in self.graph.edges.items() 
                            if e_key[0] in ['I1','I2'] or e_key[1] in ['I1','I2'] ],
            'right_edges' : [ e_info['name'] for e_key, e_info in self.graph.edges.items() 
                            if e_key[0] in ['O1','O2'] or e_key[1] in ['O1','O2'] ],
            'cuts' : [ self.graph.edges[edge[1]]['name'] for edge in self.cuts[:-1] ],
            'momentum_sink' : self.cuts[-1][0]
        },]

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
            
            all_external_pdgs = [ l.get('id') for l in self.process.get('legs') if l.get('state')==True ]
            for pdg in all_external_pdgs:
                if pdg not in alphaLoop_options['_jet_PDGs']:
                    final_state_particle_ids.append(pdg)
            if local_DEBUG: misc.sprint(final_state_particle_ids)
            return tuple(final_state_particle_ids)

        def get_njets_in_observable_process():
            #TODO do an actual correct filtering based on the coupling orderrs
            njets_in_observable_process = max( (len(self.cuts)-len(get_final_state_particle_ids())) - 
                                    sum(alphaLoop_options['perturbative_orders'].values())//2, 0)
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
                'call_signature' : self.call_signature,
            },
            subgraphs_info = self.get_subgraphs_info(),
            # The numerator specifications below are of no use in the
            # w context of using MG numerators.
            overall_numerator=1.0,
            numerator_structure={}
        )
        squared_topology.export(file_path)

class SelfEnergySuperGraph(SuperGraph):
    """ Class representing a super graph in the LTD^2 formalism corresponding to a self-energy"""

    def __init__(self, LTD2_diagram_left, LTD2_diagram_right, call_signature, **opts):
        """ Instantiate a self-nergy super graph from two LTD2Diagram instances sitting respectively
        to the left and right of the Cutkosky cut."""

        super(SelfEnergySuperGraph, self).__init__(LTD2_diagram_left, LTD2_diagram_right, call_signature, **opts)

    def should_be_considered_before_isomorphism(self):
        """ Place here all rules regarding whether this supergraph should be 
        considered without our squared LTD framework, but only if this filtering must take place
        before checking for isomorphisms between supergraphs."""

        # By default, all selection will be performed *after* checking for isomorphisms.
        is_connected = nx.is_weakly_connected(self.graph)

        return is_connected

    def should_be_considered(self):
        """ Place here all rules regarding whether this supergraph should be 
        considered without our squared LTD framework."""

        return True

    def sew_graphs(self, diag_left_of_cut, diag_right_of_cut):
        """ Combine two diagrams into one graph that corresponds to the super graph."""

        # First sew graph as normally done for any super graph built from two LTD diagram
        super(SelfEnergySuperGraph,self).sew_graphs(diag_left_of_cut, diag_right_of_cut)

        # First unpack the information about the nature of each edge stored in its PDG
        subgraph_ids = set([])
        for edge_key, edge_info in self.graph.edges.items():
            # First set which self-energy this edge belongs to.
            # 0 indicates no self-energy
            edge_info['self_energy_number'] = (abs(edge_info['pdg'])//self_energy_global_id_offset)%10
#            misc.sprint(edge_info['pdg'])
#            misc.sprint(edge_info['pdg']//self_energy_global_id_offset)
#            misc.sprint((edge_info['pdg']//self_energy_global_id_offset)%10)
            subgraph_ids.add(edge_info['self_energy_number'])
            # Then whether this edge is a bridge, 0 means this is not a bridge
            edge_info['bridge_number'] = (abs(edge_info['pdg'])//(self_energy_global_id_offset*100))%10
            # Whether this is an anchor edge
            edge_info['anchor_number'] = (abs(edge_info['pdg'])//(self_energy_global_id_offset*10))%10
            # Finally we can revert back the PDG to the physical one.
            edge_info['pdg'] =  edge_info['pdg']%self_energy_global_id_offset

        # Same for the vertex id of nodes
        for node_key, node_info in self.graph.nodes.items():

            vertex_id = abs(node_info['vertex_id'])
            # The self-energy this node belongs to
            node_info['self_energy_number'] = (vertex_id//self_energy_global_id_offset)%10
            # Type of interaction
            if node_info['self_energy_number']==0:
                node_info['interaction_type'] = 'base'
            else:
                if ((vertex_id//(self_energy_global_id_offset*10))%10)!=0:
                    node_info['interaction_type'] = 'selfEnergy_anchor'
                elif ((vertex_id//(self_energy_global_id_offset*1000))%10)!=0:
                    if vertex_id%self_energy_global_id_offset==0:
                        node_info['interaction_type'] = 'bridge_anchor'
                    else:
                        if ((vertex_id//(self_energy_global_id_offset*1000))%10)==1:
                            node_info['interaction_type'] = 'bridge_base'
                        else:
                            node_info['interaction_type'] = 'selfEnergy_bridge'
                elif ((vertex_id//(self_energy_global_id_offset*100))%10)!=0:
                    node_info['interaction_type'] = 'selfEnergy'
                else:
                    raise LTD2Error("Incorrect interaction ID in self-energy reconstruction.")
            # We can now revert the id to its original one of the base model
            node_info['vertex_id'] =  node_info['vertex_id']%self_energy_global_id_offset

        subgraph_ids = sorted(list(subgraph_ids))

        # Then identify all the edges of the self-energy subgraphs and all relevant information
        # for its use later in Rust and remove all bridges
        n_self_energies = len(subgraph_ids)-1
        self.subgraphs = [
            {
                'id' : subgraph_id, # 0 for the base graph and 1 for each self-energy
                'cuts' : [],
                'left_edge' : None,
                'right_edge' : None,
                'momentum_sink' : None,
                'nodes' : [],
                'edges' : [],
            } 
            for subgraph_id in subgraph_ids
        ]

        # First assign the infrmation each self-energy:
        for leg_number, cut_edge_key in self.cuts:
            cut_edge = self.graph.edges[cut_edge_key]
            cut_edge['leg_number'] = leg_number
            if cut_edge['anchor_number'] == 0:
                self.subgraphs[cut_edge['self_energy_number']]['cuts'].append((leg_number, cut_edge_key))
            else:
                if self.subgraphs[cut_edge['self_energy_number']]['momentum_sink'] is None:
                    self.subgraphs[cut_edge['self_energy_number']]['momentum_sink'] = leg_number
                else:
                    raise LTD2Error("Each self-energy should only have one momentum sink.")

        # Also set all the bunches that merge into a bridge, and whose momenta must therefore sum to zero.
        # Keys are subgraph_IDs
        nodes_defining_bunches = {}

        for node_key, node_info in self.graph.nodes.items():

            if 'selfEnergy' in node_info['interaction_type']:
                self.subgraphs[edge_info['self_energy_number']]['nodes'].append(node_key)

            if node_info['interaction_type']=='selfEnergy_anchor':
                for (u,v,c,edge_info) in chain(
                        self.graph.in_edges(node_key,keys=True,data=True),
                        self.graph.out_edges(node_key,keys=True,data=True)
                    ):
                    if edge_info['self_energy_number']==0:
                        if node_key.startswith('L'):
                            self.subgraphs[node_info['self_energy_number']]['left_edge'] = (u,v,c)
                            break
                        else:
                            self.subgraphs[node_info['self_energy_number']]['right_edge'] = (u,v,c)
                            break
                else:
                    raise LTD2Error("Could not find entry or exit points of a self-energy.")

            if node_info['interaction_type'] in ['selfEnergy_bridge','bridge_base']:
                if node_info['self_energy_number'] not in nodes_defining_bunches:
                    nodes_defining_bunches[node_info['self_energy_number']] = [node_key]
                else:
                    nodes_defining_bunches[node_info['self_energy_number']].append(node_key)

        def is_a_bridge_or_anchor(edge_info):
            return ( 
                (abs(edge_info['pdg'])//(self_energy_global_id_offset*10))%10!=0 or 
                (abs(edge_info['pdg'])//(self_energy_global_id_offset*100))%10!=0
            )

        # for each bunch node detected, identify all the cut edges that need to combine into zero momentum
        self.bunches = {}
        for subgraph_id in nodes_defining_bunches:
            cut_bunches = []
            for node_key in nodes_defining_bunches[subgraph_id]:
                if node_key.startswith('R'):
                    leg_number_for_this_bunch = [ int(d[1:]) for d in self.diag_right_of_cut.get_descendants(
                            'L'+node_key[1:], 
                            # Forbid moving through a bridge edge
                            veto_edge_func=is_a_bridge_or_anchor
                        ) if d[:1] in ['I','O'] ]
                else:
                    leg_number_for_this_bunch = [ int(d[1:]) for d in self.diag_left_of_cut.get_descendants(
                            node_key, 
                            # Forbid moving through a bridge edge
                            veto_edge_func=is_a_bridge_or_anchor
                        ) if d[:1] in ['I','O'] ]
                if sorted(leg_number_for_this_bunch) not in cut_bunches:
                    cut_bunches.append(sorted(leg_number_for_this_bunch))
            self.bunches[subgraph_id] = cut_bunches

        # Assign the edges of each self-energy
        for edge_key, edge_info in self.graph.edges.items():

            if edge_info['bridge_number']==0 and edge_info['anchor_number']==0:
                self.subgraphs[edge_info['self_energy_number']]['edges'].append(edge_key)

        # We can finally remove all the superfluous edges from this graph
        for edge_key, edge_info in list(self.graph.edges.items()):
            if edge_info['bridge_number']!=0 or edge_info['anchor_number']!=0:
                self.graph.remove_edge(*edge_key)

        # And also remove all remaining isolated nodes
        self.graph.remove_nodes_from(list(nx.isolates(self.graph)))

        # We can now appropriately rename the edges
        for subgraph in self.subgraphs:
            if subgraph['id']==0:
                # Skip basegraph renaming
                continue
            for edge_key in subgraph['edges']:
                self.graph.edges[edge_key]['name'] = 'SE%d_%s'%(subgraph['id'],self.graph.edges[edge_key]['name'])

            self.graph.edges[subgraph['left_edge']]['name'] = 'SEL%d_%s'%(subgraph['id'],self.graph.edges[edge_key]['name'])
            self.graph.edges[subgraph['right_edge']]['name'] = 'SER%d_%s'%(subgraph['id'],self.graph.edges[edge_key]['name'])

            for leg_number, edge_key in subgraph['cuts']:
                self.graph.edges[edge_key]['name'] = 'SEC%d_%s'%(subgraph['id'],self.graph.edges[edge_key]['name'])

        # We build here the list of edges which are actually cut in the resulting self-energy super-graph.
        # These correspond to the definition of the loop momentum basis
        actual_cuts = []
        for subgraph in self.subgraphs[1:]:
            #TODO generalise to support nest_super_graphs, using information from self.bunches
            actual_cuts.extend(subgraph['cuts'][:-1])
        actual_cuts.extend(self.subgraphs[0]['cuts'][:-n_self_energies])
        self.cuts = actual_cuts

    def get_subgraphs_info(self):
        """ Return information about subgraphs. For non self-energy supergraphs, this is actually trivial."""

#        from pprint import pformat
#        misc.sprint(pformat(self.cuts))
#        misc.sprint(pformat(self.subgraphs))
#        misc.sprint(pformat(self.bunches))
        subgraph_info = [{
            'id' : 0,
            'left_edges': [ self.graph.edges[e_key]['name'] for e_key in self.subgraphs[0]['edges'] 
                            if e_key[0] in ['I1','I2'] or e_key[1] in ['I1','I2'] ],
            'right_edges' : [ self.graph.edges[e_key]['name'] for e_key in self.subgraphs[0]['edges'] 
                            if e_key[0] in ['O1','O2'] or e_key[1] in ['O1','O2'] ],
            'cuts' : [ self.graph.edges[edge[1]]['name'] for edge in self.subgraphs[0]['cuts'][:-1] ],
            'momentum_sink' : self.subgraphs[0]['cuts'][-1][0]
        },]
        for subgraph in self.subgraphs[1:]:
            subgraph_info.append({
                'id' : subgraph['id'],
                'left_edges': [ self.graph.edges[subgraph['left_edge']]['name'] ],
                'right_edges' : [ self.graph.edges[subgraph['right_edge']]['name'] ],
                'cuts' : [ self.graph.edges[edge[1]]['name'] for edge in subgraph['cuts'] ],
                'momentum_sink' : subgraph['momentum_sink']
            })

class SuperGraphList(list):
    """ Class for storing a list of SuperGraph instances."""

    super_graph_class = SuperGraph

    def __init__(self, LTD2_diagram_list, proc_number):
        """ Instantiate a list of all super graphs from a list of LTD2 diagrams."""

        self.proc_number = proc_number

        # First build all possible super graphs 
        # TODO there are certainly optimisations to consider here. This step will be
        # very slow for large number of diagrams (100+)
        all_super_graphs = []
        logger.info("Building supergraphs from %d x %d diagramatic combinations..."%(
                                                len(LTD2_diagram_list),len(LTD2_diagram_list)))
        with progressbar.ProgressBar(
            prefix = 'Generating non-unique supergraphs : ',
            max_value=int((len(LTD2_diagram_list)*(len(LTD2_diagram_list)+1))/2)) as bar:
            i_bar = 0
            for i_diag_left, LTD2_diag_left in enumerate(LTD2_diagram_list):
                for i_diag_right, LTD2_diag_right in enumerate(LTD2_diagram_list[i_diag_left:]):
                    new_super_graph = self.super_graph_class(LTD2_diag_left,LTD2_diag_right, 
                        {   'proc_id' : proc_number, 
                            'left_diagram_id' : i_diag_left+1, 
                            'right_diagram_id' : i_diag_left+i_diag_right+1,
                        }
                    )
                    all_super_graphs.append(new_super_graph)
                    i_bar += 1
                    bar.update(i_bar)

        logger.info("Filtering %d super-graphs (before isomorphism check) to remove undesired ones ..."%len(self))
        filtered_list = []
        n_removed_super_graphs = 0
        with progressbar.ProgressBar(
            prefix = 'Removing undesired supergraphs ({variables.n_removed_graphs_so_far} removed so far) : ',
            max_value=len(all_super_graphs),
            variables = {'n_removed_graphs_so_far' : '0'}
            ) as bar:
            for i_graph, super_graph in enumerate(all_super_graphs):
                if super_graph.should_be_considered_before_isomorphism():
                    filtered_list.append(super_graph)
                else:
                    n_removed_super_graphs += 1
                bar.update(n_removed_graphs_so_far='%d'%n_removed_super_graphs)
                bar.update(i_graph+1)

        all_super_graphs = filtered_list
        logger.info("alphaLoop removed a total of %d "%n_removed_super_graphs+
                    " undesired supergraphs *before* isomorphism check.")

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


        logger.info("Filtering %d unique super-graphs to remove undesired ones ..."%len(self))
        filtered_list = []
        n_removed_super_graphs = 0
        with progressbar.ProgressBar(
            prefix = 'Removing undesired supergraphs ({variables.n_removed_graphs_so_far} removed so far) : ',
            max_value=len(self),
            variables = {'n_removed_graphs_so_far' : '0'}
            ) as bar:
            for i_graph, super_graph in enumerate(self):
                if super_graph.should_be_considered():
                    filtered_list.append(super_graph)
                else:
                    n_removed_super_graphs += 1
                bar.update(n_removed_graphs_so_far='%d'%n_removed_super_graphs)
                bar.update(i_graph+1)

        self[:] = filtered_list
        logger.info("alphaLoop removed a total of %d "%n_removed_super_graphs+
                    "supergraphs (likely self-energies that will be generated independently.)")

        # self.draw()
        # print("This process contains %d individual different super-graphs."%len(self))

    def draw(self):
        """ For debugging: use matplotlib to draw this diagram. """

        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(nrows=(len(self)//2)+1, ncols=2,figsize=(20,10))
        ax = axes.flatten()
        for a in ax:
            a.set_axis_off()

        for i_graph, super_graph in enumerate(self):
            nx.draw_networkx(super_graph.graph, ax=ax[i_graph], 
                    label="Proc #%d, Graph #%d"%(self.proc_number,i_graph+1))

        plt.show()

class SelfEnergySuperGraphList(SuperGraphList):
    """ Class for storing and generating a list of SelfEnergySuperGraph instances."""

    super_graph_class = SelfEnergySuperGraph

    def __init__(self, *args, **opts):
        super(SelfEnergySuperGraphList,self).__init__(*args, **opts)