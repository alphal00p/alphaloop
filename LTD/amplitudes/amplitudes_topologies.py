import os
import sys
import copy
import math
import mpmath
import numpy as np
from scipy.special import zeta
import networkx as nx
import matplotlib.pyplot as plt

pjoin = os.path.join
file_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, pjoin(file_path, os.pardir))

import vectors
from ltd_utils import LoopTopology, LoopLine, Propagator, TopologyCollection, TopologyGenerator
from analytical_expressions import ladder_phi, analytic_four_point_ladder, analytic_three_point_ladder, analytic_two_point_ladder

zero_lv = vectors.LorentzVector([0., 0., 0., 0.])


#############################################################################################################
# Create the collection of function to build Amplitude Topologies
#############################################################################################################

# A network is generated from the graph corresponding to a amplitude topology
# In blue it shows the notes and in red the loopline id corresponding to that edge
def plot_edge_info(edge_info, file_name=None, show_graph=True):
    # Great graph with oriented edges
    # - Internal nodes are labelled from 0 to 99
    # - External nodes are labelled from 100 to 199
    # - LoopLines nodes are labelled with negative numbers

    plt.figure(figsize=(20, 20))
    G = nx.DiGraph()   # or DiGraph, MultiGraph, MultiDiGraph, etc

    edges = []
    edge_labels = {}
    node_labels = {}

    extra_node=-1
    for k, v in edge_info.items():
        e = v['edge']
        if e[0] >= 100: # External
            G.add_edge(*e, msq = v['m_squared'])
            node_labels.update({e[0]: k, e[1]: str(e[1])})
        elif e[1] >= 100: # External
            G.add_edge(*e, msq = v['m_squared'])
            node_labels.update({e[0]: str(e[0]), e[1]: k})
        else:
            G.add_edge(e[0], extra_node, msq = v['m_squared'])
            G.add_edge(extra_node, e[1], msq = v['m_squared'])
            node_labels.update({e[0]: str(e[0]), e[1]: str(e[1]), extra_node: str(v['loopline'])})
            edge_labels.update({(e[0], extra_node): k })
            edge_labels.update({(extra_node, e[1]): k })
            extra_node -= 1

    massive_edges = []
    massless_edges = []
    for e, msq in nx.get_edge_attributes(G, 'msq').items():
        if msq is not None and msq > 1e-16:
            massive_edges += [e]
        else:
            massless_edges += [e]
    

    pos = nx.kamada_kawai_layout(G)  # positions for all nodes

    # LoopLines nodes
    ll_nodes = [key for key in list(pos.keys()) if key < 0]
    # External momenta nodes
    ext_nodes = [key for key in list(pos.keys()) if key >= 100 ]
    # Internal nodes
    nodes = [key for key in list(pos.keys()) if 0 <= key < 100]

    # Draw
    nx.draw_networkx_edges(G, pos, edgelist=massive_edges, width=4)
    nx.draw_networkx_edges(G, pos, edgelist=massless_edges, width=1)
    nx.draw_networkx_nodes(
        G, pos, nodes, node_color='#dd7070', alpha=.95, node_size=800)
    nx.draw_networkx_nodes(G, pos, ext_nodes, node_color='#dd7070',
                           node_shape="d", alpha=.95, node_size=800)
    nx.draw_networkx_nodes(
        G, pos, ll_nodes, node_color='#7070dd', alpha=.95, node_size=800)
    
    nx.draw_networkx_edge_labels(G,pos, edge_labels=edge_labels)
    nx.draw_networkx_labels(G, pos,node_labels)

    # Show
    if file_name is not None:
        plt.savefig(file_name)
    if show_graph:
        plt.show()


# From the topology generator after calling generate_loop_topology
# will identify from the flow the edges that correspond to each LoopLine
def get_edge_info(topology_generator, topology: LoopTopology):
    # Create edge infos. loopline to be filled
    edges = topology_generator.edges
    edge_info = { k: {"edge": edges[v], 'm_squared': None, "loopline": None, "signature": None} for k, v in topology_generator.edge_name_map.items()}
    
    # Get LoopLine Information
    for n, ll in enumerate(topology.loop_lines):
        for m, prop in enumerate(ll.propagators):
            if prop.name is None:
                raise Exception("EdgeInfo_Error::Propagator is missing it's name! Poor guy :( ")
            edge_info[prop.name]['m_squared'] = prop.m_squared
            edge_info[prop.name]['loopline'] = (n,m)
            edge_info[prop.name]['signature'] = list(ll.signature)
    
    return edge_info

#TODO: Ideally this should become an automatically generated from some QGRAF output
class AmplitudeTopologies(object):

    def __init__(self):
        self.build_topology = {}

    def __call__(self, build_loop_topology):
        self.build_topology.update(
            {build_loop_topology.__name__: build_loop_topology})

    def create(self, amplitude_type: str, *args, fixed_deformation=None):
        try:
            self.build_topology[amplitude_type]
        except KeyError:
            print("KeyError::Callable function to build the Amplitude's Topology:")
            for key in self.build_topology.keys():
                print("\t- %s" % key)
            raise
        return self.build_topology[amplitude_type](*args, fixed_deformation=fixed_deformation)


amplitude_topologies = AmplitudeTopologies()

##############################################################
# gg -> HH ( Extended topology )
##############################################################


@amplitude_topologies
class DiHiggsTopology(object):
    #  Diag |      Propagators      |  Externals    | routing |
    # --------------------------------------------------------
    #   1   | [k,k-p1,k+p2,k-p1-p3] | (p1,p2,p4,p3) |  (1,2)  |
    #   2   | [k,k-p1,k+p2,k-p1-p4] | (p1,p2,p3,p4) |  (1,2)  |
    #   3   | [k,k-p1,k+p3,k-p1-p4] | (p1,p3,p2,p4) |  (1,2)  |
    #
    #
    # # Sugar-Topology
    #
    #    p3 ----------------- -p2
    #           |       |
    #   p2 -----|       |----- p2+p4-p3
    #           |       |
    #    p1 ----------------- p3
    #
    # [k, k+p2, k-p1-p4, k+p3, k-p1-p3, k-p1]

    def __init__(self, qs, ms, topology_name, fixed_deformation=None):
        # self.propagator_to_loopline = {"k": (0, 0),
        #                                "k+p2": (0, 1),
        #                                "k-p1-p4": (0, 2),
        #                                "k+p3": (0, 3),
        #                                "k-p1-p3": (0, 4),
        #                                "k-p1": (0, 5)}
        self.qs = (qs[0], qs[1], qs[2], -qs[1], qs[1]+qs[3]-qs[2], qs[2])
        self.ms = ms
        self.topology_name = topology_name
        self.top_mass = 174.0

        points = len(self.qs)
        # Create the Graph for the topology
        # pi: propagators, qi: externals
        self.graph = [('p%d' % (i+1), i+1, ((i+1) % points)+1)
                      for i in range(points)]
        self.graph.extend([('q%d' % (i+1), i+101, i+1) for i in range(points)])

        self.topology = self.build_loop_topology(
            fixed_deformation=fixed_deformation)

    def build_loop_topology(self, fixed_deformation=None):
        points = len(self.qs)

        self.top_generator = TopologyGenerator(self.graph)

        return self.top_generator.create_loop_topology(
            self.topology_name,
            ext_mom={'q%d' % (i+1): self.qs[i] for i in range(points)},
            mass_map={'p%d' % (i+1): self.top_mass for i in range(points)},
            # Force the loop momentum routing on the first edge
            loop_momenta_names=('p1',),
            analytic_result=None,
            fixed_deformation=fixed_deformation,
        )

    def get_edge_info(self):
        self.edge_info = get_edge_info(self.top_generator, self.topology)


##############################################################
# qqbar -> photons
##############################################################

@amplitude_topologies
class qqbarphotonsNLO(object):
    # # Sugar-Topology
    #   p2 ----[2]--->---[3]---- p3
    #           |         |
    #           ^         v       :
    #           |         |
    #   p1 ----[1]---<---[n]---- pn
    def __init__(self, qs, ms, topology_name, fixed_deformation=None):
        self.qs = copy.deepcopy(qs)
        self.ms = copy.deepcopy(ms)
        self.topology_name = topology_name

        # Create map
        points = len(self.qs)
        self.propagator_to_loopline = {"".join(
            ["k"]+["+p%d" % (i+1) for i in range(n+1)]): (0, n) for n in range(points-1)}
        self.propagator_to_loopline.update({"k": (0, points-1)})

        # Create the Graph for the topology
        # pi: propagators, qi: externals
        self.graph = [('p%d' % (i+1), i+1, ((i+1) % points)+1)
                      for i in range(points)]

        self.graph.extend([('q%d' % (i+1), i+101, i+1) for i in range(points)])
        
        self.topology = self.build_loop_topology(fixed_deformation=fixed_deformation)

    def build_loop_topology(self, fixed_deformation=None):
        points = len(self.qs)

        self.top_generator = TopologyGenerator(self.graph)

        print(self.graph)
        return self.top_generator.create_loop_topology(
            self.topology_name,
            ext_mom={'q%d' % (i+1): self.qs[i] for i in range(points)},
            mass_map={'p%d' % (i+1): self.ms[i] for i in range(points)},
            # If not specified an arbitrary spanning tree will be used for momentum routing
            #loop_momenta_names=('p%d' % points,),
            # For triangle and box one-loop topology, the analytic result is automatically computed
            loop_momenta_names=('p%d' % points,),
            #loop_momenta_names=('p1',),
            analytic_result=None,
            fixed_deformation=fixed_deformation,
        )

    def get_edge_info(self):
        self.edge_info = get_edge_info(self.top_generator, self.topology)


@amplitude_topologies
class Nf_qqbarphotonsNNLO(object):
    # # Sugar-Topology
    #                   2
    #  p2 -->--[3]------->------[4]--<-- p3
    #           |                |        :
    #           ^  1             | 3
    #           |                |
    #          [2]               V        :
    #          / \               |
    #   n+1   ^   ^  n+2         | n-1
    #          \ /               |        :
    #  p1 -->--[1]-------<------[n+1]--<-- pn
    #                    n

    def __init__(self, qs, m_uv, topology_name, fixed_deformation=None):
        self.qs = copy.deepcopy(qs)
        self.topology_name = topology_name

        # Number of TRUE externals
        points = len(self.qs)
        # Set masses for each propagator
        self.ms = [0 for i in range(points+2)]

        # Add zero momentum insertion for UV massive propagators
        self.qs.insert(1, vectors.LorentzVector([0, 0, 0, 0]))
        self.qs += [vectors.LorentzVector([0, 0, 0, 0])] * 2 
        # Always define the first of the two identical propagtors to be the massive one 
        self.ms.insert(0, 1e2)
        self.ms.insert(points+2, 1e2) 
        self.ms.insert(points+1, 1e2)
         

        # Create the Graph for the topology
        # pi: propagators, qi: externals
        # Start with the graph for the NLO graph starting from vertex 2
        self.graph = [('p%d' % (i+1), i+2, ((i+2) % (points+2))+1)
                      for i in range(points+1)]
        # Add the bubble with two additional propagators for the UV CTs
        self.graph.extend([('p%d' % (points+2), 1, points + 3), ('p%d' % (points+3), points + 3, 2)])
        self.graph.extend([('p%d' % (points+4), 1, points + 4), ('p%d' % (points+5), points + 4, 2)])
        
        # Externals
        self.graph.extend([('q1', 101, 1)])
        self.graph.extend([('q%d' % (i+2), i+102, i+3)
                           for i in range(points+2)])

        self.topology = self.build_loop_topology(
            fixed_deformation=fixed_deformation)

    def build_loop_topology(self, fixed_deformation=None):
        points = len(self.qs)

        self.top_generator = TopologyGenerator(self.graph)

        return self.top_generator.create_loop_topology(
            self.topology_name,
            ext_mom={'q%d' % (i+1): self.qs[i] for i in range(points)},
            mass_map={'p%d' % (i+1): self.ms[i] for i in range(len(self.ms))},
            # If not specified an arbitrary spanning tree will be used for momentum routing
            #loop_momenta_names=('p%d' % (points+1), 'p%d' % (points+2)),
            loop_momenta_names=( 'p%d' % (points-2), 'p%d' % (points+2)),
            analytic_result=None,
            fixed_deformation=fixed_deformation,
        )

    def get_edge_info(self):
        self.edge_info = get_edge_info(self.top_generator, self.topology)


#############################################################################################################
# Create the collection of hard-coded topologies.
#############################################################################################################

if __name__ == "__main__":
    
    topology_collection = TopologyCollection()
    print(amplitude_topologies.build_topology)
    
    # TEST
    moms = [vectors.LorentzVector([1, 0, 0, 1]),
            vectors.LorentzVector([1, 0, 0, -1]),
            vectors.LorentzVector([-1, 0, 1.0/np.sqrt(2), 1.0/np.sqrt(2)]),
            vectors.LorentzVector([-1, 0, -1.0/np.sqrt(2), -1.0/np.sqrt(2)])
            ]
    
    topology_name = "dd2A_NLO"
    print(topology_name)
    ms = [0.]*4
    amp_top = amplitude_topologies.create(
        "qqbarphotonsNLO", moms, ms, topology_name, fixed_deformation=False)
    amp_top.get_edge_info()
    topology_collection.add_topology(amp_top.topology, topology_name)
    topology_to_graph(amp_top.edge_info, file_name="diag_"+topology_name+".pdf", show_graph=False)
    
    topology_name = "dd2A_NNLO"
    print(topology_name)
    ms = [0.]*6
    amp_top = amplitude_topologies.create(
        "Nf_qqbarphotonsNNLO", moms, ms, topology_name, fixed_deformation=False)
    topology_collection.add_topology(amp_top.topology, topology_name)
    amp_top.get_edge_info()
    topology_to_graph(amp_top.edge_info, file_name="diag_"+topology_name+".pdf", show_graph=False)
    
    topology_name = "dihiggs"
    print(topology_name)
    ms = [0.]*4
    amp_top = amplitude_topologies.create(
        "DiHiggsTopology", moms, ms, topology_name, fixed_deformation=False)
    topology_collection.add_topology(amp_top.topology, topology_name)
    amp_top.get_edge_info()
    topology_to_graph(amp_top.edge_info, file_name="diag_"+topology_name+".pdf", show_graph=False)
    # print(amp_top.propagator_to_loopline)
    #topology_collection.add_topology(amp_top.build_loop_topology(), topology_name)
    
    
    print(topology_collection)
