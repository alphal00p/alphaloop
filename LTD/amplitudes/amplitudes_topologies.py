import os
import sys
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
def topology_to_graph(amplitude_topology, file_name=None, show_graph=True):
    # Great graph with oriented edges
    # - Internal nodes are labelled from 0 to 99
    # - External nodes are labelled from 100 to 199
    # - LoopLines nodes are labelled from 200 up

    G = nx.DiGraph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
    edges = [e[1:] for e in amplitude_topology.graph]

    edges = [e[1:] for e in amplitude_topology.graph if e[1] >= 100]
    for n, ll in enumerate(amplitude_topology.ll_chains):
        for m, e in enumerate(ll):
            edges += [(e[0], (n+2)*100+m), ((n+2)*100+m, e[1])]

    G.add_edges_from(edges)
    plt.figure(figsize=(12, 12))
    pos = nx.kamada_kawai_layout(G)  # positions for all nodes
    labels = {key: str(key) for key in pos.keys()}

    # LoopLines nodes
    ll_nodes = [key for key in list(pos.keys()) if key >= 200]
    for ll_node in ll_nodes:
        labels[ll_node] = "[{},{}]".format(int(ll_node/100)-2, ll_node % 100)
    # External momenta nodes
    ext_nodes = [key for key in list(pos.keys()) if 100 <= key < 200]
    for ext_node in ext_nodes:
        labels[ext_node] = "q{}".format(ext_node % 100)

    # Internal nodes
    nodes = [key for key in list(pos.keys()) if key < 100]

    # Draw
    nx.draw_networkx_edges(G, pos, node_size=1000)
    nx.draw_networkx_nodes(
        G, pos, nodes, node_color='#dd7070', alpha=.95, node_size=800)
    nx.draw_networkx_nodes(G, pos, ext_nodes, node_color='#dd7070',
                           node_shape="d", alpha=.95, node_size=800)
    nx.draw_networkx_nodes(
        G, pos, ll_nodes, node_color='#7070dd', alpha=.95, node_size=800)
    nx.draw_networkx_labels(G, pos, labels)

    # Show
    if file_name is not None:
        plt.savefig(file_name)
    if show_graph:
        plt.show()


# From the topology generator after calling generate_loop_topology
# will identify from the flow the edges that correspond to each LoopLine
def get_looplines_edges(topology_generator, topology):
    flows = topology_generator.flows
    # [: -len(self.topology.external_kinematics)]
    edges = topology_generator.edges

    lines = {y[0] for x in flows for y in x}
    res = {l: [[0 for i in range(len(flows))], edges[l]] for l in lines}

    for n, flow in enumerate(flows):
        for l in flow:
            if l[1]:
                res[l[0]][0][n] = 1
            else:
                res[l[0]][0][n] = -1
    res = list(res.values())
    ll_chains = []
    for ll in topology.loop_lines:
        s = list(ll.signature)
        start = ll.start_node
        end = ll.end_node

        es = [e[1] for e in res if e[0] == s]
        chain = []
        node = start
        while True:
            if len(es) == 0:
                assert(node == end)
                break
            pos = np.where([e[0] == node for e in es])
            if len(pos[0]) != 1:
                raise Exception("No unique edge found for the LoopLine")

            chain += [es.pop(pos[0][0])]
            node = chain[-1][1]
        ll_chains += [chain]
    return ll_chains


# Ideally this should become an automatically generated from some QGRAF output


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
            topology_name,
            ext_mom={'q%d' % (i+1): self.qs[i] for i in range(points)},
            mass_map={'p%d' % (i+1): self.top_mass for i in range(points)},
            # Force the loop momentum routing on the first edge
            loop_momenta_names=('p1',),
            analytic_result=None,
            fixed_deformation=fixed_deformation,
        )

    def get_ll_chains(self):
        self.ll_chains = get_looplines_edges(self.top_generator, self.topology)


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
        self.qs = qs
        self.ms = ms
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

        self.topology = self.build_loop_topology(
            fixed_deformation=fixed_deformation)

    def build_loop_topology(self, fixed_deformation=None):
        points = len(self.qs)

        self.top_generator = TopologyGenerator(self.graph)

        return self.top_generator.create_loop_topology(
            topology_name,
            ext_mom={'q%d' % (i+1): self.qs[i] for i in range(points)},
            mass_map={'p%d' % (i+1): self.ms[i] for i in range(points)},
            # If not specified an arbitrary spanning tree will be used for momentum routing
            loop_momenta_names=('p%d' % points,),
            # For triangle and box one-loop topology, the analytic result is automatically computed
            analytic_result=None,
            fixed_deformation=fixed_deformation,
        )

    def get_ll_chains(self):
        self.ll_chains = get_looplines_edges(self.top_generator, self.topology)


@amplitude_topologies
class Nf_qqbarphotonsNNLO(object):
    # # Sugar-Topology
    #
    #  p2 -->--[3]------->------[4]--<-- p3
    #          / \               |
    #         ^   ^              |
    #          \ /               V        :
    #          [1]               |        :
    #           |                |        :
    #           ^                |
    #           |                |
    #  p1 -->--[1]-------<------[n+1]--<-- pn

    def __init__(self, qs, ms, topology_name, fixed_deformation=None):
        self.qs = qs
        self.ms = ms
        self.topology_name = topology_name

        # Create Map
        points = len(self.qs)
        self.propagator_to_loopline = {"".join(
            ["k"]+["+p%d" % i for i in range(2, n+3)]): (0, n+2) for n in range(points-1)}
        self.propagator_to_loopline.update({"k": (0, 0)})

        # Create the Graph for the topology
        # pi: propagators, qi: externals
        self.graph = [('p%d' % (i+1), i+1, ((i+1) % (points+1))+1)
                      for i in range(points+1)]
        self.graph.extend([('p%d' % (points+2), 2, 3)])
        #self.graph.extend([('p%d' % (points+2), 3, 2)])
        # Externals
        self.graph.extend([('q1', 101, 1)])
        self.graph.extend([('q%d' % (i+2), i+102, i+3)
                           for i in range(points-1)])

        self.topology = self.build_loop_topology(
            fixed_deformation=fixed_deformation)

    def build_loop_topology(self, fixed_deformation=None):
        points = len(self.qs)

        self.top_generator = TopologyGenerator(self.graph)

        return self.top_generator.create_loop_topology(
            topology_name,
            ext_mom={'q%d' % (i+1): self.qs[i] for i in range(points)},
            mass_map={'p%d' % (i+1): self.ms[i] for i in range(len(self.ms))},
            # If not specified an arbitrary spanning tree will be used for momentum routing
            loop_momenta_names=('p%d' % (points+1), 'p%d' % (points+2)),
            analytic_result=None,
            fixed_deformation=fixed_deformation,
        )

    def get_ll_chains(self):
        self.ll_chains = get_looplines_edges(self.top_generator, self.topology)


#############################################################################################################
# Create the collection of hard-coded topologies.
#############################################################################################################

topology_collection = TopologyCollection()
print(amplitude_topologies.build_topology)

# sys.exit(0)

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
amp_top.get_ll_chains()
topology_collection.add_topology(amp_top.topology, topology_name)
topology_to_graph(amp_top, file_name="diag_" +
                  topology_name+".pdf", show_graph=False)

topology_name = "dd2A_NNLO"
print(topology_name)
ms = [0.]*6
amp_top = amplitude_topologies.create(
    "Nf_qqbarphotonsNNLO", moms, ms, topology_name, fixed_deformation=False)
topology_collection.add_topology(amp_top.topology, topology_name)
amp_top.get_ll_chains()
topology_to_graph(amp_top, file_name="diag_" +
                  topology_name+".pdf", show_graph=False)

topology_name = "dihiggs"
print(topology_name)
ms = [0.]*4
amp_top = amplitude_topologies.create(
    "DiHiggsTopology", moms, ms, topology_name, fixed_deformation=False)
topology_collection.add_topology(amp_top.topology, topology_name)
amp_top.get_ll_chains()
topology_to_graph(amp_top, file_name="diag_" +
                  topology_name+".pdf", show_graph=False)
# print(amp_top.propagator_to_loopline)
#topology_collection.add_topology(amp_top.build_loop_topology(), topology_name)


print(topology_collection)
