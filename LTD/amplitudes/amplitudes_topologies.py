import os
import vectors
import math

import mpmath
import numpy as np
from scipy.special import zeta

zero_lv = vectors.LorentzVector([0.,0.,0.,0.])

from ltd_utils import LoopTopology, LoopLine, Propagator, TopologyCollection, TopologyGenerator
from analytical_expressions  import ladder_phi, analytic_four_point_ladder, analytic_three_point_ladder, analytic_two_point_ladder

#############################################################################################################
# Create the collection of function to build Amplitude Topologies
#############################################################################################################

# Ideally this should become an automatically generated from some QGRAF output

class AmplitudeTopologies(object):
    
    def __init__(self):
        self.types = {}
    
    def __call__(self,build_loop_topology):
        self.types.update({build_loop_topology.__name__:build_loop_topology})
    
    def call(self,amplitude_type: str, *args) -> LoopTopology:
        
        try:
            self.types[amplitude_type]
        except KeyError:
            print("KeyError::Callable function to build the Amplitude's Topology:")
            for key in self.types.keys():
                print("\t- %s" % key)
            raise 
        
        return self.types[amplitude_type](*args)

amplitude_topologies = AmplitudeTopologies()

##############################################################
# gg -> HH (aka: Box topology)
##############################################################

@amplitude_topologies
def DiHiggsTopology(qs, ms, topology_name):
    points = len(qs)
    # Create the Graph for the topology
    # pi: propagators, qi: externals
    graph = [('p%d' % (i+1), i+1, ((i+1) % points)+1) for i in range(points)]
    graph.extend([('q%d' % (i+1), i+101, i+1) for i in range(points)])
    
    mytop = TopologyGenerator(graph)
    print(graph)
    
    return mytop.create_loop_topology(
        topology_name,
        ext_mom={'q%d' % (i+1): qs[i] for i in range(points)},
        mass_map={'p%d' % (i+1): ms[i] for i in range(points)},
        # If not specified an arbitrary spanning tree will be used for momentum routing
        loop_momenta_names=('p%d' % points,),
        analytic_result=None,
        )


##############################################################
# qqbar -> photons 
##############################################################

@amplitude_topologies
def qqbarphotonsNLO(qs, ms, topology_name):
    points = len(qs)
    # Create the Graph for the topology
    # pi: propagators, qi: externals
    graph = [('p%d' % (i+1), i+1, ((i+1) % points)+1) for i in range(points)]
    graph.extend([('q%d' % (i+1), i+101, i+1) for i in range(points)])

    mytop = TopologyGenerator(graph)

    return mytop.create_loop_topology(
        topology_name,
        ext_mom={'q%d' % (i+1): qs[i] for i in range(points)},
        mass_map={'p%d' % (i+1): ms[i] for i in range(points)},
        # If not specified an arbitrary spanning tree will be used for momentum routing
        loop_momenta_names=('p%d' % points,),
        # For triangle and box one-loop topology, the analytic result is automatically computed
        analytic_result=None,
        fixed_deformation=None,
    )


@amplitude_topologies
def Nf_qqbarphotonsNNLO(qs, ms, topology_name):
    points = len(qs)
    # Create the Graph for the topology
    # pi: propagators, qi: externals

    # Propagators
    graph = [('p%d' % (i+1), i+1, ((i+1) % (points+1))+1) for i in range(points+1)]
    graph.extend([('p%d' % (points+2), 2, 3)])
    # Externals
    graph.extend([('q1', 101, 1)])
    graph.extend([('q%d' % (i+2), i+102, i+3) for i in range(points-1)])
   
    mytop = TopologyGenerator(graph)
    mytop.n_loops = 2

    return mytop.create_loop_topology(
            topology_name,
            ext_mom={'q%d' % (i+1): qs[i] for i in range(points)},
            mass_map={'p%d' % (i+1): ms[i] for i in range(len(ms))},
            # If not specified an arbitrary spanning tree will be used for momentum routing
            loop_momenta_names=('p%d' % (points+1),'p%d' % (points+2)),
            analytic_result=None,
            fixed_deformation=None,
            )


#############################################################################################################
# Create the collection of hard-coded topologies.
#############################################################################################################

topology_collection = TopologyCollection()
print(amplitude_topologies.types)

# TEST
moms = [vectors.LorentzVector([1,0,0,1]),
        vectors.LorentzVector([1,0,0,-1]),
        vectors.LorentzVector([-1,0,1.0/np.sqrt(2), 1.0/np.sqrt(2)]),
        vectors.LorentzVector([-1,0,-1.0/np.sqrt(2), -1.0/np.sqrt(2)])
        ]

topology_name = "dd2A_NNLO"
ms = [0.]*6
amp_top = amplitude_topologies.call("Nf_qqbarphotonsNNLO",moms,ms, topology_name)
amp_top = topology_collection.add_topology(amp_top, topology_name)

topology_name = "dihiggs"
ms = [0.]*4
amp_top = amplitude_topologies.call("DiHiggsTopology",moms, ms, topology_name)
topology_collection.add_topology(amp_top, topology_name)


print(topology_collection)

