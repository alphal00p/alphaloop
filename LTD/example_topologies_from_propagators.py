#!/usr/bin/env python3
import os
import vectors
import math
import sys

import mpmath
import numpy
from scipy.special import zeta

zero_lv = vectors.LorentzVector([0.,0.,0.,0.])

from ltd_utils import LoopTopology, LoopLine, Propagator, TopologyCollection, TopologyGeneratorFromPropagators 
from analytical_expressions  import ladder_phi, analytic_four_point_ladder, analytic_three_point_ladder, analytic_two_point_ladder

#############################################################################################################
# Create the collection of hard-coded topologies.
#############################################################################################################

# Add Now automatically generated topologies
def load(selected_topologies=None):

    all_topologies = TopologyCollection()

    topology_name = "GGHH_TEST_A"
    if selected_topologies is None or topology_name in selected_topologies:
        
        rescaling = 1.0
        m_h = 125.0
        m_t = 173.0
        sqrt_s = 300.0

        p1 = vectors.LorentzVector(
            [ sqrt_s/2., 0., 0.,  sqrt_s/2. ]
        )*rescaling 
        p2 = vectors.LorentzVector(
            [ sqrt_s/2., 0., 0., -sqrt_s/2. ]
        )*rescaling
        p3 = vectors.LorentzVector(
            [ -sqrt_s, 0., 0., 0., ]
        )*rescaling

        all_propagators = {
            (0, 1) : (
                Propagator( q=-p1,      m_squared = (m_t*rescaling)**2,     power=1 ), # g1
                Propagator( q=zero_lv,  m_squared = (m_t*rescaling)**2,     power=1 ), # g2
#                Propagator( q=p3,       m_squared = (m_t*rescaling)**2,     power=1 )  # g4
            ),
            (1, -1) : (
                Propagator( q=zero_lv,  m_squared = (m_t*rescaling)**2,     power=1 ), # g3
#                Propagator( q=p1+p3,    m_squared = (m_t*rescaling)**2,     power=1 ), # g7
            ),
            (1, 0) : (
                Propagator( q=zero_lv,  m_squared = (m_h*rescaling)**2,     power=1 ), # g5
                Propagator( q=p3,       m_squared = (m_h*rescaling)**2,     power=1 ), # g6
            )
        }
        factory = TopologyGeneratorFromPropagators(all_propagators)
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                # 5.69785496028533 + 3.15127029975641*I +/- ( 0.00524247904438823 + 0.00298271311291565*I )
                analytic_result= (5.69785496028533 + 3.15127029975641j)*(( ((math.pi**2)*1j) )**2 / (2.0*math.pi)**8 ) / (rescaling**2)
             ),
             entry_name = topology_name
        )


    topology_name = "GGHH_TEST_B"
    if selected_topologies is None or topology_name in selected_topologies:
        
        rescaling = 1.0
        m_h = 125.0
        m_t = 173.0
        sqrt_s = 300.0

        p1 = vectors.LorentzVector(
            [ sqrt_s/2., 0., 0.,  sqrt_s/2. ]
        )*rescaling 
        p2 = vectors.LorentzVector(
            [ sqrt_s/2., 0., 0., -sqrt_s/2. ]
        )*rescaling
        p3 = vectors.LorentzVector(
            [ -sqrt_s, 0., 0., 0., ]
        )*rescaling

        all_propagators = {
            (0, 1) : (
                Propagator( q=-p1,      m_squared = (m_t*rescaling)**2,     power=1 ), # g1
                Propagator( q=zero_lv,  m_squared = (m_t*rescaling)**2,     power=1 ), # g2
#                Propagator( q=p3,       m_squared = (m_t*rescaling)**2,     power=1 )  # g4
            ),
            (1, -1) : (
                Propagator( q=zero_lv,  m_squared = (m_t*rescaling)**2,     power=1 ), # g3
                Propagator( q=p1+p3,    m_squared = (m_t*rescaling)**2,     power=1 ), # g7
            ),
            (1, 0) : (
                Propagator( q=zero_lv,  m_squared = (m_h*rescaling)**2,     power=1 ), # g5
                Propagator( q=p3,       m_squared = (m_h*rescaling)**2,     power=1 ), # g6
            )
        }
        factory = TopologyGeneratorFromPropagators(all_propagators)
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                # 5.69785496028533 + 3.15127029975641*I +/- ( 0.00524247904438823 + 0.00298271311291565*I )
                analytic_result= (-0.614411975169297 - 0.56576230850425j)*(( ((math.pi**2)*1j) )**2 / (2.0*math.pi)**8 ) / (rescaling**4)
             ),
             entry_name = topology_name
        )

    return all_topologies

if __name__=='__main__':
    
    args = sys.argv[1:]

    # VH topologies
    VH_topologies = [ 
        'GGHH_TEST_A'
    ]

    if len(args)==0:
        print("Now processing all topologies...")
        raw_list_of_selected_topologies = [ None, ]
    else:
        raw_list_of_selected_topologies = args
        print("Now processing the following topologies: %s"%(', '.join(args)))

    list_of_selected_topologies = []
    for topo_name in raw_list_of_selected_topologies:
        if topo_name == 'VH_topologies':
            list_of_selected_topologies.append(VH_topologies)
        else:
            if topo_name is None:
                list_of_selected_topologies = [None,]
                break
            else:    
                list_of_selected_topologies.append([topo_name,])
   
    for selected_topologies in list_of_selected_topologies:
        all_topologies = load(selected_topologies)

        # Write out topoloies one by one in separate yaml files given the time it
        # takes to generate them and their size.
        for topology_name, topology in all_topologies.items():
            TopologyCollection({topology_name:topology}).export_to(os.path.join('topologies/.', '%s.yaml'%topology_name))

