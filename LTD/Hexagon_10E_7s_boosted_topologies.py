#!/usr/bin/env python3
import os
import vectors
import math

import mpmath
import numpy
from scipy.special import zeta

zero_lv = vectors.LorentzVector([0.,0.,0.,0.])

from ltd_utils import LoopTopology, LoopLine, Propagator, TopologyCollection, TopologyGenerator
from analytical_expressions  import ladder_phi, analytic_four_point_ladder, analytic_three_point_ladder, analytic_two_point_ladder

#############################################################################################################
# Create the collection of hard-coded topologies.
#############################################################################################################

hard_coded_topology_collection = TopologyCollection()

import itertools

# Add Now automatically generated topologies

scaling = 1.e-02
base_q1 = vectors.LorentzVector([42.6    ,21.2   ,-2.6   ,0])*scaling #m1^2 = 1358.6 *1e-6
base_q2 = vectors.LorentzVector([45      ,3.4    ,-40.2  ,0])*scaling #m2^2 = 397.4
base_q3 = vectors.LorentzVector([68      ,-64.4  ,6.4    ,0])*scaling #m3^2 = 435.68
base_q4 = vectors.LorentzVector([-57.4    ,-48.8  ,-4.2  ,0])*scaling #m4^2 = 895.68
base_q5 = vectors.LorentzVector([-39.5   ,65.5   ,94.0   ,0])*scaling #m5^2 = -11566
all_qs = [base_q1, base_q2, base_q3, base_q4, base_q5]

all_combinations = list(itertools.combinations(all_qs,1))
all_combinations.extend(list(itertools.combinations(all_qs,2)))
all_combinations_names = list(itertools.combinations(range(1,len(all_qs)+1),1))
all_combinations_names.extend(list(itertools.combinations(range(1,len(all_qs)+1),2)))

for ((comb, comb_elems), z_boost) in itertools.product(
                        zip(all_combinations, all_combinations_names), [None, 10., 100., 1000.]):
    
    hexagon = TopologyGenerator([
        ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),  ('p4', 4, 5),  ('p5', 5, 6), ('p6', 6, 1),
        ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6)
    ])
    boost_vecA = sum(comb)
    if z_boost:
        boost_vecA += vectors.LorentzVector([0., 0., 0., z_boost*scaling])
    if boost_vecA.square()>0.:
        boost_vecB = vectors.LorentzVector([math.sqrt(boost_vecA.square())    ,0.   ,0.   ,0])
    else:
        boost_vecB = vectors.LorentzVector([0.    ,0.   ,math.sqrt(abs(boost_vecA.square())) ,0. ])
    print("Now considering boost from %s (=%s) to %s (=its rest frame%s)."%(
        str(boost_vecA), 
        '+'.join('q%d'%i for i in comb_elems),
        str(boost_vecB),
        ', with a z-boost of %.3f'%z_boost if z_boost else ''
    ))
    q1 = base_q1.get_copy().rotoboost(boost_vecA,boost_vecB)
    q2 = base_q2.get_copy().rotoboost(boost_vecA,boost_vecB)
    q3 = base_q3.get_copy().rotoboost(boost_vecA,boost_vecB)
    q4 = base_q4.get_copy().rotoboost(boost_vecA,boost_vecB)
    q5 = base_q5.get_copy().rotoboost(boost_vecA,boost_vecB)
    
    boost_name = '_boost_%s'%('_'.join('%d'%qi for qi in comb_elems ))
    if z_boost:
        boost_name += '_z_boosted_%d'%int(z_boost)
    hard_coded_topology_collection.add_topology(hexagon.create_loop_topology(
            "Hexagon_10E_7s%s"%boost_name, 
            ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6': -q5-q4-q3-q2-q1  }, 
            mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0, 'p6': 0.0}, 
            loop_momenta_names=('p1',), # If not specified an arbitrary spanning tree will be used for momentum routing 
            analytic_result=2.11928148966e-02+6.4030325864e-03j
         ),
         entry_name = 'Hexagon_10E_7s%s'%boost_name
    )

if __name__=='__main__':

    hard_coded_topology_collection.export_to(os.path.join('.', 'Hexagon_10E_7s_boosted_topologies.yaml'))
