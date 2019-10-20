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

# The fishnet topology generation takes quite some time. It is therefore disabled by default.
# Enable it by setting the flag below to True.
_GENERATE_FISHNETS = False 

hard_coded_topology_collection = TopologyCollection()

# Add the manually crafted topologies
import manually_crafted_topologies
for name, topology in manually_crafted_topologies.hard_coded_topology_collection.items():
    topology.name='manual_%s'%topology.name
    hard_coded_topology_collection.add_topology(topology, entry_name = topology.name)

# Add Now automatically generated topologies

# triangle
triangle = TopologyGenerator([('q1', 0, 1), ('p1', 1, 2), ('p2', 2, 3),
                              ('p3', 3, 1), ('q3', 4, 3), ('q2', 5, 2)])
q1 = vectors.LorentzVector([ 0.1, 0.2, 0.5, 0.1 ])
q2 = vectors.LorentzVector([-0.3, 0.4, 0.1, 0.2 ])
hard_coded_topology_collection.add_topology(triangle.create_loop_topology(
        "Triangle_no_ellipse", 
        ext_mom={ 'q1': q1, 'q2': q2 , 'q3': -q1-q2 }, 
        mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0}, 
        loop_momenta_names=('p3',), # If not specified an arbitrary spanning tree will be used for momentum routing 
        analytic_result=None # For triangle and box one-loop topology, the analytic result is automatically computed
     ),
     entry_name = 'Triangle_no_ellipse'
)

# box
box = TopologyGenerator([
    ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),  ('p4', 4, 1),
    ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4)
])
q1 = vectors.LorentzVector([ 0.1, 0.2, 0.5, 0.1 ])
q2 = vectors.LorentzVector([-0.3, 0.4, 0.1, 0.2 ])
q3 = vectors.LorentzVector([ 0.1, 0.2, 0.5, 0.3 ])
hard_coded_topology_collection.add_topology(box.create_loop_topology(
        "Box_no_ellipse", 
        ext_mom={ 'q1': q1, 'q2': q2 , 'q3': q3, 'q4': -q1-q2-q3 }, 
        mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0}, 
        loop_momenta_names=('p1',), # If not specified an arbitrary spanning tree will be used for momentum routing 
        analytic_result=None # For triangle and box one-loop topology, the analytic result is automatically computed
     ),
     entry_name = 'Box_no_ellipse'
)

# Box with customised ellipses from mathematica

box = TopologyGenerator([
    ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),  ('p4', 4, 1),
    ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4)
])
q1 = vectors.LorentzVector([23.2, -9.4, -5.6, 0.])
q2 = vectors.LorentzVector([-50., -18.4, -17.8, 0.])
q3 = vectors.LorentzVector([30.4, 2.6, -1.8, 0.])
hard_coded_topology_collection.add_topology(box.create_loop_topology(
        "Box_3E", 
        ext_mom={ 'q1': q1, 'q2': q2 , 'q3': q3, 'q4': -q1-q2-q3 }, 
        mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0}, 
        loop_momenta_names=('p1',), # If not specified an arbitrary spanning tree will be used for momentum routing 
        analytic_result=None, # For triangle and box one-loop topology, the analytic result is automatically computed
        # For now specified by hand as the cvxpy automated implementation is not done yet
    #    fixed_deformation = [{'deformation_sources': [[0., -9.01288e-6, -5.91311e-6, 0.]], 'weight_per_source': [1], 'excluded_surface_ids': [5]},
     #                        {'deformation_sources': [[0., 18.399985727313332, 17.800014049513244, 0.]], 'excluded_surface_ids': [2]}]
     ),
     entry_name = 'Box_3E'
)

# Box with customised ellipses from mathematica

box = TopologyGenerator([
    ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),  ('p4', 4, 1),
    ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4)
])
q1 = vectors.LorentzVector([21.4,13.8,-5.2,0])
q2 = vectors.LorentzVector([5.6,-1.6,27.8,0])
q3 = vectors.LorentzVector([32.6,6.6,9.8,0])
hard_coded_topology_collection.add_topology(box.create_loop_topology(
        "Box_5E", 
        ext_mom={ 'q1': q1, 'q2': q2 , 'q3': q3, 'q4': -q1-q2-q3 }, 
        mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0}, 
        loop_momenta_names=('p1',), # If not specified an arbitrary spanning tree will be used for momentum routing 
        analytic_result=None, # For triangle and box one-loop topology, the analytic result is automatically computed
        # For now specified by hand as the cvxpy automated implementation is not done yet
      #  fixed_deformation = [{'deformation_sources': [[0., 0., 0., 0.]], 'excluded_surface_ids': []}]
     ),
     entry_name = 'Box_5E'
)
# Pentagon with customrised ellipses from mathematica

pentagon = TopologyGenerator([
    ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),  ('p4', 4, 5),  ('p5', 5, 1),
    ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q4', 105,5)
])
q1 = vectors.LorentzVector([0.123698110274668e+01, 0.652467807974028e+00, 0.633070356251368e+00, 0.682813086866660e+00])
q2 = vectors.LorentzVector([0.156570804799072e+01, 0.566351831093323e+00, 0.935202196659332e+00, 0.976187756902101e+00])
q3 = vectors.LorentzVector([-0.688161154310136e+00, +0.238451694824191e+00, +0.637562295790242e+00, +0.101098380420291e+00])
q4 = vectors.LorentzVector([-0.473302592019865e+01, -0.240729108908187e+01, -0.718450510177215e+00, -0.366510869539274e+01])
hard_coded_topology_collection.add_topology(pentagon.create_loop_topology(
        "Pentagon_pairwise_3E", 
        ext_mom={ 'q1': q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': -q4-q3-q2-q1  }, 
        mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0}, 
        loop_momenta_names=('p1',), # If not specified an arbitrary spanning tree will be used for momentum routing 
        analytic_result=(-1.52339813764031085e-3 + 2.04369604371007528e-3j), # For triangle and box one-loop topology, the analytic result is automatically computed
        # For now specified by hand as the cvxpy automated implementation is not done yet
    #    fixed_deformation = [{'deformation_sources': [[0., 0.0424834862261251, -1.5779576840628833, 0.47971132471067496]], 'excluded_surface_ids': [5]},]
     ),
     entry_name = 'Pentagon_pairwise_3E'
)


# Pentagon with customrised ellipses from mathematica

pentagon = TopologyGenerator([
    ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),  ('p4', 4, 5),  ('p5', 5, 1),
    ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q4', 105,5)
])
q1 = vectors.LorentzVector([-32.6,1.4,-2.6,0])
q2 = vectors.LorentzVector([-37,-0.4,-3,0])
q3 = vectors.LorentzVector([-32.6,-4.6,27.4, 0])
q4= vectors.LorentzVector([-27.8, -4.4, -4.2, 0])
hard_coded_topology_collection.add_topology(pentagon.create_loop_topology(
        "Pentagon_10E_1s", 
        ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': -q4-q3-q2-q1  }, 
        mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0}, 
        loop_momenta_names=('p1',), # If not specified an arbitrary spanning tree will be used for momentum routing 
        analytic_result=(9.42971e-11)/(16*(math.pi**2)), # For triangle and box one-loop topology, the analytic result is automatically computed
        # For now specified by hand as the cvxpy automated implementation is not done yet
   #     fixed_deformation = [{'deformation_sources': [[0., 0.0424834862261251, -1.5779576840628833, 0.47971132471067496]], 'excluded_surface_ids': [5]},]
     ),
     entry_name = 'Pentagon_10E_1s'
)

# Pentagon with customrised ellipses from mathematica

pentagon = TopologyGenerator([
    ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),  ('p4', 4, 5),  ('p5', 5, 1),
    ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q4', 105,5)
])
q1 = vectors.LorentzVector([-32.,19.4,-12.4, 0.])
q2 = vectors.LorentzVector([1.4,-33.2,-5.8,0.])
q3 = vectors.LorentzVector([17,8.4,-14.2,0.])
q4= vectors.LorentzVector([28.6,-4,-13.4,0])
hard_coded_topology_collection.add_topology(pentagon.create_loop_topology(
        "Pentagon_6E_4s", 
        ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': -q4-q3-q2-q1  }, 
        mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0}, 
        loop_momenta_names=('p1',), # If not specified an arbitrary spanning tree will be used for momentum routing 
        analytic_result=(-2.82935e-8 + 2.63906e-3j)/(16*(math.pi**2)), # For triangle and box one-loop topology, the analytic result is automatically computed
        # For now specified by hand as the cvxpy automated implementation is not done yet
   #     fixed_deformation = [{'deformation_sources': [[0., 0.0424834862261251, -1.5779576840628833, 0.47971132471067496]], 'excluded_surface_ids': [5]},]
     ),
     entry_name = 'Pentagon_6E_4s'
)

# Pentagon with customrised ellipses from mathematica

pentagon = TopologyGenerator([
    ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),  ('p4', 4, 5),  ('p5', 5, 1),
    ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q4', 105,5)
])
q1 = vectors.LorentzVector([42.6,21.2,-2.6,0])
q2 = vectors.LorentzVector([45,3.4,-40.2,0])
q3 = vectors.LorentzVector([68,-64.4,6.4,0])
q4= vectors.LorentzVector([-57.4,-48.8,-4.2,0])
hard_coded_topology_collection.add_topology(pentagon.create_loop_topology(
        "Pentagon_8E_5s", 
        ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': -q4-q3-q2-q1  }, 
        mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0}, 
        loop_momenta_names=('p1',), # If not specified an arbitrary spanning tree will be used for momentum routing 
        analytic_result=(-3.47721e-11 -4.37955e-13j)/(16.*(math.pi**2.)), # For triangle and box one-loop topology, the analytic result is automatically computed
        # For now specified by hand as the cvxpy automated implementation is not done yet
   #     fixed_deformation = [{'deformation_sources': [[0., 0.0424834862261251, -1.5779576840628833, 0.47971132471067496]], 'excluded_surface_ids': [5]},]
     ),
     entry_name = 'Pentagon_8E_5s'
)

# Hexagon with customised ellipses from mathematica

hexagon = TopologyGenerator([
    ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),  ('p4', 4, 5),  ('p5', 5, 6), ('p6', 6, 1),
    ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6)
])
q1 = vectors.LorentzVector([24.,-21.2,71.,0.])
q2 = vectors.LorentzVector([50.4,15.8,-18.8,0.])
q3 = vectors.LorentzVector([-0.2,46.2,8.6,0.])
q4 = vectors.LorentzVector([-33.2,2.6,-70.8,0.])
q5 = vectors.LorentzVector([-80.,-5.6,-40,0.])
hard_coded_topology_collection.add_topology(hexagon.create_loop_topology(
        "Hexagon_6E_4s", 
        ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6': -q5-q4-q3-q2-q1  }, 
        mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0, 'p6': 0.0}, 
        loop_momenta_names=('p1',), # If not specified an arbitrary spanning tree will be used for momentum routing 
        analytic_result=-2.7216841734816e-15-1.2089609024538e-14j,
     ),
     entry_name = 'Hexagon_6E_4s'
)

hexagon = TopologyGenerator([
    ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),  ('p4', 4, 5),  ('p5', 5, 6), ('p6', 6, 1),
    ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6)
])
q1 = vectors.LorentzVector([-80.,29.,-70.,0.])
q2 = vectors.LorentzVector([83.5,14.0,70.,0.])
q3 = vectors.LorentzVector([88.5,6.5,-6.,0.])
q4 = vectors.LorentzVector([36.5,-71.,97.5,0.])
q5 = vectors.LorentzVector([12.5,-83.5,-57.5,0.])
hard_coded_topology_collection.add_topology(hexagon.create_loop_topology(
        "Hexagon_10E_4s", 
        ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6': -q5-q4-q3-q2-q1  }, 
        mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0, 'p6': 0.0}, 
        loop_momenta_names=('p1',), # If not specified an arbitrary spanning tree will be used for momentum routing 
        analytic_result=-3.0193937848736e-17-7.73337287748906e-17j,
     ),
     entry_name = 'Hexagon_10E_4s'
)

hexagon = TopologyGenerator([
    ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),  ('p4', 4, 5),  ('p5', 5, 6), ('p6', 6, 1),
    ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6)
])
q1 = vectors.LorentzVector([-100.   ,-80.5  ,-56.5  ,0.])
q2 = vectors.LorentzVector([-64.    ,83.    ,-51.   ,0.])
q3 = vectors.LorentzVector([100.    ,-87.0  ,-47.5  ,0.])
q4 = vectors.LorentzVector([100.    ,92.5   ,-17.5  ,0.])
q5 = vectors.LorentzVector([100.    ,-39.5  ,-88.5  ,0.])
hard_coded_topology_collection.add_topology(hexagon.create_loop_topology(
        "Hexagon_9E_4s", 
        ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6': -q5-q4-q3-q2-q1  }, 
        mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0, 'p6': 0.0}, 
        loop_momenta_names=('p1',), # If not specified an arbitrary spanning tree will be used for momentum routing 
        analytic_result=2.83771892273698e-17+8.3141576190839e-18j,
     ),
     entry_name = 'Hexagon_9E_4s'
)

hexagon = TopologyGenerator([
    ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),  ('p4', 4, 5),  ('p5', 5, 6), ('p6', 6, 1),
    ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6)
])
q1 = vectors.LorentzVector([42.6    ,21.2   ,-2.6   ,0])*1e-03 #m1^2 = 1358.6 *1e-6
q2 = vectors.LorentzVector([45      ,3.4    ,-40.2  ,0])*1e-03 #m2^2 = 397.4
q3 = vectors.LorentzVector([68      ,-64.4  ,6.4    ,0])*1e-03 #m3^2 = 435.68
q4 = vectors.LorentzVector([-57.4    ,-48.8  ,-4.2  ,0])*1e-03 #m4^2 = 895.68
q5 = vectors.LorentzVector([-39.5   ,65.5   ,94.0   ,0])*1e-03 #m5^2 = -11566
hard_coded_topology_collection.add_topology(hexagon.create_loop_topology(
        "Hexagon_10E_7s", 
        ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6': -q5-q4-q3-q2-q1  }, 
        mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0, 'p6': 0.0}, 
        loop_momenta_names=('p1',), # If not specified an arbitrary spanning tree will be used for momentum routing 
        analytic_result=2.11928148966e-02+6.4030325864e-03j
     ),
     entry_name = 'Hexagon_10E_7s'
)

hexagon = TopologyGenerator([
    ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),  ('p4', 4, 5),  ('p5', 5, 6), ('p6', 6, 1),
    ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6)
])
q1 = vectors.LorentzVector([42.6    ,21.2   ,-2.6   ,0])*1e-02
q2 = vectors.LorentzVector([45      ,3.4    ,-40.2  ,0])*1e-02
q3 = vectors.LorentzVector([68      ,-64.4  ,6.4    ,0])*1e-02
q4 = vectors.LorentzVector([-57.4    ,-32.5  ,-4.2   ,0])*1e-02
q5 = vectors.LorentzVector([-39.5   ,65.5   ,14.5   ,0.])*1e-02
hard_coded_topology_collection.add_topology(hexagon.create_loop_topology(
        "Hexagon_10E_5s", 
        ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6': -q5-q4-q3-q2-q1  }, 
        mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0, 'p6': 0.0}, 
        loop_momenta_names=('p1',), # If not specified an arbitrary spanning tree will be used for momentum routing 
        analytic_result=-2.81475384+2.4732705j
     ),
     entry_name = 'Hexagon_10E_5s'
)

# Hexagon with a disconnected ellipsoid
hexagon = TopologyGenerator([
    ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),  ('p4', 4, 5),  ('p5', 5, 6), ('p6', 6, 1),
    ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6)
])
q1 = vectors.LorentzVector([-12.    ,-16.   ,-53.5  ,0])*1e-02
q2 = vectors.LorentzVector([47.     ,-3.5    ,-14.0  ,0])*1e-02
q3 = vectors.LorentzVector([-4.     ,-12.  ,14.5    ,0])*1e-02
q4 = vectors.LorentzVector([31.    ,10.5  ,-5.5   ,0])*1e-02
q5 = vectors.LorentzVector([3.5   ,39.5   ,88.5   ,0.])*1e-02
hard_coded_topology_collection.add_topology(hexagon.create_loop_topology(
        "Hexagon_6E_2s", 
        ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6': -q5-q4-q3-q2-q1  }, 
        mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0, 'p6': 0.0}, 
        loop_momenta_names=('p1',), # If not specified an arbitrary spanning tree will be used for momentum routing 
        analytic_result=-20.7013988797142+4.234325146404j
     ),
     entry_name = 'Hexagon_6E_2s'
)

# double triangle
doubletriangle = TopologyGenerator([
    ('q', 0, 1), ('p1', 1, 2), ('p2', 3, 1), ('p3', 2, 3),
    ('p4', 3, 4), ('p5', 4, 2), ('-q', 4, 5)])
q = vectors.LorentzVector([1.0, 1.3, 0.5, 2.1])
hard_coded_topology_collection.add_topology(doubletriangle.create_loop_topology(
        "DoubleTriangle_no_ellipse", 
        ext_mom={'q': q , '-q' : -q}, 
        mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0}, 
        loop_momenta_names=('p1', 'p5'), 
        analytic_result= (-(6.*zeta(3))/((16.*(math.pi**2))**2*q.square())),
        contour_closure = [0,0]
    ),
    entry_name = 'DoubleTriangle_no_ellipse'
)


# double triangle
doubletriangle = TopologyGenerator([
    ('q', 0, 1), ('p1', 1, 2), ('p2', 3, 1), ('p3', 2, 3),
    ('p4', 3, 4), ('p5', 4, 2), ('-q', 4, 5)])
q = vectors.LorentzVector([1.0, 0, 0., 0])
hard_coded_topology_collection.add_topology(doubletriangle.create_loop_topology(
        "DoubleTriangle_massive_physical", 
        ext_mom={'q': q , '-q' : -q}, 
        mass_map={'p1': 0., 'p2': 0., 'p3': 0., 'p4': 0., 'p5':0.}, 
        loop_momenta_names=('p1', 'p5'), 
        analytic_result= 0.,
        contour_closure = [0,0],
    #    fixed_deformation = [{'deformation_sources': [[0., 0., 0., 0.],[0.,0.,0.,0.]], 'weight_per_source': [1,1], 'excluded_loop_lines': [0, 1, 2], 'excluded_surface_ids': []},]
    ),
    entry_name = 'DoubleTriangle_massive_physical'
)


# trianglebox
trianglebox = TopologyGenerator([
    ('q1', 0, 1), ('p1', 2, 1), ('p2', 1, 3), ('p3', 2, 3),
    ('p4', 4, 3), ('p5', 6, 4), ('p6', 2,6), ('q2', 5, 4), ('q3', 7, 6)])
q1 = vectors.LorentzVector([1, 0, 0., 0.])
q2 = vectors.LorentzVector([-0.5, 0., 0., 0.2])
q3=-q1-q2
hard_coded_topology_collection.add_topology(trianglebox.create_loop_topology(
        "TriangleBox_physical", 
        ext_mom={'q1': q1 , 'q2' : q2, 'q3' : -q1-q2}, 
        mass_map={'p1': 0., 'p2': 0., 'p3': 0., 'p4': 0., 'p5':0., 'p6':0.}, 
        loop_momenta_names=('p1', 'p6'), 
        analytic_result= analytic_three_point_ladder(q2.square(),q3.square(),q1.square(),2),
        contour_closure = [1,0],
    #    fixed_deformation = [{'deformation_sources': [[0., 0., 0., 0.],[0.,0.,0.,0.]], 'weight_per_source': [1,1], 'excluded_loop_lines': [0, 1, 2], 'excluded_surface_ids': []},]
        
    ),
    entry_name = 'TriangleBox_physical'
)

# PRL_6p_2L
PRL_6p_2L = TopologyGenerator([
        ('q1', 101, 1), ('q2', 102, 2), ('q3', 103, 3), ('q4', 104, 4),
        ('p1', 1, 6), ('p2', 6, 8), ('p3', 7, 2), ('p4', 2, 1),
        ('p5', 7, 3), ('p6', 3, 4), ('p7', 4, 6), ('p8', 8, 9), ('p9', 9, 7),
        ('q5', 8, 108), ('q6', 9, 109)
])
q1 = vectors.LorentzVector([ 0.2, 0.3, 0.5, 0.6 ])
q2 = vectors.LorentzVector([-0.1, 0.7, 0.2, 0.1])
q3 = vectors.LorentzVector([ 0.1, 0.5, -0.3, -0.4])
q4 = vectors.LorentzVector([-0.3, 0.4, 0.5, 0.2])
q5 = vectors.LorentzVector([-0.2, 0.3, 0.2, -0.5])
q6 = -q1-q2-q3-q4-q5
hard_coded_topology_collection.add_topology(PRL_6p_2L.create_loop_topology(
        "PRL_6p_2L", 
        ext_mom={ 'q1': q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6': q6 }, 
        mass_map={}, # no masses 
        loop_momenta_names=('p4', 'p6'), 
        analytic_result = 0.
    ),
    entry_name = 'PRL_6p_2L'
)



# double box
doublebox = TopologyGenerator([
        ('q1', 101, 1), ('q2', 102, 2), ('q3', 103, 3), ('q4', 104, 4),
        ('p1', 1, 6), ('p2', 6, 7), ('p3', 7, 2), ('p4', 2, 1),
        ('p5', 7, 3), ('p6', 3, 4), ('p7', 4, 6),
])
q1 = vectors.LorentzVector([  1.2,  2.2,   1.0, 0.4 ])
q2 = vectors.LorentzVector([  2.0, -5.2,   2.1, 0.0 ])
q3 = vectors.LorentzVector([ -1.6,  0.1, -12.5, 2.4 ])
q4 = -q1-q2-q3
hard_coded_topology_collection.add_topology(doublebox.create_loop_topology(
        "DoubleBox_no_ellipse", 
        ext_mom={ 'q1': q1, 'q2': q2 , 'q3': q3, 'q4': q4 }, 
        mass_map={}, # no masses 
        loop_momenta_names=('p4', 'p6'), 
        analytic_result = analytic_four_point_ladder( 
            q1.square(), q2.square(), q3.square(), q4.square(),
            (q1+q2).square(), (q2+q3).square(), 2)
    ),
    entry_name = 'DoubleBox_no_ellipse'
)


# double box
doublebox = TopologyGenerator([
        ('q1', 101, 1), ('q2', 102, 2), ('q3', 103, 3), ('q4', 104, 4),
        ('p1', 6, 1), ('p2', 1, 2), ('p3', 2, 7), ('p4', 7, 6),
        ('p5', 6, 3), ('p6', 3, 4), ('p7', 4, 7),
])
q1 = vectors.LorentzVector([  1,  0., 0., 0.3])
q2 = vectors.LorentzVector([  1, 0.,   0., -0.3 ])
q3 = vectors.LorentzVector([ -1.3,  0., 0., 0.2 ])
q4 = -q1-q2-q3
hard_coded_topology_collection.add_topology(doublebox.create_loop_topology(
        "DoubleBox_physical", 
        ext_mom={ 'q1': q1, 'q2': q2 , 'q3': q3, 'q4': q4 }, 
        mass_map={}, # no masses 
        loop_momenta_names=('p5', 'p1'), 
        analytic_result = analytic_four_point_ladder( 
            q1.square(), q2.square(), q3.square(), q4.square(),
            (q1+q2).square(), (q2+q3).square(), 2),
    #    fixed_deformation = [{'deformation_sources': [[0., 0., 0., 0.],[0.,0.,0.,0.]], 'weight_per_source': [1,1], 'excluded_loop_lines': [0, 1, 2], 'excluded_surface_ids': []},]
        
    ),
    entry_name = 'DoubleBox_physical'
)




# Fishnets generation. This is typically pretty slow and thus disabled by default.

if _GENERATE_FISHNETS:
    # FISHNET_1x2
    print("DOING 1x2...")
    FISHNET_1x2 = TopologyGenerator([
        ('q1', 101, 1), ('q2', 102, 4), ('q3', 103, 6), ('q4', 104, 3),
        ('p1', 1, 2), ('p2', 2, 3),
        ('p3', 1, 4), ('p4', 2, 5), ('p5', 3, 6),
         ('p6', 4, 5), ('p7', 5, 6)
    ])
    q1 = vectors.LorentzVector([  1.2,  2.2,   1.0, 0.4 ])
    q2 = vectors.LorentzVector([  2.0, -5.2,   2.1, 0.0 ])
    q3 = vectors.LorentzVector([ -1.6,  0.1, -12.5, 2.4 ])
    q4 = -q1-q2-q3
    hard_coded_topology_collection.add_topology(FISHNET_1x2.create_loop_topology(
        "FISHNET_1x2", 
        ext_mom={ 'q1': q1, 'q2': q2, 'q3': q3, 'q4': q4 }, 
        mass_map={}, # no masses 
        loop_momenta_names=('p1','p2'), 
        analytic_result = analytic_four_point_ladder( 
            q1.square(), q2.square(), q3.square(), q4.square(),
            (q1+q2).square(), (q2+q3).square(), 2)
    ),
    entry_name = 'FISHNET_1x2'
    )
    print("DONE")

    # FISHNET_1x4
    print("DOING 1x4...")
    FISHNET_1x4 = TopologyGenerator([
            ('q1', 101, 1), ('q2', 102, 6), ('q3', 103, 10), ('q4', 104, 5),
            ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4), ('p4', 4, 5),
            ('p5', 1, 6), ('p6', 2, 7), ('p7', 3, 8), ('p8', 4, 9), ('p9', 5, 10),
            ('p10', 6, 7), ('p11', 7, 8), ('p12', 8, 9), ('p13', 9, 10)
    ])
    q1 = vectors.LorentzVector([  1.2,  2.2,   1.0, 0.4 ])
    q2 = vectors.LorentzVector([  2.0, -5.2,   2.1, 0.0 ])
    q3 = vectors.LorentzVector([ -1.6,  0.1, -12.5, 2.4 ])
    q4 = -q1-q2-q3
    hard_coded_topology_collection.add_topology(FISHNET_1x4.create_loop_topology(
            "FISHNET_1x4", 
            ext_mom={ 'q1': q1, 'q2': q2, 'q3': q3, 'q4': q4 }, 
            mass_map={}, # no masses 
            loop_momenta_names=('p1','p2','p3','p4'), 
            analytic_result = analytic_four_point_ladder( 
                q1.square(), q2.square(), q3.square(), q4.square(),
                (q1+q2).square(), (q2+q3).square(), 4)
        ),
        entry_name = 'FISHNET_1x4'
    )
    print("DONE")

    # FISHNET_1x5
    print("DOING 1x5...")
    FISHNET_1x5 = TopologyGenerator([
            ('q1', 101, 1), ('q2', 102, 7), ('q3', 103, 12), ('q4', 104, 6),
            ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4), ('p4', 4, 5), ('p5', 5, 6),
            ('p6', 1, 7), ('p7', 2, 8), ('p8', 3, 9), ('p9', 4, 10), ('p10', 5, 11), ('p11', 6, 12),
            ('p12', 7, 8), ('p13', 8, 9), ('p14', 9, 10), ('p15', 10, 11), ('p16', 11, 12)
    ])
    q1 = vectors.LorentzVector([  1.2,  2.2,   1.0, 0.4 ])
    q2 = vectors.LorentzVector([  2.0, -5.2,   2.1, 0.0 ])
    q3 = vectors.LorentzVector([ -1.6,  0.1, -12.5, 2.4 ])
    q4 = -q1-q2-q3
    hard_coded_topology_collection.add_topology(FISHNET_1x5.create_loop_topology(
            "FISHNET_1x5", 
            ext_mom={ 'q1': q1, 'q2': q2, 'q3': q3, 'q4': q4 }, 
            mass_map={}, # no masses 
            loop_momenta_names=('p1','p2','p3','p4','p5'), 
            analytic_result = analytic_four_point_ladder( 
                q1.square(), q2.square(), q3.square(), q4.square(),
                (q1+q2).square(), (q2+q3).square(), 5)
        ),
        entry_name = 'FISHNET_1x5'
    )
    print("DONE")

    # FISHNET_1x6
    print("DOING 1x6...")
    FISHNET_1x6 = TopologyGenerator([
            ('q1', 101, 1), ('q2', 102, 8), ('q3', 103, 14), ('q4', 104, 7),
            ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4), ('p4', 4, 5), ('p5', 5, 6), ('p6', 6, 7),
            ('p7', 1, 8), ('p8', 2, 9), ('p9', 3, 10), ('p10', 4, 11), ('p11', 5, 12), ('p12', 6, 13), ('p13', 7, 14),
            ('p14', 8, 9), ('p15', 9, 10), ('p16', 10, 11), ('p17', 11, 12), ('p18', 12, 13), ('p19', 13, 14)
    ])
    q1 = vectors.LorentzVector([  1.2,  2.2,   1.0, 0.4 ])
    q2 = vectors.LorentzVector([  2.0, -5.2,   2.1, 0.0 ])
    q3 = vectors.LorentzVector([ -1.6,  0.1, -12.5, 2.4 ])
    q4 = -q1-q2-q3
    hard_coded_topology_collection.add_topology(FISHNET_1x6.create_loop_topology(
            "FISHNET_1x6", 
            ext_mom={ 'q1': q1, 'q2': q2, 'q3': q3, 'q4': q4 }, 
            mass_map={}, # no masses 
            loop_momenta_names=('p1','p2','p3','p4','p5','p6'), 
            analytic_result = analytic_four_point_ladder( 
                q1.square(), q2.square(), q3.square(), q4.square(),
                (q1+q2).square(), (q2+q3).square(), 6)
        ),
        entry_name = 'FISHNET_1x6'
    )
    print("DONE")

    # FISHNET_2x2
    print("DOING 2x2...")
    FISHNET_2x2 = TopologyGenerator([
            ('q1', 101, 1), ('q2', 102, 7), ('q3', 103, 9), ('q4', 104, 3),
            ('p1', 1, 2), ('p2', 2, 3), 
            ('p3', 1, 4), ('p4', 2, 5), ('p5', 3, 6),
            ('p6', 4, 5), ('p7', 5, 6), 
            ('p8', 4, 7), ('p9', 5, 8), ('p10', 6, 9),
            ('p11', 7, 8), ('p12', 8, 9),
    ])
    q1 = vectors.LorentzVector([  1.2,  2.2,   1.0, 0.4 ])
    q2 = vectors.LorentzVector([  2.0, -5.2,   2.1, 0.0 ])
    q3 = vectors.LorentzVector([ -1.6,  0.1, -12.5, 2.4 ])
    q4 = -q1-q2-q3
    hard_coded_topology_collection.add_topology(FISHNET_2x2.create_loop_topology(
            "FISHNET_2x2", 
            ext_mom={ 'q1': q1, 'q2': q2 , 'q3': q3, 'q4': q4 }, 
            mass_map={}, # no masses 
            loop_momenta_names=('p1', 'p2', 'p11', 'p12'), 
            analytic_result = 2.6918653677981387e-14
        ),
        entry_name = 'FISHNET_2x2'
    )
    print("DONE")

    # FISHNET_2x3
    FISHNET_2x3 = TopologyGenerator([
            ('q1', 101, 1), ('q2', 102, 9), ('q3', 103, 12), ('q4', 104, 4),
            ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),
            ('p4', 1, 5), ('p5', 2, 6), ('p6', 3, 7), ('p7', 4, 8),
            ('p8', 5, 6), ('p9', 6, 7), ('p10', 7, 8),
            ('p11', 5, 9), ('p12', 6, 10), ('p13', 7, 11), ('p14', 8, 12),
            ('p15', 9, 10), ('p16', 10, 11), ('p17', 11, 12)
    ])
    q1 = vectors.LorentzVector([  1.2,  2.2,   1.0, 0.4 ])
    q2 = vectors.LorentzVector([  2.0, -5.2,   2.1, 0.0 ])
    q3 = vectors.LorentzVector([ -1.6,  0.1, -12.5, 2.4 ])
    q4 = -q1-q2-q3
    print("DOING 2x3...")
    hard_coded_topology_collection.add_topology(FISHNET_2x3.create_loop_topology(
            "FISHNET_2x3", 
            ext_mom={ 'q1': q1, 'q2': q2 , 'q3': q3, 'q4': q4 }, 
            mass_map={}, # no masses 
            loop_momenta_names=('p1', 'p2', 'p3', 'p15', 'p16', 'p17'), 
            analytic_result = 8.4044862640909e-19
        ),
        entry_name = 'FISHNET_2x3'
    )
    print("Done 2x3")

# PRL two-loops
# PRL_6p_2L
PRL_6p_2L = TopologyGenerator([
        ('q1', 101, 1), ('q2', 102, 2), ('q3', 103, 3), ('q4', 104, 4),
        ('p1', 1, 6), ('p2', 6, 8), ('p3', 7, 2), ('p4', 2, 1),
        ('p5', 7, 3), ('p6', 3, 4), ('p7', 4, 6), ('p8', 8, 9), ('p9', 9, 7),
        ('q5', 108, 8), ('q6', 109, 9)
])
q1 = vectors.LorentzVector([ 0.2, 0.3, 0.5, 0.6 ])
q2 = vectors.LorentzVector([-0.1, 0.7, 0.2, 0.1])
q3 = vectors.LorentzVector([ 0.1, 0.5, -0.3, -0.4])
q4 = vectors.LorentzVector([-0.3, 0.4, 0.5, 0.2])
q5 = vectors.LorentzVector([-0.2, 0.3, 0.2, -0.5])
q6 = -q1-q2-q3-q4-q5
hard_coded_topology_collection.add_topology(PRL_6p_2L.create_loop_topology(
        "PRL_6p_2L",
        ext_mom={ 'q1': q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6': q6 }, 
        mass_map={}, # no masses 
        loop_momenta_names=('p4', 'p6'),
        analytic_result = 0.
    ),
    entry_name = 'PRL_6p_2L'
)

# PRL_6p_2L
PRL_6p_2L = TopologyGenerator([
        ('q1', 101, 1), ('q2', 102, 2), ('q3', 103, 3), ('q4', 104, 4),
        ('p1', 1, 6), ('p2', 6, 8), ('p3', 7, 2), ('p4', 2, 1),
        ('p5', 7, 3), ('p6', 3, 4), ('p7', 4, 6), ('p8', 8, 9), ('p9', 9, 7),
        ('q5', 108, 8), ('q6', 109, 9)
])
q1 = vectors.LorentzVector([ 0.2, 0.3, 0.5, 0.6 ])
q2 = vectors.LorentzVector([-0.1, 0.7, 0.2, 0.1])
q3 = vectors.LorentzVector([ 0.1, 0.5, -0.3, -0.4])
q4 = vectors.LorentzVector([-0.3, 0.4, 0.5, 0.2])
q5 = vectors.LorentzVector([-0.2, 0.3, 0.2, -0.5])
q6 = -q1-q2-q3-q4-q5
hard_coded_topology_collection.add_topology(PRL_6p_2L.create_loop_topology(
        "PRL_6p_2L_massive",
        ext_mom={ 'q1': q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6': q6 }, 
        mass_map={ 'p1':1., 'p2':1., 'p3':1., 'p4':1., 'p5':1., 'p6':1.,
                   'p7':1., 'p8':1., 'p9':1. }, # no masses 
        loop_momenta_names=('p4', 'p6'),
        analytic_result = 0.
    ),
    entry_name = 'PRL_6p_2L_massive'
)

# PRL_8p_2L
PRL_8p_2L = TopologyGenerator([
        ('q1', 101, 1), ('q2', 102, 2), ('q3', 103, 3), ('q4', 104, 4),
        ('p1', 1, 6), ('p2', 6, 8), ('p3', 7, 2), ('p4', 2, 1),
        ('p5', 7, 3), ('p6', 3, 4), ('p7', 4, 6), ('p8', 8, 9), ('p9', 9, 10), ('p10', 10, 11), ('p11', 11, 7),
        ('q5', 108, 8), ('q6', 109, 9), ('q7', 110, 10), ('q8', 111, 11),
])
q1 = vectors.LorentzVector([ 0.15, 0.09, 0.23, 0. ])
q2 = vectors.LorentzVector([-0.2, 0.79, 0.12, 0.11])
q3 = vectors.LorentzVector([-0.23, 0.14, -0.47, -0.22])
q4 = vectors.LorentzVector([0.11, -0.59, 0.54, 0.12])
q5 = vectors.LorentzVector([-0.15, 0.21, 0.10, -0.32])
q6 = vectors.LorentzVector([0.32,0.84,0.27,0.49])
q7 = vectors.LorentzVector([0.11,-0.3,-0.12,-0.1])
q8 = -q1-q2-q3-q4-q5-q6-q7
hard_coded_topology_collection.add_topology(PRL_8p_2L.create_loop_topology(
        "PRL_8p_2L_massive",
        ext_mom={ 'q1': q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6': q6, 'q7': q7, 'q8': q8},
        mass_map={ 'p1':1., 'p2':1., 'p3':1., 'p4':1., 'p5':1., 'p6':1.,
                   'p7':1., 'p8':1., 'p9':1., 'p10':1., 'p11':1. }, # no masses 
        loop_momenta_names=('p4', 'p6'),
        analytic_result = 0.
    ),
    entry_name = 'PRL_8p_2L_massive'
)

# mercedes three loop 6p
PRL_mercedes_6p = TopologyGenerator([
    ('q1', 0, 5), ('p0', 1, 13), ('p1', 5, 2), ('p2', 10, 3), ('p10',13,5),
    ('p3', 3, 1), ('p4', 1, 8), ('p5', 2, 7), ('p6', 3, 4), ('p7', 2 ,10), ('p8', 7,4), ('p9',8,4),
    ('q2', 6, 3), ('q3',110,10), ('q4',107,7), ('q5',108,8), ('q6', 113,13),
])
q1 = vectors.LorentzVector([ 0.2, 0.3, 0.5, 0.6 ])
q2 = vectors.LorentzVector([-0.1, 0.7, 0.2, 0.1])
q3 = vectors.LorentzVector([ 0.1, 0.5, -0.3, -0.4])
q4 = vectors.LorentzVector([-0.3, 0.4, 0.5, 0.2])
q5 = vectors.LorentzVector([-0.2, 0.3, 0.2, -0.5])
q6 = -q1-q2-q3-q4-q5
hard_coded_topology_collection.add_topology(PRL_mercedes_6p.create_loop_topology(
    'PRL_Mercedes_6p_massive',
    ext_mom={'q1': q1, 'q2': q2, 'q3': q3, 'q4': q4, 'q5': q5, 'q6': q6},
    mass_map={ 'p1':1., 'p2':1., 'p3':1., 'p4':1., 'p5':1., 'p6':1.,
               'p7':1., 'p8':1., 'p9':1., 'p10':1., 'p0':1. }, # no masses
    loop_momenta_names=('p3', 'p4', 'p6',),
    analytic_result = 0.,
    contour_closure = [0,0,1],
    ),
    entry_name = 'PRL_Mercedes_6p_massive'
)

# non planar four loop
PRL_non_planar_four_loop = TopologyGenerator([
    ('q', 0, 1), ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),
    ('p4', 4, 1), ('p5', 5, 1), ('p6', 4, 6), ('p7', 5, 6), 
    ('p8', 6, 2),('p9', 5, 3), ('-q', 7, 3)
])
q = vectors.LorentzVector([ 5., 3., 4., 1.])
hard_coded_topology_collection.add_topology(PRL_non_planar_four_loop.create_loop_topology(
    'PRL_non_planar_four_loop',
    ext_mom={'q': q, '-q': -q},
    loop_momenta_names=('p4', 'p5', 'p6', 'p7',),
    analytic_result = 36. * zeta(3)**2 /(16.*math.pi**2)**4 / q.square(),
    ),
    entry_name = 'PRL_non_planar_four_loop'
)

# PRL physical BoxBox
PRL_physical_BoxBox = TopologyGenerator([
        ('q1', 101, 1), ('q2', 102, 2), ('q3', 103, 3), ('q4', 104, 4),
        ('p1', 1, 6), ('p2', 6, 7), ('p3', 7, 2), ('p4', 2, 1),
        ('p5', 7, 3), ('p6', 3, 4), ('p7', 4, 6)
    ])
q1 = vectors.LorentzVector([ 1.5, 0.2, 0.1, -0.4 ])
q2 = vectors.LorentzVector([ 0.4, -0.3, 0.1, 0.2 ])
q3 = vectors.LorentzVector([ -0.7, 0.2, 0.2, 0.3 ])
q4 = -q3-q2-q1
hard_coded_topology_collection.add_topology(PRL_physical_BoxBox.create_loop_topology(
    'PRL_physical_BoxBox',
    ext_mom={'q1': q1, 'q2': q2, 'q3': q3, 'q4': q4},
    mass_map={},
    loop_momenta_names=('p4','p2'),
    analytic_result = analytic_four_point_ladder(
                    q1.square(), q2.square(), q3.square(), q4.square(),
                    (q1+q2).square(), (q2+q3).square(), 2),
    ),
    entry_name = 'PRL_physical_BoxBox'
)


# PRL physical BoxBox
PRL_physical_BoxBox_1ellipse = TopologyGenerator([
        ('q1', 101, 1), ('q2', 102, 2), ('q3', 103, 3), ('q4', 104, 4),
        ('p1', 1, 6), ('p2', 6, 7), ('p3', 7, 2), ('p4', 2, 1),
        ('p5', 7, 3), ('p6', 3, 4), ('p7', 4, 6)
    ])
q1 = vectors.LorentzVector([1.082,0.2891,0.5276,0.119])
q2 = vectors.LorentzVector([-0.3978,0.4261,0.1091,0.2143])
q3 = vectors.LorentzVector([0.1182,0.2192,0.5019,0.3210])
q4 = -q3-q2-q1
hard_coded_topology_collection.add_topology(PRL_physical_BoxBox_1ellipse.create_loop_topology(
    'PRL_physical_BoxBox_1ellipse',
    ext_mom={'q1': q1, 'q2': q2, 'q3': q3, 'q4': q4},
    mass_map={},
    loop_momenta_names=('p4','p2'),
    analytic_result = analytic_four_point_ladder(
                    q1.square(), q2.square(), q3.square(), q4.square(),
                    (q1+q2).square(), (q2+q3).square(), 2),
    ),
    entry_name = 'PRL_physical_BoxBox_1ellipse'
)

# PRL physical BoxBox few ellipse
PRL_physical_BoxBox_fewellipse = TopologyGenerator([
        ('q1', 101, 1), ('q2', 102, 2), ('q3', 103, 3), ('q4', 104, 4),
        ('p1', 1, 6), ('p2', 6, 7), ('p3', 7, 2), ('p4', 2, 1),
        ('p5', 7, 3), ('p6', 3, 4), ('p7', 4, 6)
    ])
q1 = vectors.LorentzVector([0.812,-0.2891,0.4372,-0.012])
q2 = vectors.LorentzVector([0.609,0.4333,-0.291,0.2813])
q3 = vectors.LorentzVector([-0.7182,0.1229,1.5019,0.3640])
q4 = -q3-q2-q1
hard_coded_topology_collection.add_topology(PRL_physical_BoxBox_fewellipse.create_loop_topology(
    'PRL_physical_BoxBox_fewellipse',
    ext_mom={'q1': q1, 'q2': q2, 'q3': q3, 'q4': q4},
    mass_map={},
    loop_momenta_names=('p4','p2'),
    analytic_result = analytic_four_point_ladder(
                    q1.square(), q2.square(), q3.square(), q4.square(),
                    (q1+q2).square(), (q2+q3).square(), 2),
    ),
    entry_name = 'PRL_physical_BoxBox_fewellipse'
)

# PRL BoxBoxBoxBox
PRL_BoxBoxBoxBox = TopologyGenerator([
        ('q1', 101, 1), ('q2', 102, 2), ('q3', 103, 3), ('q4', 104, 4),
        ('p1', 1, 6), ('p2', 6, 7), ('p3', 7, 2), ('p4', 2, 1),
        ('p5', 7, 9), ('p6', 9, 8), ('p7', 8, 6),
        ('p8', 8,10), ('p9',10,11), ('p10',11, 9),
        ('p11',11, 3), ('p12', 3, 4), ('p13', 4, 10)
])
q1 = vectors.LorentzVector([ 0.1, 0.2, 0.5, 0.1 ])
q2 = vectors.LorentzVector([-0.3, 0.4, 0.1, 0.2 ])
q3 = vectors.LorentzVector([ 0.1, 0.2, 0.5, 0.3 ])
q4 = -q3-q2-q1
hard_coded_topology_collection.add_topology(PRL_BoxBoxBoxBox.create_loop_topology(
    'PRL_BoxBoxBoxBox',
    ext_mom={'q1': q1, 'q2': q2, 'q3': q3, 'q4': q4},
    mass_map={},
    loop_momenta_names=('p4','p2','p6','p9'),
    analytic_result = analytic_four_point_ladder(
                    q1.square(), q2.square(), q3.square(), q4.square(),
                    (q1+q2).square(), (q2+q3).square(), 4),
    ),
    entry_name = 'PRL_BoxBoxBoxBox'
)

# mercedes three loop
PRL_mercedes = TopologyGenerator([
    ('q', 0, 5), ('p0', 1, 5), ('p1', 5, 2), ('p2', 2, 3),
    ('p3', 3, 1), ('p4', 1, 4), ('p5', 2, 4), ('p6', 3, 4),
    ('-q', 6, 3)
])
q = vectors.LorentzVector([ 5., 3., 4., 1.])
hard_coded_topology_collection.add_topology(PRL_mercedes.create_loop_topology(
    'PRL_Mercedes',
    ext_mom={'q': q, '-q': -q},
    loop_momenta_names=('p3', 'p4', 'p6',),
    analytic_result = 20. * zeta(5) * -1j/(16.*math.pi**2)**3 / q.square(),
    contour_closure = [0,0,1],
    ),
    entry_name = 'PRL_Mercedes'
)

# triangle box triangle
TriangleBoxTriangle = TopologyGenerator([
    ('q', 0, 1), ('p0', 1, 2), ('p1', 2, 3), ('p2', 3, 4),
    ('p3', 4, 5), ('p4', 5, 6), ('p5', 6, 1), ('p6', 6, 2),
    ('p7', 5, 3), ('-q', 7, 4)
])
q = vectors.LorentzVector([ 1., 1., 1., 0.])
hard_coded_topology_collection.add_topology(TriangleBoxTriangle.create_loop_topology(
    'TriangleBoxTriangle',
    ext_mom={'q': q, '-q': -q},
    loop_momenta_names=('p5', 'p6', 'p7',),
    analytic_result = analytic_two_point_ladder(q.square(),3),
    contour_closure = [0,0,1],
    ),
    entry_name = 'TriangleBoxTriangle'
)

# triangle box triangle
TriangleBoxTriangle = TopologyGenerator([
    ('q', 0, 1), ('p0', 1, 2), ('p1', 2, 3), ('p2', 3, 4),
    ('p3', 4, 5), ('p4', 5, 6), ('p5', 6, 1), ('p6', 6, 2),
    ('p7', 5, 3), ('-q', 7, 4)
])
q = vectors.LorentzVector([ 1., 0., 0., 0.])
hard_coded_topology_collection.add_topology(TriangleBoxTriangle.create_loop_topology(
    'TriangleBoxTriangle_physical',
    ext_mom={'q': q, '-q': -q},
    loop_momenta_names=('p5', 'p6', 'p7',),
    analytic_result = analytic_two_point_ladder(q.square(),3),
    contour_closure = [0,0,1],
    #fixed_deformation = [{'deformation_sources': [[0.,0.,0.,0.],[0.,0.,0.,0.],[0.,0.,0.,0.]], 'weight_per_source': [1,1], 'excluded_loop_lines': [0,1,2,3,4,5], 'excluded_surface_ids':[]}],
    fixed_deformation = [], # this takes a long time to generate
    ),
    entry_name = 'TriangleBoxTriangle_physical'
)

two_loop_6pt = TopologyGenerator([
    ('q1', 101, 1), ('q2', 102, 2), ('q3', 103, 3), ('q4', 104, 4), 
    ('q5', 105, 5), ('q6', 106, 6), ('p1', 1, 2), ('p2', 2, 3),
    ('p3', 3, 4), ('p4', 4, 5), ('p5', 5, 6), ('p6', 6, 1),
    ('p7', 6, 7), ('p8', 7, 8), ('p9', 8, 2)
])
q1 = vectors.LorentzVector([ 1., 3., 0., 0.])
q2 = vectors.LorentzVector([-1., 0., 2., 0.])
q3 = vectors.LorentzVector([ 2.,-4., 0., 1.])
q4 = vectors.LorentzVector([-3., 5., 1., 2.])
q5 = vectors.LorentzVector([ 2.,-2., 0.,-2.])
q6 = vectors.LorentzVector([-1.,-2.,-3.,-1.])
hard_coded_topology_collection.add_topology(two_loop_6pt.create_loop_topology(
    'two_loop_6pt_no_ellipse',
    ext_mom={'q1': q1, 'q2': q2, 'q3': q3, 'q4': q4, 'q5': q5, 'q6': q6},
    loop_momenta_names=('p1', 'p3',),
    analytic_result = None,
    contour_closure = [0,0],
    ),
    entry_name = 'two_loop_6pt_no_ellipse'
)

box = TopologyGenerator([
    ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),  ('p4', 4, 1),
    ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4)
])
q1 = vectors.LorentzVector([14.,-6.6,-40.,0.])
q2 = vectors.LorentzVector([-43., 12., 33., 0.])
q3 = vectors.LorentzVector([-28., -50., 10., 0.])
hard_coded_topology_collection.add_topology(box.create_loop_topology(
        "Box_RodrigoFail2", 
        ext_mom={ 'q1': q1, 'q2': q2 , 'q3': q3, 'q4': -q1-q2-q3 }, 
        mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0}, 
        loop_momenta_names=('p4',), # If not specified an arbitrary spanning tree will be used for momentum routing 
        analytic_result=None, # For triangle and box one-loop topology, the analytic result is automatically computed
        # For now specified by hand as the cvxpy automated implementation is not done yet
        #fixed_deformation =  [{'deformation_sources': [[0., 6.2, 11.8, 0.]], 'weight_per_source': [1], 'excluded_loop_lines': [0], 'excluded_surface_ids': []}]
     ),
     entry_name = 'Box_RodrigoFail2'
)


# Example printout
# ----------------
#hard_coded_topology_collection['non_planar_four_loop_no_ellipse'].print_topology()
#hard_coded_topology_collection['TriangleBox'].print_topology()
#hard_coded_topology_collection['DoubleBox'].print_topology()
#hard_coded_topology_collection['TriangleBoxBox'].print_topology()

# Example of yaml export and import
# ---------------------------------
#hard_coded_topology_collection['DoubleTriangle'].export_to('ExampleDoubleTriangleExport.yaml')
#test = LoopTopology.import_from('ExampleDoubleTriangleExport.yaml')
#test.print_topology()

# Example of a yaml export and import of a complete TopologyCollection
# --------------------------------------------------------------------
#hard_coded_topology_collection.export_to('ExampleTopologyCollectionExport.yaml')
#test = TopologyCollection.import_from('ExampleTopologyCollectionExport.yaml')
#test['DoubleTriange'].print_topology()


#print(hard_coded_topology_collection['manual_Box_1_cutgroup'].export_to('',format='mathematica'))
#print(hard_coded_topology_collection['manual_Box_no_ellipse'].export_to('',format='mathematica'))
#print(hard_coded_topology_collection['Hexagon_10E_7s'].export_to('',format='mathematica'))
