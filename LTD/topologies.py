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
        analytic_result= (-(6.*zeta(3))/((16.*(math.pi**2))**2*q.square()))
    ),
    entry_name = 'DoubleTriangle_no_ellipse'
)

# double box
doublebox = TopologyGenerator([
    ('q1', 101, 1), ('q1', 102, 2), ('q1', 103, 3), ('q1', 104, 4),
    ('p1', 1, 6), ('p2', 6, 7), ('p3', 7, 2), ('p4', 2, 1),
    ('p5', 6, 3), ('p6', 3, 4), ('p7', 4, 7)
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

# Example printout
# ----------------
#hard_coded_topology_collection['DoubleTriangle'].print_topology()
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

