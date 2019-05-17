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
                                ('p3',3, 2), ('q2', 4, 3), ('q3', 5, 2)])
q1 = vectors.LorentzVector([0.1, 0.2, 0.3, 0.4])
q2 = vectors.LorentzVector([0.2, 0.5, 0.6, 0.7])
#hard_coded_topology_collection.add_topology(triangle.create_loop_topology(
#        "Triangle_no_ellipse", 
#        ext_mom={ 'q1': q1, 'q2': q2 }, 
#        mass_map={'p1': 1.0, 'p2': 2.0, 'p3': 3.0}, 
#        loop_momenta_names=('p1',), 
#        analytic_result=None
#    ),
#    entry_name = 'Triangle_no_ellipse'
#)

# double-triangle
doubletriangle = TopologyGenerator([
    ('q', 0, 1), ('p1', 1, 2), ('p2', 3, 1), ('p3', 2, 3),
    ('p4', 3, 4), ('p5', 4, 2), ('-q', 4, 5)])
q = vectors.LorentzVector([0.1, 0.2, 0.3, 0.4])
hard_coded_topology_collection.add_topology(doubletriangle.create_loop_topology(
        "DoubleTriangle_no_ellipse", 
        ext_mom={'q': q }, 
        mass_map={'p1': 1.0, 'p2': 2.0, 'p3': 3.0}, 
        loop_momenta_names=('p1', 'p5'), 
        analytic_result= (-(6.*zeta(3))/((16.*(math.pi**2))**2*q.square()))
    ),
    entry_name = 'DoubleTriangle_no_ellipse'
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

