#!/usr/bin/env python2
import ltd

l1 = ltd.LTD("LTD/topologies.yaml", "TriangleBox", "LTD/hyperparameters.yaml")

# evaluate with unit hypercube input
print l1.evaluate([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])

# deform a loop momentum. The result is kappa and the jacobian
print l1.deform(loop_momenta=[[2.0, 3.0, 4.0], [5.0, 6.0, 7.0]])

# evaluate a cut given a complex loop momentum
print l1.evaluate_cut(loop_momenta=[[(0.2,0.1), (0.3,0.1), (0.4,0.1)],[(0.2,0.1), (0.3,0.1), (0.4,0.1)]], 
    cut_structure_index=0, cut_index=1)