import sys
import os
from pathlib import Path

if True:
    pjoin = os.path.join
    root_path = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, pjoin(root_path, os.path.pardir))
    sys.path.insert(0, pjoin(root_path, os.path.pardir, os.path.pardir))
    from LTD.squared_topologies import SquaredTopologyGenerator

ST_FORCER_2P_3L = []
ST_FORCER_2P_3L.append([(0, 3), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4), (4, 5)])   # + 2*ep^-1*z3
ST_FORCER_2P_3L.append([(0, 1), (1, 2), (1, 5), (2, 3), (2, 4), (3, 4), (3, 5), (4, 5), (5, 6)])   # - 2*ep^-1*z3
ST_FORCER_2P_3L.append([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 5), (4, 6), (5, 6), (2, 7)])   # - 2*ep^-1*z3

ST_FORCER_2P_4L = []
ST_FORCER_2P_4L.append([(0, 3), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5), (3, 5), (4, 5), (4, 6)]) # +5*ep^-1*z5
ST_FORCER_2P_4L.append([(0, 3), (1, 2), (1, 3), (1, 5), (2, 4), (2, 5), (3, 4), (3, 5), (4, 5), (4, 6)]) # +5*ep^-1*z5
ST_FORCER_2P_4L.append([(0, 5), (1, 2), (1, 3), (1, 5), (2, 4), (2, 6), (3, 4), (3, 6), (4, 5), (5, 6), (6, 7)])  # -5*ep^-1*z5
ST_FORCER_2P_4L.append([(0, 5), (1, 2), (1, 3), (1, 5), (2, 4), (2, 5), (3, 4), (3, 6), (4, 6), (5, 6), (6, 7)])  # -5*ep^-1*z5
ST_FORCER_2P_4L.append([(0, 1), (1, 5), (1, 6), (2, 3), (2, 4), (2, 6), (3, 5), (3, 6), (4, 5), (4, 6), (5, 7)])  # -5*ep^-1*z5
ST_FORCER_2P_4L.append([(0, 1), (1, 2), (1, 5), (2, 3), (2, 6), (3, 4), (3, 6), (4, 5), (4, 6), (5, 6), (5, 7)])  # -5*ep^-1*z5
ST_FORCER_2P_4L.append([(0, 1), (1, 2), (1, 7), (2, 3), (2, 4), (3, 5), (3, 6), (4, 5), (4, 7), (5, 6), (6, 7), (7, 8)])  # +5/4*ep^-1*z3
ST_FORCER_2P_4L.append([(0, 1), (1, 2), (1, 3), (2, 4), (2, 5), (3, 6), (3, 7), (4, 6), (4, 7), (5, 6), (5, 7), (7, 8)])  # -5*ep^-1*z5
ST_FORCER_2P_4L.append([(0, 1), (1, 2), (1, 3), (2, 4), (2, 5), (3, 4), (3, 7), (4, 6), (5, 6), (5, 7), (6, 7), (7, 8)])  # -5*ep^-1*z5
ST_FORCER_2P_4L.append([(0, 1), (1, 3), (1, 7), (2, 4), (2, 7), (3, 5), (3, 6), (4, 5), (4, 6), (5, 7), (6, 7), (2, 8)])  # -5*ep^-1*z5
ST_FORCER_2P_4L.append([(0, 1), (1, 3), (1, 7), (2, 4), (2, 7), (3, 4), (3, 5), (4, 6), (5, 6), (5, 7), (6, 7), (2, 8)])  # -5*ep^-1*z5
ST_FORCER_2P_4L.append([(0, 1), (1, 3), (1, 4), (2, 3), (2, 7), (3, 5), (4, 6), (4, 7), (5, 6), (5, 7), (6, 7), (2, 8)])  # -5*ep^-1*z5
ST_FORCER_2P_4L.append([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 7), (4, 6), (4, 7), (5, 6), (5, 7), (6, 7), (2, 8)])  # -5*ep^-1*z5
ST_FORCER_2P_4L.append([(0, 1), (1, 2), (1, 3), (2, 7), (3, 4), (3, 7), (4, 5), (4, 6), (5, 6), (5, 7), (6, 7), (2, 8)])  # +5/4*ep^-1*z3
ST_FORCER_2P_4L.append([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 7), (4, 8), (5, 8), (6, 7), (6, 8), (2, 9)]) # -5*ep^-1*z5
ST_FORCER_2P_4L.append([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 6), (4, 8), (5, 8), (6, 7), (7, 8), (2, 9)]) # -10*ep^-1*z5
ST_FORCER_2P_4L.append([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 6), (4, 7), (5, 8), (6, 8), (7, 8), (2, 9)]) # -10*ep^-1*z5
ST_FORCER_2P_4L.append([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 7), (4, 8), (5, 7), (5, 8), (6, 7), (6, 8), (2, 9)]) # -5/2*ep^-1*z5
ST_FORCER_2P_4L.append([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 6), (4, 7), (5, 7), (5, 8), (6, 8), (7, 8), (2, 9)]) # -5/2*ep^-1*z5
ST_FORCER_2P_4L.append([(0, 1), (1, 2), (1, 3), (2, 4), (3, 5), (3, 6), (4, 5), (4, 7), (5, 8), (6, 7), (6, 8), (7, 8), (2, 9)]) # +5/12*ep^-1*z3

if __name__ == '__main__':
    for n, edges_raw in enumerate(ST_FORCER_2P_3L):
        n_loops = 3
        name = "STF_3L_%d"%n
        edges = [('p%d'%i , *e) for i, e in enumerate(edges_raw)]
        edges[0] = ('q1',*edges[0][1:])
        edges[-1] = ('q2',*edges[-1][1:])
        incoming_momenta_names = ['q1']
        n_jets=0
        print(edges)
        external_momenta = {'q1': [1.0, 0.0, 0.0 ,0.0], 'q2':[1.0, 0.0, 0.0, 0.0]}
        stf_3L = SquaredTopologyGenerator(edges, 
                                          name, 
                                          incoming_momenta_names,
                                          n_jets,
                                          external_momenta,
                                          overall_numerator=1.0,
        )
        stf_3L.export('%s.yaml'%name)
    for n, edges_raw in enumerate(ST_FORCER_2P_4L):
        n_loops = 4
        name = "STF_4L_%d"%n
        edges = [('p%d'%i , *e) for i, e in enumerate(edges_raw)]
        edges[0] = ('q1',*edges[0][1:])
        edges[-1] = ('q2',*edges[-1][1:])
        incoming_momenta_names = ['q1']
        n_jets=0
        print(edges)
        external_momenta = {'q1': [1.0, 0.0, 0.0 ,0.0], 'q2':[1.0, 0.0, 0.0, 0.0]}
        stf_4L = SquaredTopologyGenerator(edges, 
                                          name, 
                                          incoming_momenta_names,
                                          n_jets,
                                          external_momenta,
                                          overall_numerator=1.0,
        )
        stf_4L.export('%s.yaml'%name)
