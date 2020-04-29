#!/usr/bin/env python3
from ltd_utils import TopologyGenerator
import copy
import math
from itertools import combinations_with_replacement
import vectors
from squared_topologies import SquaredTopologyGenerator

if __name__ == "__main__":

        ttHNNLO_4Box = SquaredTopologyGenerator(
            [
                ('q1', 101, 1), ('q2', 102, 1), ('q3', 2, 103), ('q4', 2, 104), 
                ('p5', 1, 3), ('p6', 4, 2), ('p7', 3, 5), ('p8', 6, 3), ('p9',7, 4), 
                ('p10', 4, 8), ('p11', 7, 5), ('p12', 5, 9), ('p13',8, 6), ('p14', 9, 6),
                ('p15', 10, 7), ('p16', 8, 10), ('p17', 9, 10)
            ],
        "ttHNNLO_4Box", ['q1', 'q2'], 
        3, 
        {'q1': [1000., 0., 0., 1000.], 'q2': [1000., 0., 0., -1000.], 'q3': [1000., 0., 0., 1000.], 'q4': [1000., 0., 0., -1000.]},
        # This is the basis generated from QGRAF
        loop_momenta_names=(
            'p7', # k1
            'p9', # -k2
            'p11', # -k3, Higgs
            'p13' # -k4
        ),
        # This is the basis that puts the gluon on loop momenta:
#        loop_momenta_names=(
#            'p7', # k1
#            'p14', # k1+k4-p1-p2
#            'p11', # -k3
#            'p16' # -k2+k4-p1-p2
#        ),
        final_state_particle_ids=(),#(6,6,25),
        particle_ids={
            'q1' : -11,
            'q2' : 11,
            'q3' : -11,
            'q4' : 11,
            'p5' : 22,
            'p6' : 22,
            'p7' : 6,
            'p8' : 6,
            'p9' : 6,
            'p10' : 6,
            'p11' : 25,
            'p12' : 6,
            'p13' : 6,
            'p14' : 21,
            'p15' : 6,
            'p16' : 21,
            'p17' : 6,
        },
        overall_numerator=1.0
#        ,numerator_structure={('p2', 'p5', 'p7'):
#            { (): # uv structure
#            [
#            [[0,4],[0.,+2.954161761482786e8]],
#            [[0,4,4],[0.,+1.477080880741393e8]],
#            [[3,3,7,7],[0.,-7.385404403706966e7]]
#            ]
#            },
#            ('p3', 'p6', 'p7'):
#            { (): # uv structure
#            [
#            [[0,4],[0.,+2.954161761482786e8]],
#            [[0,4,4],[0.,+1.477080880741393e8]],            
#            ]
#            },
#            }
        )
        ttHNNLO_4Box.export('ttHNNLO_4Box.yaml')

