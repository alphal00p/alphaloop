#!/usr/bin/env python3
from ltd_utils import TopologyGenerator
import copy
import math
from itertools import combinations_with_replacement
import vectors
from squared_topologies import SquaredTopologyGenerator
print("Now loading coefficients")
from ttHNNLO_4Box_coefficients_from_mathematica import ttHNNLO_4Box_coeffs
print("Coefficient loading done.")

# Swap real and imaginary part of all coefficients and multiply by minus, 
# emulating a multiplication by i which is in the overall numerator in mathematica.
for cut in ttHNNLO_4Box_coeffs:
    ttHNNLO_4Box_coeffs[cut] =[[[],[1.0,0.0]],]
#    for i_coef, coef in enumerate(ttHNNLO_4Box_coeffs[cut]):
#        ttHNNLO_4Box_coeffs[cut][i_coef][1] = [ -coef[1][1], coef[1][0] ]

# N[(2/3)^2 (I 1/3 ge^4 gs^4 yt^2) /. {yt -> 0.99366614581500623, 
#    gs -> Sqrt[0.118 4 \[Pi]], ge -> Sqrt[1/132.507 4 \[Pi]]}, 
#  16] // FullForm

# And divide by two in order to compensate for the incoorrect symmetry
# factor of 4 (it should only be 2 due to the presence of two identical
# gluons in the final state). This incorrect factor comes from the fact
# that my implementation below does not differentiate top and antitop.

overall_numerator = 0.0028926976458677253 / 2.

if __name__ == "__main__":

        ttHNNLO_4Box = SquaredTopologyGenerator(
            [
                ('q1', 101, 1), ('q2', 102, 1), ('q3', 2, 103), ('q4', 2, 104), 
                ('p5', 1, 3), ('p6', 4, 2), ('p7', 3, 5), ('p8', 6, 3), ('p9',7, 4), 
                ('p10', 4, 8), ('p11', 7, 5), ('p12', 5, 9), ('p13',8, 6), ('p14', 9, 6),
                ('p15', 10, 7), ('p16', 8, 10), ('p17', 9, 10)
            ],
        "ttHNNLO_4Box_NNLO", ['q1', 'q2'], 
        0, 
        {'q1': [1000., 0., 0., 1000.], 'q2': [1000., 0., 0., -1000.], 'q3': [1000., 0., 0., 1000.], 'q4': [1000., 0., 0., -1000.]},
        # This is the basis generated from QGRAF
#        loop_momenta_names=(
#            'p7', # k1
#            'p9', # -k2
#            'p11', # -k3, Higgs
#            'p13' # -k4
#        ),
        # This is the basis that puts the gluon on loop momenta:
        loop_momenta_names=(
            'p7', # k1
            'p14', # k1+k4-p1-p2
            'p11', # -k3
            'p16' # -k2+k4-p1-p2
        ),
        final_state_particle_ids=(6,6,25),
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
        overall_numerator=overall_numerator
        ,numerator_structure={
            ('p10', 'p11', 'p12', 'p14', 'p16'):
            {
                (): ttHNNLO_4Box_coeffs[(10,11,12,14,16)]
            },
            ('p10', 'p11', 'p15'):
            {
                (): ttHNNLO_4Box_coeffs[(10,11,15)]
            },
            ('p10', 'p11', 'p16', 'p17'):
            {
                (): ttHNNLO_4Box_coeffs[(10,11,16,17)]
            },
            ('p11', 'p12', 'p13', 'p14'):
            {
                (): ttHNNLO_4Box_coeffs[(11,12,13,14)]
            },
            ('p11', 'p12', 'p8'):
            {
                (): ttHNNLO_4Box_coeffs[(11,12,8)]
            },
            ('p11', 'p13', 'p15', 'p16'):
            {
                (): ttHNNLO_4Box_coeffs[(11,13,15,16)]
            },
            ('p11', 'p13', 'p17'):
            {
                (): ttHNNLO_4Box_coeffs[(11,13,17)]
            },
            ('p11', 'p14', 'p15', 'p16', 'p8'):
            {
                (): ttHNNLO_4Box_coeffs[(11,14,15,16,8)]
            },
            ('p11', 'p14', 'p17', 'p8'):
            {
                (): ttHNNLO_4Box_coeffs[(11,14,17,8)]
            },
        }
        )
        ttHNNLO_4Box.export('ttHNNLO_4Box_no_numerator.yaml')

        ttHNLO_4Box = SquaredTopologyGenerator(
            [
                ('q1', 101, 1), ('q2', 102, 1), ('q3', 2, 103), ('q4', 2, 104), 
                ('p5', 1, 3), ('p6', 4, 2), ('p7', 3, 5), ('p8', 6, 3), ('p9',7, 4), 
                ('p10', 4, 8), ('p11', 7, 5), ('p12', 5, 9), ('p13',8, 6), ('p14', 9, 6),
                ('p15', 10, 7), ('p16', 8, 10), ('p17', 9, 10)
            ],
        "ttHNLO_4Box", ['q1', 'q2'], 
        1, 
        {'q1': [1000., 0., 0., 1000.], 'q2': [1000., 0., 0., -1000.], 'q3': [1000., 0., 0., 1000.], 'q4': [1000., 0., 0., -1000.]},
        # This is the basis generated from QGRAF
#        loop_momenta_names=(
#            'p7', # k1
#            'p9', # -k2
#            'p11', # -k3, Higgs
#            'p13' # -k4
#        ),
        # This is the basis that puts the gluon on loop momenta:
        loop_momenta_names=(
            'p7', # k1
            'p14', # k1+k4-p1-p2
            'p11', # -k3
            'p16' # -k2+k4-p1-p2
        ),
        final_state_particle_ids=(6,6,25),
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
        overall_numerator=overall_numerator
        ,numerator_structure={
            ('p10', 'p11', 'p12', 'p14', 'p16'):
            {
                (): ttHNNLO_4Box_coeffs[(10,11,12,14,16)]
            },
            ('p10', 'p11', 'p15'):
            {
                (): ttHNNLO_4Box_coeffs[(10,11,15)]
            },
            ('p10', 'p11', 'p16', 'p17'):
            {
                (): ttHNNLO_4Box_coeffs[(10,11,16,17)]
            },
            ('p11', 'p12', 'p13', 'p14'):
            {
                (): ttHNNLO_4Box_coeffs[(11,12,13,14)]
            },
            ('p11', 'p12', 'p8'):
            {
                (): ttHNNLO_4Box_coeffs[(11,12,8)]
            },
            ('p11', 'p13', 'p15', 'p16'):
            {
                (): ttHNNLO_4Box_coeffs[(11,13,15,16)]
            },
            ('p11', 'p13', 'p17'):
            {
                (): ttHNNLO_4Box_coeffs[(11,13,17)]
            },
            ('p11', 'p14', 'p15', 'p16', 'p8'):
            {
                (): ttHNNLO_4Box_coeffs[(11,14,15,16,8)]
            },
            ('p11', 'p14', 'p17', 'p8'):
            {
                (): ttHNNLO_4Box_coeffs[(11,14,17,8)]
            },
        }
        )
        ttHNLO_4Box.export('ttHNLO_4Box_no_numerator.yaml')

        ttHLO_4Box = SquaredTopologyGenerator(
            [
                ('q1', 101, 1), ('q2', 102, 1), ('q3', 2, 103), ('q4', 2, 104), 
                ('p5', 1, 3), ('p6', 4, 2), ('p7', 3, 5), ('p8', 6, 3), ('p9',7, 4), 
                ('p10', 4, 8), ('p11', 7, 5), ('p12', 5, 9), ('p13',8, 6), ('p14', 9, 6),
                ('p15', 10, 7), ('p16', 8, 10), ('p17', 9, 10)
            ],
        "ttHLO_4Box", ['q1', 'q2'], 
        2, 
        {'q1': [1000., 0., 0., 1000.], 'q2': [1000., 0., 0., -1000.], 'q3': [1000., 0., 0., 1000.], 'q4': [1000., 0., 0., -1000.]},
        # This is the basis generated from QGRAF
#        loop_momenta_names=(
#            'p7', # k1
#            'p9', # -k2
#            'p11', # -k3, Higgs
#            'p13' # -k4
#        ),
        # This is the basis that puts the gluon on loop momenta:
        loop_momenta_names=(
            'p7', # k1
            'p14', # k1+k4-p1-p2
            'p11', # -k3
            'p16' # -k2+k4-p1-p2
        ),
        final_state_particle_ids=(6,6,25),
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
        overall_numerator=overall_numerator
        ,numerator_structure={
            ('p10', 'p11', 'p12', 'p14', 'p16'):
            {
                (): ttHNNLO_4Box_coeffs[(10,11,12,14,16)]
            },
            ('p10', 'p11', 'p15'):
            {
                (): ttHNNLO_4Box_coeffs[(10,11,15)]
            },
            ('p10', 'p11', 'p16', 'p17'):
            {
                (): ttHNNLO_4Box_coeffs[(10,11,16,17)]
            },
            ('p11', 'p12', 'p13', 'p14'):
            {
                (): ttHNNLO_4Box_coeffs[(11,12,13,14)]
            },
            ('p11', 'p12', 'p8'):
            {
                (): ttHNNLO_4Box_coeffs[(11,12,8)]
            },
            ('p11', 'p13', 'p15', 'p16'):
            {
                (): ttHNNLO_4Box_coeffs[(11,13,15,16)]
            },
            ('p11', 'p13', 'p17'):
            {
                (): ttHNNLO_4Box_coeffs[(11,13,17)]
            },
            ('p11', 'p14', 'p15', 'p16', 'p8'):
            {
                (): ttHNNLO_4Box_coeffs[(11,14,15,16,8)]
            },
            ('p11', 'p14', 'p17', 'p8'):
            {
                (): ttHNNLO_4Box_coeffs[(11,14,17,8)]
            },
        }
        )
        ttHLO_4Box.export('ttHLO_4Box_no_numerator.yaml')

