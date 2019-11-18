#!/usr/bin/env python3
import os
import vectors
import math
import sys

import mpmath
import numpy
from scipy.special import zeta

zero_lv = vectors.LorentzVector([0.,0.,0.,0.])

from ltd_utils import LoopTopology, LoopLine, Propagator, TopologyCollection, TopologyGenerator
from analytical_expressions  import ladder_phi, analytic_four_point_ladder, analytic_three_point_ladder, analytic_two_point_ladder

#############################################################################################################
# Create the collection of hard-coded topologies.
#############################################################################################################

# Add Now automatically generated topologies
def load(selected_topologies=None):

    all_topologies = TopologyCollection()

    topology_name = "T1_Pentabox_physical"
    if selected_topologies is None or topology_name in selected_topologies:
        
        rescaling = 1.0e-3
        # This topology correspond to a two-loop pentabox with physical kinematics borrowed from
        # the process `t t~ > w+ z w-` at 1 TeV E_com.
        q1 = vectors.LorentzVector(
            # mass: 173.0 incoming top quark
            [0.5980260048915123e+03,   0.0000000000000000e+00,   0.0000000000000000e+00,   0.5724562014045295e+03]
        )*rescaling 
        q2 = vectors.LorentzVector(
            # mass: 173.0 incoming anti-top quark
            [0.5980260048915123e+03,   0.0000000000000000e+00,   0.0000000000000000e+00,   -0.5724562014045295e+03]
        )*rescaling
        q3 = vectors.LorentzVector(
            # mass: 80.419 outgoing W+
            [-0.5394473213122507e+03,   -0.1971081698462961e+03,   -0.4416135519343869e+03,  0.2250822886064787e+03]
        )*rescaling
        q4 = vectors.LorentzVector(
            # mass: 91.188 outgoing Z
            [-0.2255538754188549e+03,  0.1757868459829899e+03,   0.3716353112335996e+02,   -0.1013763093935658e+03]            
        )*rescaling
        # mass: 80.419 outgoing W-
        q5 = -q4-q3-q2-q1

        factory = TopologyGenerator([
            ('p1', 1, 6), ('p2', 6, 7), ('p3', 7, 2), ('p4', 2, 1),  
            ('p5', 7, 3), ('p6', 3, 4), ('p7', 4, 5), ('p8', 5, 6),
            ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5),
        ])
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5}, 
                mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0, 'p6': 0.0, 'p7': 0.0, 'p8': 0.0}, 
                loop_momenta_names=('p4','p5'), # If not specified an arbitrary spanning tree will be used for momentum routing 
                analytic_result=None
             ),
             entry_name = topology_name
        )

    topology_name = "T1_Pentabox_physical_bis"
    if selected_topologies is None or topology_name in selected_topologies:

        rescaling = 1.0e0
        # This topology correspond to a two-loop pentabox with physical kinematics
        q1 = vectors.LorentzVector(
            [0.149500000000000E+01,    0.000000000000000E+00,    0.000000000000000E+00,    0.149165176901313E+01]
        )*rescaling
        q2 = vectors.LorentzVector(
            [0.150500000000000E+01,    0.000000000000000E+00,    0.000000000000000E+00,   -0.149165176901313E+01]
        )*rescaling
        q3 = vectors.LorentzVector(
            [-0.126041949101381e+01,    -0.452362952912639e+00,    -0.101350243653045e+01,   0.516563513332600e+00]
        )*rescaling
        q4 = vectors.LorentzVector(
            [-0.105098730574850e+01,   0.489324061520790e-01,   0.928212188578101e+00,    -0.283905035967510e+00]
        )*rescaling
        q5 = -q4-q3-q2-q1

        factory = TopologyGenerator([
            ('p1', 1, 6), ('p2', 6, 7), ('p3', 7, 2), ('p4', 2, 1),
            ('p5', 7, 3), ('p6', 3, 4), ('p7', 4, 5), ('p8', 5, 6),
            ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5),
        ])
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name,
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5},
                mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0, 'p6': 0.0, 'p7': 0.0, 'p8': 0.0},
                loop_momenta_names=('p4','p5'), # If not specified an arbitrary spanning tree will be used for momentum routing
                analytic_result=None
             ),
             entry_name = topology_name
        )

    topology_name = "T1_Pentabox_physical_bis_massive"
    if selected_topologies is None or topology_name in selected_topologies:

        rescaling = 1.0e0
        # This topology correspond to a two-loop pentabox with physical kinematics
        q1 = vectors.LorentzVector(
            [0.149500000000000E+01,    0.000000000000000E+00,    0.000000000000000E+00,    0.149165176901313E+01]
        )*rescaling
        q2 = vectors.LorentzVector(
            [0.150500000000000E+01,    0.000000000000000E+00,    0.000000000000000E+00,   -0.149165176901313E+01]
        )*rescaling
        q3 = vectors.LorentzVector(
            [-0.126041949101381e+01,    -0.452362952912639e+00,    -0.101350243653045e+01,   0.516563513332600e+00]
        )*rescaling
        q4 = vectors.LorentzVector(
            [-0.105098730574850e+01,   0.489324061520790e-01,   0.928212188578101e+00,    -0.283905035967510e+00]
        )*rescaling
        q5 = -q4-q3-q2-q1

        factory = TopologyGenerator([
            ('p1', 1, 6), ('p2', 6, 7), ('p3', 7, 2), ('p4', 2, 1),
            ('p5', 7, 3), ('p6', 3, 4), ('p7', 4, 5), ('p8', 5, 6),
            ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5),
        ])
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name,
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5},
                mass_map={'p1': 0.35, 'p2': 0.35, 'p3': 0.35, 'p4': 0.35, 'p5': 0.35, 'p6': 0.35, 'p7': 0.35, 'p8': 0.35},
                loop_momenta_names=('p4','p5'), # If not specified an arbitrary spanning tree will be used for momentum routing
                analytic_result=None
             ),
             entry_name = topology_name
        )

    # Setting up Weinzierl's kinematic for all six 6-point 2-loop topologies of Fig.8 of ref. https://arxiv.org/pdf/1211.0509.pdf
    rescaling = 1.0e-2
    q1 = vectors.LorentzVector(
        [-12.0588,1.00017,2.55373,-2.65288]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [-18.6089,8.9195,-6.43508,-9.61832]
    )*rescaling
    q3 = vectors.LorentzVector(
        [-13.8389,4.73227,6.55009,1.55854]
    )*rescaling
    q4 = vectors.LorentzVector(
        [-25.5377,-20.892,-4.32472,8.89684] 
    )*rescaling
    q6 = vectors.LorentzVector(
        [90.0,0.0,0.0,0.0]            
    )*rescaling
    q5 = -q6-q4-q3-q2-q1

    topology_name = "T2_6P_2L_Weinzierl_A"
    if selected_topologies is None or topology_name in selected_topologies:
        factory = TopologyGenerator([
            ('p1', 1, 7), ('p2', 7, 8), ('p3', 8, 1),
            ('p4', 7, 2), ('p5', 2, 3), ('p6', 3, 4), ('p7', 4, 5), ('p8', 5, 6), ('p9', 6, 8),
            ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6)
        ])
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0, 'p6': 0.0, 'p7': 0.0, 'p8': 0.0, 'p9': 0.0}, 
                loop_momenta_names=('p1','p4'), 
                # 86.0768134710165628 - 0.0552231059543505382*I +/- ( 0.0858487113390352118 + 0.0297582753855115178*I )
                analytic_result=complex(-86.07 , 0.0j)
             ),
             entry_name = topology_name
        )

    topology_name = "T2_6P_2L_Weinzierl_B"
    if selected_topologies is None or topology_name in selected_topologies:
        factory = TopologyGenerator([
            ('p1', 2, 7), ('p2', 7, 8), ('p3', 8, 1), ('p4', 1, 2),
            ('p5', 7, 3), ('p6', 3, 4), ('p7', 4, 5), ('p8', 5, 6), ('p9', 6, 8), 
            ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6)
        ])
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0, 'p6': 0.0, 'p7': 0.0, 'p8': 0.0, 'p9': 0.0}, 
                loop_momenta_names=('p1','p5'), 
                # 118.097115126810919 - 0.0587716403920794492*I +/- ( 0.341363419341972997 + 0.0760029534962764215*I )
                analytic_result=complex(-118.09,0.0j)
             ),
             entry_name = topology_name
        )

    topology_name = "T2_6P_2L_Weinzierl_C"
    if selected_topologies is None or topology_name in selected_topologies:
        factory = TopologyGenerator([
            ('p1', 3, 7), ('p2', 7, 8), ('p3', 8, 1), ('p4', 1, 2), ('p5', 2, 3),
            ('p6', 7, 4), ('p7', 4, 5), ('p8', 5, 6), ('p9', 6, 8),
            ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6)
        ])
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0, 'p6': 0.0, 'p7': 0.0, 'p8': 0.0, 'p9': 0.0}, 
                loop_momenta_names=('p1','p6'), 
                # 75.487251662971792 - 0.125562683853285491*I +/- ( 0.187605106060742998 + 0.0363941907232131749*I )
                analytic_result=complex(-75.48,0.0j)
             ),
             entry_name = topology_name
        )

    topology_name = "T2_6P_2L_Weinzierl_D"
    if selected_topologies is None or topology_name in selected_topologies:
        factory = TopologyGenerator([
            ('p1', 1, 7), ('p2', 7, 6), ('p3', 6, 8), ('p4', 8, 1),
            ('p5', 7, 2), ('p6', 2, 3), ('p7', 3, 4), ('p8', 4, 5), ('p9', 5, 8),
            ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6)
        ])
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0, 'p6': 0.0, 'p7': 0.0, 'p8': 0.0, 'p9': 0.0}, 
                loop_momenta_names=('p1','p5'),
                # eps^0: 18.3271249575150019 - 0.0102603219383815086*I +/- ( 0.0171281153680031097 + 0.00439019843081463525*I )
                # eps^0: 18.2273146863270607 - 0.0351072120390151432*I +/- ( 0.0506710500921387608 + 0.0125415378798960689*I )
                analytic_result=complex(-18.22,0.0j)
             ),
             entry_name = topology_name
        )

    topology_name = "T2_6P_2L_Weinzierl_E"
    if selected_topologies is None or topology_name in selected_topologies:
        factory = TopologyGenerator([
            ('p1', 1, 7), ('p2', 7, 6), ('p3', 6, 5), ('p4', 5, 8), ('p5', 8, 1),
            ('p6', 7, 2), ('p7', 2, 3), ('p8', 3, 4), ('p9', 4, 8),
            ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6)
        ])
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0, 'p6': 0.0, 'p7': 0.0, 'p8': 0.0, 'p9': 0.0}, 
                loop_momenta_names=('p1','p6'),
                # diff. seed : 45.9727585576039985 - 0.00210426401936783294*I +/- ( 0.0443921298216177022 + 0.0027359584173322363*I )
                # higher-stat : 45.9306553282958365 - 0.0227134030759781315*I +/- ( 0.0448969823894826872 + 0.00814578717986984305*I )
                # low-stat : 45.4459598582354531 - 0.188225725121877622*I +/- ( 0.115356519900322018 + 0.0285153450900536506*I )
                analytic_result=complex(-45.45,0.0j)
             ),
             entry_name = topology_name
        )

    topology_name = "T2_6P_2L_Weinzierl_F"
    if selected_topologies is None or topology_name in selected_topologies:
        factory = TopologyGenerator([
            ('p1', 2, 7), ('p2', 7, 6), ('p3', 6, 5), ('p4', 5, 8), ('p5', 8, 1), ('p6', 1, 2),
            ('p7', 7, 3), ('p8', 3, 4), ('p9', 4, 8),
            ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6)
        ])
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0, 'p6': 0.0, 'p7': 0.0, 'p8': 0.0, 'p9': 0.0}, 
                loop_momenta_names=('p1','p7'),
                # 102.713355708622871 - 0.00969731708704666883*I +/- ( 0.0288577968836020429 + 0.00604760342011912482*I )
                analytic_result=complex(-102.71,0.0j)
             ),
             entry_name = topology_name
        )



    #2>N 3LOOP kinematics
    q1 = vectors.LorentzVector(
        [3, 0., 0., 1.21]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [3.893, 0., 0., -1.21]
    )*rescaling
    q4 = vectors.LorentzVector(
        [-2.527, 0., 0.9736, -0.4657] 
    )*rescaling
    q3 = -q1-q2-q4
    topology_name = "2to2_TripleBox_v1"
    if selected_topologies is None or topology_name in selected_topologies:
        factory = TopologyGenerator([
            ('p1', 2, 5), ('p2', 5, 6), ('p3', 6, 1), ('p4', 1, 2),
            ('p5', 5, 7), ('p6', 7, 8), ('p7', 8, 6),
            ('p8', 7, 3), ('p9', 3, 4), ('p10', 4, 8),
            ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4)
        ])
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4}, 
                mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 
                          'p5': 0.0, 'p6': 0.0, 'p7': 0.0, 'p8': 0.0, 
                          'p9': 0.0, 'p10': 0.0}, 
                loop_momenta_names=('p1','p5','p8'), 
                analytic_result=analytic_four_point_ladder(q1.square(),q2.square(),q3.square(),q4.square(), (q1+q2).square(), (q2+q3).square(), 3)
             ),
             entry_name = topology_name
        )


    #DoubleBox 2to2 kinematics
    rescaling = 1
    q1 = vectors.LorentzVector(
        [2.,0.,0.,0.45]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [1.43,0.,0.,-0.45]
    )*rescaling
    q4 = vectors.LorentzVector(
        [-1.012, 0., 0.5281, -0.278] 
    )*rescaling
    q3 = -q1-q2-q4

    topology_name = "2to2_DoubleBox"
    if selected_topologies is None or topology_name in selected_topologies:
        factory = TopologyGenerator([
            ('p1', 2, 5), ('p2', 5, 6), ('p3', 6, 1), ('p4', 1, 2),
            ('p5', 5, 3), ('p6', 3, 4), ('p7', 4, 6),
            ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4)
        ])
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4,}, 
                mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0, 'p6': 0.0, 'p7': 0.0}, 
                loop_momenta_names=('p1','p5'), 
                analytic_result=analytic_four_point_ladder(q1.square(),q2.square(),q3.square(),q4.square(),(q1+q2).square(),(q2+q3).square(),2)
             ),
             entry_name = topology_name
        )





    # Setting up Weinzierl's kinematic for the four-point double and triple box of Weinzierl's paper.
    # See Eq. 35 of ref. https://arxiv.org/pdf/1211.0509.pdf
    rescaling = 1.0e-2
    q1 = vectors.LorentzVector(
        [-19.6586,7.15252,0.206016,-8.96383]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [-26.874,-7.04203,0.0501295,12.9055]
    )*rescaling
    q4 = vectors.LorentzVector(
        [90.0, 0., 0., 0.] 
    )*rescaling
    q3 = -q1-q2-q4

    topology_name = "T3_DoubleBox_Weinzierl"
    if selected_topologies is None or topology_name in selected_topologies:
        factory = TopologyGenerator([
            ('p1', 2, 5), ('p2', 5, 6), ('p3', 6, 1), ('p4', 1, 2),
            ('p5', 5, 3), ('p6', 3, 4), ('p7', 4, 6),
            ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4)
        ])
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4,}, 
                mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0, 'p6': 0.0, 'p7': 0.0}, 
                loop_momenta_names=('p1','p5'), 
                analytic_result=complex(-5.897e-2,0.0j)
             ),
             entry_name = topology_name
        )

    topology_name = "T4_TripleBox_Weinzierl"
    if selected_topologies is None or topology_name in selected_topologies:
        factory = TopologyGenerator([
            ('p1', 2, 5), ('p2', 5, 6), ('p3', 6, 1), ('p4', 1, 2),
            ('p5', 5, 7), ('p6', 7, 8), ('p7', 8, 6),
            ('p8', 7, 3), ('p9', 3, 4), ('p10', 4, 8),
            ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4)
        ])
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4}, 
                mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 
                          'p5': 0.0, 'p6': 0.0, 'p7': 0.0, 'p8': 0.0, 
                          'p9': 0.0, 'p10': 0.0}, 
                loop_momenta_names=('p1','p5','p8'), 
                analytic_result=complex(0.0,-6.744e-3j)
             ),
             entry_name = topology_name
        )

    topology_name = "T5_3L_5P"
    if selected_topologies is None or topology_name in selected_topologies:

        rescaling = 1.0e0
        # This topology correspond to a two-loop pentabox with physical kinematics
        q1 = vectors.LorentzVector(
            [0.149500000000000E+01,    0.000000000000000E+00,    0.000000000000000E+00,    0.149165176901313E+01]
        )*rescaling
        q2 = vectors.LorentzVector(
            [0.150500000000000E+01,    0.000000000000000E+00,    0.000000000000000E+00,   -0.149165176901313E+01]
        )*rescaling
        q3 = vectors.LorentzVector(
            [-0.126041949101381e+01,    -0.452362952912639e+00,    -0.101350243653045e+01,   0.516563513332600e+00]
        )*rescaling
        q4 = vectors.LorentzVector(
            [-0.105098730574850e+01,   0.489324061520790e-01,   0.928212188578101e+00,    -0.283905035967510e+00]
        )*rescaling
        q5 = -q4-q3-q2-q1

        factory = TopologyGenerator([
            ('p1', 1, 6), ('p2', 6, 7), ('p3', 7, 2), ('p4', 2, 1),
            ('p5', 7, 8), ('p6', 8, 9), ('p7', 9, 6),
            ('p8', 8, 3), ('p9', 3, 4), ('p10', 4, 5), ('p11', 5, 9),
            ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5),
        ])
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name,
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5},
                mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0, 'p6': 0.0, 
                          'p7': 0.0, 'p8': 0.0, 'p9': 0.0, 'p10': 0.0, 'p11': 0.0},
                loop_momenta_names=('p3','p5','p8'), # If not specified an arbitrary spanning tree will be used for momentum routing
                analytic_result=None
             ),
             entry_name = topology_name
        )

    topology_name = "T5_3L_5P_massive"
    if selected_topologies is None or topology_name in selected_topologies:

        rescaling = 1.0e0
        # This topology correspond to a two-loop pentabox with physical kinematics
        q1 = vectors.LorentzVector(
            [0.149500000000000E+01,    0.000000000000000E+00,    0.000000000000000E+00,    0.149165176901313E+01]
        )*rescaling
        q2 = vectors.LorentzVector(
            [0.150500000000000E+01,    0.000000000000000E+00,    0.000000000000000E+00,   -0.149165176901313E+01]
        )*rescaling
        q3 = vectors.LorentzVector(
            [-0.126041949101381e+01,    -0.452362952912639e+00,    -0.101350243653045e+01,   0.516563513332600e+00]
        )*rescaling
        q4 = vectors.LorentzVector(
            [-0.105098730574850e+01,   0.489324061520790e-01,   0.928212188578101e+00,    -0.283905035967510e+00]
        )*rescaling
        q5 = -q4-q3-q2-q1

        factory = TopologyGenerator([
            ('p1', 1, 6), ('p2', 6, 7), ('p3', 7, 2), ('p4', 2, 1),
            ('p5', 7, 8), ('p6', 8, 9), ('p7', 9, 6),
            ('p8', 8, 3), ('p9', 3, 4), ('p10', 4, 5), ('p11', 5, 9),
            ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5),
        ])
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name,
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5},
                mass_map={'p1': 0.35, 'p2': 0.35, 'p3': 0.35, 'p4': 0.35, 'p5': 0.35, 'p6': 0.35, 
                          'p7': 0.35, 'p8': 0.35, 'p9': 0.35, 'p10': 0.35, 'p11': 0.35},
                loop_momenta_names=('p3','p5','p8'), # If not specified an arbitrary spanning tree will be used for momentum routing
                analytic_result=None
             ),
             entry_name = topology_name
        )



    # Setting up One Loop topologies with uniformly generated PS points

    topology_name = "Pentagon_1s"
    # 8 ellipsoids 1 source with a point from MadGraph's FLatInvertiblePhaseSpace
    if selected_topologies is None or topology_name in selected_topologies:
        pentagon = TopologyGenerator([
            ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),  ('p4', 4, 5),  ('p5', 5, 1),
            ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5)
        ])
        q1 = vectors.LorentzVector([5.1000000000000000e+01, 0.0000000000000000e+00, 0.0000000000000000e+00, 4.8744230427815765e+01])    #m = 15
        q2 = vectors.LorentzVector([4.9000000000000000e+01, 0.0000000000000000e+00, 0.0000000000000000e+00, -4.8744230427815765e+01])   #m = 5
        q3 = vectors.LorentzVector([-2.7232566433392297e+01, -5.5485178014592513e+00, 1.6209175964465906e+01, 6.9346405320676556e+00])  #m = 20 
        q4 = vectors.LorentzVector([-3.2743923031870672e+01, 1.2995389184547482e+01, -5.1743016921802543e-01, 1.7368423875527190e+00])  #m = 30
        q5 = vectors.LorentzVector([-4.0023510534737028e+01, -7.4468713830882312e+00, -1.5691745795247881e+01, -8.6714829196203738e+00])#m = 35
        all_topologies.add_topology(pentagon.create_loop_topology(
                "Pentagon_1s", 
                ext_mom={ 'q1': q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5  }, 
                mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0}, 
                loop_momenta_names=('p1',), # If not specified an arbitrary spanning tree will be used for momentum routing 
                analytic_result=-3.44342331910881e-012-2.56487078110481e-012j,
             ),
             entry_name = 'Pentagon_1s'
        )

    topology_name = "Pentagon_2s"
    # 8 ellipsoids 2 sources with a point from MadGraph's FLatInvertiblePhaseSpace
    if selected_topologies is None or topology_name in selected_topologies:
        pentagon = TopologyGenerator([
            ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),  ('p4', 4, 5),  ('p5', 5, 1),
            ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5)
        ])
        q1 = vectors.LorentzVector([ 4.812500e+01,0.000000e+00,0.000000e+00,4.786455e+01])      # m = 5
        q2 = vectors.LorentzVector([ 5.187500e+01,0.000000e+00,0.000000e+00,-4.786455e+01])     # m = 20
        q3 = vectors.LorentzVector([ -3.319087e+01,4.234357e+00,-3.141040e+00,-1.318475e+01])   # m = 30
        q4 = vectors.LorentzVector([ -2.090788e+01,-1.143512e+01,-4.165504e-01,9.011335e+00])   # m = 15
        q5 = vectors.LorentzVector([ -4.590125e+01,7.200760e+00,3.557590e+00,4.173417e+00])     # m = 45
        all_topologies.add_topology(pentagon.create_loop_topology(
                "Pentagon_2s", 
                ext_mom={ 'q1': q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5  }, 
                mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0}, 
                loop_momenta_names=('p1',), # If not specified an arbitrary spanning tree will be used for momentum routing 
                analytic_result=-8.39048452581577133E-013-1.71340504218085295E-012j,
             ),
             entry_name = 'Pentagon_2s' #8E 2s with a point from MadGraph's FLatInvertiblePhaseSpace
        )

    topology_name = "Pentagon_3s"
    # 8 ellipsoids 3 sources with a point from MadGraph's FLatInvertiblePhaseSpace
    if selected_topologies is None or topology_name in selected_topologies:
        pentagon = TopologyGenerator([
            ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),  ('p4', 4, 5),  ('p5', 5, 1),
            ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5)
        ])
        q1 = vectors.LorentzVector([ 4.809394e+01,0.000000e+00,0.000000e+00,4.805087e+01]) #m=2
        q2 = vectors.LorentzVector([ 5.190606e+01,0.000000e+00,0.000000e+00,-4.805087e+01]) #m=19
        q3 = vectors.LorentzVector([ -3.068412e+01,6.865600e+00,-2.846056e+01,5.413178e-01]) #m=9
        q4 = vectors.LorentzVector([ -3.607736e+01,-9.406076e+00,2.938819e+01,1.081775e+01]) #m=15
        q5 = vectors.LorentzVector([ -3.323852e+01,2.540476e+00,-9.276298e-01,-1.135907e+01]) #m=31
        all_topologies.add_topology(pentagon.create_loop_topology(
                "Pentagon_3s", 
                ext_mom={ 'q1': q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5  }, 
                mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0}, 
                loop_momenta_names=('p1',), # If not specified an arbitrary spanning tree will be used for momentum routing 
                analytic_result=-3.48997234169132E-012-3.90012989047481E-012j,
             ),
             entry_name = 'Pentagon_3s' #8E 3s with a point from MadGraph's FLatInvertiblePhaseSpace
        )

    topology_name = "Hexagon_1s"
    # 12 ellipsoids 1 source with a point from MadGraph's FLatInvertiblePhaseSpace
    if selected_topologies is None or topology_name in selected_topologies:
        hexagon = TopologyGenerator([
            ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),  ('p4', 4, 5),  ('p5', 5, 6), ('p6', 6, 1),
            ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6)
        ])
        q1 = vectors.LorentzVector([ 5.226362e+01,0.000000e+00,0.000000e+00,2.682755e+01]) #m=44.852746
        q2 = vectors.LorentzVector([ 4.773638e+01,0.000000e+00,0.000000e+00,-2.682755e+01]) #m=39.484736
        q3 = vectors.LorentzVector([ -2.285915e+01,-3.498550e+00,-2.596035e+00,-2.026016e+00]) #m=22.348531
        q4 = vectors.LorentzVector([ -3.149256e+00,-6.632673e-01,1.218285e+00,2.711957e+00]) #m=0.799353
        q5 = vectors.LorentzVector([ -2.363008e+01,-2.345272e+00,5.121577e-01,9.109623e+00]) #m=21.671014
        q6 = vectors.LorentzVector([ -5.036151e+01,6.507089e+00,8.655920e-01,-9.795563e+00]) #m=48.961591
        all_topologies.add_topology(hexagon.create_loop_topology(
                "Hexagon_1s", 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6': q6  }, 
                mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0, 'p6': 0.0}, 
                loop_momenta_names=('p1',), # If not specified an arbitrary spanning tree will be used for momentum routing 
                analytic_result=3.03979043760586194E-015-1.17682533261426436E-013j,
             ),
             entry_name = 'Hexagon_1s'
        )


    topology_name = "Decagon_phys"
    # 12 ellipsoids 1 source with a point from MadGraph's FLatInvertiblePhaseSpace
    if selected_topologies is None or topology_name in selected_topologies:
        Decagon_phys = TopologyGenerator([
            ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),  ('p4', 4, 5),  ('p5', 5, 6), ('p6', 6, 7),('p7', 7, 8),('p8', 8, 9),('p9', 9, 10),('p10', 10, 1),
            ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6),('q7', 107,7),('q8', 108,8),('q9', 109,9),('q10', 110,10),
        ])
        q1 = vectors.LorentzVector([0.8249090909090908E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.8248484759423774E+00]) #m=44.852746
        q2 = vectors.LorentzVector([0.8250909090909090E+00,0.0000000000000000E+00,0.0000000000000000E+00,-0.8248484759423774E+00]) #m=39.484736
        q3 = vectors.LorentzVector([-0.1085133239146588E+00,0.1167616774783045E-01,-0.6413747334200047E-01,0.8139528909609182E-01]) #m=22.348531
        q4 = vectors.LorentzVector([-0.1871911023039417E+00,0.7334738260816014E-01,0.1671281314565328E+00,-0.1135156033965376E-01]) #m=0.799353
        q5 = vectors.LorentzVector([-0.1046974473941988E+00,0.8070338217617294E-01,0.4122487429560785E-01,-0.1578066334153676E-01]) #m=21.671014
        q6 = vectors.LorentzVector([-0.5167577589181960E+00,-0.2202964129865497E+00,-0.4438512772975853E+00,0.1338062609473788E+00]) #m=48.961591
        q7 = vectors.LorentzVector([-0.9935413704930492E-01,0.3830140390088815E-01,-0.3044895751950545E-02,-0.5911831880130597E-01]) #m=22.348531
        q8 = vectors.LorentzVector([-0.2749069123381626E+00,0.3654135698965530E-01,0.7048580448475607E-01,-0.2507394884073777E+00]) #m=0.799353
        q9 = vectors.LorentzVector([-0.1942267082273213E+00,-0.7883355524646646E-01,0.1183158353754920E+00,0.9700849372824423E-01]) #m=21.671014
        q10 = vectors.LorentzVector([-0.1643526098542161E+00,0.5856027481030913E-01,0.1138790007791476E+00,0.2477998711815934E-01]) #m=48.961591

        all_topologies.add_topology(Decagon_phys.create_loop_topology(
                "Decagon_phys", 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6': q6, 'q7': q7, 'q8': q8, 'q9': q9, 'q10': q10}, 
                mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0, 'p6': 0.0, 'p7': 0., 'p8': 0., 'p9': 0., 'p10': 0.}, 
                loop_momenta_names=('p1',), # If not specified an arbitrary spanning tree will be used for momentum routing 
                analytic_result=232.45130851452109 - 4342.3577516872392j,
             ),
             entry_name = 'Decagon_phys'
        )


    topology_name = "Decagon_phys_massive"
    # 12 ellipsoids 1 source with a point from MadGraph's FLatInvertiblePhaseSpace
    if selected_topologies is None or topology_name in selected_topologies:
        Decagon_phys_massive = TopologyGenerator([
            ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),  ('p4', 4, 5),  ('p5', 5, 6), ('p6', 6, 7),('p7', 7, 8),('p8', 8, 9),('p9', 9, 10),('p10', 10, 1),
            ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6),('q7', 107,7),('q8', 108,8),('q9', 109,9),('q10', 110,10),
        ])
        q1 = vectors.LorentzVector([0.8249090909090908E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.8248484759423774E+00]) #m=44.852746
        q2 = vectors.LorentzVector([0.8250909090909090E+00,0.0000000000000000E+00,0.0000000000000000E+00,-0.8248484759423774E+00]) #m=39.484736
        q3 = vectors.LorentzVector([-0.1085133239146588E+00,0.1167616774783045E-01,-0.6413747334200047E-01,0.8139528909609182E-01]) #m=22.348531
        q4 = vectors.LorentzVector([-0.1871911023039417E+00,0.7334738260816014E-01,0.1671281314565328E+00,-0.1135156033965376E-01]) #m=0.799353
        q5 = vectors.LorentzVector([-0.1046974473941988E+00,0.8070338217617294E-01,0.4122487429560785E-01,-0.1578066334153676E-01]) #m=21.671014
        q6 = vectors.LorentzVector([-0.5167577589181960E+00,-0.2202964129865497E+00,-0.4438512772975853E+00,0.1338062609473788E+00]) #m=48.961591
        q7 = vectors.LorentzVector([-0.9935413704930492E-01,0.3830140390088815E-01,-0.3044895751950545E-02,-0.5911831880130597E-01]) #m=22.348531
        q8 = vectors.LorentzVector([-0.2749069123381626E+00,0.3654135698965530E-01,0.7048580448475607E-01,-0.2507394884073777E+00]) #m=0.799353
        q9 = vectors.LorentzVector([-0.1942267082273213E+00,-0.7883355524646646E-01,0.1183158353754920E+00,0.9700849372824423E-01]) #m=21.671014
        q10 = vectors.LorentzVector([-0.1643526098542161E+00,0.5856027481030913E-01,0.1138790007791476E+00,0.2477998711815934E-01]) #m=48.961591

        all_topologies.add_topology(Decagon_phys_massive.create_loop_topology(
                "Decagon_phys_massive", 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6': q6, 'q7': q7, 'q8': q8, 'q9': q9, 'q10': q10}, 
                mass_map={'p1': 0.045, 'p2': 0.045, 'p3': 0.045, 'p4': 0.045, 'p5': 0.045, 'p6': 0.045, 'p7': 0.045, 'p8': 0.045, 'p9': 0.045, 'p10': 0.045}, 
                loop_momenta_names=('p1',), # If not specified an arbitrary spanning tree will be used for momentum routing 
                analytic_result=-3000.1965833281529 +2125.2823773318205j,
             ),
             entry_name = 'Decagon_phys_massive'
        )

    topology_name = "Decagon_phys_ez"
    # 12 ellipsoids 1 source with a point from MadGraph's FLatInvertiblePhaseSpace
    if selected_topologies is None or topology_name in selected_topologies:
        Decagon_phys_ez = TopologyGenerator([
            ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),  ('p4', 4, 5),  ('p5', 5, 6), ('p6', 6, 7),('p7', 7, 8),('p8', 8, 9),('p9', 9, 10),('p10', 10, 1),
            ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6),('q7', 107,7),('q8', 108,8),('q9', 109,9),('q10', 110,10),
        ])
        q1 = vectors.LorentzVector([0.9999999999999999E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.9949874371066199E+00]) #m=44.852746
        q2 = vectors.LorentzVector([0.9999999999999999E+00,0.0000000000000000E+00,0.0000000000000000E+00,-0.9949874371066199E+00]) #m=39.484736
        q3 = vectors.LorentzVector([-0.1571649989067209E+00,0.1357543946677973E-01,-0.7457021907451697E-01,0.9463522997178132E-01]) #m=22.348531
        q4 = vectors.LorentzVector([-0.2349560134737885E+00,0.8527823290555461E-01,0.1943135693820209E+00,-0.1319803068717482E-01]) #m=0.799353
        q5 = vectors.LorentzVector([-0.1464178649067739E+00,0.9383077591537656E-01,0.4793060511053515E-01,-0.1834758154947442E-01]) #m=21.671014
        q6 = vectors.LorentzVector([ -0.6050718384523494E+00,-0.2561303232221291E+00,-0.5160491247931777E+00,0.1555714884367599E+00]) #m=48.961591
        q7 = vectors.LorentzVector([-0.1293060184016942E+00,0.4453159644317389E-01,-0.3540185346424128E-02,-0.6873463756240042E-01]) #m=22.348531
        q8 = vectors.LorentzVector([-0.3217266170955511E+00,0.4248525634099662E-01,0.8195118404562783E-01,-0.2915253377922865E+00]) #m=0.799353
        q9 = vectors.LorentzVector([-0.2237079654788912E+00,-0.9165679872989961E-01,0.1375613553856241E+00,0.1127881136013963E+00]) #m=21.671014
        q10 = vectors.LorentzVector([-0.1816486832842308E+00,0.6808582088014734E-01,0.1324028152903109E+00,0.2881075558139869E-01]) #m=48.961591

        all_topologies.add_topology(Decagon_phys_ez.create_loop_topology(
                "Decagon_phys_ez", 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6': q6, 'q7': q7, 'q8': q8, 'q9': q9, 'q10': q10}, 
                mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0, 'p6': 0.0, 'p7': 0.0, 'p8': 0.0, 'p9': 0.0, 'p10': 0.0}, 
                loop_momenta_names=('p1',), # If not specified an arbitrary spanning tree will be used for momentum routing 
                analytic_result=2.2076906500021636 - 85.761425296521352j,
             ),
             entry_name = 'Decagon_phys_ez'
        )

    topology_name = "Decagon_phys_massive_ez"
    # 12 ellipsoids 1 source with a point from MadGraph's FLatInvertiblePhaseSpace
    if selected_topologies is None or topology_name in selected_topologies:
        Decagon_phys_massive_ez = TopologyGenerator([
            ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),  ('p4', 4, 5),  ('p5', 5, 6), ('p6', 6, 7),('p7', 7, 8),('p8', 8, 9),('p9', 9, 10),('p10', 10, 1),
            ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6),('q7', 107,7),('q8', 108,8),('q9', 109,9),('q10', 110,10),
        ])
        q1 = vectors.LorentzVector([0.9999999999999999E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.9949874371066199E+00]) #m=44.852746
        q2 = vectors.LorentzVector([0.9999999999999999E+00,0.0000000000000000E+00,0.0000000000000000E+00,-0.9949874371066199E+00]) #m=39.484736
        q3 = vectors.LorentzVector([-0.1571649989067209E+00,0.1357543946677973E-01,-0.7457021907451697E-01,0.9463522997178132E-01]) #m=22.348531
        q4 = vectors.LorentzVector([-0.2349560134737885E+00,0.8527823290555461E-01,0.1943135693820209E+00,-0.1319803068717482E-01]) #m=0.799353
        q5 = vectors.LorentzVector([-0.1464178649067739E+00,0.9383077591537656E-01,0.4793060511053515E-01,-0.1834758154947442E-01]) #m=21.671014
        q6 = vectors.LorentzVector([ -0.6050718384523494E+00,-0.2561303232221291E+00,-0.5160491247931777E+00,0.1555714884367599E+00]) #m=48.961591
        q7 = vectors.LorentzVector([-0.1293060184016942E+00,0.4453159644317389E-01,-0.3540185346424128E-02,-0.6873463756240042E-01]) #m=22.348531
        q8 = vectors.LorentzVector([-0.3217266170955511E+00,0.4248525634099662E-01,0.8195118404562783E-01,-0.2915253377922865E+00]) #m=0.799353
        q9 = vectors.LorentzVector([-0.2237079654788912E+00,-0.9165679872989961E-01,0.1375613553856241E+00,0.1127881136013963E+00]) #m=21.671014
        q10 = vectors.LorentzVector([-0.1816486832842308E+00,0.6808582088014734E-01,0.1324028152903109E+00,0.2881075558139869E-01]) #m=48.961591

        all_topologies.add_topology(Decagon_phys_massive_ez.create_loop_topology(
                "Decagon_phys_massive_ez", 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6': q6, 'q7': q7, 'q8': q8, 'q9': q9, 'q10': q10}, 
                mass_map={'p1': 0.1, 'p2': 0.1, 'p3': 0.1, 'p4': 0.1, 'p5': 0.1, 'p6': 0.1, 'p7': 0.1, 'p8': 0.1, 'p9': 0.1, 'p10': 0.1}, 
                loop_momenta_names=('p1',), # If not specified an arbitrary spanning tree will be used for momentum routing 
                analytic_result=-9.5855441923909925 + 67.414157164112297j,
             ),
             entry_name = 'Decagon_phys_massive_ez'
        )

    topology_name = "Hexagon_2s"
    # 12 ellipsoids 2 sources with a point from MadGraph's FLatInvertiblePhaseSpace
    if selected_topologies is None or topology_name in selected_topologies:
        hexagon = TopologyGenerator([
            ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),  ('p4', 4, 5),  ('p5', 5, 6), ('p6', 6, 1),
            ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6)
        ])
        q1 = vectors.LorentzVector([ 4.944317e+01,0.000000e+00,0.000000e+00,4.859427e+01]) #m=9
        q2 = vectors.LorentzVector([ 5.055683e+01,0.000000e+00,0.000000e+00,-4.859427e+01]) #m=13
        q3 = vectors.LorentzVector([ -1.720780e+01,1.103588e+01,8.354782e+00,7.675268e-01]) #m=10
        q4 = vectors.LorentzVector([ -1.744801e+01,9.007121e+00,-9.498643e+00,-5.609868e+00]) #m=10
        q5 = vectors.LorentzVector([ -3.656425e+01,-7.735831e+00,1.165555e+01,1.083451e+01]) #m=31
        q6 = vectors.LorentzVector([ -2.877995e+01,-1.230717e+01,-1.051169e+01,-5.992167e+00]) #m=23
        all_topologies.add_topology(hexagon.create_loop_topology(
                "Hexagon_2s", 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6': q6  }, 
                mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0, 'p6': 0.0}, 
                loop_momenta_names=('p1',), # If not specified an arbitrary spanning tree will be used for momentum routing 
                analytic_result=1.36918026543054139E-015-2.25900866711100386E-015j,
             ),
             entry_name = 'Hexagon_2s'
        )

    topology_name = "Hexagon_4s"
    # 12 ellipsoids 4 sources with a point from MadGraph's FLatInvertiblePhaseSpace
    if selected_topologies is None or topology_name in selected_topologies:
        hexagon = TopologyGenerator([
            ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),  ('p4', 4, 5),  ('p5', 5, 6), ('p6', 6, 1),
            ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6)
        ])
        q1 = vectors.LorentzVector([ 4.2745528739922953e+01,0.0000000000000000e+00,0.0000000000000000e+00,3.9734245060394706e+01]) #m=15.759759
        q2 = vectors.LorentzVector([ 5.7254471260077047e+01,0.0000000000000000e+00,0.0000000000000000e+00,-3.9734245060394706e+01]) #m=41.222133
        q3 = vectors.LorentzVector([ -3.4238985949421256e+01,-1.7923525330384898e+01,7.8375531169491950e+00,1.8894060234548277e+01]) #m=20.800064
        q4 = vectors.LorentzVector([ -2.6647614743187443e+01,1.0053514550024701e+01,6.9280420963147478e+00,-2.1196780134595919e+01]) #m=10.569814
        q5 = vectors.LorentzVector([ -1.7206142231165188e+01,-8.3577938813141195e+00,-1.4264811689133330e+01,-2.4369559632230517e+00]) #m=4.095730
        q6 = vectors.LorentzVector([ -2.1907257076226109e+01,1.6227804661674316e+01,-5.0078352413061289e-01,4.7396758632706932e+00]) #m=13.923755
        all_topologies.add_topology(hexagon.create_loop_topology(
                "Hexagon_4s", 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6': q6  }, 
                mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0, 'p6': 0.0}, 
                loop_momenta_names=('p1',), # If not specified an arbitrary spanning tree will be used for momentum routing 
                analytic_result=7.93962056654444536E-017-1.15281776541741777E-014j,
             ),
             entry_name = 'Hexagon_4s'
        )

    topology_name = "Hexagon_3s"        
    # 12 ellipsoids 4 sources with a point from MadGraph's FLatInvertiblePhaseSpace
    if selected_topologies is None or topology_name in selected_topologies:
        hexagon = TopologyGenerator([
            ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4),  ('p4', 4, 5),  ('p5', 5, 6), ('p6', 6, 1),
            ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6)
        ])
        q1 = vectors.LorentzVector([ 5.4189671909701467e+01,0.0000000000000000e+00,0.0000000000000000e+00,4.4892435038200482e+01]) #m=30.351109
        q2 = vectors.LorentzVector([ 4.5810328090298533e+01,0.0000000000000000e+00,0.0000000000000000e+00,-4.4892435038200482e+01]) #m=9.124442
        q3 = vectors.LorentzVector([ -1.1052241235704336e+01,4.7965011605217480e+00,2.9270203579361931e-01,-9.9496941946437296e+00]) #m=0.252039
        q4 = vectors.LorentzVector([ -3.6743682919549229e+01,2.0343318097720404e+00,1.7505904174804598e+01,-1.2330041917636034e+01]) #m=29.790487
        q5 = vectors.LorentzVector([ -2.7515446756514880e+01,-1.0275942918391802e+01,-1.1985662984602909e+01,1.3882072879012041e+01]) #m=17.752091
        q6 = vectors.LorentzVector([ -2.4688629088231536e+01,3.4451099480980139e+00,-5.8129432259953102e+00,8.3976632332677230e+00]) #m=22.211451
        all_topologies.add_topology(hexagon.create_loop_topology(
                "Hexagon_3s", 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6': q6  }, 
                mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0, 'p5': 0.0, 'p6': 0.0}, 
                loop_momenta_names=('p1',), # If not specified an arbitrary spanning tree will be used for momentum routing 
                analytic_result=1.29770369586567909E-015-2.16589510654016577E-015j,
             ),
             entry_name = 'Hexagon_3s'
        )

    topology_name = "T4_Quadruple_Box_Weinzierl"

    if selected_topologies is None or topology_name in selected_topologies:
        rescaling = 1.0e-2
        q1 = vectors.LorentzVector([-19.6586,7.15252,0.206016,-8.96383])*rescaling
        q2 = vectors.LorentzVector([-26.874,-7.04203,0.0501295,12.9055])*rescaling
        q4 = vectors.LorentzVector([90.0, 0., 0., 0.])*rescaling
        q3 = -q1-q2-q4
        analytic_result = analytic_four_point_ladder(q1.square(), q2.square(), q3.square(), q4.square(), (q1+q2).square(), (q1+q3).square(), 4)
        analytic_result = complex(analytic_result.real, 0.0)
        factory = TopologyGenerator([
                        ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4), ('p4', 4, 1),
                        ('p5', 1, 5), ('p6', 5, 6), ('p7', 6, 4),
                        ('p8', 5, 7), ('p9', 7, 8), ('p10', 8, 6),
                        ('p11', 7, 9), ('p12', 9, 10), ('p13', 10, 8),
                        ('q1', 101,2), ('q2', 102,3), ('q3', 103,9), ('q4', 104,10)])
        all_topologies.add_topology(factory.create_loop_topology(
            topology_name,
            ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4},
            mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0,
                'p5': 0.0, 'p6': 0.0, 'p7': 0.0, 'p8': 0.0,
                'p9': 0.0, 'p10': 0.0, 'p11': 0.0, 'p12': 0.0, 'p13': 0.0},
            loop_momenta_names=('p1','p5','p8','p11'),
            analytic_result=analytic_result,
            ),
            entry_name = topology_name
        )

    topology_name = "Evil_box"
    if selected_topologies is None or topology_name in selected_topologies:
        q1 = vectors.LorentzVector([1.1180339887498949e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 2.5000000000000000e+00])
        q2 = vectors.LorentzVector([-1.1180339887498949e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, -1.5000000000000000e+00])
        q3 = vectors.LorentzVector([-4.6733820729745602e+01, -1.5580614736124998e+01, -4.4068633339876151e+01, -5.0000000000000000e-01])
        q4 = -q1-q2-q3
        analytic_result = analytic_four_point_ladder(q1.square(), q2.square(), q3.square(), q4.square(), (q1+q2).square(), (q2+q3).square(), 1)
        factory = TopologyGenerator([
            ('q1', 101, 1), ('q2', 102, 2), ('q3', 103, 3), ('q4', 104, 4),
            ('p1', 1, 4), ('p2', 4, 3), ('p3', 3, 2), ('p4', 2, 1),
        ])
        all_topologies.add_topology(factory.create_loop_topology(
            topology_name,
            ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4},
            mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0},
            loop_momenta_names=('p4',),
            analytic_result=analytic_result,
            ),
            entry_name = topology_name
        )



    topology_name = "Evil_box_v2"
    if selected_topologies is None or topology_name in selected_topologies:
        q1 = vectors.LorentzVector([1.1180339887498949e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 2.5000000000000000e+00])
        q2 = vectors.LorentzVector([-1.1180339887498949e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, -1.5000000000000000e+00])
        q3 = vectors.LorentzVector([-33.31741286474686, -11.109555446651418, -31.42256796918065, -0.49999999999999994])
        q4 = -q1-q2-q3
        factory = TopologyGenerator([
            ('q1', 101, 1), ('q2', 102, 2), ('q3', 103, 3), ('q4', 104, 4),
            ('p1', 1, 4), ('p2', 4, 3), ('p3', 3, 2), ('p4', 2, 1),
        ])
        all_topologies.add_topology(factory.create_loop_topology(
            topology_name,
            ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4},
            mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0},
            loop_momenta_names=('p4',),
            analytic_result=-0.0018314315398-0.00044180461139j,
            ),
            entry_name = topology_name
        )

    topology_name = "Evil_box_v3"
    if selected_topologies is None or topology_name in selected_topologies:
        q1 = vectors.LorentzVector([1.1180339887498949e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 2.5000000000000000e+00])
        q2 = vectors.LorentzVector([-1.1180339887498949e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, -1.5000000000000000e+00])
        q3 = vectors.LorentzVector([-1.118033988749895, -0.4714045207910317, -1.3333333333333333, -0.5])
        q4 = -q1-q2-q3
        factory = TopologyGenerator([
            ('q1', 101, 1), ('q2', 102, 2), ('q3', 103, 3), ('q4', 104, 4),
            ('p1', 1, 4), ('p2', 4, 3), ('p3', 3, 2), ('p4', 2, 1),
        ])
        all_topologies.add_topology(factory.create_loop_topology(
            topology_name,
            ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4},
            mass_map={'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'p4': 0.0},
            loop_momenta_names=('p4',),
            analytic_result=0.007968390283228252,
            ),
            entry_name = topology_name
        )

    ###############################################################################
    # PAPER: Now a series of topologies with specific external_kinematics for the paper
    ###############################################################################

    # 6P PS1
    rescaling = 1.0
    q1 = vectors.LorentzVector(
        [0.900000000000000E+01, 0.000000000000000E+00, 0.000000000000000E+00, 0.894427190999916E+01]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [0.900000000000000E+01, 0.000000000000000E+00, 0.000000000000000E+00, -0.894427190999916E+01]
    )*rescaling
    q3 = -vectors.LorentzVector(
        [0.183442509122858E+01, -0.383828222192743E+00, 0.696085529916260E+00, -0.131653190094982E+01]
    )*rescaling
    q4 = -vectors.LorentzVector(
        [0.578920098524940E+01, -0.180358221330469E+01, -0.524375913342836E+01, 0.132850645389857E+01] 
    )*rescaling
    q5 = -vectors.LorentzVector(
        [0.282869851482228E+01, -0.183886113889963E+01, -0.169694775511281E+01, 0.860519213045309E+00]            
    )*rescaling
    q6 = -q5-q4-q3-q2-q1

    # 1L
    factory = TopologyGenerator([
        ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4), ('p4', 4, 5), ('p5', 5, 6), ('p6', 6, 1),
        ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6)
    ])

    mass=0.
    topology_name = "1L_6P_PS1"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass}, 
                loop_momenta_names=('p1',), 
                analytic_result=complex(-6.04399581694419651E-009, +6.96338774260588204E-008),
            ),
            entry_name = topology_name
        )


    mass=0.25
    topology_name = "1L_6P_PS1_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass}, 
                loop_momenta_names=('p1',), 
                analytic_result=complex(2.19936376909663045E-008, +6.37930575804718533E-008)
            ),
            entry_name = topology_name
        )

    # 2L
    factory = TopologyGenerator([
        ('p1', 1, 7), ('p2', 7, 8), ('p3', 8, 1),
        ('p4', 7, 2), ('p5', 2, 3), ('p6', 3, 4), ('p7', 4, 5), ('p8', 5, 6), ('p9', 6, 8),
        ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6)
    ])

    mass=0.
    topology_name = "2L_6P_A_PS1"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass}, 
                loop_momenta_names=('p1','p4'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )


    mass=0.25
    topology_name = "2L_6P_A_PS1_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass}, 
                loop_momenta_names=('p1','p4'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )


    rescaling = 1.0
    q1 = vectors.LorentzVector(
        [0.314761904761905E+02, 0.000000000000000E+00, 0.000000000000000E+00, 0.314603014431430E+02]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [0.315238095238095E+02, 0.000000000000000E+00, 0.000000000000000E+00, -0.314603014431430E+02]
    )*rescaling
    q3 = -vectors.LorentzVector(
        [0.608892674669345E+01, -0.132242436194504E+01, 0.239826153871624E+01, -0.453591934732614E+01]
    )*rescaling
    q4 = -vectors.LorentzVector(
        [0.200491020189906E+02, -0.621398042077058E+01, -0.180666100752103E+02, 0.457717592937829E+01] 
    )*rescaling
    q5 = -vectors.LorentzVector(
        [0.103976683309230E+02, -0.633552883220179E+01, -0.584658689110735E+01, 0.296479389856186E+01]            
    )*rescaling
    q6 = -q5-q4-q3-q2-q1

    # 1L
    factory = TopologyGenerator([
        ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4), ('p4', 4, 5), ('p5', 5, 6), ('p6', 6, 1),
        ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6)
    ])

    mass=0.
    topology_name = "1L_6P_PS2"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass}, 
                loop_momenta_names=('p1',), 
                analytic_result=complex(-4.06595380412787312E-013, +2.51955969166405794E-012),
            ),
            entry_name = topology_name
        )


    mass=1.0
    topology_name = "1L_6P_PS2_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass}, 
                loop_momenta_names=('p1',), 
                analytic_result=complex(1.27623847914047270E-012, + 2.26088552166129903E-012)
            ),
            entry_name = topology_name
        )

    # 2L
    factory = TopologyGenerator([
        ('p1', 1, 7), ('p2', 7, 8), ('p3', 8, 1),
        ('p4', 7, 2), ('p5', 2, 3), ('p6', 3, 4), ('p7', 4, 5), ('p8', 5, 6), ('p9', 6, 8),
        ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6)
    ])

    mass=0.
    topology_name = "2L_6P_A_PS2"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass}, 
                loop_momenta_names=('p1','p4'), 
                analytic_result=None,
            ),
            entry_name = topology_name
        )

    mass=1.
    topology_name = "2L_6P_A_PS2_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass}, 
                loop_momenta_names=('p1','p4'), 
                analytic_result=None,
            ),
            entry_name = topology_name
        )

    factory = TopologyGenerator([
        ('p1', 2, 7), ('p2', 7, 8), ('p3', 8, 1), ('p4', 1, 2),
        ('p5', 7, 3), ('p6', 3, 4), ('p7', 4, 5), ('p8', 5, 6), ('p9', 6, 8), 
        ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6)
    ])

    rescaling = 1.0
    q1 = vectors.LorentzVector(
        [0.900000000000000E+01, 0.000000000000000E+00, 0.000000000000000E+00, 0.894427190999916E+01]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [0.900000000000000E+01, 0.000000000000000E+00, 0.000000000000000E+00, -0.894427190999916E+01]
    )*rescaling
    q3 = -vectors.LorentzVector(
        [0.183442509122858E+01, -0.383828222192743E+00, 0.696085529916260E+00, -0.131653190094982E+01]
    )*rescaling
    q4 = -vectors.LorentzVector(
        [0.578920098524940E+01, -0.180358221330469E+01, -0.524375913342836E+01, 0.132850645389857E+01] 
    )*rescaling
    q5 = -vectors.LorentzVector(
        [0.282869851482228E+01, -0.183886113889963E+01, -0.169694775511281E+01, 0.860519213045309E+00]            
    )*rescaling
    q6 = -q5-q4-q3-q2-q1

    mass=0.
    topology_name = "2L_6P_B_PS1"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass}, 
                loop_momenta_names=('p1','p5'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )


    mass=0.25
    topology_name = "2L_6P_B_PS1_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass}, 
                loop_momenta_names=('p1','p5'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )

    rescaling = 1.0
    q1 = vectors.LorentzVector(
        [0.314761904761905E+02, 0.000000000000000E+00, 0.000000000000000E+00, 0.314603014431430E+02]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [0.315238095238095E+02, 0.000000000000000E+00, 0.000000000000000E+00, -0.314603014431430E+02]
    )*rescaling
    q3 = -vectors.LorentzVector(
        [0.608892674669345E+01, -0.132242436194504E+01, 0.239826153871624E+01, -0.453591934732614E+01]
    )*rescaling
    q4 = -vectors.LorentzVector(
        [0.200491020189906E+02, -0.621398042077058E+01, -0.180666100752103E+02, 0.457717592937829E+01] 
    )*rescaling
    q5 = -vectors.LorentzVector(
        [0.103976683309230E+02, -0.633552883220179E+01, -0.584658689110735E+01, 0.296479389856186E+01]            
    )*rescaling
    q6 = -q5-q4-q3-q2-q1

    mass=0.
    topology_name = "2L_6P_B_PS2"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass}, 
                loop_momenta_names=('p1','p5'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )

    mass=1.
    topology_name = "2L_6P_B_PS2_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass}, 
                loop_momenta_names=('p1','p5'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )

    factory = TopologyGenerator([
        ('p1', 3, 7), ('p2', 7, 8), ('p3', 8, 1), ('p4', 1, 2), ('p5', 2, 3),
        ('p6', 7, 4), ('p7', 4, 5), ('p8', 5, 6), ('p9', 6, 8),
        ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6)
    ])

    rescaling = 1.0
    q1 = vectors.LorentzVector(
        [0.900000000000000E+01, 0.000000000000000E+00, 0.000000000000000E+00, 0.894427190999916E+01]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [0.900000000000000E+01, 0.000000000000000E+00, 0.000000000000000E+00, -0.894427190999916E+01]
    )*rescaling
    q3 = -vectors.LorentzVector(
        [0.183442509122858E+01, -0.383828222192743E+00, 0.696085529916260E+00, -0.131653190094982E+01]
    )*rescaling
    q4 = -vectors.LorentzVector(
        [0.578920098524940E+01, -0.180358221330469E+01, -0.524375913342836E+01, 0.132850645389857E+01] 
    )*rescaling
    q5 = -vectors.LorentzVector(
        [0.282869851482228E+01, -0.183886113889963E+01, -0.169694775511281E+01, 0.860519213045309E+00]            
    )*rescaling
    q6 = -q5-q4-q3-q2-q1

    mass=0.
    topology_name = "2L_6P_C_PS1"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass}, 
                loop_momenta_names=('p1','p6'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )

    mass=0.25
    topology_name = "2L_6P_C_PS1_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass}, 
                loop_momenta_names=('p1','p6'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )

    rescaling = 1.0
    q1 = vectors.LorentzVector(
        [0.314761904761905E+02, 0.000000000000000E+00, 0.000000000000000E+00, 0.314603014431430E+02]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [0.315238095238095E+02, 0.000000000000000E+00, 0.000000000000000E+00, -0.314603014431430E+02]
    )*rescaling
    q3 = -vectors.LorentzVector(
        [0.608892674669345E+01, -0.132242436194504E+01, 0.239826153871624E+01, -0.453591934732614E+01]
    )*rescaling
    q4 = -vectors.LorentzVector(
        [0.200491020189906E+02, -0.621398042077058E+01, -0.180666100752103E+02, 0.457717592937829E+01] 
    )*rescaling
    q5 = -vectors.LorentzVector(
        [0.103976683309230E+02, -0.633552883220179E+01, -0.584658689110735E+01, 0.296479389856186E+01]            
    )*rescaling
    q6 = -q5-q4-q3-q2-q1

    mass=0.
    topology_name = "2L_6P_C_PS2"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass}, 
                loop_momenta_names=('p1','p6'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )


    mass=1.
    topology_name = "2L_6P_C_PS2_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass}, 
                loop_momenta_names=('p1','p6'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )



    factory = TopologyGenerator([
        ('p1', 1, 7), ('p2', 7, 6), ('p3', 6, 8), ('p4', 8, 1),
        ('p5', 7, 2), ('p6', 2, 3), ('p7', 3, 4), ('p8', 4, 5), ('p9', 5, 8),
        ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6)
    ])
    rescaling = 1.0
    q1 = vectors.LorentzVector(
        [0.900000000000000E+01, 0.000000000000000E+00, 0.000000000000000E+00, 0.894427190999916E+01]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [0.900000000000000E+01, 0.000000000000000E+00, 0.000000000000000E+00, -0.894427190999916E+01]
    )*rescaling
    q3 = -vectors.LorentzVector(
        [0.183442509122858E+01, -0.383828222192743E+00, 0.696085529916260E+00, -0.131653190094982E+01]
    )*rescaling
    q4 = -vectors.LorentzVector(
        [0.578920098524940E+01, -0.180358221330469E+01, -0.524375913342836E+01, 0.132850645389857E+01] 
    )*rescaling
    q5 = -vectors.LorentzVector(
        [0.282869851482228E+01, -0.183886113889963E+01, -0.169694775511281E+01, 0.860519213045309E+00]            
    )*rescaling
    q6 = -q5-q4-q3-q2-q1

    mass=0.
    topology_name = "2L_6P_D_PS1"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass}, 
                loop_momenta_names=('p1','p5'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )


    mass=0.25
    topology_name = "2L_6P_D_PS1_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass}, 
                loop_momenta_names=('p1','p5'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )



    rescaling = 1.0
    q1 = vectors.LorentzVector(
        [0.314761904761905E+02, 0.000000000000000E+00, 0.000000000000000E+00, 0.314603014431430E+02]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [0.315238095238095E+02, 0.000000000000000E+00, 0.000000000000000E+00, -0.314603014431430E+02]
    )*rescaling
    q3 = -vectors.LorentzVector(
        [0.608892674669345E+01, -0.132242436194504E+01, 0.239826153871624E+01, -0.453591934732614E+01]
    )*rescaling
    q4 = -vectors.LorentzVector(
        [0.200491020189906E+02, -0.621398042077058E+01, -0.180666100752103E+02, 0.457717592937829E+01] 
    )*rescaling
    q5 = -vectors.LorentzVector(
        [0.103976683309230E+02, -0.633552883220179E+01, -0.584658689110735E+01, 0.296479389856186E+01]            
    )*rescaling
    q6 = -q5-q4-q3-q2-q1

    mass=0.
    topology_name = "2L_6P_D_PS2"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass}, 
                loop_momenta_names=('p1','p5'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )

    mass=1.
    topology_name = "2L_6P_D_PS2_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass}, 
                loop_momenta_names=('p1','p5'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )

    factory = TopologyGenerator([
        ('p1', 1, 7), ('p2', 7, 6), ('p3', 6, 5), ('p4', 5, 8), ('p5', 8, 1),
        ('p6', 7, 2), ('p7', 2, 3), ('p8', 3, 4), ('p9', 4, 8),
        ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6)
    ])

    rescaling = 1.0
    q1 = vectors.LorentzVector(
        [0.900000000000000E+01, 0.000000000000000E+00, 0.000000000000000E+00, 0.894427190999916E+01]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [0.900000000000000E+01, 0.000000000000000E+00, 0.000000000000000E+00, -0.894427190999916E+01]
    )*rescaling
    q3 = -vectors.LorentzVector(
        [0.183442509122858E+01, -0.383828222192743E+00, 0.696085529916260E+00, -0.131653190094982E+01]
    )*rescaling
    q4 = -vectors.LorentzVector(
        [0.578920098524940E+01, -0.180358221330469E+01, -0.524375913342836E+01, 0.132850645389857E+01] 
    )*rescaling
    q5 = -vectors.LorentzVector(
        [0.282869851482228E+01, -0.183886113889963E+01, -0.169694775511281E+01, 0.860519213045309E+00]            
    )*rescaling
    q6 = -q5-q4-q3-q2-q1

    mass=0.
    topology_name = "2L_6P_E_PS1"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass}, 
                loop_momenta_names=('p1','p6'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )


    mass=0.25
    topology_name = "2L_6P_E_PS1_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass}, 
                loop_momenta_names=('p1','p6'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )


    rescaling = 1.0
    q1 = vectors.LorentzVector(
        [0.314761904761905E+02, 0.000000000000000E+00, 0.000000000000000E+00, 0.314603014431430E+02]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [0.315238095238095E+02, 0.000000000000000E+00, 0.000000000000000E+00, -0.314603014431430E+02]
    )*rescaling
    q3 = -vectors.LorentzVector(
        [0.608892674669345E+01, -0.132242436194504E+01, 0.239826153871624E+01, -0.453591934732614E+01]
    )*rescaling
    q4 = -vectors.LorentzVector(
        [0.200491020189906E+02, -0.621398042077058E+01, -0.180666100752103E+02, 0.457717592937829E+01] 
    )*rescaling
    q5 = -vectors.LorentzVector(
        [0.103976683309230E+02, -0.633552883220179E+01, -0.584658689110735E+01, 0.296479389856186E+01]            
    )*rescaling
    q6 = -q5-q4-q3-q2-q1

    mass=0.
    topology_name = "2L_6P_E_PS2"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass}, 
                loop_momenta_names=('p1','p6'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )


    mass=1.
    topology_name = "2L_6P_E_PS2_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass}, 
                loop_momenta_names=('p1','p6'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )

    factory = TopologyGenerator([
        ('p1', 2, 7), ('p2', 7, 6), ('p3', 6, 5), ('p4', 5, 8), ('p5', 8, 1), ('p6', 1, 2),
        ('p7', 7, 3), ('p8', 3, 4), ('p9', 4, 8),
        ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6)
    ])

    rescaling = 1.0
    q1 = vectors.LorentzVector(
        [0.900000000000000E+01, 0.000000000000000E+00, 0.000000000000000E+00, 0.894427190999916E+01]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [0.900000000000000E+01, 0.000000000000000E+00, 0.000000000000000E+00, -0.894427190999916E+01]
    )*rescaling
    q3 = -vectors.LorentzVector(
        [0.183442509122858E+01, -0.383828222192743E+00, 0.696085529916260E+00, -0.131653190094982E+01]
    )*rescaling
    q4 = -vectors.LorentzVector(
        [0.578920098524940E+01, -0.180358221330469E+01, -0.524375913342836E+01, 0.132850645389857E+01] 
    )*rescaling
    q5 = -vectors.LorentzVector(
        [0.282869851482228E+01, -0.183886113889963E+01, -0.169694775511281E+01, 0.860519213045309E+00]            
    )*rescaling
    q6 = -q5-q4-q3-q2-q1

    mass=0.
    topology_name = "2L_6P_F_PS1"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass}, 
                loop_momenta_names=('p1','p7'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )


    mass=0.25
    topology_name = "2L_6P_F_PS1_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass}, 
                loop_momenta_names=('p1','p7'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )

    rescaling = 1.0
    q1 = vectors.LorentzVector(
        [0.314761904761905E+02, 0.000000000000000E+00, 0.000000000000000E+00, 0.314603014431430E+02]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [0.315238095238095E+02, 0.000000000000000E+00, 0.000000000000000E+00, -0.314603014431430E+02]
    )*rescaling
    q3 = -vectors.LorentzVector(
        [0.608892674669345E+01, -0.132242436194504E+01, 0.239826153871624E+01, -0.453591934732614E+01]
    )*rescaling
    q4 = -vectors.LorentzVector(
        [0.200491020189906E+02, -0.621398042077058E+01, -0.180666100752103E+02, 0.457717592937829E+01] 
    )*rescaling
    q5 = -vectors.LorentzVector(
        [0.103976683309230E+02, -0.633552883220179E+01, -0.584658689110735E+01, 0.296479389856186E+01]            
    )*rescaling
    q6 = -q5-q4-q3-q2-q1

    mass=0.
    topology_name = "2L_6P_F_PS2"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass}, 
                loop_momenta_names=('p1','p7'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )

    mass=1.
    topology_name = "2L_6P_F_PS2_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass}, 
                loop_momenta_names=('p1','p7'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )

    rescaling = 1.0
    q1 = vectors.LorentzVector(
        [0.120000000000000E+02, 0.000000000000000E+00, 0.000000000000000E+00, 0.119582607431014E+02]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [0.120000000000000E+02, 0.000000000000000E+00, 0.000000000000000E+00, -0.119582607431014E+02]
    )*rescaling
    q3 = -vectors.LorentzVector(
        [0.219782618596980E+01, -0.188136703758383E+00, 0.827510477484041E+00, -0.176359602349507E+01]
    )*rescaling
    q4 = -vectors.LorentzVector(
        [0.427024821284048E+01, -0.132396937288338E+01, -0.393263718110345E+01, -0.128412229700889E+00] 
    )*rescaling
    q5 = -vectors.LorentzVector(
        [0.214335972215852E+01, -0.150816191066566E+01, -0.114171719555026E+01, 0.126176048335524E+00]            
    )*rescaling
    q6 = -vectors.LorentzVector(
        [0.878601318026432E+01, 0.435187618706639E+01, 0.668261360838255E+01, -0.354934875960129E+01]            
    )*rescaling
    q7 = -vectors.LorentzVector(
        [0.159371582791415E+01, -0.712394954938346E+00, -0.181410093648943E+00, 0.999756943584162E+00]            
    )*rescaling
    q8 = -q7-q6-q5-q4-q3-q2-q1

    # 1L
    factory = TopologyGenerator([
        ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4), ('p4', 4, 5), ('p5', 5, 6), ('p6', 6, 7), ('p7', 7, 8), ('p8', 8, 1),
        ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6), ('q7', 107,7), ('q8', 108,8)
    ])

    mass=0.
    topology_name = "1L_8P_PS1"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6, 'q7': q7, 'q8':q8}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass}, 
                loop_momenta_names=('p1',), 
                analytic_result=complex(4.20915144753003114E-012, -1.95288751055433100E-012),
            ),
            entry_name = topology_name
        )


    mass=0.25
    topology_name = "1L_8P_PS1_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6, 'q7': q7, 'q8':q8}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass}, 
                loop_momenta_names=('p1',), 
                analytic_result=complex(-1.14718008873370900E-012, -2.70586855351436128E-012),
            ),
            entry_name = topology_name
        )

    # 2L
    factory = TopologyGenerator([
        ('p1', 10, 1), ('p2', 1, 2), ('p3', 2, 11), ('p4', 11, 3),  
        ('p5', 3, 4), ('p6', 4, 10), ('p7', 10, 5), ('p8', 5, 6),
        ('p9', 6, 7), ('p10', 7, 8), ('p11', 8, 11), 
        ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), 
        ('q6', 106, 6), ('q7', 107, 7), ('q8', 108, 8)
    ])

    mass=0.
    topology_name = "2L_8P_PS1"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6, 'q7':q7, 'q8':q8}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass, 'p10': mass, 'p11': mass}, 
                loop_momenta_names=('p1','p4'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )

    mass=0.25
    topology_name = "2L_8P_PS1_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6, 'q7': q7, 'q8':q8}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass, 'p10': mass, 'p11': mass}, 
                loop_momenta_names=('p1','p4'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )

    q1 = vectors.LorentzVector(
        [0.539861111111111E+02, 0.000000000000000E+00, 0.000000000000000E+00, 0.539768486751610E+02]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [0.540138888888889E+02, 0.000000000000000E+00, 0.000000000000000E+00, -0.539768486751610E+02]
    )*rescaling
    q3 = -vectors.LorentzVector(
        [0.909124865343181E+01, -0.824969594628605E+00, 0.362859011305769E+01, -0.773327621631893E+01]
    )*rescaling
    q4 = -vectors.LorentzVector(
        [0.186384260049275E+02, -0.580553637343942E+01, -0.172444081155096E+02, -0.563080903223156E+00] 
    )*rescaling
    q5 = -vectors.LorentzVector(
        [0.970074302406015E+01, -0.661321100679014E+01, -0.500637011905572E+01, 0.553275365028595E+00]            
    )*rescaling
    q6 = -vectors.LorentzVector(
        [0.387432693086229E+02, 0.190827492041572E+02, 0.293029107528482E+02, -0.155637084572518E+02]            
    )*rescaling
    q7 = -vectors.LorentzVector(
        [0.886619512838972E+01, -0.312381457445814E+01, -0.795473761522395E+00, 0.438388184761108E+01]            
    )*rescaling
    q8 = -q7-q6-q5-q4-q3-q2-q1

    # 1L
    factory = TopologyGenerator([
        ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4), ('p4', 4, 5), ('p5', 5, 6), ('p6', 6, 7), ('p7', 7, 8), ('p8', 8, 1),
        ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), ('q6', 106,6), ('q7', 107,7), ('q8', 108,8)
    ])

    mass=0.
    topology_name = "1L_8P_PS2"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6, 'q7': q7, 'q8':q8}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass}, 
                loop_momenta_names=('p1',), 
                analytic_result=complex(1.27379153209995714E-019, -8.25671397891757203E-020),
            ),
            entry_name = topology_name
        )

    mass=1.0
    topology_name = "1L_8P_PS2_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6, 'q7': q7, 'q8':q8}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass}, 
                loop_momenta_names=('p1',), 
                analytic_result=complex(-5.68050041076084999E-021, -4.11546542907940264E-020),
            ),
            entry_name = topology_name
        )

    # 2L
    factory = TopologyGenerator([
        ('p1', 10, 1), ('p2', 1, 2), ('p3', 2, 11), ('p4', 11, 3),  
        ('p5', 3, 4), ('p6', 4, 10), ('p7', 10, 5), ('p8', 5, 6),
        ('p9', 6, 7), ('p10', 7, 8), ('p11', 8, 11), 
        ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5), 
        ('q6', 106, 6), ('q7', 107, 7), ('q8', 108, 8)
    ])

    mass=0.
    topology_name = "2L_8P_PS2"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6, 'q7':q7, 'q8':q8}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass, 'p10': mass, 'p11': mass}, 
                loop_momenta_names=('p1','p4'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )

    mass=1.0
    topology_name = "2L_8P_PS2_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5, 'q6':q6, 'q7':q7, 'q8':q8}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 'p9': mass, 'p10': mass, 'p11': mass}, 
                loop_momenta_names=('p1','p4'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )

    rescaling = 1.0

    q1 = vectors.LorentzVector(
        [0.750000000000000E+01, 0.000000000000000E+00, 0.000000000000000E+00, 0.743303437365925E+01]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [0.750000000000000E+01, 0.000000000000000E+00, 0.000000000000000E+00, -0.743303437365925E+01]
    )*rescaling
    q3 = -vectors.LorentzVector(
        [0.679092057543442E+01, 0.248201303108541E+01, 0.556085824073040E+01, -0.283426696023532E+01]
    )*rescaling
    q4 = -vectors.LorentzVector(
        [0.542550551344167E+01, -0.268481025976676E+00, -0.509288997436551E+01, 0.155772260819559E+01]            
    )*rescaling
    q5 = -q4-q3-q2-q1

    # 1L
    factory = TopologyGenerator([
        ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4), ('p4', 4, 5), ('p5', 5, 1),
        ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5),
    ])

    mass=0.
    topology_name = "1L_5P_PS1"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5 }, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass }, 
                loop_momenta_names=('p1',), 
                analytic_result=complex(-1.51075539128628032E-007, -1.80678677132336518E-006),
            ),
            entry_name = topology_name
        )

    mass=0.25
    topology_name = "1L_5P_PS1_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5 }, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass },  
                loop_momenta_names=('p1',), 
                analytic_result=complex(-4.83026877573232252E-007, -3.27695172838049434E-006),
            ),
            entry_name = topology_name
        )

    # 2L
    factory = TopologyGenerator([
        ('p1', 1, 6), ('p2', 6, 7), ('p3', 7, 2), ('p4', 2, 1),  
        ('p5', 7, 3), ('p6', 3, 4), ('p7', 4, 5), ('p8', 5, 6),
        ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5),
    ])

    topology_name = "2L_5P_Planar_PS1"
    mass=0.
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass}, 
                loop_momenta_names=('p4','p5'), # If not specified an arbitrary spanning tree will be used for momentum routing 
                analytic_result=None
            ),
            entry_name = topology_name
        )

    topology_name = "2L_5P_Planar_PS1_massive"
    mass=0.25
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass}, 
                loop_momenta_names=('p4','p5'), # If not specified an arbitrary spanning tree will be used for momentum routing 
                analytic_result=None
            ),
            entry_name = topology_name
        )


    q1 = vectors.LorentzVector(
        [0.224666666666667E+02, 0.000000000000000E+00, 0.000000000000000E+00, 0.224444004400009E+02]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [0.225333333333333E+02, 0.000000000000000E+00, 0.000000000000000E+00, -0.224444004400009E+02]
    )*rescaling
    q3 = -vectors.LorentzVector(
        [0.198399446430962E+02, 0.724692418432910E+01, 0.162364651456933E+02, -0.827542705123993E+01]
    )*rescaling
    q4 = -vectors.LorentzVector(
        [0.160754655359130E+02, -0.783904683744953E+00, -0.148701022360130E+02, 0.454820240684721E+01]            
    )*rescaling
    q5 = -q4-q3-q2-q1

    # 1L
    factory = TopologyGenerator([
        ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4), ('p4', 4, 5), ('p5', 5, 1),
        ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5),
    ])

    mass=0.
    topology_name = "1L_5P_PS2"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5 }, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass }, 
                loop_momenta_names=('p1',), 
                analytic_result=complex(-6.62401397282335361E-010, -1.23531216614665428E-009),
            ),
            entry_name = topology_name
        )

    mass=1.0
    topology_name = "1L_5P_PS2_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5 }, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass },  
                loop_momenta_names=('p1',), 
                analytic_result=complex(-1.21496673675206204E-009, -1.53816973650855041E-009),
            ),
            entry_name = topology_name
        )

    # 2L
    factory = TopologyGenerator([
        ('p1', 1, 6), ('p2', 6, 7), ('p3', 7, 2), ('p4', 2, 1),  
        ('p5', 7, 3), ('p6', 3, 4), ('p7', 4, 5), ('p8', 5, 6),
        ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5),
    ])

    topology_name = "2L_5P_Planar_PS2"
    mass=0.
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass}, 
                loop_momenta_names=('p4','p5'), # If not specified an arbitrary spanning tree will be used for momentum routing 
                analytic_result=None
            ),
            entry_name = topology_name
        )

    topology_name = "2L_5P_Planar_PS2_massive"
    mass=1.
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass}, 
                loop_momenta_names=('p4','p5'), # If not specified an arbitrary spanning tree will be used for momentum routing 
                analytic_result=None
            ),
            entry_name = topology_name
        )

    rescaling = 1.0e0
    factory = TopologyGenerator([
        ('p1', 1, 6), ('p2', 6, 7), ('p3', 7, 2), ('p4', 2, 1),
        ('p5', 7, 8), ('p6', 8, 9), ('p7', 9, 6),
        ('p8', 8, 3), ('p9', 3, 4), ('p10', 4, 5), ('p11', 5, 9),
        ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5),
    ])
    q1 = vectors.LorentzVector(
        [0.149500000000000E+01,    0.000000000000000E+00,    0.000000000000000E+00,    0.149165176901313E+01]
    )*rescaling
    q2 = vectors.LorentzVector(
        [0.150500000000000E+01,    0.000000000000000E+00,    0.000000000000000E+00,   -0.149165176901313E+01]
    )*rescaling
    q3 = vectors.LorentzVector(
        [-0.126041949101381e+01,    -0.452362952912639e+00,    -0.101350243653045e+01,   0.516563513332600e+00]
    )*rescaling
    q4 = vectors.LorentzVector(
        [-0.105098730574850e+01,   0.489324061520790e-01,   0.928212188578101e+00,    -0.283905035967510e+00]
    )*rescaling
    q5 = -q4-q3-q2-q1

    topology_name = "3L_5P_Planar_PS1"
    mass=0.
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name,
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5},
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 
                        'p7': mass, 'p8': mass, 'p9': mass, 'p10': mass, 'p11': mass},
                loop_momenta_names=('p3','p5','p8'), # If not specified an arbitrary spanning tree will be used for momentum routing
                analytic_result=None
            ),
            entry_name = topology_name
        )

    topology_name = "3L_5P_Planar_PS1_massive"
    mass=0.25
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name,
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5},
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 
                        'p7': mass, 'p8': mass, 'p9': mass, 'p10': mass, 'p11': mass},
                loop_momenta_names=('p3','p5','p8'), # If not specified an arbitrary spanning tree will be used for momentum routing
                analytic_result=None
            ),
            entry_name = topology_name
        )

    q1 = vectors.LorentzVector(
        [0.224666666666667E+02, 0.000000000000000E+00, 0.000000000000000E+00, 0.224444004400009E+02]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [0.225333333333333E+02, 0.000000000000000E+00, 0.000000000000000E+00, -0.224444004400009E+02]
    )*rescaling
    q3 = -vectors.LorentzVector(
        [0.198399446430962E+02, 0.724692418432910E+01, 0.162364651456933E+02, -0.827542705123993E+01]
    )*rescaling
    q4 = -vectors.LorentzVector(
        [0.160754655359130E+02, -0.783904683744953E+00, -0.148701022360130E+02, 0.454820240684721E+01]            
    )*rescaling
    q5 = -q4-q3-q2-q1


    topology_name = "3L_5P_Planar_PS2"
    mass=0.
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name,
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5},
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 
                        'p7': mass, 'p8': mass, 'p9': mass, 'p10': mass, 'p11': mass},
                loop_momenta_names=('p3','p5','p8'), # If not specified an arbitrary spanning tree will be used for momentum routing
                analytic_result=None
            ),
            entry_name = topology_name
        )



    topology_name = "3L_5P_Planar_PS2_massive"
    mass=1.
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name,
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5},
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 
                        'p7': mass, 'p8': mass, 'p9': mass, 'p10': mass, 'p11': mass},
                loop_momenta_names=('p3','p5','p8'), # If not specified an arbitrary spanning tree will be used for momentum routing
                analytic_result=None
            ),
            entry_name = topology_name
        )

    mass=0
    rescaling = 1

    q1 = vectors.LorentzVector(
        [6.0,0.,0.,5.91607978309962]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [6.00000000000000, 0.000000000000000,  0.000000000000000,  -5.91607978309962  ]
    )*rescaling
    q3 = -vectors.LorentzVector(
        [0.600000000000000E+01, 0.131247383330590E+01, 0.526330888118183E+01, -0.236114210884473E+01] 
    )*rescaling
    q4 = -q1-q2-q3

    # 1L
    factory = TopologyGenerator([
        ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4), ('p4', 4, 1),
        ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4)
    ])

    mass=0.
    topology_name = "1L_4P_PS1"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4 }, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass }, 
                loop_momenta_names=('p1',), 
                analytic_result=complex(5.71928373468859384E-005, -7.24005031692044203E-005),
            ),
            entry_name = topology_name
        )

    mass=0.25
    topology_name = "1L_4P_PS1_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4 }, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass }, 
                loop_momenta_names=('p1',), 
                analytic_result=complex(3.02700896981466728E-005, -1.08125136760749764E-004),
            ),
            entry_name = topology_name
        )

    # 2L
    factory = TopologyGenerator([
        ('p1', 2, 5), ('p2', 5, 6), ('p3', 6, 1), ('p4', 1, 2),
        ('p5', 5, 3), ('p6', 3, 4), ('p7', 4, 6),
        ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4)
    ])
    mass=0.
    topology_name = "2L_4P_Ladder_PS1"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4,}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass}, 
                loop_momenta_names=('p1','p5'), 
                # Complex[3.110526761263476`*^-8,9.538852709003075`*^-8]
                analytic_result=analytic_four_point_ladder(q1.square(),q2.square(),q3.square(),q4.square(),(q1+q2).square(),(q2+q3).square(),2)
            ),
            entry_name = topology_name
        )

    mass=0.25
    topology_name = "2L_4P_Ladder_PS1_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4,}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass}, 
                loop_momenta_names=('p1','p5'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )


    q1 = vectors.LorentzVector(
        [0.149500000000000E+02, 0.000000000000000E+00, 0.000000000000000E+00, 0.149165176901313E+02 ]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [0.150500000000000E+02, 0.000000000000000E+00, 0.000000000000000E+00, -0.149165176901313E+02 ]
    )*rescaling
    q3 = -vectors.LorentzVector(
        [0.148833333333333E+02, 0.323407440276709E+01, 0.129693500125724E+02, -0.581810399699641E+01] 
    )*rescaling
    q4 = -q1-q2-q3

    # 1L
    factory = TopologyGenerator([
        ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4), ('p4', 4, 1),
        ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4)
    ])
    mass=0.
    topology_name = "1L_4P_PS2"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4 }, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass }, 
                loop_momenta_names=('p1',), 
                analytic_result=complex(1.55382259557521411E-006, -2.06993924914093002E-006),
            ),
            entry_name = topology_name
        )

    mass=1.0
    topology_name = "1L_4P_PS2_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4 }, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass }, 
                loop_momenta_names=('p1',), 
                analytic_result=complex(-1.79348897991097843E-007, -2.29184464605969399E-006),
            ),
            entry_name = topology_name
        )

    # 2L
    factory = TopologyGenerator([
        ('p1', 2, 5), ('p2', 5, 6), ('p3', 6, 1), ('p4', 1, 2),
        ('p5', 5, 3), ('p6', 3, 4), ('p7', 4, 6),
        ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4)
    ])
    mass=0.
    topology_name = "2L_4P_Ladder_PS2"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4,}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass}, 
                loop_momenta_names=('p1','p5'), 
                # Complex[1.7037231361181443`*^-10,4.564972190075811`*^-10]
                analytic_result=analytic_four_point_ladder(q1.square(),q2.square(),q3.square(),q4.square(),(q1+q2).square(),(q2+q3).square(),2)
            ),
            entry_name = topology_name
        )

    mass=1.0
    topology_name = "2L_4P_Ladder_PS2_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4,}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 'p5': mass, 'p6': mass, 'p7': mass}, 
                loop_momenta_names=('p1','p5'), 
                analytic_result=None,
            ),
            entry_name = topology_name
        )

    factory = TopologyGenerator([
        ('p1', 2, 5), ('p2', 5, 6), ('p3', 6, 1), ('p4', 1, 2),
        ('p5', 5, 7), ('p6', 7, 8), ('p7', 8, 6),
        ('p8', 7, 3), ('p9', 3, 4), ('p10', 4, 8),
        ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4)
    ])
    q1 = vectors.LorentzVector(
        [6.0,0.,0.,5.91607978309962]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [6.00000000000000, 0.000000000000000,  0.000000000000000,  -5.91607978309962  ]
    )*rescaling
    q3 = -vectors.LorentzVector(
        [0.600000000000000E+01, 0.131247383330590E+01, 0.526330888118183E+01, -0.236114210884473E+01] 
    )*rescaling
    q4 = -q1-q2-q3

    mass=0.
    topology_name = "3L_4P_Ladder_PS1"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 
                        'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 
                        'p9': mass, 'p10': mass}, 
                loop_momenta_names=('p1','p5','p8'), 
                # Complex[-5.3030947706862534`*^-11,-1.0780299639363004`*^-11]
                analytic_result=analytic_four_point_ladder(q1.square(),q2.square(),q3.square(),q4.square(), (q1+q2).square(), (q2+q3).square(), 3)
            ),
            entry_name = topology_name
        )

    mass=0.25
    topology_name = "3L_4P_Ladder_PS1_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 
                        'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 
                        'p9': mass, 'p10': mass}, 
                loop_momenta_names=('p1','p5','p8'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )

    q1 = vectors.LorentzVector(
        [0.149500000000000E+02, 0.000000000000000E+00, 0.000000000000000E+00, 0.149165176901313E+02 ]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [0.150500000000000E+02, 0.000000000000000E+00, 0.000000000000000E+00, -0.149165176901313E+02 ]
    )*rescaling
    q3 = -vectors.LorentzVector(
        [0.148833333333333E+02, 0.323407440276709E+01, 0.129693500125724E+02, -0.581810399699641E+01] 
    )*rescaling
    q4 = -q1-q2-q3

    mass=0.
    topology_name = "3L_4P_Ladder_PS2"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 
                        'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 
                        'p9': mass, 'p10': mass}, 
                loop_momenta_names=('p1','p5','p8'), 
                # Complex[-4.470474117170373`*^-14,-6.638346523199431`*^-15]
                analytic_result=analytic_four_point_ladder(q1.square(),q2.square(),q3.square(),q4.square(), (q1+q2).square(), (q2+q3).square(), 3)
            ),
            entry_name = topology_name
        )

    mass=1.
    topology_name = "3L_4P_Ladder_PS2_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
                topology_name, 
                ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4}, 
                mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass, 
                        'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass, 
                        'p9': mass, 'p10': mass}, 
                loop_momenta_names=('p1','p5','p8'), 
                analytic_result=None
            ),
            entry_name = topology_name
        )

    factory = TopologyGenerator([
                    ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4), ('p4', 4, 1),
                    ('p5', 1, 5), ('p6', 5, 6), ('p7', 6, 4),
                    ('p8', 5, 7), ('p9', 7, 8), ('p10', 8, 6),
                    ('p11', 7, 9), ('p12', 9, 10), ('p13', 10, 8),
                    ('q1', 101,2), ('q2', 102,3), ('q3', 103,9), ('q4', 104,10)])
     
    q1 = vectors.LorentzVector(
        [6.0,0.,0.,5.91607978309962]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [6.00000000000000, 0.000000000000000,  0.000000000000000,  -5.91607978309962  ]
    )*rescaling
    q3 = -vectors.LorentzVector(
        [0.600000000000000E+01, 0.131247383330590E+01, 0.526330888118183E+01, -0.236114210884473E+01] 
    )*rescaling
    q4 = -q1-q2-q3

    mass=0.
    topology_name = "4L_4P_Ladder_PS1"
    if selected_topologies is None or topology_name in selected_topologies:     
        all_topologies.add_topology(factory.create_loop_topology(
            topology_name,
            ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4},
            mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass,
                'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass,
                'p9': mass, 'p10': mass, 'p11': mass, 'p12': mass, 'p13': mass},
            loop_momenta_names=('p1','p5','p8','p11'),
            # Complex[1.1020376702707358`*^-14,-1.4741193204679985`*^-14]
            analytic_result=analytic_four_point_ladder(q1.square(),q2.square(),q3.square(),q4.square(), (q1+q2).square(), (q2+q3).square(), 4),
            ),
            entry_name = topology_name
        )

    mass=0.25
    topology_name = "4L_4P_Ladder_PS1_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
            topology_name,
            ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4},
            mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass,
                'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass,
                'p9': mass, 'p10': mass, 'p11': mass, 'p12': mass, 'p13': mass},
            loop_momenta_names=('p1','p5','p8','p11'),
            analytic_result=None,
            ),
            entry_name = topology_name
        )

    q1 = vectors.LorentzVector(
        [0.149500000000000E+02, 0.000000000000000E+00, 0.000000000000000E+00, 0.149165176901313E+02 ]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [0.150500000000000E+02, 0.000000000000000E+00, 0.000000000000000E+00, -0.149165176901313E+02 ]
    )*rescaling
    q3 = -vectors.LorentzVector(
        [0.148833333333333E+02, 0.323407440276709E+01, 0.129693500125724E+02, -0.581810399699641E+01] 
    )*rescaling
    q4 = -q1-q2-q3

    mass=0.
    topology_name = "4L_4P_Ladder_PS2"
    if selected_topologies is None or topology_name in selected_topologies:                
        all_topologies.add_topology(factory.create_loop_topology(
            topology_name,
            ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4},
            mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass,
                'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass,
                'p9': mass, 'p10': mass, 'p11': mass, 'p12': mass, 'p13': mass},
            loop_momenta_names=('p1','p5','p8','p11'),
            # Complex[1.4474899563458067`*^-18,-2.188966109775833`*^-18]
            analytic_result=analytic_four_point_ladder(q1.square(),q2.square(),q3.square(),q4.square(), (q1+q2).square(), (q2+q3).square(), 4),
            ),
            entry_name = topology_name
        )

    mass=1.
    topology_name = "4L_4P_Ladder_PS2_massive"
    if selected_topologies is None or topology_name in selected_topologies:       
        all_topologies.add_topology(factory.create_loop_topology(
            topology_name,
            ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4},
            mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass,
                'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass,
                'p9': mass, 'p10': mass, 'p11': mass, 'p12': mass, 'p13': mass},
            loop_momenta_names=('p1','p5','p8','p11'),
            analytic_result=None,
            ),
            entry_name = topology_name
        )

    factory = TopologyGenerator([
        ('q1', 101, 1), ('q2', 102, 7), ('q3', 103, 9), ('q4', 104, 3),
        ('p1', 1, 2), ('p2', 2, 3),
        ('p3', 1, 4), ('p4', 2, 5), ('p5', 3, 6),
        ('p6', 4, 5), ('p7', 5, 6),
        ('p8', 4, 7), ('p9', 5, 8), ('p10', 6, 9),
        ('p11', 7, 8), ('p12', 8, 9),])

    q1 = vectors.LorentzVector(
        [6.0,0.,0.,5.91607978309962]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [6.00000000000000, 0.000000000000000,  0.000000000000000,  -5.91607978309962  ]
    )*rescaling
    q3 = -vectors.LorentzVector(
        [0.600000000000000E+01, 0.131247383330590E+01, 0.526330888118183E+01, -0.236114210884473E+01] 
    )*rescaling
    q4 = -q1-q2-q3

    mass=0.
    topology_name = "FISHNET_2x2_PS1"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
            topology_name,
            ext_mom={ 'q1': q1, 'q2': q2 , 'q3': q3, 'q4': q4 },
            mass_map={}, # no masses 
            loop_momenta_names=('p1', 'p2', 'p11', 'p12'),
            # Complex[-2.2872050115559266`*^-11,2.5813298062484836`*^-11]
            analytic_result = complex(-2.2872050115559266e-11, 2.5813298062484836e-11)
            ),
            entry_name = topology_name
        )
                                    

    mass=0.25
    topology_name = "FISHNET_2x2_PS1_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
            topology_name,
            ext_mom={ 'q1': q1, 'q2': q2 , 'q3': q3, 'q4': q4 },
            mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass,
                'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass,
                'p9': mass, 'p10': mass, 'p11': mass, 'p12': mass},
            loop_momenta_names=('p1', 'p2', 'p11', 'p12'),
            analytic_result = None
            ),
            entry_name = topology_name
        )
                                    

    q1 = vectors.LorentzVector(
        [0.149500000000000E+02, 0.000000000000000E+00, 0.000000000000000E+00, 0.149165176901313E+02 ]
    )*rescaling 
    q2 = vectors.LorentzVector(
        [0.150500000000000E+02, 0.000000000000000E+00, 0.000000000000000E+00, -0.149165176901313E+02 ]
    )*rescaling
    q3 = -vectors.LorentzVector(
        [0.148833333333333E+02, 0.323407440276709E+01, 0.129693500125724E+02, -0.581810399699641E+01] 
    )*rescaling
    q4 = -q1-q2-q3

    mass=0.
    topology_name = "FISHNET_2x2_PS2"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
            topology_name,
            ext_mom={ 'q1': q1, 'q2': q2 , 'q3': q3, 'q4': q4 },
            mass_map={}, # no masses 
            loop_momenta_names=('p1', 'p2', 'p11', 'p12'),
            analytic_result = complex(-1.8967404778606168e-14, 2.4601018170911512e-14)
            ),
            entry_name = topology_name
        )

    mass=1.
    topology_name = "FISHNET_2x2_PS2_massive"
    if selected_topologies is None or topology_name in selected_topologies:
        all_topologies.add_topology(factory.create_loop_topology(
            topology_name,
            ext_mom={ 'q1': q1, 'q2': q2 , 'q3': q3, 'q4': q4 },
            mass_map={'p1': mass, 'p2': mass, 'p3': mass, 'p4': mass,
                'p5': mass, 'p6': mass, 'p7': mass, 'p8': mass,
                'p9': mass, 'p10': mass, 'p11': mass, 'p12': mass},
            loop_momenta_names=('p1', 'p2', 'p11', 'p12'),
            analytic_result = None
            ),
            entry_name = topology_name
        )

    return all_topologies

if __name__=='__main__':
    
    args = sys.argv[1:]

    # PS1 / PS2 paper topologies
    PS1PS2_paper_topologies = [
'1L_6P_PS1',
'1L_6P_PS1_massive',
'1L_6P_PS2',
'1L_6P_PS2_massive',

'2L_6P_A_PS1',
'2L_6P_A_PS1_massive',
'2L_6P_A_PS2',
'2L_6P_A_PS2_massive',

'2L_6P_B_PS1',
'2L_6P_B_PS1_massive',
'2L_6P_B_PS2',
'2L_6P_B_PS2_massive',

'2L_6P_C_PS1',
'2L_6P_C_PS1_massive',
'2L_6P_C_PS2',
'2L_6P_C_PS2_massive',

'2L_6P_D_PS1',
'2L_6P_D_PS1_massive',
'2L_6P_D_PS2',
'2L_6P_D_PS2_massive',

'2L_6P_E_PS1',
'2L_6P_E_PS1_massive',
'2L_6P_E_PS2',
'2L_6P_E_PS2_massive',

'2L_6P_F_PS1',
'2L_6P_F_PS1_massive',
'2L_6P_F_PS2',
'2L_6P_F_PS2_massive',

'1L_8P_PS1',
'1L_8P_PS1_massive',
'1L_8P_PS2',
'1L_8P_PS2_massive',

'2L_8P_PS1',
'2L_8P_PS1_massive',
'2L_8P_PS2',
'2L_8P_PS2_massive',

'1L_5P_PS1',
'1L_5P_PS1_massive',
'1L_5P_PS2',
'1L_5P_PS2_massive',

'2L_5P_Planar_PS1',
'2L_5P_Planar_PS1_massive',
'2L_5P_Planar_PS2',
'2L_5P_Planar_PS2_massive',

'3L_5P_Planar_PS1',
'3L_5P_Planar_PS1_massive',
'3L_5P_Planar_PS2',
'3L_5P_Planar_PS2_massive',

'1L_4P_PS1',
'1L_4P_PS1_massive',
'1L_4P_PS2',
'1L_4P_PS2_massive',

'2L_4P_Ladder_PS1',
'2L_4P_Ladder_PS1_massive',
'2L_4P_Ladder_PS2',
'2L_4P_Ladder_PS2_massive',

'3L_4P_Ladder_PS1',
'3L_4P_Ladder_PS1_massive',
'3L_4P_Ladder_PS2',
'3L_4P_Ladder_PS2_massive',

'4L_4P_Ladder_PS1',
'4L_4P_Ladder_PS1_massive',
'4L_4P_Ladder_PS2',
'4L_4P_Ladder_PS2_massive',

'FISHNET_2x2_PS1',
'FISHNET_2x2_PS1_massive',
'FISHNET_2x2_PS2',
'FISHNET_2x2_PS2_massive',
    ]


    if len(args)==0:
        print("Now processing all topologies...")
        raw_list_of_selected_topologies = [None,]
    else:
        raw_list_of_selected_topologies = args
        print("Now processing the following topologies: %s"%(', '.join(args)))

    list_of_selected_topologies = []
    for topo_name in raw_list_of_selected_topologies:
        if topo_name == 'PS1PS2_paper':
            list_of_selected_topologies.append(PS1PS2_paper_topologies)
        elif topo_name == 'PS1PS2_paper_individual':
            list_of_selected_topologies.extend([[t,] for t in PS1PS2_paper_topologies])
        else:
            list_of_selected_topologies.append([topo_name,])
   
    for selected_topologies in list_of_selected_topologies:
        all_topologies = load(selected_topologies)

        # Write out topoloies one by one in separate yaml files given the time it
        # takes to generate them and their size.
        for topology_name, topology in all_topologies.items():
            TopologyCollection({topology_name:topology}).export_to(os.path.join('topologies/.', '%s.yaml'%topology_name))

