#! /usr/bin/env python
from pySecDec.loop_integral import loop_package
import pySecDec as psd

repl_rules = [
                        ('p1*p1', 'p1sq'),
                        ('p2*p2', 'p2sq'),
                        ('p3*p3', 'p3sq'),
                        ('p4*p4', 'p4sq'),
                        #('p5*p5', 'p5sq'),
                        ('p5*p5', 'p1sq+p2sq+p3sq+p4sq+2*p1p2+2*p1p3+2*p1p4+2*p2p3+2*p2p4+2*p3p4'),

                        ('p1*p2', 'p1p2'),
                        ('p1*p3', 'p1p3'),
                        ('p1*p4', 'p1p4'),
                        #('p1*p5', 'p1p5'),
                        ('p1*p5', '-p1sq-p1p2-p1p3-p1p4'),

                        ('p2*p3', 'p2p3'),
                        ('p2*p4', 'p2p4'),
                        #('p2*p4', 'p2p5'),
                        ('p2*p5', '-p1p2-p2sq-p2p3-p2p4'),

                        ('p3*p4', 'p3p4'),
                        #('p3*p5', 'p3p5'),
                        ('p3*p5', '-p1p3-p2p3-p3sq-p3p4'),

                        #('p4*p5', 'p4p5'),
                        ('p4*p5', '-p1p4-p2p4-p3p4-p4sq'),

                    ]

repl_rules.append(('m**2', 'msq'))

li = psd.loop_integral.LoopIntegralFromGraph(
#internal_lines = [['m',[3,4]],['m',[4,5]],['m',[3,5]],[0,[1,2]],[0,[4,1]],[0,[2,5]]],
#external_lines = [['p1',1],['p2',2],['p3',3]],

#internal_lines = [ [0,[2,1]],[0,[1,6]],[0,[6,8]],
#                   [0,[8,9]],[0,[9,7]],[0,[7,2]],
#                   [0,[7,3]],[0,[3,4]],[0,[4,6]] ],

internal_lines = [ 
    ['m',[1,6]], ['m',[6,7]], ['m',[7,2]], ['m',[2,1]], 
    ['m',[7,3]], ['m',[3,4]], ['m',[4,5]], ['m',[5,6]], 
],

external_lines = [['p1',1],['p2',2],['p3',3],['p4',4],['p5',5],],

replacement_rules = repl_rules 
)

Mandelstam_symbols = [
'p1sq',
'p2sq',
'p3sq',
'p4sq',
#'p5sq',
'p1p2',
'p1p3',
'p1p4',
#'p1p5',
'p2p3',
'p2p4',
#'p2p5',
'p3p4',
#'p3p5',
#'p4p5',
]
mass_symbols = ['msq']


loop_package(

name = 'Pentabox_physical',

loop_integral = li,

real_parameters = Mandelstam_symbols,
complex_parameters = mass_symbols,

# the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
requested_order = 0,

# the optimization level to use in FORM (can be 0, 1, 2, 3)
form_optimization_level = 2,

# the WorkSpace parameter for FORM
form_work_space = '1000M',

# the method to be used for the sector decomposition
# valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
decomposition_method = 'iterative',
# if you choose ``geometric[_ku]`` and 'normaliz' is not in your
# $PATH, you can set the path to the 'normaliz' command-line
# executable here
#normaliz_executable='/path/to/normaliz',

contour_deformation = True 
 
)
