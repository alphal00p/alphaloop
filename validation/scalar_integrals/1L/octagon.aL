import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_integrand_type PF
set_alphaLoop_option FORM_compile_optimization 3
set_FORM_option optimisation_strategy 'greedy'
set_FORM_option extra-options {'OPTIMITERATIONS':1000}
!rm -rf TEST_SCALAR_INTEGRAL_1L_OCTAGON_0 
output_scalar_integral TEST_SCALAR_INTEGRAL_1L_OCTAGON_0 \
 --topology=\
 [('q1', 1, 101), ('q2', 2, 102), ('q3', 3, 103), ('q4', 4, 104),\ 
  ('p1', 101, 102), ('p2', 102, 103), ('p3', 103, 104), ('p4', 104, 105), ('p5', 105, 106), ('p6', 106, 107), ('p7', 107, 108), ('p8', 108, 101),\
  ('pout5', 105, 5), ('pout6', 106, 6), ('pout7', 107, 7), ('pout8', 108, 8)] \
 --name=scalar_integral_1L_octagon \
 --analytical_result=5.10252154860455440E-007-1.54755519924903470E-006j \
 --externals=[('q1','q2','q3','q4'),('pout5','pout6','pout7','pout8')] \
 --default_kinematics={'q1':[0.411227272727273E+01,    	0.000000000000000E+00,    0.000000000000000E+00,    0.398883278459640E+01],\
                       'q2':[0.413772727272727E+01,    	0.000000000000000E+00,    0.000000000000000E+00,   -0.398883278459640E+01],\
		       'q3':[0.327181222232223E+01, 	0.102374858888E+00, 	0.000005647382777E+00, 		2.99919222229222E+00],\
		       'q4':[1.65477777377717273, 	0.83617211728, 		0.7629384777172838,		0.0722222],\
	                   'pout5':[0.130849361300178E+01,   -0.130202499338920E+00,    0.236126659032492E+00,   -0.446594945478000E+00],\
                       'pout6':[0.233055937421053E+01,   -0.611812520178803E+00,   -0.177879193250363E+01,    0.450656962370573E+00],\
                       'pout7':[0.166303088102142E+01,   -0.623779863956214E+00,   -0.575639898757384E+00,    0.291905977177980E+00]} \
 --numerator='(k1.p1+k1.p2)^0'
 --lmb=['p1',]
