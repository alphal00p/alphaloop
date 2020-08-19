import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_integrand_type PF
set_alphaLoop_option FORM_compile_optimization 3
set_FORM_option optimisation_strategy 'CSEgreedy'
set_FORM_option FORM_setup {\
'MaxTermSize':'10M',\
'Workspace':'100G',\
'LargeSize':'20G',\
'SmallSize':'1G',\
'TermsInSmall':'10M',\
'SortIOSize':'2M',\
'ScratchSize':'10G'\
}
set_FORM_option extra-options {'OPTIMITERATIONS':1000}
!rm -rf TEST_SCALAR_INTEGRAL_2x2_FISHNET_PF
output_scalar_integral TEST_SCALAR_INTEGRAL_2x2_FISHNET_PF \
 --topology=[ \
            ('q1', 101, 1), ('q2', 102, 7), ('pout3', 9, 103), ('pout4', 3, 104), \
            ('p1', 1, 2), ('p2', 2, 3), \
            ('p3', 1, 4), ('p4', 2, 5), ('p5', 3, 6), \
            ('p6', 4, 5), ('p7', 5, 6), \
            ('p8', 4, 7), ('p9', 5, 8), ('p10', 6, 9), \
            ('p11', 7, 8), ('p12', 8, 9), \
    ] \
 --name=scalar_integral_4L_2x2_fishnet \
 --analytical_result=2.6918653677981387e-14 \
 --externals=[('q1','q2'),('pout3','pout4')] \
 --default_kinematics={'q1':[  1.2,  2.2,   1.0, 0.4 ],\
                       'q2':[  2.0, -5.2,   2.1, 0.0 ],\
	                   'pout3':[ 1.6,  -0.1, 12.5, -2.4 ] } \
 --lmb=['p3','p4','p8','p9']
