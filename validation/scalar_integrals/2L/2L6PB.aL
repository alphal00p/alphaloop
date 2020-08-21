import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_integrand_type LTD
set_alphaLoop_option FORM_compile_optimization 3
set_FORM_option extra-options {'OPTIMITERATIONS':0}
!rm -rf TEST_SCALAR_INTEGRAL_2L_2L6PB
output_scalar_integral TEST_SCALAR_INTEGRAL_2L_2L6PB \
--topology=[ \
        ('q1', 101, 1), ('q2', 102, 2), ('q3', 103, 3), \
        ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4), ('p5', 4, 5), ('p6', 5, 1), ('p7', 4, 6), ('p8', 6, 7), ('p9', 7, 8), ('p10', 8, 5), \
        ('pout4', 6, 106), ('pout5', 7, 107), ('pout6', 8, 108) \
    ] \
 --externals=[('q1', 'q2', 'q3'),('pout4','pout5', 'pout6')] \
 --default_kinematics={'q1':[3., 4., 2., 1.], 'q2':[1., 1.2, 6.4, 8.9], 'q3':[3./2., 4./2., 2./2., 1./2.], 'pout4':[0.56, 0.21, 0.5, 0.11], 'pout5':[0.43, 0.11, 0.1,0.14]} \
 --name=scalar_integral_2L_2L6PB \
 --analytical_result=2.6918653677981387e-14 \
 --numerator='((k1.p1+k2.p1)^2+(k1.k2))^0'
 --lmb=['p1','p5']
