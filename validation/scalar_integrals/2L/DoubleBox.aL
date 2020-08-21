import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_integrand_type LTD
set_alphaLoop_option FORM_compile_optimization 3
set_FORM_option extra-options {'OPTIMITERATIONS':0}
!rm -rf TEST_SCALAR_INTEGRAL_2L_DOUBLEBOX
output_scalar_integral TEST_SCALAR_INTEGRAL_2L_DOUBLEBOX \
--topology=[ \
        ('q1', 101, 1), ('q2', 102, 2), \
        ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4), ('p5', 4, 1), ('p6', 3, 5), ('p7', 5, 6), ('p8', 6, 4), \
        ('pout2A', 5, 105),('pout2B', 6, 106) \
    ] \
 --externals=[('q1', 'q2'),('pout2A','pout2B')] \
 --default_kinematics={'q1':[3., 4., 2., 1.], 'q2':[1., 1.2, 6.4, 8.9], 'pout2A':[3./2., 4./2., 2./2., 1./2.]} \
 --name=scalar_integral_2L_doublebox \
 --analytical_result=2.6918653677981387e-14 \
 --numerator='((k1.p1+k2.p1)^2+(k1.k2))^2'
 --lmb=['p1','p3']
