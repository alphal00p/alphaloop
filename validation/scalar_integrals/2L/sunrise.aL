import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_integrand_type PF
set_alphaLoop_option FORM_compile_optimization 3
set_FORM_option extra-options {'OPTIMITERATIONS':1000}
!rm -rf TEST_SCALAR_INTEGRAL_2L_SUNRISE
output_scalar_integral TEST_SCALAR_INTEGRAL_2L_SUNRISE \
--topology=[ \
        ('q1', 101, 1), \
        ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4), ('p5', 4, 1), ('p6', 2, 4), \
        ('pout2A', 3, 103),('pout2B', 3, 104) \
    ] \
 --externals=[('q1',),('pout2A','pout2B')] \
 --default_kinematics={'q1':[3., 4., 2., 1.],'pout2A':[3./2., 4./2., 2./2., 1./2.]} \
 --name=scalar_integral_2L_sunrise \
 --analytical_result=2.6918653677981387e-14 \
 --numerator='((k1.p1+k2.p1)^2+(k1.k2))^3'
 --lmb=['p5','p6']
