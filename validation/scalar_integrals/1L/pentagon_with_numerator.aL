import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_integrand_type PF
set_alphaLoop_option FORM_compile_optimization 3
set_FORM_option extra-options {'OPTIMITERATIONS':1000}
!rm -rf TEST_SCALAR_INTEGRAL_1L_PENTAGON_NUMERATOR_0 
output_scalar_integral TEST_SCALAR_INTEGRAL_1L_PENTAGON_NUMERATOR_0 \
 --topology=\
 [('q1', 1, 101), ('q2', 2, 102),\ 
  ('p1', 101, 102), ('p2', 102, 103), ('p3', 103, 104), ('p4', 104, 105), ('p5', 105, 101),\
  ('pout3', 103, 3), ('pout4', 104, 4), ('pout5', 105, 5)] \
 --name=scalar_integral_1L_pentagon \
 --analytical_result=-1.51075539128628032E-007-1.80678677132336518E-006j \
 --externals=[('q1','q2'),('pout3','pout4','pout5')] \
 --default_kinematics={'q1':[0.750000000000000E+01, 0.000000000000000E+00, 0.000000000000000E+00, 0.743303437365925E+01],\
                       'q2':[0.750000000000000E+01, 0.000000000000000E+00, 0.000000000000000E+00, -0.743303437365925E+01],\
	                   'pout3':[0.679092057543442E+01, 0.248201303108541E+01, 0.556085824073040E+01, -0.283426696023532E+01],\
                       'pout4':[0.542550551344167E+01, -0.268481025976676E+00, -0.509288997436551E+01, 0.155772260819559E+01]} \
 --lmb=['p1',] \
 --numerator='(k1.k1+k1.p1+k1.p2)^2'
