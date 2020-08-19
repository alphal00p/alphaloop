import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_integrand_type PF
set_alphaLoop_option FORM_compile_optimization 3
set_FORM_option extra-options {'OPTIMITERATIONS':1000}
!rm -rf TEST_SCALAR_INTEGRAL_1L_HEXAGON_0 
output_scalar_integral TEST_SCALAR_INTEGRAL_1L_HEXAGON_0 \
 --topology=\
 [('q1', 1, 101), ('q2', 2, 102),\ 
  ('p1', 101, 102), ('p2', 102, 103), ('p3', 103, 104), ('p4', 104, 105), ('p5', 105, 106), ('p6', 106, 101),\
  ('pout3', 103, 3), ('pout4', 104, 4), ('pout5', 105, 5), ('pout6', 106, 6)] \
 --name=scalar_integral_1L_hexagon \
 --analytical_result=1337.0+3.0j \
 --externals=[('q1','q2'),('pout3','pout4','pout5','pout6')] \
 --default_kinematics={'q1':[0.5,0.,0.,1.0],'q2':[0.5,0.,0.,2.0],\
	                   'pout3':[0.5,0.,0.,3.0],'pout4':[0.5,0.,0.,4.0],'pout5':[0.5,0.,0.,5.0]} \
 --lmb=['p1',]
