import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_integrand_type LTD
set_alphaLoop_option FORM_compile_optimization 3
set_FORM_option extra-options {'OPTIMITERATIONS':0}
!rm -rf TEST_SCALAR_INTEGRAL_1L_BOX_RANK0
output_scalar_integral TEST_SCALAR_INTEGRAL_1L_BOX_RANK0 \
 --topology=\
 [('q1', 1, 101), ('q2', 2, 102),\
  ('p1', 101, 102), ('p2', 102, 103), ('p3', 103, 104), ('p4', 104, 101), \
  ('pout3', 103, 3), ('pout4', 104, 4)] \
 --name=scalar_integral_1L_box_rank0 \
 --analytical_result=3.15521080399766E-10-2.432412361669758E-10j \
 --externals=[('q1','q2'),('pout3','pout4')] \
 --default_kinematics={'q1':[62.80274, -49.71968, -5.53340, -79.44048],\
                       'q2':[48.59375, -1.65947, 34.91140, 71.89564],\
                           'pout3':[76.75934, -19.14334,-17.10279,30.22959]} \
 --lmb=['p1',]
