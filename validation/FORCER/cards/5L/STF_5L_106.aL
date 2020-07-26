import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_5L_106 
output_LU_scalar TEST_SCALAR_STF_5L_106 \
 --topology=\
 [('q1', 0, 1), ('p1', 1, 2), ('p2', 1, 9), ('p3', 2, 3), ('p4', 2, 4), ('p5', 3, 5), ('p6', 3, 6), ('p7', 4, 5), ('p8', 4, 7), ('p9', 5, 8), ('p10', 6, 8), ('p11', 6, 9), ('p12', 7, 8), ('p13', 7, 9), ('q2', 9, 10)] \
 --name=STF_5L_106 \
 --analytical_result=3.31741444780746e-10 \
 --externals=['q1'] \
 --lmb=None
