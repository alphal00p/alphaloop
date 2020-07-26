import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_5L_138 
output_LU_scalar TEST_SCALAR_STF_5L_138 \
 --topology=\
 [('q1', 0, 1), ('p1', 1, 2), ('p2', 1, 3), ('p3', 2, 9), ('p4', 3, 4), ('p5', 3, 5), ('p6', 4, 6), ('p7', 4, 7), ('p8', 5, 6), ('p9', 5, 9), ('p10', 6, 8), ('p11', 7, 8), ('p12', 7, 9), ('p13', 8, 9), ('q2', 2, 10)] \
 --name=STF_5L_138 \
 --analytical_result=1.65870722390373e-10 \
 --externals=['q1'] \
 --lmb=None
