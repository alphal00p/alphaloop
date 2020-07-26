import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_5L_136 
output_LU_scalar TEST_SCALAR_STF_5L_136 \
 --topology=\
 [('q1', 0, 1), ('p1', 1, 2), ('p2', 1, 3), ('p3', 2, 9), ('p4', 3, 4), ('p5', 3, 9), ('p6', 4, 5), ('p7', 4, 6), ('p8', 5, 7), ('p9', 5, 8), ('p10', 6, 7), ('p11', 6, 9), ('p12', 7, 8), ('p13', 8, 9), ('q2', 2, 10)] \
 --name=STF_5L_136 \
 --analytical_result=-4.80713448676683e-11 \
 --externals=['q1'] \
 --lmb=None
