import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_5L_101 
output_LU_scalar TEST_SCALAR_STF_5L_101 \
 --topology=\
 [('q1', 0, 1), ('p1', 1, 2), ('p2', 1, 9), ('p3', 2, 3), ('p4', 2, 9), ('p5', 3, 4), ('p6', 3, 5), ('p7', 4, 6), ('p8', 4, 7), ('p9', 5, 6), ('p10', 5, 8), ('p11', 6, 7), ('p12', 7, 8), ('p13', 8, 9), ('q2', 9, 10)] \
 --name=STF_5L_101 \
 --analytical_result=-1.60237816225561e-11 \
 --externals=['q1'] \
 --lmb=None
