import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_5L_108 
output_LU_scalar TEST_SCALAR_STF_5L_108 \
 --topology=\
 [('q1', 0, 1), ('p1', 1, 2), ('p2', 1, 3), ('p3', 2, 4), ('p4', 2, 5), ('p5', 3, 4), ('p6', 3, 9), ('p7', 4, 6), ('p8', 5, 6), ('p9', 5, 7), ('p10', 6, 8), ('p11', 7, 8), ('p12', 7, 9), ('p13', 8, 9), ('q2', 9, 10)] \
 --name=STF_5L_108 \
 --analytical_result=-1.12909443847047e-9 \
 --externals=['q1'] \
 --lmb=None
