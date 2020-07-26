import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_5L_126 
output_LU_scalar TEST_SCALAR_STF_5L_126 \
 --topology=\
 [('q1', 0, 1), ('p1', 1, 3), ('p2', 1, 4), ('p3', 2, 3), ('p4', 2, 9), ('p5', 3, 5), ('p6', 4, 5), ('p7', 4, 6), ('p8', 5, 7), ('p9', 6, 8), ('p10', 6, 9), ('p11', 7, 8), ('p12', 7, 9), ('p13', 8, 9), ('q2', 2, 10)] \
 --name=STF_5L_126 \
 --analytical_result=-1.12909443847047e-9 \
 --externals=['q1'] \
 --lmb=None
