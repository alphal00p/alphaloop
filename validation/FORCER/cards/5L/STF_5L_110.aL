import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_5L_110 
output_LU_scalar TEST_SCALAR_STF_5L_110 \
 --topology=\
 [('q1', 0, 1), ('p1', 1, 3), ('p2', 1, 9), ('p3', 2, 4), ('p4', 2, 9), ('p5', 3, 5), ('p6', 3, 6), ('p7', 4, 7), ('p8', 4, 8), ('p9', 5, 7), ('p10', 5, 8), ('p11', 6, 7), ('p12', 6, 9), ('p13', 8, 9), ('q2', 2, 10)] \
 --name=STF_5L_110 \
 --analytical_result=-1.12909443847047e-9 \
 --externals=['q1'] \
 --lmb=None
