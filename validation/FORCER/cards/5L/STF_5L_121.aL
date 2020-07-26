import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_5L_121 
output_LU_scalar TEST_SCALAR_STF_5L_121 \
 --topology=\
 [('q1', 0, 1), ('p1', 1, 3), ('p2', 1, 4), ('p3', 2, 5), ('p4', 2, 6), ('p5', 3, 5), ('p6', 3, 7), ('p7', 4, 6), ('p8', 4, 8), ('p9', 5, 9), ('p10', 6, 9), ('p11', 7, 8), ('p12', 7, 9), ('p13', 8, 9), ('q2', 2, 10)] \
 --name=STF_5L_121 \
 --analytical_result=-3.55664748118199e-9 \
 --externals=['q1'] \
 --lmb=None
