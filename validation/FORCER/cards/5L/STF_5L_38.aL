import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_5L_38 
output_LU_scalar TEST_SCALAR_STF_5L_38 \
 --topology=\
 [('q1', 0, 1), ('p1', 1, 7), ('p2', 1, 8), ('p3', 2, 3), ('p4', 2, 4), ('p5', 2, 5), ('p6', 3, 6), ('p7', 3, 8), ('p8', 4, 6), ('p9', 4, 8), ('p10', 5, 7), ('p11', 5, 8), ('p12', 6, 7), ('q2', 7, 9)] \
 --name=STF_5L_38 \
 --analytical_result=4.97612167171119e-10 \
 --externals=['q1'] \
 --lmb=None
