import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_5L_57 
output_LU_scalar TEST_SCALAR_STF_5L_57 \
 --topology=\
 [('q1', 0, 1), ('p1', 1, 2), ('p2', 1, 7), ('p3', 2, 3), ('p4', 2, 4), ('p5', 3, 5), ('p6', 3, 7), ('p7', 4, 6), ('p8', 4, 8), ('p9', 5, 6), ('p10', 5, 8), ('p11', 6, 8), ('p12', 7, 8), ('q2', 7, 9)] \
 --name=STF_5L_57 \
 --analytical_result=4.97612167171119e-10 \
 --externals=['q1'] \
 --lmb=None
