import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_5L_37 
output_LU_scalar TEST_SCALAR_STF_5L_37 \
 --topology=\
 [('q1', 0, 7), ('p1', 1, 2), ('p2', 1, 3), ('p3', 1, 4), ('p4', 2, 3), ('p5', 2, 5), ('p6', 3, 7), ('p7', 4, 5), ('p8', 4, 6), ('p9', 5, 8), ('p10', 6, 7), ('p11', 6, 8), ('p12', 7, 8), ('q2', 8, 9)] \
 --name=STF_5L_37 \
 --analytical_result=2.48806083585559e-10 \
 --externals=['q1'] \
 --lmb=None
