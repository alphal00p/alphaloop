import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_4L_17 
output_LU_scalar TEST_SCALAR_STF_4L_17 \
 --topology=\
 [('q1', 0, 1), ('p1', 1, 3), ('p2', 1, 4), ('p3', 2, 3), ('p4', 2, 5), ('p5', 3, 6), ('p6', 4, 7), ('p7', 4, 8), ('p8', 5, 7), ('p9', 5, 8), ('p10', 6, 7), ('p11', 6, 8), ('q2', 2, 9)] \
 --name=STF_4L_17 \
 --analytical_result=-5.23865091749086e-8 \
 --externals=['q1'] \
 --lmb=None
