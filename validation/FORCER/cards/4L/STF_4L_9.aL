import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_4L_9 
output_LU_scalar TEST_SCALAR_STF_4L_9 \
 --topology=\
 [('q1', 0, 1), ('p1', 1, 3), ('p2', 1, 7), ('p3', 2, 4), ('p4', 2, 7), ('p5', 3, 5), ('p6', 3, 6), ('p7', 4, 5), ('p8', 4, 6), ('p9', 5, 7), ('p10', 6, 7), ('q2', 2, 8)] \
 --name=STF_4L_9 \
 --analytical_result=-1.04773018349817e-7 \
 --externals=['q1'] \
 --lmb=None
