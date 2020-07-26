import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_4L_3 
output_LU_scalar TEST_SCALAR_STF_4L_3 \
 --topology=\
 [('q1', 0, 5), ('p1', 1, 2), ('p2', 1, 3), ('p3', 1, 5), ('p4', 2, 4), ('p5', 2, 5), ('p6', 3, 4), ('p7', 3, 6), ('p8', 4, 6), ('p9', 5, 6), ('q2', 6, 7)] \
 --name=STF_4L_3 \
 --analytical_result=-1.04773018349817e-7 \
 --externals=['q1'] \
 --lmb=None
