import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_4L_4 
output_LU_scalar TEST_SCALAR_STF_4L_4 \
 --topology=\
 [('q1', 0, 1), ('p1', 1, 5), ('p2', 1, 6), ('p3', 2, 3), ('p4', 2, 4), ('p5', 2, 6), ('p6', 3, 5), ('p7', 3, 6), ('p8', 4, 5), ('p9', 4, 6), ('q2', 5, 7)] \
 --name=STF_4L_4 \
 --analytical_result=-1.04773018349817e-7 \
 --externals=['q1'] \
 --lmb=None
