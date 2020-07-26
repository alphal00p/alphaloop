import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_4L_7 
output_LU_scalar TEST_SCALAR_STF_4L_7 \
 --topology=\
 [('q1', 0, 1), ('p1', 1, 2), ('p2', 1, 3), ('p3', 2, 4), ('p4', 2, 5), ('p5', 3, 6), ('p6', 3, 7), ('p7', 4, 6), ('p8', 4, 7), ('p9', 5, 6), ('p10', 5, 7), ('q2', 7, 8)] \
 --name=STF_4L_7 \
 --analytical_result=-1.04773018349817e-7 \
 --externals=['q1'] \
 --lmb=None
