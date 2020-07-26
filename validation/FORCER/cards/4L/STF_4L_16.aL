import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_4L_16 
output_LU_scalar TEST_SCALAR_STF_4L_16 \
 --topology=\
 [('q1', 0, 1), ('p1', 1, 3), ('p2', 1, 4), ('p3', 2, 5), ('p4', 2, 6), ('p5', 3, 5), ('p6', 3, 7), ('p7', 4, 6), ('p8', 4, 7), ('p9', 5, 8), ('p10', 6, 8), ('p11', 7, 8), ('q2', 2, 9)] \
 --name=STF_4L_16 \
 --analytical_result=-2.09546036699635e-7 \
 --externals=['q1'] \
 --lmb=None
