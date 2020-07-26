import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_4L_19 
output_LU_scalar TEST_SCALAR_STF_4L_19 \
 --topology=\
 [('q1', 0, 1), ('p1', 1, 2), ('p2', 1, 3), ('p3', 2, 4), ('p4', 3, 5), ('p5', 3, 6), ('p6', 4, 5), ('p7', 4, 7), ('p8', 5, 8), ('p9', 6, 7), ('p10', 6, 8), ('p11', 7, 8), ('q2', 2, 9)] \
 --name=STF_4L_19 \
 --analytical_result=1.01214966799408e-8 \
 --externals=['q1'] \
 --lmb=None
