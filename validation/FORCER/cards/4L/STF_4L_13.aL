import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_4L_13 
output_LU_scalar TEST_SCALAR_STF_4L_13 \
 --topology=\
 [('q1', 0, 1), ('p1', 1, 2), ('p2', 1, 3), ('p3', 2, 7), ('p4', 3, 4), ('p5', 3, 7), ('p6', 4, 5), ('p7', 4, 6), ('p8', 5, 6), ('p9', 5, 7), ('p10', 6, 7), ('q2', 2, 8)] \
 --name=STF_4L_13 \
 --analytical_result=3.03644900398223e-8 \
 --externals=['q1'] \
 --lmb=None
