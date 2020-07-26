import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_4L_6 
output_LU_scalar TEST_SCALAR_STF_4L_6 \
 --topology=\
 [('q1', 0, 1), ('p1', 1, 2), ('p2', 1, 7), ('p3', 2, 3), ('p4', 2, 4), ('p5', 3, 5), ('p6', 3, 6), ('p7', 4, 5), ('p8', 4, 7), ('p9', 5, 6), ('p10', 6, 7), ('q2', 7, 8)] \
 --name=STF_4L_6 \
 --analytical_result=3.03644900398223e-8 \
 --externals=['q1'] \
 --lmb=None
