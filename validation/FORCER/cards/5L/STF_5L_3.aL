import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_5L_3 
output_LU_scalar TEST_SCALAR_STF_5L_3 \
 --topology=\
 [('q1', 0, 3), ('p1', 1, 3), ('p2', 1, 4), ('p3', 1, 5), ('p4', 2, 3), ('p5', 2, 5), ('p6', 2, 6), ('p7', 3, 6), ('p8', 4, 5), ('p9', 4, 6), ('p10', 5, 6), ('q2', 4, 7)] \
 --name=STF_5L_3 \
 --analytical_result=1.77832374059099e-9 \
 --externals=['q1'] \
 --lmb=None
