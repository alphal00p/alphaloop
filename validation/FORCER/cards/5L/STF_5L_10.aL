import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_5L_10 
output_LU_scalar TEST_SCALAR_STF_5L_10 \
 --topology=\
 [('q1', 0, 5), ('p1', 1, 2), ('p2', 1, 3), ('p3', 1, 5), ('p4', 2, 4), ('p5', 2, 7), ('p6', 3, 6), ('p7', 3, 7), ('p8', 4, 6), ('p9', 4, 7), ('p10', 5, 6), ('p11', 5, 7), ('q2', 6, 8)] \
 --name=STF_5L_10 \
 --analytical_result=-1.77832374059099e-9 \
 --externals=['q1'] \
 --lmb=None
