import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_5L_70 
output_LU_scalar TEST_SCALAR_STF_5L_70 \
 --topology=\
 [('q1', 0, 1), ('p1', 1, 2), ('p2', 1, 3), ('p3', 2, 4), ('p4', 2, 5), ('p5', 3, 6), ('p6', 3, 7), ('p7', 4, 6), ('p8', 4, 8), ('p9', 5, 7), ('p10', 5, 8), ('p11', 6, 8), ('p12', 7, 8), ('q2', 7, 9)] \
 --name=STF_5L_70 \
 --analytical_result=-1.77832374059099e-9 \
 --externals=['q1'] \
 --lmb=None
