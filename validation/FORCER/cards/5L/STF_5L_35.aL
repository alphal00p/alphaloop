import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_5L_35 
output_LU_scalar TEST_SCALAR_STF_5L_35 \
 --topology=\
 [('q1', 0, 7), ('p1', 1, 2), ('p2', 1, 3), ('p3', 1, 4), ('p4', 2, 5), ('p5', 2, 6), ('p6', 3, 5), ('p7', 3, 7), ('p8', 4, 6), ('p9', 4, 7), ('p10', 5, 8), ('p11', 6, 8), ('p12', 7, 8), ('q2', 8, 9)] \
 --name=STF_5L_35 \
 --analytical_result=9.95224334342237e-10 \
 --externals=['q1'] \
 --lmb=None
