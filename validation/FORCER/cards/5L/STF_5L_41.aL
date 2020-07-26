import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_5L_41 
output_LU_scalar TEST_SCALAR_STF_5L_41 \
 --topology=\
 [('q1', 0, 1), ('p1', 1, 2), ('p2', 1, 8), ('p3', 2, 3), ('p4', 2, 8), ('p5', 3, 7), ('p6', 3, 8), ('p7', 4, 5), ('p8', 4, 6), ('p9', 4, 7), ('p10', 5, 6), ('p11', 5, 7), ('p12', 6, 8), ('q2', 7, 9)] \
 --name=STF_5L_41 \
 --analytical_result=-2.77365561323262e-9 \
 --externals=['q1'] \
 --lmb=None
