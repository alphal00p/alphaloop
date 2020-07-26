import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_5L_13 
output_LU_scalar TEST_SCALAR_STF_5L_13 \
 --topology=\
 [('q1', 0, 5), ('p1', 1, 2), ('p2', 1, 3), ('p3', 1, 5), ('p4', 2, 3), ('p5', 2, 6), ('p6', 3, 7), ('p7', 4, 5), ('p8', 4, 6), ('p9', 4, 7), ('p10', 5, 7), ('p11', 6, 7), ('q2', 6, 8)] \
 --name=STF_5L_13 \
 --analytical_result=-2.77365561323262e-9 \
 --externals=['q1'] \
 --lmb=None
