import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_5L_15 
output_LU_scalar TEST_SCALAR_STF_5L_15 \
 --topology=\
 [('q1', 0, 1), ('p1', 1, 6), ('p2', 1, 7), ('p3', 2, 3), ('p4', 2, 4), ('p5', 2, 5), ('p6', 3, 4), ('p7', 3, 6), ('p8', 4, 7), ('p9', 5, 6), ('p10', 5, 7), ('p11', 6, 7), ('q2', 5, 8)] \
 --name=STF_5L_15 \
 --analytical_result=-2.77365561323262e-9 \
 --externals=['q1'] \
 --lmb=None
