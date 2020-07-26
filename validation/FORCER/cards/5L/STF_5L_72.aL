import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_5L_72 
output_LU_scalar TEST_SCALAR_STF_5L_72 \
 --topology=\
 [('q1', 0, 1), ('p1', 1, 2), ('p2', 1, 3), ('p3', 2, 3), ('p4', 2, 8), ('p5', 3, 8), ('p6', 4, 5), ('p7', 4, 6), ('p8', 4, 7), ('p9', 5, 6), ('p10', 5, 7), ('p11', 6, 8), ('p12', 7, 8), ('q2', 7, 9)] \
 --name=STF_5L_72 \
 --analytical_result=-2.77365561323262e-9 \
 --externals=['q1'] \
 --lmb=None
