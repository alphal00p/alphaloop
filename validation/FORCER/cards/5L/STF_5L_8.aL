import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_5L_8 
output_LU_scalar TEST_SCALAR_STF_5L_8 \
 --topology=\
 [('q1', 0, 5), ('p1', 1, 2), ('p2', 1, 5), ('p3', 1, 7), ('p4', 2, 6), ('p5', 2, 7), ('p6', 3, 4), ('p7', 3, 5), ('p8', 3, 7), ('p9', 4, 6), ('p10', 4, 7), ('p11', 5, 6), ('q2', 6, 8)] \
 --name=STF_5L_8 \
 --analytical_result=-1.66419336793957e-9 \
 --externals=['q1'] \
 --lmb=None
