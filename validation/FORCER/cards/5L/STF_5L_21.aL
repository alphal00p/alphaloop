import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_5L_21 
output_LU_scalar TEST_SCALAR_STF_5L_21 \
 --topology=\
 [('q1', 0, 1), ('p1', 1, 2), ('p2', 1, 5), ('p3', 2, 6), ('p4', 2, 7), ('p5', 3, 5), ('p6', 3, 6), ('p7', 3, 7), ('p8', 4, 5), ('p9', 4, 6), ('p10', 4, 7), ('p11', 6, 7), ('q2', 5, 8)] \
 --name=STF_5L_21 \
 --analytical_result=-1.66419336793957e-9 \
 --externals=['q1'] \
 --lmb=None
