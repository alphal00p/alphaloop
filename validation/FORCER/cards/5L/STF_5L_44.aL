import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_5L_44 
output_LU_scalar TEST_SCALAR_STF_5L_44 \
 --topology=\
 [('q1', 0, 1), ('p1', 1, 2), ('p2', 1, 8), ('p3', 2, 3), ('p4', 2, 7), ('p5', 3, 4), ('p6', 3, 8), ('p7', 4, 5), ('p8', 4, 6), ('p9', 5, 7), ('p10', 5, 8), ('p11', 6, 7), ('p12', 6, 8), ('q2', 7, 9)] \
 --name=STF_5L_44 \
 --analytical_result=-1.66419336793957e-9 \
 --externals=['q1'] \
 --lmb=None
