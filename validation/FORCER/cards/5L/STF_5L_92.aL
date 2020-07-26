import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_5L_92 
output_LU_scalar TEST_SCALAR_STF_5L_92 \
 --topology=\
 [('q1', 0, 1), ('p1', 1, 3), ('p2', 1, 4), ('p3', 2, 3), ('p4', 2, 7), ('p5', 3, 8), ('p6', 4, 5), ('p7', 4, 6), ('p8', 5, 7), ('p9', 5, 8), ('p10', 6, 7), ('p11', 6, 8), ('p12', 7, 8), ('q2', 2, 9)] \
 --name=STF_5L_92 \
 --analytical_result=-1.66419336793957e-9 \
 --externals=['q1'] \
 --lmb=None
