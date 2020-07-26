import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_5L_7 
output_LU_scalar TEST_SCALAR_STF_5L_7 \
 --topology=\
 [('q1', 0, 3), ('p1', 1, 2), ('p2', 1, 3), ('p3', 1, 4), ('p4', 2, 5), ('p5', 2, 6), ('p6', 3, 5), ('p7', 3, 6), ('p8', 4, 5), ('p9', 4, 6), ('p10', 5, 6), ('q2', 4, 7)] \
 --name=STF_5L_7 \
 --analytical_result=1.66419336793957e-9 \
 --externals=['q1'] \
 --lmb=None
