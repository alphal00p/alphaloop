import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_5L_4 
output_LU_scalar TEST_SCALAR_STF_5L_4 \
 --topology=\
 [('q1', 0, 3), ('p1', 1, 2), ('p2', 1, 5), ('p3', 1, 6), ('p4', 2, 5), ('p5', 2, 6), ('p6', 3, 4), ('p7', 3, 5), ('p8', 3, 6), ('p9', 4, 5), ('p10', 4, 6), ('q2', 4, 7)] \
 --name=STF_5L_4 \
 --analytical_result=1.66419336793957e-9 \
 --externals=['q1'] \
 --lmb=None
