import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_3L_1 
output_LU_scalar TEST_SCALAR_STF_3L_1 \
 --topology=\
 [('q1', 0, 1), ('p1', 1, 2), ('p2', 1, 5), ('p3', 2, 3), ('p4', 2, 4), ('p5', 3, 4), ('p6', 3, 5), ('p7', 4, 5), ('q2', 5, 6)] \
 --name=STF_3L_1 \
 --analytical_result=-5.75396168705020e-6 \
 --externals=['q1'] \
 --lmb=None
