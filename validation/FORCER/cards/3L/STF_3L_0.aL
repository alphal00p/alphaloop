import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_3L_0 
output_LU_scalar TEST_SCALAR_STF_3L_0 \
 --topology=\
 [('q1', 0, 3), ('p1', 1, 2), ('p2', 1, 3), ('p3', 1, 4), ('p4', 2, 3), ('p5', 2, 4), ('p6', 3, 4), ('q2', 4, 5)] \
 --name=STF_3L_0 \
 --analytical_result=5.75396168705020e-6 \
 --externals=['q1'] \
 --lmb=None
