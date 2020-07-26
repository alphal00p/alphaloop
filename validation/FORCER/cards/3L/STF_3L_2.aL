import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 0
!rm -rf TEST_SCALAR_STF_3L_2 
output_LU_scalar TEST_SCALAR_STF_3L_2 \
 --topology=\
 [('q1', 0, 1), ('p1', 1, 3), ('p2', 1, 4), ('p3', 2, 3), ('p4', 2, 5), ('p5', 3, 6), ('p6', 4, 5), ('p7', 4, 6), ('p8', 5, 6), ('q2', 2, 7)] \
 --name=STF_3L_2 \
 --analytical_result=-5.75396168705020e-6 \
 --externals=['q1'] \
 --lmb=('p2', 'p4', 'p6')
