import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
!rm -rf TEST_SCALAR_mercedes 
output_LU_scalar TEST_SCALAR_mercedes \
 --topology=[\
	('q1', 0, 1), ('p1', 1, 2), ('p2', 2, 7), ('p3', 7, 3), ('p4', 3, 6),\
    ('p5', 6, 5), ('p6', 5, 1), ('p7', 2, 4), ('p8', 3, 4), ('p9', 4, 5),\
	('p10', 7, 4), ('q2', 6, 8) ] \
 --name="Mercedes" \
 --analytical_result=5.75396e-6 \
 --externals=('q1',) \
 --lmb=('p1', 'p2', 'p3', 'p4')
