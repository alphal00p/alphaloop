import model aL_sm-no_widths
# And finally let us also enable the FORM processing of numerators
set_alphaLoop_option FORM_processing_output_format c 
!rm -rf TEST_QGRAF_epem_a_ttx
qgraf_generate epem_ttx TEST_QGRAF_epem_a_ttx
