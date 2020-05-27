import model ./PLUGIN/alphaloop/models/aL_sm-no_widths
# And finally let us also enable the FORM processing of numerators
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option n_rust_inputs_to_generate -1
!rm -rf TEST_QGRAF_epem_a_httxgg_LO
qgraf_generate epem_a_ttxhgg_LO TEST_QGRAF_epem_a_httxgg_LO
