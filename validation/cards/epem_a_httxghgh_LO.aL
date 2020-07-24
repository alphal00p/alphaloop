import model aL_sm-no_widths

set_alphaLoop_option n_rust_inputs_to_generate -1
set_alphaLoop_option differentiate_particle_from_antiparticle_in_graph_isomorphism False
set_alphaLoop_option consider_edge_orientation_in_graph_isomorphism False
set_alphaLoop_option consider_vertex_id_in_graph_isomorphism False 
set_alphaLoop_option apply_graph_isomorphisms True

# Memory limit (max 1T)
set_alphaLoop_option virtual_memory 100G
set_alphaLoop_option FORM_compile_cores 10

# FORM options
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_optimization 3 
set_alphaLoop_option FORM_compile_arg all 

# QGRAF model 
set_alphaLoop_option qgraf_template_model epem

# Process
qgraf_generate e+ e- > a > h t t~ gh gh~ /d u s c b QCD^2==1 QED^2==3 []
!rm -rf TEST_QGRAF_epem_a_httxghgh_LO
output qgraf TEST_QGRAF_epem_a_httxghgh_LO

