import model aL_sm-no_widths
set_alphaLoop_option n_rust_inputs_to_generate -1
set_alphaLoop_option include_self_energies_from_squared_amplitudes True
set_alphaLoop_option differentiate_particle_from_antiparticle_in_graph_isomorphism False
set_alphaLoop_option consider_vertex_id_in_graph_isomorphism False 
set_alphaLoop_option consider_edge_orientation_in_graph_isomorphism False
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option apply_graph_isomorphisms True
set_alphaLoop_option qgraf_template_model epem
qgraf_generate e+ e- > a > d d~ /t u s c b QCD^2==0 QED^2==2 []
!rm -rf TEST_QGRAF_epem_a_ddx_LO
output qgraf TEST_QGRAF_epem_a_ddx_LO
