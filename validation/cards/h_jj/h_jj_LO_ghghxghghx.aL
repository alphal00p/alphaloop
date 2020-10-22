import model aL_sm-no_widths
set_alphaLoop_option n_rust_inputs_to_generate -1
set_alphaLoop_option include_self_energies_from_squared_amplitudes True
set_alphaLoop_option differentiate_particle_from_antiparticle_in_graph_isomorphism False
set_alphaLoop_option consider_vertex_id_in_graph_isomorphism False
set_alphaLoop_option consider_edge_orientation_in_graph_isomorphism False
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option apply_graph_isomorphisms True
set_alphaLoop_option qgraf_template_model no_s 
set_alphaLoop_option loop_induced true
set_alphaLoop_option qgraf_cut_filter True
set_alphaLoop_option n_jets 0
set_alphaLoop_option final_state_pdgs (82,-82,82,-82)
set_FORM_option number_of_lmbs None
set_FORM_option generate_integrated_UV_CTs True 
qgraf_define j = d d~
qgraf_generate h > gh gh~ gh gh~ / d u s c b QCD^2==4 QED^2==1 []
!rm -rf h_jj_LO_ghghxghghx
output qgraf h_jj_LO_ghghxghghx
