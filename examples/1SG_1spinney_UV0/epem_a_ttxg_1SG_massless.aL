import model aL_sm_no_masses-no_widths
set_alphaLoop_option perturbative_orders {'QCD':0}
# Force including *all* cutkosky cuts!
set_alphaLoop_option n_jets 0
set_alphaLoop_option final_state_pdgs []
set_alphaLoop_option n_rust_inputs_to_generate -1
set_alphaLoop_option include_self_energies_from_squared_amplitudes True
set_alphaLoop_option differentiate_particle_from_antiparticle_in_graph_isomorphism False
set_alphaLoop_option consider_vertex_id_in_graph_isomorphism False 
set_alphaLoop_option consider_edge_orientation_in_graph_isomorphism False
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option apply_graph_isomorphisms True
set_FORM_option number_of_lmbs -1
set_FORM_option reference_lmb None
set_alphaLoop_option supergraphs_selection [(1,2),]
add process e+ e- > a > t t~ g @0
!rm -rf single_SG_paper_epem_a_ttxg_1SG_massless
output alphaLoop single_SG_paper_epem_a_ttxg_1SG_massless
