import model sm-no_widths
set_alphaLoop_option perturbative_orders {'QCD':0}
set_alphaLoop_option n_rust_inputs_to_generate -1
set_alphaLoop_option include_self_energies_from_squared_amplitudes True
set_alphaLoop_option differentiate_particle_from_antiparticle_in_graph_isomorphism False 
set_alphaLoop_option FORM_processing_output_format c
generate e+ e- > a > h t t~ g g
!rm -rf TEST_epem_a_httxgg
output alphaLoop TEST_epem_a_httxgg
